! Conformal Cubic Atmospheric Model
    
! Copyright 2015-2020 Commonwealth Scientific Industrial Research Organisation (CSIRO)
    
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

! This module solves for the atmosphere Helmholtz equation and the ocean free
! surface equation using a SOR, conjugate gradient or multi-grid approach.

! Design notes:

! The solution to the Helmholtz equation is currently the limiting factor on
! the model scaling with increasing cores.  However, the mass-flux,
! split-explicit dynamical core should avoid the use of an implicit solution
! and hence this bottleneck should be avoided with the next generation of
! CCAM.
    
! Notes on the geometric multigrid method:
    
! Note that the solver can be greatly optimised for grids that permit many
! local subdivisions.  For example if the number of grid points on a processor
! is M x M, where M = 2^N and N is an integer greater than 1, then the
! multi-grid solver will be considerably more efficent as more processors can
! be applied in parallel.  Example grids that typically work well include
! C48, C96, C192, C384, C768 and C1536.  An alternate family of grids includes
! C32, C64, C128, C256, C512 and C1024.  In cases were a non-optimal grid
! size is specified, then we use a coloured SOR solver for the coarse grid
! so as to solve for the sparse matrix in a reasonable amount of walltime.
    
module helmsolve

implicit none

private
public helmsor
public helmsol
public mghelm,mgmlo,mgsor_init,mgzz_init

! Congujate gradient

integer, dimension(:,:), allocatable, private, save :: ileft
integer, dimension(:), allocatable, private, save   :: nleft, nright
real, dimension(:,:,:), allocatable, private, save  :: pleft
real, dimension(:,:), allocatable, private, save    :: ppinv

! Geometric multigrid

integer, save :: mg_maxsize, mg_minsize, gmax ! grid sizes for automatic arrays
integer, save :: mg_maxlevel_decomp           ! maximum level for shared memory
integer, save :: itr_mg    = 20               ! maximum number of iterations for atmosphere MG solver
integer, save :: itr_mgice = 20               ! maximum number of iterations for ocean/ice MG solver
integer, save :: itrbgn    = 2                ! number of iterations relaxing the solution after MG restriction
integer, save :: itrend    = 2                ! number of iterations relaxing the solution after MG interpolation
real, parameter :: dfac = 0.25                ! adjustment for grid spacing after MG restriction
real, parameter :: dfaci = 1.                 ! seaice adjustment for grid spacing after MG restriction
logical, save :: sorfirst = .true.            ! first call to mgsor_init
logical, save :: zzfirst  = .true.            ! first call to mgzz_init

contains

! ****************************************************************************
! This is the SOR solver from JLM
subroutine helmsor(zz,zzn,zze,zzw,zzs,helm,s,irhs)
!     Solve Helmholtz equation - experimental jlm version
!     e.g. use  -2540, -2335, -2635, -2425 (in decr. order)

use cc_mpi
use diag_m
use indices_m
use newmpar_m
use parm_m
use parmdyn_m
use parmgeom_m
use sumdd_m
use vecs_m

implicit none

integer, parameter :: ntest=0 
integer, parameter :: itmax=400 ! maximum number of iterations allowed
!     Arguments
real, dimension(ifull), intent(in) :: zz,zzn,zze,zzw,zzs
real, dimension(ifull,kl), intent(in) :: helm         ! Helmholtz coefficients
real, dimension(ifull+iextra,kl), intent(inout) :: s  ! Solution
real, dimension(ifull,kl), intent(in) :: irhs         ! RHS
real, dimension(ifull,kl) :: s_new
real, dimension(ifull,kl) :: rhs
real, dimension(ifull_maxcolour,kl,maxcolour) :: helmc,rhsc
real, dimension(ifull,kl) :: sb, sa, snew
real, dimension(ifull,kl) :: dsol
real, dimension(ifull_maxcolour,maxcolour) :: zznc, zzec
real, dimension(ifull_maxcolour,maxcolour) :: zzwc, zzsc
real, dimension(kl) ::  dsolmax, dsolmax_g, smax, smax_g
real, dimension(kl) ::  smin, smin_g, savg
real, dimension(:), allocatable, save :: accel
real, dimension(ifull) :: aa, bb, cc
real, save ::  axel
real, save :: dtsave = 0.
integer, save :: meth, nx_max
integer iq, iter, k, nx, klim, klimnew
integer isc, iec
integer, dimension(kl) :: iters

call START_LOG(helm_begin)
      
rhs=irhs ! allows subroutine to modify rhs

s(ifull+1:ifull+iextra,:)=0. ! For IBM compiler
s_new=0.                     ! For IBM compiler

if (abs(dt-dtsave)>=1.e-20) then ! dt/=dtsave
  dtsave=dt
  if (.not.allocated(accel)) then
    allocate(accel(kl))
  end if

  if(precon==-1)precon=-2325  ! i.e. 2, 3, .25
  nx_max=abs(precon)/1000
  meth=abs(precon)/100-10*nx_max
  axel=-.01*real(precon)-10*nx_max-meth
  if(myid==0)write(6,*)'in helmsor nx_max,meth,axel: ',nx_max,meth,axel
 
  if (nx_max/=maxcolour) then
    if (myid==0) then
      write(6,*) "WARN: mismatched number of colours"
      write(6,*) "changing nx_max ",nx_max,maxcolour
    end if
    nx_max=maxcolour
  end if
  
  if (il_g<=200) then
    ! usual
    do k=1,kl
      call optmx(il_g,schmidt,dt,bam(k),accel(k))
      ! MJT - not sure about the following line
      accel(k)=1.+.55*(accel(k)-1.) ! for uniform-dec 22/4/08
      if(myid==0)write(6,*)'k,accel ',k,accel(k)
    enddo
  else
    ! large grid
    accel=1. ! gauss-seidel
  end if
endif ! dt/=dtsave

if ( precon>=-2899 ) then  ! e.g. not -2900 or -3900
  klim = kl
  iter = 1
  sa = s(1:ifull,:)
  sb = s(1:ifull,:)
  do while ( iter<itmax .and. klim>1)
    call bounds(s, klim=klim)
    do k = 1,klim        
      do nx = 1,nx_max
        dsol(iqx(:,nx),k)=                                                       &
             ( zzn(iqx(:,nx))*s(iqn(:,nx),k) + zzw(iqx(:,nx))*s(iqw(:,nx),k)     &
             + zze(iqx(:,nx))*s(iqe(:,nx),k) + zzs(iqx(:,nx))*s(iqs(:,nx),k)     &
             + (zz(iqx(:,nx))-helm(iqx(:,nx),k))*s(iqx(:,nx),k)                  &
             - rhs(iqx(:,nx),k) )/(helm(iqx(:,nx),k)-zz(iqx(:,nx)))
        snew(iqx(:,nx),k) = s(iqx(:,nx),k) + dsol(iqx(:,nx),k)

!       following are jlm methods for improving guess
        if(iter>=3)then
          select case(meth)
            case(3)
              aa(iqx(:,nx))=(sb(iqx(:,nx),k)-3.*sa(iqx(:,nx),k)+3.*s(iqx(:,nx),k)+19.*snew(iqx(:,nx),k))/20.  
              bb(iqx(:,nx))=(9.*sb(iqx(:,nx),k)-17.*sa(iqx(:,nx),k)-13.*s(iqx(:,nx),k)+21.*snew(iqx(:,nx),k))/20.
              snew(iqx(:,nx),k)=aa(iqx(:,nx))+axel*bb(iqx(:,nx))
            case(4)   ! oscill
              aa(iqx(:,nx))=(7.*snew(iqx(:,nx),k)+3.*s(iqx(:,nx),k)-3.*sa(iqx(:,nx),k)+sb(iqx(:,nx),k))/8. ! oscill
              bb(iqx(:,nx))=snew(iqx(:,nx),k)-.5*s(iqx(:,nx),k)-sa(iqx(:,nx),k)+.5*sb(iqx(:,nx),k)         ! oscill
              snew(iqx(:,nx),k)=aa(iqx(:,nx))+axel*bb(iqx(:,nx))
            case(5)   ! wqls
              aa(iqx(:,nx))=(2.*sb(iqx(:,nx),k)-6.*sa(iqx(:,nx),k)+6*s(iqx(:,nx),k)+68.*snew(iqx(:,nx),k))/70.  
              bb(iqx(:,nx))=(5.*sb(iqx(:,nx),k)-8.*sa(iqx(:,nx),k)-13.*s(iqx(:,nx),k)+16.*snew(iqx(:,nx),k))/14.  
              snew(iqx(:,nx),k)=aa(iqx(:,nx))+axel*bb(iqx(:,nx))          
            case(6)   ! wqls again
              aa(iqx(:,nx))=(2.*sb(iqx(:,nx),k)-6.*sa(iqx(:,nx),k)+6*s(iqx(:,nx),k)+68.*snew(iqx(:,nx),k))/70.  
              bb(iqx(:,nx))=(5.*sb(iqx(:,nx),k)-8.*sa(iqx(:,nx),k)-13*s(iqx(:,nx),k)+16.*snew(iqx(:,nx),k))/14.  
              cc(iqx(:,nx))=(3.*sb(iqx(:,nx),k)-2.*sa(iqx(:,nx),k)-5*s(iqx(:,nx),k)+4.*snew(iqx(:,nx),k))/14.  
              snew(iqx(:,nx),k)=aa(iqx(:,nx))+axel*(bb(iqx(:,nx))+axel*cc(iqx(:,nx)))          
          end select
        end if

        sb(iqx(:,nx),k)=sa(iqx(:,nx),k)
        sa(iqx(:,nx),k)=s(iqx(:,nx),k)
        s(iqx(:,nx),k)=snew(iqx(:,nx),k)
      enddo  ! nx loop
      iters(k)=iter
    enddo ! k loop   
    if(ntest>0.and.meth>=1.and.mydiag) write (6,"('Iter ,s',i4,4f14.5)") iter,(s(iq,1),iq=1,4)
    if(iter==1)then
      do k=1,klim
        smax(k) = maxval(s(1:ifull,k))
        smin(k) = minval(s(1:ifull,k))
      end do
      call ccmpi_allreduce(smax(1:klim),smax_g(1:klim),"max",comm_world)
      call ccmpi_allreduce(smin(1:klim),smin_g(1:klim),"min",comm_world)
    endif
    do k=1,klim
      dsolmax(k) = maxval(abs(dsol(1:ifull,k)))
    enddo
    call ccmpi_allreduce(dsolmax(1:klim),dsolmax_g(1:klim),"max",comm_world)
    if(myid==0.and.ntest>0)then
      write(6,*)'smin_g ',smin_g(:)
      write(6,*)'smax_g ',smax_g(:)
      write(6,*)'dsolmax_g ',dsolmax_g(:)
    endif  ! (myid==0)
    klimnew=klim
    do k=klim,1,-1
      if(dsolmax_g(k)<restol*(smax_g(k)-smin_g(k)))then
        klimnew=k
      endif
    enddo
    klim=klimnew
    iter = iter + 1
  enddo   ! while( iter<itmax .and. klim>1)

  if(myid==0.and.(diag.or.ktau<6))then
    do k=1,kl
      write(6,*)'helmjlm ktau,k,Iterations ',ktau,k,iters(k)
    enddo
  endif

else ! e.g. -2900 or -3900

  klim = kl
  iter = 1
  do k = 1,klim
    smax(k) = maxval(s(1:ifull,k))
    smin(k) = minval(s(1:ifull,k))
  enddo
  if(ntest>0.and.diag)write(6,*)' before smax call myid ',myid
  call ccmpi_allreduce(smax(1:klim),smax_g(1:klim),"max",comm_world)
  if(ntest>0.and.diag)write(6,*)' before smin call myid ',myid
  call ccmpi_allreduce(smin(1:klim),smin_g(1:klim),"min",comm_world)
  if((ntest>0.or.nmaxpr==1).and.myid==0)then
    write(6,*)'ktau,smin_g ',ktau,smin_g(:)
    write(6,*)'ktau,smax_g ',ktau,smax_g(:)
  endif  ! (myid==0)

  ! JLM suggestion
  do k=1,kl
    savg(k)=0.5*(smax_g(k)+smin_g(k))
    s(1:ifull,k)=s(1:ifull,k)-savg(k)
    rhs(:,k)=rhs(:,k)+(helm(:,k)-zz-zzn-zzs-zze-zzw)*savg(k)
  end do

  do nx=1,nx_max
    zznc(1:ifull_colour(nx),nx) =zzn(iqx(1:ifull_colour(nx),nx))
    zzwc(1:ifull_colour(nx),nx) =zzw(iqx(1:ifull_colour(nx),nx))
    zzec(1:ifull_colour(nx),nx) =zze(iqx(1:ifull_colour(nx),nx))
    zzsc(1:ifull_colour(nx),nx) =zzs(iqx(1:ifull_colour(nx),nx))
    do k=1,kl
      helmc(1:ifull_colour(nx),k,nx)=helm(iqx(1:ifull_colour(nx),nx),k)-zz(iqx(1:ifull_colour(nx),nx))
      rhsc(1:ifull_colour(nx),k,nx) =rhs(iqx(1:ifull_colour(nx),nx),k)
    end do
  end do
     
  call bounds(s)     
  do while ( iter<itmax .and. klim>0)
    do nx=1,nx_max
      isc = 1
      iec = ifull_colour_border(nx)
      do k=1,klim
        dsol(iqx(isc:iec,nx),k)=                                                                &
            ( zznc(isc:iec,nx)*s(iqn(isc:iec,nx),k) + zzwc(isc:iec,nx)*s(iqw(isc:iec,nx),k)     &
            + zzec(isc:iec,nx)*s(iqe(isc:iec,nx),k) + zzsc(isc:iec,nx)*s(iqs(isc:iec,nx),k)     &
            -helmc(isc:iec,k,nx)*s(iqx(isc:iec,nx),k)-rhsc(isc:iec,k,nx) )/helmc(isc:iec,k,nx)
        s_new(iqx(isc:iec,nx),k) = s(iqx(isc:iec,nx),k) + accel(k)*dsol(iqx(isc:iec,nx),k)
      end do ! k loop
      call bounds_colour_send(s_new, nx, klim=klim)
      isc = ifull_colour_border(nx) + 1
      iec = ifull_colour(nx)
      do k=1,klim
        dsol(iqx(isc:iec,nx),k)=                                                               &
            ( zznc(isc:iec,nx)*s(iqn(isc:iec,nx),k) + zzwc(isc:iec,nx)*s(iqw(isc:iec,nx),k)    &
            + zzec(isc:iec,nx)*s(iqe(isc:iec,nx),k) + zzsc(isc:iec,nx)*s(iqs(isc:iec,nx),k)    &
            -helmc(isc:iec,k,nx)*s(iqx(isc:iec,nx),k)-rhsc(isc:iec,k,nx) )/helmc(isc:iec,k,nx)
        s(iqx(isc:iec,nx),k) = s(iqx(isc:iec,nx),k) + accel(k)*dsol(iqx(isc:iec,nx),k)
      end do ! k loop
      iec = ifull_colour_border(nx)
      s(iqx(1:iec,nx),1:klim) = s_new(iqx(1:iec,nx),1:klim)      
      call bounds_colour_recv(s, nx, klim=klim)
    end do  ! nx loop  
    do k=1,klim
      iters(k)=iter
    end do
      
       
    do k=1,klim
      dsolmax(k) = maxval(abs(dsol(1:ifull,k)))
    enddo  ! k loop
    klimnew=klim
    call ccmpi_allreduce(dsolmax(1:klim),dsolmax_g(1:klim),"max",comm_world)
    do k=klim,1,-1
      if(dsolmax_g(k)<restol*(smax_g(k)-smin_g(k)))then
        klimnew=k-1
      endif
    enddo
    klim=klimnew
     iter = iter + 1
  enddo   ! while( iter<itmax .and. klim>1)
      
  do k=1,kl
    s(1:ifull,k)=s(1:ifull,k)+savg(k)
  end do

  if(myid==0)then
    if(nmaxpr==1) then
      write(6,*)'helmjlm ktau,k,Iterations ',ktau,1,iters(1)
    end if
    if(diag.or.ktau<6.or.iters(1)>itmax-5)then
      do k=1,kl
        write(6,*)'helmjlm ktau,k,Iterations ',ktau,k,iters(k)
      enddo
    endif
  endif
      
end if ! precon>=-2899 ..else..
      
call END_LOG(helm_end)
      
return
end subroutine helmsor


! ****************************************************************************
! This is the conjugate gradient solver from MRD
subroutine helmsol(zz,zzn,zze,zzw,zzs,helm,s,rhs)

!     Solve Helmholtz equation using preconditioned conjugate gradient method.
!     Each mode is solved separately. Highest numbered modes 
!     converge fastest.
!     Note that the ILU preconditioner is expensive to set up. This code 
!     works best with ntbar=0 so the matrix doesn't change in time.

use cc_mpi
use indices_m
use newmpar_m
use parm_m
use parmdyn_m
use sumdd_m

implicit none

integer, parameter :: itmax=600 ! maximum number of iterations allowed
!     Arguments
real, intent(in), dimension(ifull) :: zz,zzn,zze,zzw,zzs
!     WHY are helm and rhs ifull+iextra?????????
!     Not just for printa call ?????
real, intent(in) :: helm(ifull,kl)             ! Helmholtz coefficients
real, intent(inout) :: s(ifull+iextra,kl)      ! Solution
real, intent(in) :: rhs(ifull,kl)              ! RHS
real, dimension(ifull,kl) :: fac, invfac, v, sx
real, dimension(ifull+iextra,kl) :: d, r, h

real, dimension(kl) :: alpha, beta, smag, &
                       gamma_1, sigma, delta
real, dimension(kl) :: gsmag, ggamma_0, &
                       ggamma_1, gsigma, gdelta
integer :: iq, iter, k, klim
logical, save :: ilustart = .true.
real, save :: factest
complex, dimension(3*kl) :: local_sum, global_sum
!     Temporary array for the drpdr_local function
real, dimension(ifull) :: tmparr, tmparr2 

call START_LOG(helm_begin)
      
klim = kl ! All modes at first

do k=1,klim
   fac(1:ifull,k) = zz(1:ifull) - helm(1:ifull,k)
   invfac(1:ifull,k) = 1.0/fac(1:ifull,k)
end do
if ( precon /= 0 ) then
   if ( ilustart) then
      call iludecomp(precon,fac(:,1:precon),zzn,zze,zzs,zzw)
      ilustart = .false.
!     Save this as a check on whether array has changed.
!     This will catch both changes in the time step and the reference
!     temperature profile
      factest = fac(1,1) 
   else
      if ( abs(factest-fac(1,1))>=1.e-20 ) then ! facttest /= fac(1,1)
         call iludecomp(precon,fac(:,1:precon),zzn,zze,zzs,zzw)
         factest = fac(1,1) 
      end if
   end if
end if

! Use D'Azevedo method
smag = 0.
local_sum = (0.,0.)
do k=1,klim
   if ( k <= precon) then
      do iq=1,ifull
         r(iq,k) = - ( zze(iq)*s(ie(iq),k) + zzw(iq)*s(iw(iq),k) + &
                       zzn(iq)*s(in(iq),k) + zzs(iq)*s(is(iq),k) - &
                       rhs(iq,k) ) - s(iq,k)*fac(iq,k)
         tmparr(iq) = s(iq,k)*s(iq,k)
      end do
      call ilusolve(h,r,k)
      d(1:ifull,k) = h(1:ifull,k)
   else
      do iq=1,ifull
         r(iq,k) = - ( zze(iq)*s(ie(iq),k) + zzw(iq)*s(iw(iq),k) + &
                       zzn(iq)*s(in(iq),k) + zzs(iq)*s(is(iq),k) - &
                       rhs(iq,k) ) * invfac(iq,k) - s(iq,k)
         d(iq,k) = r(iq,k)
         tmparr(iq) = s(iq,k)*s(iq,k)
      end do
   end if
   call drpdr_local(tmparr, local_sum(2*klim+k))
   smag(k)=real(local_sum(2*klim+k))
end do

call bounds(d, klim=klim)
sigma = 0.
gamma_1 = 0.
do k=1,klim
   if ( k <= precon ) then
      do iq=1,ifull
         v(iq,k) = zze(iq)*d(ie(iq),k) + zzw(iq)*d(iw(iq),k) + &
                   zzn(iq)*d(in(iq),k) + zzs(iq)*d(is(iq),k) + &
                    d(iq,k)*fac(iq,k) 
         tmparr(iq) = d(iq,k)*v(iq,k)
         tmparr2(iq) = r(iq,k)*h(iq,k)
      end do
   else
      do iq=1,ifull
         v(iq,k) = ( zze(iq)*d(ie(iq),k) + zzw(iq)*d(iw(iq),k) + &
                     zzn(iq)*d(in(iq),k) + zzs(iq)*d(is(iq),k) ) *  &
                     invfac(iq,k) + d(iq,k)
         tmparr(iq) = d(iq,k)*v(iq,k)
         tmparr2(iq) = r(iq,k)*r(iq,k)
      end do
   end if
   call drpdr_local(tmparr, local_sum(k))
   call drpdr_local(tmparr2, local_sum(klim+k))
   sigma(k)=real(local_sum(k))
   gamma_1(k)=real(local_sum(klim+k))
end do
call ccmpi_allreduce(local_sum(1:2*klim),global_sum(1:2*klim),"sumdr",comm_world)
gsigma(1:klim) = real(global_sum(1:klim))
ggamma_1(1:klim) = real(global_sum(klim+1:2*klim))
local_sum(1:2*klim) = (0.,0.) ! Still need the later part.
alpha(1:klim) = ggamma_1(1:klim) / gsigma(1:klim)
do k=1,klim
   s(1:ifull,k) = s(1:ifull,k) + alpha(k) * d(1:ifull,k)
   r(1:ifull,k) = r(1:ifull,k) - alpha(k) * v(1:ifull,k)
   if ( k <= precon) then
      call ilusolve(h,r,k)
   else
      ! Need this so that a single bounds call does both h and r
      h(1:ifull,k) = r(1:ifull,k)
   end if
end do
ggamma_0(1:klim) = ggamma_1(1:klim)
      
select case (helmmeth)
   case(0) ! D'Azevedo Method 1
         
      do iter=2,itmax
         call bounds(h, klim=klim)
         gamma_1(1:klim) = 0.0
         delta(1:klim) = 0.0
         do k=1,klim
            if ( k <= precon ) then
               do iq=1,ifull
                  sx(iq,k) = zze(iq)*h(ie(iq),k) + zzw(iq)*h(iw(iq),k) + &
                             zzn(iq)*h(in(iq),k) + zzs(iq)*h(is(iq),k) + &
                             h(iq,k)*fac(iq,k)
                  tmparr(iq) = r(iq,k)*h(iq,k)
                  tmparr2(iq) = sx(iq,k)*h(iq,k)
               end do
            else
               do iq=1,ifull
                  ! Here h is a copy of r
                  sx(iq,k) = ( zze(iq)*h(ie(iq),k) + zzw(iq)*h(iw(iq),k) + &
                               zzn(iq)*h(in(iq),k) + zzs(iq)*h(is(iq),k) ) *&
                               invfac(iq,k) + h(iq,k)
                  tmparr(iq) = h(iq,k)*h(iq,k)
                  tmparr2(iq) = sx(iq,k)*h(iq,k)
               end do
            end if
            call drpdr_local(tmparr, local_sum(k))
            call drpdr_local(tmparr2, local_sum(klim+k))
            gamma_1(k)=real(local_sum(k))
            delta(k)=real(local_sum(klim+k))
         end do
            
         call ccmpi_allreduce(local_sum(1:3*klim),global_sum(1:3*klim),"sumdr",comm_world)
         ggamma_1(1:klim) = real(global_sum(1:klim))
         gdelta(1:klim) = real(global_sum(klim+1:2*klim))
         gsmag(1:klim) = real(global_sum(2*klim+1:3*klim))
         local_sum = (0.,0.)
         if ( (diag .or. ktau<6 .or. itmax-iter.lt.50) .and. myid == 0 ) then
            write(6,'("Iterations",i4,i3,6g13.6/(10x,6g13.6))')   &
     &                 iter, klim, sqrt(abs(ggamma_1(1:klim)))
         end if
         !  Check which modes have converged
         do k=klim,1,-1
            if ( sqrt(abs(ggamma_1(k))) >= restol*sqrt(gsmag(k)) ) exit
         end do
!        Now k is the lowest mode yet to converge
         klim = k
         if ( klim == 0 ) exit

         beta(1:klim) = ggamma_1(1:klim) / ggamma_0(1:klim)
         ggamma_0(1:klim) = ggamma_1(1:klim)
         gsigma(1:klim) = gdelta(1:klim) - beta(1:klim)**2*gsigma(1:klim)
         alpha(1:klim) = ggamma_1(1:klim) / gsigma(1:klim)
         smag(1:klim) = 0.
         do k=1,klim
            if ( k <= precon) then
               do iq=1,ifull
                  d(iq,k) = h(iq,k) + beta(k) * d(iq,k)
                  v(iq,k) = sx(iq,k) + beta(k) * v(iq,k)
                  s(iq,k) = s(iq,k) + alpha(k) * d(iq,k)
                  r(iq,k) = r(iq,k) - alpha(k) * v(iq,k)
                  tmparr(iq) = s(iq,k)*s(iq,k)
               end do
               call ilusolve(h,r,k)
            else
              do iq=1,ifull
                  ! Use h in place of r
                  d(iq,k) = h(iq,k) + beta(k) * d(iq,k)
                  v(iq,k) = sx(iq,k) + beta(k) * v(iq,k)
                  s(iq,k) = s(iq,k) + alpha(k) * d(iq,k)
                  h(iq,k) = h(iq,k) - alpha(k) * v(iq,k)
                  tmparr(iq) = s(iq,k)*s(iq,k)
               end do
            end if
            call drpdr_local(tmparr, local_sum(2*klim+k))
            smag(k)=real(local_sum(2*klim+k))
         end do

      end do
      
      
   case(1) ! D'Azevedo Method standard
      
      local_sum(klim+1:2*klim)=local_sum(2*klim+1:3*klim)
      do iter=2,itmax
         gamma_1(1:klim) = 0.0
         do k=1,klim
            if ( k <= precon ) then
               do iq=1,ifull
                  tmparr(iq) = r(iq,k)*h(iq,k)
               end do
            else
               do iq=1,ifull
                  ! Here h is a copy of r
                  tmparr(iq) = h(iq,k)*h(iq,k)
               end do
            end if
            call drpdr_local(tmparr, local_sum(k))
            gamma_1(k)=real(local_sum(k))
         end do
            
         call ccmpi_allreduce(local_sum(1:2*klim),global_sum(1:2*klim),"sumdr",comm_world) 
         ggamma_1(1:klim) = real(global_sum(1:klim))
         gsmag(1:klim) = real(global_sum(klim+1:2*klim))
         local_sum = (0.,0.)
         if ( (diag .or. ktau<6 .or. itmax-iter<50) .and. myid == 0 ) then
            write(6,'("Iterations",i4,i3,6g13.6/(10x,6g13.6))')   &
     &                 iter, klim, sqrt(abs(ggamma_1(1:klim)))
         end if
         !  Check which modes have converged
         do k=klim,1,-1
            if ( sqrt(abs(ggamma_1(k))) >= restol*sqrt(gsmag(k)) ) exit
         end do
!        Now k is the lowest mode yet to converge
         klim = k
         if ( klim == 0 ) exit

         beta(1:klim) = ggamma_1(1:klim) / ggamma_0(1:klim)
         ggamma_0(1:klim) = ggamma_1(1:klim)
         
         sigma(1:klim)=0.
         do k=1,klim
            if ( k <= precon) then
               do iq=1,ifull
                  d(iq,k) = h(iq,k) + beta(k) * d(iq,k)
               end do
            else
               do iq=1,ifull
                  ! Use h in place of r
                  d(iq,k) = h(iq,k) + beta(k) * d(iq,k)
               end do
            end if
         end do
         call bounds(d, klim=klim)
         do k=1,klim
            if ( k <= precon) then
               do iq=1,ifull
                  v(iq,k) = zze(iq)*d(ie(iq),k) + zzw(iq)*d(iw(iq),k) + &
                            zzn(iq)*d(in(iq),k) + zzs(iq)*d(is(iq),k) + &
                            d(iq,k)*fac(iq,k)
                  tmparr(iq) = d(iq,k)*v(iq,k)
               end do
            else
               do iq=1,ifull
                  v(iq,k) = ( zze(iq)*d(ie(iq),k) + zzw(iq)*d(iw(iq),k) + &
                              zzn(iq)*d(in(iq),k) + zzs(iq)*d(is(iq),k) ) *  &
                              invfac(iq,k) + d(iq,k)
                  tmparr(iq) = d(iq,k)*v(iq,k)
               end do
            end if
            call drpdr_local(tmparr, local_sum(k))
            sigma(k)=real(local_sum(k))
         end do     

         call ccmpi_allreduce(local_sum(1:klim),global_sum(1:klim),"sumdr",comm_world)
         gsigma(1:klim) = real(global_sum(1:klim))
         local_sum(1:klim) = (0.,0.) ! Still need the later part.
         
         alpha(1:klim) = ggamma_1(1:klim) / gsigma(1:klim)
         smag(1:klim) = 0.
         do k=1,klim
            if ( k <= precon) then
               do iq=1,ifull
                  s(iq,k) = s(iq,k) + alpha(k) * d(iq,k)
                  r(iq,k) = r(iq,k) - alpha(k) * v(iq,k)
                  tmparr(iq) = s(iq,k)*s(iq,k)
               end do
               call ilusolve(h,r,k)
            else
               do iq=1,ifull
                  ! Use h in place of r
                  s(iq,k) = s(iq,k) + alpha(k) * d(iq,k)
                  h(iq,k) = h(iq,k) - alpha(k) * v(iq,k)
                  tmparr(iq) = s(iq,k)*s(iq,k)
               end do
            end if
            call drpdr_local(tmparr, local_sum(klim+k))
            smag(k)=real(local_sum(klim+k))
         end do

      end do      
      
end select

call END_LOG(helm_end)

end subroutine helmsol

! Routines for working with the sparse form of the incomplete LU
! decomposition.
! Values outside the processors local region are just ignored (no attempt
! to rebalance sum)
subroutine iludecomp(ilumax,fac,zzn,zze,zzs,zzw)

use cc_mpi
use indices_m
use newmpar_m

implicit none

integer, intent(in) :: ilumax
real, dimension(ifull), intent(in)  :: zzn, zze, zzs, zzw
real, dimension(ifull,ilumax), intent(in)  :: fac
real, dimension(ifull,4,ilumax) :: pp4
real, dimension(ifull,ilumax) :: pp
integer, dimension(ifull,4) :: ii4
integer :: i, j, k, ir, itmp, jtmp, idir, jdir, ni, nj
real, dimension(ilumax) :: e, dval
integer, dimension(4) :: rlist, list2
      
if (.not.allocated(ileft)) then
   allocate(ileft(ifull,4))
   allocate(nleft(ifull),nright(ifull))
   allocate(pleft(ifull,4,kl))
   allocate(ppinv(ifull,kl))
end if


! Need thse arrays so can index over the various cases
ii4(:,1) = in
ii4(:,2) = ie
ii4(:,3) = is
ii4(:,4) = iw
pp = fac
do k=1,ilumax
   pp4(:,1,k) = zzn
   pp4(:,2,k) = zze
   pp4(:,3,k) = zzs
   pp4(:,4,k) = zzw
end do

do ir=1,ifull
   dval = 1.0/pp(ir,:)
   ! i = r+1 .. n
   ! Look at all the neighbours. If i_x(r) > r, then this is a point
   ! below too by matrix symmetry.
   ! Use <= ifull check to get only points in this processor's region
   ni = count(ii4(ir,:)>ir .and. ii4(ir,:)<=ifull)
   rlist(1:ni) = pack(ii4(ir,:), ii4(ir,:)>ir .and. ii4(ir,:)<=ifull)
   do itmp=1,ni
      i = rlist(itmp)
      ! By the symmetry a(i,r) is non zero. To get value need to know
      ! for which direction x is i_x(i) = r. Direction index idir
      ! Use sum trick to get a scalar.
      idir = sum(pack( (/1, 2, 3, 4 /), ii4(i,:) == ir ))
      e = pp4(i,idir,:)*dval
      pp4(i,idir,:) = e
            
      ! j = r+1 .. n

      ! Point on diagonal certainly meets a[i,j] != 0, but is there
      ! another point in the same row?
      j = i
      if ( any ( ii4(ir,:) == j ) ) then
         jdir = sum(pack( (/1, 2, 3, 4 /), ii4(ir,:) == j ))
         pp(j,:) = pp(j,:) - e*pp4(ir,jdir,:)
      end if

      ! Loop checking other points for which j>r
      nj = count(ii4(i,:)>ir .and. ii4(ir,:)<=ifull)
      list2(1:nj) = pack(ii4(i,:), ii4(i,:)>ir .and. ii4(ir,:)<=ifull)
      do jtmp=1,nj ! points for which a(i,j) is non zero
         j = list2(jtmp)
         ! Check if a(r,j) is non-zero
         if ( any ( ii4(ir,:) == j ) ) then
            jdir = sum(pack( (/1, 2, 3, 4 /), ii4(ir,:) == j ))
            pp(j,:) = pp(j,:) - e*pp4(ir,jdir,:)
         end if
      end do
   end do
end do

nleft = 0
do i=1,ifull
   do idir=1,4
      ! Pack indices and values that are to the left of the diagonal
      if ( ii4(i,idir) < i ) then
         nleft(i) = nleft(i) + 1
         ileft(i,nleft(i)) = ii4(i,idir)
         pleft(i,nleft(i),1:ilumax) = pp4(i,idir,:)
      end if
   end do
   nright(i) = nleft(i)
   do idir=1,4
      ! Pack indices and values that are to the right of the diagonal
      if ( ii4(i,idir) > i .and. ii4(i,idir)<= ifull ) then
         nright(i) = nright(i) + 1
         ileft(i,nright(i)) = ii4(i,idir)
         pleft(i,nright(i),1:ilumax) = pp4(i,idir,:)
      end if
   end do
end do
ppinv(:,1:ilumax) = 1./pp
      
end subroutine iludecomp
   
subroutine ilusolve(x,rhs,k)

use cc_mpi
use indices_m
use newmpar_m

implicit none

real, dimension(:,:), intent(in)  :: rhs
real, dimension(:,:), intent(out) :: x
integer, intent(in) :: k
real, dimension(ifull) :: y
integer :: i, itmp
real :: tmpsum

call START_LOG(precon_begin)

! Solve Cx=rhs, where C is in sparse LU form

! ifull = size(rhs)
! First solve L y = rhs
y(1) = rhs(1,k)
do i=2,ifull
   ! Sum over points to left of diagonal
   tmpsum = 0.0
   do itmp=1,nleft(i)
      tmpsum = tmpsum+pleft(i,itmp,k)*y(ileft(i,itmp))
   end do
   y(i) = rhs(i,k) - tmpsum
end do
      
! Now solve Ux = y by backsubstitution
x(ifull,k) = y(ifull)*ppinv(ifull,k)
do i=ifull-1,1,-1
   tmpsum = 0.0
   do itmp=nleft(i)+1,nright(i)
      tmpsum = tmpsum+pleft(i,itmp,k)*x(ileft(i,itmp),k)
   end do
   x(i,k) = (y(i) -tmpsum) * ppinv(i,k)
end do

call END_LOG(precon_end)
      
end subroutine ilusolve


! ****************************************************************************    
! This version of the geometric multigrid solver is for the atmosphere
subroutine mghelm(izz,izzn,izze,izzw,izzs,ihelm,iv,jrhs)

use cc_mpi
use indices_m
use newmpar_m
use parm_m
use parmdyn_m

implicit none

integer, dimension(kl) :: iters
integer itr, ng, ng4, g, k, jj, i, iq
integer knew, klim, ir, ic, nc, n, iq_a, iq_c
integer isc, iec, klimc, itrc
real, dimension(ifull+iextra,kl), intent(inout) :: iv
real, dimension(ifull+iextra) :: vdum
real, dimension(ifull,kl), intent(in) :: ihelm, jrhs
real, dimension(ifull,kl) :: iv_new, iv_old, irhs
real, dimension(ifull), intent(in) :: izz, izzn, izze, izzw, izzs
real, dimension(ifull) :: dummy2, iv_n, iv_s, iv_e, iv_w
real, dimension(ifull_maxcolour,kl,maxcolour) :: rhelmc, rhsc
real, dimension(ifull_maxcolour,maxcolour) :: zznc, zzec, zzwc, zzsc
real, dimension(ifull_maxcolour) :: xdum
real, dimension(mg_ifull_maxcolour,3,kl) :: helmc_c, rhsc_c
real, dimension(mg_ifull_maxcolour,3) :: zznc_c, zzec_c, zzwc_c, zzsc_c
real, dimension(mg_maxsize,2*kl,2:gmax+1) :: rhs
real, dimension(mg_maxsize,kl,gmax+1) :: v, helm
real, dimension(mg_maxsize,2*kl) :: w
real, dimension(mg_maxsize) :: v_n, v_e, v_w, v_s
real, dimension(mg_minsize) :: vsavc
real, dimension(2*kl,2) :: smaxmin_g
real, dimension(kl) :: dsolmax_g, savg, sdif, dsolmaxc, sdifc

if ( sorfirst .or. zzfirst ) then
  write(6,*) "ERROR: mghelm requires mgsor_init and mgzz_init to be called first"
  call ccmpi_abort(-1)
end if

! zz*(DIV^2 v) - helm*v = rhs
! zz*(DIV^2 vd+v0) - helm*(vd+v0) = rhs
! zz*(DIV^2 vd) - helm*vd = rhs + helm*v0 - zz*(DIV^2 v0)

! MJT notes - to remove an additional syncronisation step we update rhs and helm
! parameters during the first iteration of the solution.  Effectively the mgsetup
! stage also becomes the first iteration of the solution.

call START_LOG(helm_begin)

call START_LOG(mgsetup_begin)

ng  = 0
ng4 = 0
klim = kl

do k = 1,kl

  iv(ifull+1:ifull+iextra,k) = 0. ! for IBM compiler
  iv_new(:,k) = 0.                ! for IBM compiler
  vdum(:) = 0.

  ! determine max/min for convergence calculations
  smaxmin_g(k,1) = maxval(iv(1:ifull,k))
  smaxmin_g(k,2) = minval(iv(1:ifull,k))
  smaxmin_g(k+kl,1) = 0.
  smaxmin_g(k+kl,2) = 0. 
  
  ! pack colour arrays at fine level
  ! note that the packing reorders the calculation to update the border points first
  dummy2(1:ifull) = 1./(ihelm(1:ifull,k)-izz(1:ifull))
  do nc = 1,maxcolour
    do iq = 1,ifull_colour(nc)  
      rhelmc(iq,k,nc) = dummy2(iqx(iq,nc))
      rhsc(iq,k,nc)   = jrhs(iqx(iq,nc),k)
    end do  
  end do
  
end do

do nc = 1,maxcolour
  do iq = 1,ifull_colour(nc)
    zznc(iq,nc) = izzn(iqx(iq,nc))
    zzwc(iq,nc) = izzw(iqx(iq,nc))
    zzec(iq,nc) = izze(iqx(iq,nc))
    zzsc(iq,nc) = izzs(iqx(iq,nc))
  end do  
end do

! solver assumes boundaries have been updated
call bounds(iv)

! Before sending convegence testing data in smaxmin_g and ihelm weights, we perform one iteration of the solver
! that can be updated with the smaxmin_g and ihelm arrays


! update on model grid using colours
do i = 1,itrbgn
  do nc = 1,maxcolour
    ! first calculate border points and send out halo
    isc = 1
    iec = ifull_colour_border(nc)
    do k = 1,kl
      do iq = isc,iec  
        iv_new(iqx(iq,nc),k) = ( zznc(iq,nc)*iv(iqn(iq,nc),k)     &
                               + zzwc(iq,nc)*iv(iqw(iq,nc),k)     &
                               + zzec(iq,nc)*iv(iqe(iq,nc),k)     &
                               + zzsc(iq,nc)*iv(iqs(iq,nc),k)     &
                               - rhsc(iq,k,nc) )*rhelmc(iq,k,nc)
      end do  
    end do
    call bounds_colour_send(iv_new,nc)
    ! calcuate non-border points while waiting for halo to update
    isc = ifull_colour_border(nc) + 1
    iec = ifull_colour(nc)
    do k = 1,kl
      do iq = isc,iec  
        xdum(iq) = ( zznc(iq,nc)*iv(iqn(iq,nc),k)                 &
                   + zzwc(iq,nc)*iv(iqw(iq,nc),k)                 &
                   + zzec(iq,nc)*iv(iqe(iq,nc),k)                 &
                   + zzsc(iq,nc)*iv(iqs(iq,nc),k)                 &
                   - rhsc(iq,k,nc) )*rhelmc(iq,k,nc)
      end do
      do iq = 1,isc-1
        iv(iqx(iq,nc),k) = iv_new(iqx(iq,nc),k)
      end do  
      do iq = isc,iec
        iv(iqx(iq,nc),k) = xdum(iq)
      end do
    end do
    call bounds_colour_recv(iv,nc)
  end do
end do

! calculate residual
do k = 1,kl
  call unpack_nsew(iv(:,k),iv_n,iv_s,iv_e,iv_w)
  w(1:ifull,k) = -izzn(:)*iv_n(:) - izzw(:)*iv_w(:) - izze(:)*iv_e(:) - izzs(:)*iv_s(:) &
                + jrhs(:,k) + iv(1:ifull,k)*(ihelm(:,k)-izz(:))
  ! also include ihelm weights for upscaled grid
  w(1:ifull,k+kl) = ihelm(1:ifull,k)
end do

! For when the inital grid cannot be upscaled - note helm and smaxmin_g are also included
call mgcollect(1,w(:,1:2*kl),smaxmin_g(1:2*kl,1:2))

if ( mg_maxlevel_local>0 ) then

  ! restriction
  ! (since this always operates within a panel, then ine = ien is always true)
  ng4 = mg(1)%ifull_fine
  do k = 1,2*kl
    do iq = 1,ng4
      rhs(iq,k,2) = 0.25*(w(mg(1)%fine(iq)  ,k) + w(mg(1)%fine_n(iq) ,k)  &
                        + w(mg(1)%fine_e(iq),k) + w(mg(1)%fine_ne(iq),k))
    end do  
  end do  

  ! merge grids if insufficent points on this processor - note helm and smaxmin_g are also included
  call mgcollect(2,rhs(:,1:2*kl,2),smaxmin_g(1:2*kl,1:2))

  if ( 2<=mg_maxlevel_local ) then
    helm(1:mg(2)%ifull,1:kl,2) = rhs(1:mg(2)%ifull,kl+1:2*kl,2)
  end if  
  
  
  ! upscale grid
  do g = 2,gmax
  
    ng = mg(g)%ifull
                
    ! update scalar field
    ! assume zero for first guess of residual (also avoids additional bounds call)
    do k = 1,kl
      v(1:ng,k,g) = -rhs(1:ng,k,g)/(helm(1:ng,k,g)-mg(g)%zz(1:ng))
    end do
    call mgbounds(g,v(:,1:kl,g))
    
    ! MJT notes - some groups use colours for the following smoothing, due to better convegence.
    do i = 2,itrbgn
      do k = 1,kl
        ! post smoothing
        call mgunpack_nsew(g,v(:,k,g),v_n,v_s,v_e,v_w)  
        v(1:ng,k,g) = (mg(g)%zze(1:ng)*v_e(1:ng) + mg(g)%zzw(1:ng)*v_w(1:ng) &
                     + mg(g)%zzn(1:ng)*v_n(1:ng) + mg(g)%zzs(1:ng)*v_s(1:ng) &
                     - rhs(1:ng,k,g))/(helm(1:ng,k,g)-mg(g)%zz(1:ng))
      end do
      call mgbounds(g,v(:,1:kl,g))
    end do
  
    ng4 = mg(g)%ifull_fine
    do k = 1,kl
      ! residual
      call mgunpack_nsew(g,v(:,k,g),v_n,v_s,v_e,v_w)  
      w(1:ng,k) = -mg(g)%zze(1:ng)*v_e(1:ng) - mg(g)%zzw(1:ng)*v_w(1:ng)   &
                 - mg(g)%zzn(1:ng)*v_n(1:ng) - mg(g)%zzs(1:ng)*v_s(1:ng)   &
                 + rhs(1:ng,k,g)+(helm(1:ng,k,g)-mg(g)%zz(1:ng))*v(1:ng,k,g)
      
      do iq = 1,ng4
        ! restriction
        rhs(iq,k,g+1) = 0.25*(w(mg(g)%fine(iq)  ,k) + w(mg(g)%fine_n(iq) ,k)  &
                            + w(mg(g)%fine_e(iq),k) + w(mg(g)%fine_ne(iq),k))
        ! restriction helm weights
        rhs(iq,k+kl,g+1) = 0.25*(helm(mg(g)%fine(iq)  ,k,g) + helm(mg(g)%fine_n(iq) ,k,g)   &
                               + helm(mg(g)%fine_e(iq),k,g) + helm(mg(g)%fine_ne(iq),k,g))
      end do  
    end do

    ! merge grids if insufficent points on this processor - note helm and smaxmin_g are also included
    call mgcollect(g+1,rhs(:,1:2*kl,g+1),smaxmin_g(1:2*kl,1:2))

    if ( g+1<=mg_maxlevel_local ) then
      helm(1:mg(g+1)%ifull,1:kl,g+1) = rhs(1:mg(g+1)%ifull,kl+1:2*kl,g+1)
    end if  

  end do


  ! solve for coarse grid
  ! MJT notes - use iterative method in case the coarse grid is still large (e.g., C25)
  if ( mg_maxlevel_local==mg_maxlevel ) then
    g = mg_maxlevel  
    ng = mg(g)%ifull
    do k = 1,kl
      do nc = 1,3
        do iq = 1,mg_ifull_maxcolour
          helmc_c(iq,nc,k) = helm(col_iq(iq,nc),k,g) - mg(g)%zz(col_iq(iq,nc))
          rhsc_c(iq,nc,k) = rhs(col_iq(iq,nc),k,g)
        end do  
      end do
      v(1:ng,k,g) = -rhs(1:ng,k,g)/(helm(1:ng,k,g)-mg(g)%zz(1:ng))
    end do
    do nc = 1,3
      do iq = 1,mg_ifull_maxcolour  
        zznc_c(iq,nc) = mg(g)%zzn(col_iq(iq,nc))
        zzec_c(iq,nc) = mg(g)%zze(col_iq(iq,nc))
        zzsc_c(iq,nc) = mg(g)%zzs(col_iq(iq,nc))
        zzwc_c(iq,nc) = mg(g)%zzw(col_iq(iq,nc))
      end do  
    end do
    sdifc(1:kl) = max( maxval(v(1:ng,1:kl,g),dim=1) - minval(v(1:ng,1:kl,g),dim=1), 1.e-20 )
    klimc = kl
    do itrc = 1,itr_mg
      do k = 1,klimc
        vsavc(1:ng) = v(1:ng,k,g)
        do nc = 1,3
          do iq = 1,mg_ifull_maxcolour  
            v(col_iq(iq,nc),k,g) = ( zznc_c(iq,nc)*v(col_iqn(iq,nc),k,g) &
                                   + zzec_c(iq,nc)*v(col_iqe(iq,nc),k,g) &
                                   + zzsc_c(iq,nc)*v(col_iqs(iq,nc),k,g) &
                                   + zzwc_c(iq,nc)*v(col_iqw(iq,nc),k,g) &
                                   - rhsc_c(iq,nc,k) )/helmc_c(iq,nc,k)
          end do  
        end do
        dsolmaxc(k) = maxval( abs( v(1:ng,k,g) - vsavc(1:ng) ) )
      end do
      ! test for convergence
      knew = klimc
      do k = klimc,1,-1
        if ( dsolmaxc(k)>=restol*sdifc(k) ) exit
        knew = k - 1
      end do
      klimc = knew
      if ( klimc<1 ) exit
    end do
  end if

  ! downscale grid
  do g = gmax,2,-1

    ! send coarse grid solution to processors and also bcast global smaxmin_g
    call mgbcastxn(g+1,v(:,1:kl,g+1),smaxmin_g(:,1:2))

    ng = mg(g)%ifull
    
    do k = 1,kl
      ! interpolation
      do iq = 1,ng  
        w(iq,k) = 0.5*v(mg(g+1)%coarse_a(iq),k,g+1)  + 0.25*v(mg(g+1)%coarse_b(iq),k,g+1) &
                + 0.25*v(mg(g+1)%coarse_c(iq),k,g+1)
      end do  

      ! extension
      ! No mgbounds as the v halo has already been updated and
      ! the coarse interpolation also updates the w halo
      v(1:ng,k,g) = v(1:ng,k,g) + w(1:ng,k)
    end do

    do i = 1,itrend
      call mgbounds(g,v(:,1:kl,g))  
      do k = 1,kl
        ! post smoothing
        call mgunpack_nsew(g,v(:,k,g),v_n,v_s,v_e,v_w)  
        v(1:ng,k,g) = (mg(g)%zze(1:ng)*v_e(1:ng) + mg(g)%zzw(1:ng)*v_w(1:ng) &
                     + mg(g)%zzn(1:ng)*v_n(1:ng) + mg(g)%zzs(1:ng)*v_s(1:ng) &
                     - rhs(1:ng,k,g))/(helm(1:ng,k,g)-mg(g)%zz(1:ng))
      end do
    end do
    
    call mgbounds(g,v(:,1:kl,g)) ! for next mgbcast

  end do

 
  ! broadcast coarse solution to fine grid, as well as global smaxmin_g
  call mgbcastxn(2,v(:,1:kl,2),smaxmin_g(:,1:2))

  ! interpolation
  ng = mg(1)%ifull
  do k = 1,kl
    do iq = 1,ng 
      w(iq,k) = 0.5*v(mg(2)%coarse_a(iq),k,2)  &
              + 0.25*v(mg(2)%coarse_b(iq),k,2) &
              + 0.25*v(mg(2)%coarse_c(iq),k,2)
    end do  
  end do
  
end if


! multi-grid solver bounds indices do not match standard iextra indices, so we need to remap the halo
if ( mg(1)%merge_len>1 ) then
  call mgbcastxn(1,w(:,1:kl),smaxmin_g(:,:),nobounds=.true.)
  ir = mod(mg(1)%merge_pos-1, mg(1)%merge_row) + 1   ! index for proc row
  ic = (mg(1)%merge_pos-1)/mg(1)%merge_row + 1       ! index for proc col
  do k = 1,kl
    do n = 1,npan
      do jj = 1,jpan
        iq_a = 1 + (jj-1)*ipan + (n-1)*ipan*jpan
        iq_c = 1 + (ir-1)*ipan + (jj-1+(ic-1)*jpan)*ipan*mg(1)%merge_row + (n-1)*ipan*jpan*mg(1)%merge_len
        vdum(iq_a:iq_a+ipan-1) = w(iq_c:iq_c+ipan-1,k)
      end do
    end do  
    ! extension
    iv(1:ifull,k) = iv(1:ifull,k) + vdum(1:ifull)
  end do
else
  ! remap mg halo to normal halo
  do k = 1,kl
    ! extension
    iv(1:ifull,k) = iv(1:ifull,k) + w(1:ifull,k)
  end do
end if

call bounds(iv)

! post smoothing
do i = 1,itrend
  do nc = 1,maxcolour
    ! update boundary grid points and send halo
    isc = 1
    iec = ifull_colour_border(nc)
    do k = 1,kl
      do iq = isc,iec  
        iv_new(iqx(iq,nc),k) = ( zznc(iq,nc)*iv(iqn(iq,nc),k)      &
                               + zzwc(iq,nc)*iv(iqw(iq,nc),k)      &
                               + zzec(iq,nc)*iv(iqe(iq,nc),k)      &
                               + zzsc(iq,nc)*iv(iqs(iq,nc),k)      &
                               - rhsc(iq,k,nc) )*rhelmc(iq,k,nc)
      end do  
    end do
    call bounds_colour_send(iv_new,nc)
    ! calculate non-boundary grid points while waiting for the halo to be updated
    isc = ifull_colour_border(nc) + 1
    iec = ifull_colour(nc)
    do k = 1,kl
      do iq = isc,iec  
        xdum(iq) = ( zznc(iq,nc)*iv(iqn(iq,nc),k)      &
                   + zzwc(iq,nc)*iv(iqw(iq,nc),k)      &
                   + zzec(iq,nc)*iv(iqe(iq,nc),k)      &
                   + zzsc(iq,nc)*iv(iqs(iq,nc),k)      &
                   - rhsc(iq,k,nc) )*rhelmc(iq,k,nc)
      end do  
      do iq = 1,isc-1
        iv(iqx(iq,nc),k) = iv_new(iqx(iq,nc),k)
      end do  
      do iq = isc,iec
        iv(iqx(iq,nc),k) = xdum(iq)
      end do  
    end do
    call bounds_colour_recv(iv,nc)
  end do
end do

do k = 1,kl

  ! remove offsets
  savg(k) = 0.5*(smaxmin_g(k,1)+smaxmin_g(k,2))
  sdif(k) = smaxmin_g(k,1) - smaxmin_g(k,2)

  iv(:,k) = iv(:,k) - savg(k)
  irhs(:,k) = jrhs(:,k) + (ihelm(:,k)-izz(:)-izzn(:)-izzs(:)-izze(:)-izzw(:))*savg(k)
  ! re-pack colour arrays at fine level to remove offsets
  do nc = 1,maxcolour
    do iq = 1,ifull_colour(nc)  
      rhsc(iq,k,nc) = irhs(iqx(iq,nc),k)
    end do  
  end do
end do

call END_LOG(mgsetup_end)


! Main loop
iters = 0
do itr = 2,itr_mg

    
  call START_LOG(mgfine_begin)
  
  ! update on model grid using colours
  do i = 1,itrbgn
    
    ! store previous solution to test for convergence  
    iv_old(1:ifull,1:klim) = iv(1:ifull,1:klim)
    
    do nc = 1,maxcolour
      ! MJT notes - We calculate the halo of the fine grid first, so that the message can
      ! be sent while the non-halo points are being updated, thereby overlapping communication
      ! and computation.

      ! update boundary grid points and send halo
      isc = 1
      iec = ifull_colour_border(nc)
      do k = 1,klim
        do iq = isc,iec  
          iv_new(iqx(iq,nc),k) = ( zznc(iq,nc)*iv(iqn(iq,nc),k)      &
                                 + zzwc(iq,nc)*iv(iqw(iq,nc),k)      &
                                 + zzec(iq,nc)*iv(iqe(iq,nc),k)      &
                                 + zzsc(iq,nc)*iv(iqs(iq,nc),k)      &
                                 - rhsc(iq,k,nc) )*rhelmc(iq,k,nc)
        end do  
      end do
      call bounds_colour_send(iv_new,nc,klim=klim)
      ! calculate non-boundary grid points while waiting for halo to update
      isc = ifull_colour_border(nc) + 1
      iec = ifull_colour(nc)
      do k = 1,klim
        do iq = isc,iec  
          xdum(iq) = ( zznc(iq,nc)*iv(iqn(iq,nc),k)      &
                     + zzwc(iq,nc)*iv(iqw(iq,nc),k)      &
                     + zzec(iq,nc)*iv(iqe(iq,nc),k)      &
                     + zzsc(iq,nc)*iv(iqs(iq,nc),k)      &
                     - rhsc(iq,k,nc) )*rhelmc(iq,k,nc)
        end do  
        do iq = 1,isc-1
          iv(iqx(iq,nc),k) = iv_new(iqx(iq,nc),k)
        end do  
        do iq = isc,iec
          iv(iqx(iq,nc),k) = xdum(iq)
        end do
      end do
      call bounds_colour_recv(iv,nc,klim=klim)
    end do
  end do
  
  do k = 1,klim
    ! test for convergence
    dsolmax_g(k) = maxval(abs(iv(1:ifull,k)-iv_old(1:ifull,k))) ! cannot vectorise with -fp-precise
  
    ! residual
    call unpack_nsew(iv(:,k),iv_n,iv_s,iv_e,iv_w)
    w(1:ifull,k)=-izzn(:)*iv_n(:)-izzw(:)*iv_w(:)-izze(:)*iv_e(:)-izzs(:)*iv_s(:) &
                 +irhs(:,k)+iv(1:ifull,k)*(ihelm(:,k)-izz(:))
  end do
  
  ! For when the inital grid cannot be upscaled
  call mgcollect(1,w(:,1:klim),dsolmax_g(1:klim),klim=klim)

  call END_LOG(mgfine_end)

  
  if ( mg_maxlevel_local>0 ) then
    
    call START_LOG(mgfine_begin)
    
    ! restriction
    ! (since this always operates within a panel, then ine = ien is always true)
    ng4 = mg(1)%ifull_fine
    do k = 1,klim
      do iq = 1,ng4
        rhs(iq,k,2) = 0.25*(w(mg(1)%fine(iq)  ,k) + w(mg(1)%fine_n(iq) ,k)  &
                          + w(mg(1)%fine_e(iq),k) + w(mg(1)%fine_ne(iq),k))
      end do  
    end do  
                             
    ! merge grids if insufficent points on this processor
    call mgcollect(2,rhs(:,1:klim,2),dsolmax_g(:),klim=klim)
 
    call END_LOG(mgfine_end)
    
    
    call START_LOG(mgup_begin)
  
    ! upscale grid
    do g = 2,gmax
  
      ng = mg(g)%ifull
                
      ! update scalar field
      ! assume zero for first guess of residual (also avoids additional bounds call)
      do k = 1,klim
        v(1:ng,k,g) = -rhs(1:ng,k,g)/( helm(1:ng,k,g) - mg(g)%zz(1:ng) )
      end do
      call mgbounds(g,v(:,1:klim,g),klim=klim)
      
      ! MJT notes - As the grid size per processor becomes small with continual upscaling, this part
      ! of the code becomes dominated by communications.  We then neglect to decompose this iteration
      ! into colours so that the number of messages is reduced.
      do i = 2,itrbgn
        do k = 1,klim
          ! post smoothing
          call mgunpack_nsew(g,v(:,k,g),v_n,v_s,v_e,v_w)   
          v(1:ng,k,g) = ( mg(g)%zze(1:ng)*v_e(1:ng) + mg(g)%zzw(1:ng)*v_w(1:ng) &
                        + mg(g)%zzn(1:ng)*v_n(1:ng) + mg(g)%zzs(1:ng)*v_s(1:ng) &
                        - rhs(1:ng,k,g) )/( helm(1:ng,k,g) - mg(g)%zz(1:ng) )
        end do
        call mgbounds(g,v(:,1:klim,g),klim=klim)
      end do
    
      ng4 = mg(g)%ifull_fine
      do k = 1,klim
        ! residual
        call mgunpack_nsew(g,v(:,k,g),v_n,v_s,v_e,v_w)     
        w(1:ng,k) = -mg(g)%zze(1:ng)*v_e(1:ng)-mg(g)%zzw(1:ng)*v_w(1:ng)   &
                    -mg(g)%zzn(1:ng)*v_n(1:ng)-mg(g)%zzs(1:ng)*v_s(1:ng)   &
                    +rhs(1:ng,k,g)+(helm(1:ng,k,g)-mg(g)%zz(1:ng))*v(1:ng,k,g)
        ! restriction
        rhs(1:ng4,k,g+1) = 0.25*( w(mg(g)%fine(1:ng4)  ,k) + w(mg(g)%fine_n(1:ng4) ,k)   &
                                + w(mg(g)%fine_e(1:ng4),k) + w(mg(g)%fine_ne(1:ng4),k) )
      end do

      ! merge grids if insufficent points on this processor
      call mgcollect(g+1,rhs(:,1:klim,g+1),dsolmax_g(1:klim),klim=klim)

    end do
  
    call END_LOG(mgup_end)
  
  
    if ( mg_maxlevel_local==mg_maxlevel ) then
         
      call START_LOG(mgcoarse_begin)
         
      g = mg_maxlevel
      ng = mg(g)%ifull
      do k = 1,klim
        do nc = 1,3
          do iq = 1,mg_ifull_maxcolour  
            helmc_c(iq,nc,k) = helm(col_iq(iq,nc),k,g) - mg(g)%zz(col_iq(iq,nc))
            rhsc_c(iq,nc,k) = rhs(col_iq(iq,nc),k,g)
          end do  
        end do
        v(1:ng,k,g) = -rhs(1:ng,k,g)/(helm(1:ng,k,g)-mg(g)%zz(1:ng))
        sdifc(k) = max( maxval(v(1:ng,k,g)) - minval(v(1:ng,k,g)), 1.e-20 )        
      end do
      ! usually it takes 6 iterations for the following to converge with 35 eigenvectors
      klimc = klim
      do itrc = 1,itr_mg
        do k = 1,klimc
          vsavc(1:ng) = v(1:ng,k,g)
          do nc = 1,3
            do iq = 1,mg_ifull_maxcolour  
              v(col_iq(iq,nc),k,g) = ( zznc_c(iq,nc)*v(col_iqn(iq,nc),k,g) &
                                     + zzec_c(iq,nc)*v(col_iqe(iq,nc),k,g) &
                                     + zzsc_c(iq,nc)*v(col_iqs(iq,nc),k,g) &
                                     + zzwc_c(iq,nc)*v(col_iqw(iq,nc),k,g) &
                                     - rhsc_c(iq,nc,k) )/helmc_c(iq,nc,k)
            end do  
          end do
          dsolmaxc(k) = maxval( abs( v(1:ng,k,g) - vsavc(1:ng) ) )
        end do
        ! test for convergence
        knew = klimc
        do k = klimc,1,-1
          if ( dsolmaxc(k)>=restol*sdifc(k) ) exit
          knew = k - 1
        end do
        klimc = knew
        if ( klimc<1 ) exit
      end do
      
      call END_LOG(mgcoarse_end)
      
    end if
    

    call START_LOG(mgdown_begin)

    ! downscale grid
    do g = gmax,2,-1

      call mgbcast(g+1,v(:,1:klim,g+1),dsolmax_g(:),klim=klim)

      ng = mg(g)%ifull 
      
      do k = 1,klim
        ! interpolation
        do iq = 1,ng
          w(iq,k) = 0.5*v(mg(g+1)%coarse_a(iq),k,g+1)  &
                  + 0.25*v(mg(g+1)%coarse_b(iq),k,g+1) &
                  + 0.25*v(mg(g+1)%coarse_c(iq),k,g+1)
        end do  

        ! extension
        ! No mgbounds as the v halo has already been updated and
        ! the coarse interpolation also updates the w halo
        v(1:ng,k,g) = v(1:ng,k,g) + w(1:ng,k)
      end do
    
      do i = 1,itrend
        call mgbounds(g,v(:,1:klim,g),klim=klim)
        do k = 1,klim
          ! post smoothing - all iterations except final iteration
          call mgunpack_nsew(g,v(:,k,g),v_n,v_s,v_e,v_w)  
          v(1:ng,k,g) = (mg(g)%zze(1:ng)*v_e(1:ng)+mg(g)%zzw(1:ng)*v_w(1:ng) &
                       + mg(g)%zzn(1:ng)*v_n(1:ng)+mg(g)%zzs(1:ng)*v_s(1:ng) &
                       - rhs(1:ng,k,g))/(helm(1:ng,k,g)-mg(g)%zz(1:ng))
        end do
      end do

      call mgbounds(g,v(:,1:klim,g),klim=klim) ! for next mgbcast
      
    end do
  
    
    ! fine grid
    call mgbcast(2,v(:,1:klim,2),dsolmax_g(:),klim=klim)

    ! interpolation
    ng = mg(1)%ifull
    do k = 1,klim
      do iq = 1,ng 
        w(iq,k) = 0.5*v(mg(2)%coarse_a(iq),k,2)  &
                + 0.25*v(mg(2)%coarse_b(iq),k,2) &
                + 0.25*v(mg(2)%coarse_c(iq),k,2)
      end do  
    end do

    
    call END_LOG(mgdown_end)
    
  end if
  
  
  call START_LOG(mgfine_begin)

  ! multi-grid solver bounds indices do not match standard iextra indicies, so we need to remap the halo
  if ( mg(1)%merge_len>1 ) then
    call mgbcast(1,w(:,1:klim),dsolmax_g(:),klim=klim,nobounds=.true.)
    ir = mod(mg(1)%merge_pos-1, mg(1)%merge_row) + 1   ! index for proc row
    ic = (mg(1)%merge_pos-1)/mg(1)%merge_row + 1       ! index for proc col
    do k = 1,klim
      do n = 1,npan
        do jj = 1,jpan
          iq_a = 1 + (jj-1)*ipan + (n-1)*ipan*jpan
          iq_c = 1 + (ir-1)*ipan + (jj-1+(ic-1)*jpan)*ipan*mg(1)%merge_row + (n-1)*ipan*jpan*mg(1)%merge_len
          vdum(iq_a:iq_a+ipan-1) = w(iq_c:iq_c+ipan-1,k)
        end do
      end do
      ! extension
      iv(1:ifull,k) = iv(1:ifull,k) + vdum(1:ifull)
    end do
  else
    ! remap mg halo to normal halo
    do k = 1,klim
      ! extension
      iv(1:ifull,k) = iv(1:ifull,k) + w(1:ifull,k)
    end do
  end if
  
  call bounds(iv(:,1:klim),klim=klim)
  
  do i = 1,itrend
    ! post smoothing
    do nc = 1,maxcolour
      isc = 1
      iec = ifull_colour_border(nc)
      do k = 1,klim
        do iq = isc,iec  
          iv_new(iqx(iq,nc),k) = ( zznc(iq,nc)*iv(iqn(iq,nc),k)      &
                                 + zzwc(iq,nc)*iv(iqw(iq,nc),k)      &
                                 + zzec(iq,nc)*iv(iqe(iq,nc),k)      &
                                 + zzsc(iq,nc)*iv(iqs(iq,nc),k)      &
                                 - rhsc(iq,k,nc) )*rhelmc(iq,k,nc)
        end do  
      end do
      call bounds_colour_send(iv_new,nc,klim=klim)
      isc = ifull_colour_border(nc) + 1
      iec = ifull_colour(nc)
      do k = 1,klim
        do iq = isc,iec  
          xdum(iq) = ( zznc(iq,nc)*iv(iqn(iq,nc),k)      &
                     + zzwc(iq,nc)*iv(iqw(iq,nc),k)      &
                     + zzec(iq,nc)*iv(iqe(iq,nc),k)      &
                     + zzsc(iq,nc)*iv(iqs(iq,nc),k)      &
                     - rhsc(iq,k,nc) )*rhelmc(iq,k,nc)
        end do  
        do iq = 1,isc-1
          iv(iqx(iq,nc),k) = iv_new(iqx(iq,nc),k)
        end do  
        do iq = isc,iec
          iv(iqx(iq,nc),k) = xdum(iq)
        end do  
      end do
      call bounds_colour_recv(iv,nc,klim=klim)
    end do
  end do

  call END_LOG(mgfine_end)

  ! test for convergence.  Test lags by one iteration, due to combining
  ! multi-grid and convergence communications.
  knew = klim
  do k = klim,1,-1
    iters(k) = itr
    if ( dsolmax_g(k)>=restol*sdif(k) ) exit
    knew = k - 1
  end do
  klim = knew
  if ( klim<1 ) exit
 
end do


! Combine offsets back into solution
do k = 1,kl
  iv(1:ifull,k) = iv(1:ifull,k) + savg(k)
end do

! Display convergence diagnostics
if ( myid==0 ) then
  if ( ktau<6 .or. iters(1)>=itr_mg ) then
    do k = 1,kl
      write(6,*) "mg ktau,k,iter ",ktau,k,iters(k),dsolmax_g(k)
    end do
  end if
end if

call END_LOG(helm_end)

return
end subroutine mghelm


! This version of the geometric multigrid solver is for the ocean and ice
subroutine mgmlo(neta,ipice,iyy,iyyn,iyys,iyye,iyyw,                   &
                 izz,izzn,izzs,izze,izzw,                              &
                 ihh,irhs,tol,itol,totits,maxglobseta,maxglobip,ipmax, &
                 ee,dd)

use cc_mpi
use indices_m
use newmpar_m

implicit none

integer, intent(out) :: totits
integer itr, itrc, g, ng, ng4, n, i, j, ir, ic, iq
integer iq_a, iq_c
integer nc, isc, iec, k
real, intent(in) :: tol, itol
real, intent(out) :: maxglobseta, maxglobip
real, dimension(ifull+iextra), intent(inout) :: neta, ipice
real, dimension(ifull+iextra), intent(in) :: ee, dd
real, dimension(ifull+iextra), intent(in) :: ipmax
real, dimension(ifull,2), intent(in) :: izz, izzn, izzs, izze, izzw
real, dimension(ifull,2), intent(in) :: irhs
real, dimension(ifull), intent(in) :: ihh
real, dimension(ifull), intent(in) :: iyy, iyyn, iyys, iyye, iyyw
real, dimension(ifull+iextra,2) :: dumc
real, dimension(ifull+iextra,2) :: vduma
real, dimension(ifull,2) :: dumc_n, dumc_s, dumc_e, dumc_w
real, dimension(ifull_maxcolour,maxcolour) :: rhsc, rhscice, ddc, eec, ipmaxc
real, dimension(ifull_maxcolour,maxcolour) :: zzhhc, zznc, zzsc, zzec, zzwc
real, dimension(ifull_maxcolour,maxcolour) :: zzcice, zzncice, zzscice, zzecice, zzwcice
real, dimension(ifull_maxcolour,maxcolour) :: yyc, yync, yysc, yyec, yywc
real, dimension(mg_maxsize,2,gmax+1) :: v
real, dimension(mg_maxsize,2:gmax+1) :: zz, zzn, zzs, zze, zzw
real, dimension(mg_maxsize,2:gmax+1) :: yyn, yys, yye, yyw, yyz
real, dimension(mg_maxsize,2:gmax+1) :: hh
real, dimension(mg_maxsize,2:gmax+1) :: zzi, zzin, zzis, zzie, zziw
real, dimension(mg_maxsize,2:gmax+1) :: rhs, rhsi
real, dimension(mg_maxsize,18) :: w
real, dimension(mg_maxsize) :: bu, cu
real, dimension(mg_maxsize,2) :: ws
real, dimension(mg_maxsize) :: v_n, v_s, v_e, v_w
real, dimension(mg_ifull_maxcolour,3) :: zzhhcu, zzncu, zzscu, zzecu, zzwcu, rhscu
real, dimension(mg_ifull_maxcolour,3) :: zzicu, zzincu, zziscu, zziecu, zziwcu, rhsicu
real, dimension(mg_ifull_maxcolour,3) :: yyzcu, yyncu, yyscu, yyecu, yywcu
real, dimension(2) :: dsolmax
real, dimension(8) :: dsolmax_g

if ( sorfirst ) then
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

! zz*(d2ipice/dx2 + d2ipice/dy2) + zz*(d2ipice/dxdy) = rhs

call START_LOG(waterhelm_begin)

call START_LOG(mgsetup_begin)

ng = 0
ng4 = 0
dumc = 0.
dsolmax = 0.
dsolmax_g = 0.

! pack colour arrays
do nc = 1,maxcolour
  do iq = 1,ifull_colour(nc)  
    yyc(iq,nc)      = iyy(iqx(iq,nc))
    yync(iq,nc)     = iyyn(iqx(iq,nc))
    yysc(iq,nc)     = iyys(iqx(iq,nc))
    yyec(iq,nc)     = iyye(iqx(iq,nc))
    yywc(iq,nc)     = iyyw(iqx(iq,nc))
    zzhhc(iq,nc)    = izz(iqx(iq,nc),1) + ihh(iqx(iq,nc))
    zznc(iq,nc)     = izzn(iqx(iq,nc),1)
    zzsc(iq,nc)     = izzs(iqx(iq,nc),1)
    zzec(iq,nc)     = izze(iqx(iq,nc),1)
    zzwc(iq,nc)     = izzw(iqx(iq,nc),1)
    rhsc(iq,nc)     = irhs(iqx(iq,nc),1)
    zzcice(iq,nc)   = izz(iqx(iq,nc),2)
    zzncice(iq,nc)  = izzn(iqx(iq,nc),2)
    zzscice(iq,nc)  = izzs(iqx(iq,nc),2)
    zzecice(iq,nc)  = izze(iqx(iq,nc),2)
    zzwcice(iq,nc)  = izzw(iqx(iq,nc),2)
    rhscice(iq,nc)  = irhs(iqx(iq,nc),2)
    ddc(iq,nc)      = dd(iqx(iq,nc))
    eec(iq,nc)      = ee(iqx(iq,nc))
    ipmaxc(iq,nc)   = ipmax(iqx(iq,nc))
  end do  
end do

! solver requires bounds to be updated
dumc(1:ifull,1) = neta(1:ifull)
dumc(1:ifull,2) = ipice(1:ifull)
call bounds(dumc(:,1:2))

do i = 1,itrbgn
  do nc = 1,maxcolour

    do k = 1,2    
      do iq = 1,ifull_colour(nc)  
        dumc_n(iq,k) = dumc(iqn(iq,nc),k)
        dumc_s(iq,k) = dumc(iqs(iq,nc),k)
        dumc_e(iq,k) = dumc(iqe(iq,nc),k)
        dumc_w(iq,k) = dumc(iqw(iq,nc),k)
      end do  
    end do  
      
    ! update halo
    isc = 1
    iec = ifull_colour_border(nc)
  
    ! ocean
    bu(isc:iec)=yync(isc:iec,nc)*dumc_n(isc:iec,1)+yysc(isc:iec,nc)*dumc_s(isc:iec,1)      &
               +yyec(isc:iec,nc)*dumc_e(isc:iec,1)+yywc(isc:iec,nc)*dumc_w(isc:iec,1)      &
               +zzhhc(isc:iec,nc)
    cu(isc:iec)=zznc(isc:iec,nc)*dumc_n(isc:iec,1)+zzsc(isc:iec,nc)*dumc_s(isc:iec,1)      &
               +zzec(isc:iec,nc)*dumc_e(isc:iec,1)+zzwc(isc:iec,nc)*dumc_w(isc:iec,1)      &
               -rhsc(isc:iec,nc)        
    do iq = isc,iec
      dumc(iqx(iq,nc),1) = eec(iq,nc)*max( -ddc(iq,nc),    &
         -2.*cu(iq)/(bu(iq)+sqrt(max(bu(iq)**2-4.*yyc(iq,nc)*cu(iq),0.1))) )
    end do  
    
    ! sea-ice (cavitating fluid)
    do iq = isc,iec
      dumc(iqx(iq,nc),2) = max(0.,min(ipmaxc(iq,nc),                       &
         ( -zzncice(iq,nc)*dumc_n(iq,2) - zzscice(iq,nc)*dumc_s(iq,2)      &
           -zzecice(iq,nc)*dumc_e(iq,2) - zzwcice(iq,nc)*dumc_w(iq,2)      &
          + rhscice(iq,nc) ) / zzcice(iq,nc) ))
    end do  

    call bounds_colour_send(dumc(:,1:2),nc)
    
    ! update interior
    isc = ifull_colour_border(nc) + 1
    iec = ifull_colour(nc)

    ! ocean
    bu(isc:iec)=yync(isc:iec,nc)*dumc_n(isc:iec,1)+yysc(isc:iec,nc)*dumc_s(isc:iec,1)      &
               +yyec(isc:iec,nc)*dumc_e(isc:iec,1)+yywc(isc:iec,nc)*dumc_w(isc:iec,1)      &
               +zzhhc(isc:iec,nc)
    cu(isc:iec)=zznc(isc:iec,nc)*dumc_n(isc:iec,1)+zzsc(isc:iec,nc)*dumc_s(isc:iec,1)      &
               +zzec(isc:iec,nc)*dumc_e(isc:iec,1)+zzwc(isc:iec,nc)*dumc_w(isc:iec,1)      &
               -rhsc(isc:iec,nc)   
    do iq = isc,iec
      dumc(iqx(iq,nc),1) = eec(iq,nc)*max( -ddc(iq,nc),     &
         -2.*cu(iq)/(bu(iq)+sqrt(max(bu(iq)**2-4.*yyc(iq,nc)*cu(iq),0.1))) )
    end do  
    
    ! sea-ice (cavitating fluid)
    do iq = isc,iec
      dumc(iqx(iq,nc),2) = max(0.,min(ipmaxc(iq,nc),                       &
         ( -zzncice(iq,nc)*dumc_n(iq,2) - zzscice(iq,nc)*dumc_s(iq,2)      &
           -zzecice(iq,nc)*dumc_e(iq,2) - zzwcice(iq,nc)*dumc_w(iq,2)      &
          + rhscice(iq,nc) ) / zzcice(iq,nc) ))
    end do  

    call bounds_colour_recv(dumc(:,1:2),nc)
    
  end do
end do

! test for convergence
neta(1:ifull+iextra)  = dumc(1:ifull+iextra,1)
ipice(1:ifull+iextra) = dumc(1:ifull+iextra,2)  

! residual - ocean
call unpack_nsew(dumc(:,1),dumc_n(:,1),dumc_s(:,1),dumc_e(:,1),dumc_w(:,1))
w(1:ifull,1)=(-neta(1:ifull)*(iyy(1:ifull)*neta(1:ifull)                             &
              +iyyn(1:ifull)*dumc_n(1:ifull,1)+iyys(1:ifull)*dumc_s(1:ifull,1)       &
              +iyye(1:ifull)*dumc_e(1:ifull,1)+iyyw(1:ifull)*dumc_w(1:ifull,1))      &
             -(izz(1:ifull,1)*neta(1:ifull)                                              &
              +izzn(1:ifull,1)*dumc_n(1:ifull,1)+izzs(1:ifull,1)*dumc_s(1:ifull,1)       &
              +izze(1:ifull,1)*dumc_e(1:ifull,1)+izzw(1:ifull,1)*dumc_w(1:ifull,1))      &
             -ihh(1:ifull)*neta(1:ifull)+irhs(1:ifull,1))*ee(1:ifull)

! upscale ocean fields
w(1:ifull,2)  = iyy(1:ifull)
w(1:ifull,3)  = iyyn(1:ifull)
w(1:ifull,4)  = iyys(1:ifull)
w(1:ifull,5)  = iyye(1:ifull)
w(1:ifull,6)  = iyyw(1:ifull)
w(1:ifull,7)  = izz(:,1) +  iyy*neta(1:ifull)
w(1:ifull,8)  = izzn(:,1) + iyyn*neta(1:ifull)
w(1:ifull,9)  = izzs(:,1) + iyys*neta(1:ifull)
w(1:ifull,10) = izze(:,1) + iyye*neta(1:ifull)
w(1:ifull,11) = izzw(:,1) + iyyw*neta(1:ifull)
w(1:ifull,12) = iyyn*dumc_n(1:ifull,1) + iyys*dumc_s(1:ifull,1)                         &
              + iyye*dumc_e(1:ifull,1) + iyyw*dumc_w(1:ifull,1)                         &
              + iyy*neta(1:ifull) + ihh(:)

! residual - ice
call unpack_nsew(dumc(:,2),dumc_n(:,2),dumc_s(:,2),dumc_e(:,2),dumc_w(:,2))
w(1:ifull,13) =(- izzn(1:ifull,2)*dumc_n(1:ifull,2) - izzs(1:ifull,2)*dumc_s(1:ifull,2)     &
                - izze(1:ifull,2)*dumc_e(1:ifull,2) - izzw(1:ifull,2)*dumc_w(1:ifull,2)     &
                - izz(1:ifull,2)*ipice(1:ifull) + irhs(1:ifull,2))*ee(1:ifull)
where ( ipice(1:ifull)>=ipmax(1:ifull) )
  w(1:ifull,13) = 0. ! improves convergence
end where

! update ice fields
w(1:ifull,14) = izz(1:ifull,2)
w(1:ifull,15) = izzn(1:ifull,2)
w(1:ifull,16) = izzs(1:ifull,2)
w(1:ifull,17) = izze(1:ifull,2)
w(1:ifull,18) = izzw(1:ifull,2)

call mgcollect(1,w(:,1:18))
  
if ( mg_maxlevel_local>0 ) then
  
  ! restriction
  ! (since this always operates within a panel, then ine = ien is always true)
  ng4 = mg(1)%ifull_fine
  do iq = 1,ng4
    rhs(iq,2)=0.25*(w(mg(1)%fine(iq)  ,1)+w(mg(1)%fine_n(iq) ,1)           &
                   +w(mg(1)%fine_e(iq),1)+w(mg(1)%fine_ne(iq),1))
    yyz(iq,2) =0.25*dfac*(w(mg(1)%fine(iq)  ,2)+w(mg(1)%fine_n(iq) ,2)     &
                         +w(mg(1)%fine_e(iq),2)+w(mg(1)%fine_ne(iq),2))
    yyn(iq,2) =0.25*dfac*(w(mg(1)%fine(iq)  ,3)+w(mg(1)%fine_n(iq) ,3)     &
                         +w(mg(1)%fine_e(iq),3)+w(mg(1)%fine_ne(iq),3))
    yys(iq,2) =0.25*dfac*(w(mg(1)%fine(iq)  ,4)+w(mg(1)%fine_n(iq) ,4)     &
                         +w(mg(1)%fine_e(iq),4)+w(mg(1)%fine_ne(iq),4))
    yye(iq,2) =0.25*dfac*(w(mg(1)%fine(iq)  ,5)+w(mg(1)%fine_n(iq) ,5)     &
                         +w(mg(1)%fine_e(iq),5)+w(mg(1)%fine_ne(iq),5))
    yyw(iq,2) =0.25*dfac*(w(mg(1)%fine(iq)  ,6)+w(mg(1)%fine_n(iq) ,6)     &
                         +w(mg(1)%fine_e(iq),6)+w(mg(1)%fine_ne(iq),6))
    zz(iq,2) =0.25*dfac*(w(mg(1)%fine(iq)  ,7)+w(mg(1)%fine_n(iq) ,7)    &
                        +w(mg(1)%fine_e(iq),7)+w(mg(1)%fine_ne(iq),7))
    zzn(iq,2)=0.25*dfac*(w(mg(1)%fine(iq)  ,8)+w(mg(1)%fine_n(iq) ,8)    &
                        +w(mg(1)%fine_e(iq),8)+w(mg(1)%fine_ne(iq),8))
    zzs(iq,2)=0.25*dfac*(w(mg(1)%fine(iq)  ,9)+w(mg(1)%fine_n(iq) ,9)    &
                        +w(mg(1)%fine_e(iq),9)+w(mg(1)%fine_ne(iq),9))
    zze(iq,2)=0.25*dfac*(w(mg(1)%fine(iq)  ,10)+w(mg(1)%fine_n(iq) ,10)    &
                        +w(mg(1)%fine_e(iq),10)+w(mg(1)%fine_ne(iq),10))
    zzw(iq,2)=0.25*dfac*(w(mg(1)%fine(iq)  ,11)+w(mg(1)%fine_n(iq) ,11)    &
                        +w(mg(1)%fine_e(iq),11)+w(mg(1)%fine_ne(iq),11))
    hh(iq,2)    =0.25*(w(mg(1)%fine(iq)  ,12)+w(mg(1)%fine_n(iq) ,12)      &
                      +w(mg(1)%fine_e(iq),12)+w(mg(1)%fine_ne(iq),12))
    
    rhsi(iq,2)=0.25*(w(mg(1)%fine(iq)  ,13)+w(mg(1)%fine_n(iq) ,13)        &
                    +w(mg(1)%fine_e(iq),13)+w(mg(1)%fine_ne(iq),13))
    ! special treatment of ice
    zzi(iq,2) =0.25*dfaci*(w(mg(1)%fine(iq)  ,14)+w(mg(1)%fine_n(iq) ,14)  &
                          +w(mg(1)%fine_e(iq),14)+w(mg(1)%fine_ne(iq),14))
    zzin(iq,2)=0.25*dfaci*(w(mg(1)%fine(iq)  ,15)+w(mg(1)%fine_n(iq) ,15)  &
                          +w(mg(1)%fine_e(iq),15)+w(mg(1)%fine_ne(iq),15))
    zzis(iq,2)=0.25*dfaci*(w(mg(1)%fine(iq)  ,16)+w(mg(1)%fine_n(iq) ,16)  &
                          +w(mg(1)%fine_e(iq),16)+w(mg(1)%fine_ne(iq),16))
    zzie(iq,2)=0.25*dfaci*(w(mg(1)%fine(iq)  ,17)+w(mg(1)%fine_n(iq) ,17)  &
                          +w(mg(1)%fine_e(iq),17)+w(mg(1)%fine_ne(iq),17))
    zziw(iq,2)=0.25*dfaci*(w(mg(1)%fine(iq)  ,18)+w(mg(1)%fine_n(iq) ,18)  &
                          +w(mg(1)%fine_e(iq),18)+w(mg(1)%fine_ne(iq),18))
  end do  

  ! merge grids if insufficent points on this processor
  if ( mg(2)%merge_len>1 ) then
    w(1:ng4,1)  = rhs(1:ng4,2)
    w(1:ng4,2)  = yyz(1:ng4,2)
    w(1:ng4,3)  = yyn(1:ng4,2)
    w(1:ng4,4)  = yys(1:ng4,2)
    w(1:ng4,5)  = yye(1:ng4,2)
    w(1:ng4,6)  = yyw(1:ng4,2)
    w(1:ng4,7)  = zz(1:ng4,2)
    w(1:ng4,8)  = zzn(1:ng4,2)
    w(1:ng4,9)  = zzs(1:ng4,2)
    w(1:ng4,10) = zze(1:ng4,2)
    w(1:ng4,11) = zzw(1:ng4,2)
    w(1:ng4,12) = hh(1:ng4,2)
    w(1:ng4,13) = rhsi(1:ng4,2)
    w(1:ng4,14) = zzi(1:ng4,2)
    w(1:ng4,15) = zzin(1:ng4,2)
    w(1:ng4,16) = zzis(1:ng4,2)
    w(1:ng4,17) = zzie(1:ng4,2)
    w(1:ng4,18) = zziw(1:ng4,2)
    call mgcollect(2,w(:,1:18))
    if ( 2<=mg_maxlevel_local ) then
      ng = mg(2)%ifull
      rhs(1:ng,2)   = w(1:ng,1)
      yyz(1:ng,2)   = w(1:ng,2)
      yyn(1:ng,2)   = w(1:ng,3)
      yys(1:ng,2)   = w(1:ng,4)
      yye(1:ng,2)   = w(1:ng,5)
      yyw(1:ng,2)   = w(1:ng,6)
      zz(1:ng,2)    = w(1:ng,7)
      zzn(1:ng,2)   = w(1:ng,8)
      zzs(1:ng,2)   = w(1:ng,9)
      zze(1:ng,2)   = w(1:ng,10)
      zzw(1:ng,2)   = w(1:ng,11)
      hh(1:ng,2)    = w(1:ng,12)
      rhsi(1:ng,2)  = w(1:ng,13)
      zzi(1:ng,2)   = w(1:ng,14)
      zzin(1:ng,2)  = w(1:ng,15)
      zzis(1:ng,2)  = w(1:ng,16)
      zzie(1:ng,2)  = w(1:ng,17)
      zziw(1:ng,2)  = w(1:ng,18)
    end if  
  end if
    
  
  ! upscale grid
  do g = 2,gmax
  
    ng = mg(g)%ifull

    ! update
    ! possibly use colours here, although v is reset to zero every iteration
    ! assume zero for first guess of residual (also avoids additional bounds call)
    bu(1:ng) = zz(1:ng,g) + hh(1:ng,g)
    cu(1:ng) = -rhs(1:ng,g)    
    v(:,1:2,g) = 0.
    v(1:ng,1,g) = -2.*cu(1:ng)/(bu(1:ng)+sqrt(max(bu(1:ng)**2-4.*yyz(1:ng,g)*cu(1:ng),0.1)))
    v(1:ng,2,g) = rhsi(1:ng,g) / zzi(1:ng,g)
    call mgbounds(g,v(:,1:2,g))
    
    do i = 2,itrbgn

      ! ocean - post smoothing
      call mgunpack_nsew(g,v(:,1,g),v_n,v_s,v_e,v_w)
      bu(1:ng) = yyn(1:ng,g)*v_n(1:ng)+yys(1:ng,g)*v_s(1:ng)       &
               + yye(1:ng,g)*v_e(1:ng)+yyw(1:ng,g)*v_w(1:ng)       &
               + zz(1:ng,g) + hh(1:ng,g)
      cu(1:ng) = zzn(1:ng,g)*v_n(1:ng) + zzs(1:ng,g)*v_s(1:ng)     &
               + zze(1:ng,g)*v_e(1:ng) + zzw(1:ng,g)*v_w(1:ng)     &
               - rhs(1:ng,g)
      v(1:ng,1,g) = -2.*cu(1:ng)/(bu(1:ng)+sqrt(max(bu(1:ng)**2-4.*yyz(1:ng,g)*cu(1:ng),0.1)))
      
      ! ice - post smoothing
      call mgunpack_nsew(g,v(:,2,g),v_n,v_s,v_e,v_w)
      v(1:ng,2,g) = ( - zzin(1:ng,g)*v_n(1:ng) - zzis(1:ng,g)*v_s(1:ng)     &
                      - zzie(1:ng,g)*v_e(1:ng) - zziw(1:ng,g)*v_w(1:ng)     &
                      + rhsi(1:ng,g) ) / zzi(1:ng,g)
      
      call mgbounds(g,v(:,1:2,g))
    end do
  
    ! restriction
    ! (calculate finer grid before mgcollect as the messages sent/recv are shorter)

    ng4 = mg(g)%ifull_fine

    ! ocean residual
    call mgunpack_nsew(g,v(:,1,g),v_n,v_s,v_e,v_w)
    ws(1:ng,1)=-v(1:ng,1,g)*(yyz(1:ng,g)*v(1:ng,1,g)+yyn(1:ng,g)*v_n(1:ng)+yys(1:ng,g)*v_s(1:ng)      &
                                                    +yye(1:ng,g)*v_e(1:ng)+yyw(1:ng,g)*v_w(1:ng))     &
               -(zz(1:ng,g)*v(1:ng,1,g)+zzn(1:ng,g)*v_n(1:ng)+zzs(1:ng,g)*v_s(1:ng)                   &
                                       +zze(1:ng,g)*v_e(1:ng)+zzw(1:ng,g)*v_w(1:ng))                  &
                -hh(1:ng,g)*v(1:ng,1,g)+rhs(1:ng,g)
    do iq = 1,ng4
      w(iq,1)=0.25*(ws(mg(g)%fine(iq)  ,1)+ws(mg(g)%fine_n(iq) ,1) &
                       +ws(mg(g)%fine_e(iq),1)+ws(mg(g)%fine_ne(iq),1))
      w(iq,2) = 0.25*dfac*(yyz(mg(g)%fine(iq)  ,g)+yyz(mg(g)%fine_n(iq) ,g) &
                              +yyz(mg(g)%fine_e(iq),g)+yyz(mg(g)%fine_ne(iq),g))
      w(iq,3) = 0.25*dfac*(yyn(mg(g)%fine(iq)  ,g)+yyn(mg(g)%fine_n(iq) ,g) &
                              +yyn(mg(g)%fine_e(iq),g)+yyn(mg(g)%fine_ne(iq),g))
      w(iq,4) = 0.25*dfac*(yys(mg(g)%fine(iq)  ,g)+yys(mg(g)%fine_n(iq) ,g) &
                              +yys(mg(g)%fine_e(iq),g)+yys(mg(g)%fine_ne(iq),g))
      w(iq,5) = 0.25*dfac*(yye(mg(g)%fine(iq)  ,g)+yye(mg(g)%fine_n(iq) ,g) &
                              +yye(mg(g)%fine_e(iq),g)+yye(mg(g)%fine_ne(iq),g))
      w(iq,6) = 0.25*dfac*(yyw(mg(g)%fine(iq)  ,g)+yyw(mg(g)%fine_n(iq) ,g) &
                              +yyw(mg(g)%fine_e(iq),g)+yyw(mg(g)%fine_ne(iq),g))
    end do
    
    ws(1:ng,1) = zz(1:ng,g) + yyz(1:ng,g)*v(1:ng,1,g)
    do iq = 1,ng4
      w(iq,7) = 0.25*dfac*(ws(mg(g)%fine(iq)  ,1)+ws(mg(g)%fine_n(iq) ,1) &
                          +ws(mg(g)%fine_e(iq),1)+ws(mg(g)%fine_ne(iq),1))
    end do
    
    ws(1:ng,1) = zzn(1:ng,g) + yyn(1:ng,g)*v(1:ng,1,g)
    do iq = 1,ng4
      w(iq,8) = 0.25*dfac*(ws(mg(g)%fine(iq)  ,1)+ws(mg(g)%fine_n(iq) ,1) &
                          +ws(mg(g)%fine_e(iq),1)+ws(mg(g)%fine_ne(iq),1))
    end do  
    
    ws(1:ng,1) = zzs(1:ng,g) + yys(1:ng,g)*v(1:ng,1,g)
    do iq = 1,ng4
      w(iq,9) = 0.25*dfac*(ws(mg(g)%fine(iq)  ,1)+ws(mg(g)%fine_n(iq) ,1) &
                          +ws(mg(g)%fine_e(iq),1)+ws(mg(g)%fine_ne(iq),1))
    end do
    
    ws(1:ng,1) = zze(1:ng,g) + yye(1:ng,g)*v(1:ng,1,g)
    do iq = 1,ng4
      w(iq,10) = 0.25*dfac*(ws(mg(g)%fine(iq)  ,1)+ws(mg(g)%fine_n(iq) ,1) &
                           +ws(mg(g)%fine_e(iq),1)+ws(mg(g)%fine_ne(iq),1))
    end do  
    
    ws(1:ng,1) = zzw(1:ng,g) + yyw(1:ng,g)*v(1:ng,1,g)
    do iq = 1,ng4
      w(iq,11) = 0.25*dfac*(ws(mg(g)%fine(iq)  ,1)+ws(mg(g)%fine_n(iq) ,1) &
                           +ws(mg(g)%fine_e(iq),1)+ws(mg(g)%fine_ne(iq),1))
    end do 
    
    ws(1:ng,1) = yyn(1:ng,g)*v_n(1:ng)+yys(1:ng,g)*v_s(1:ng)     &
                +yye(1:ng,g)*v_e(1:ng)+yyw(1:ng,g)*v_w(1:ng)     &
                +yyz(1:ng,g)*v(1:ng,1,g) + hh(1:ng,g)
    do iq = 1,ng4
      w(iq,12) = 0.25*(ws(mg(g)%fine(iq)  ,1)+ws(mg(g)%fine_n(iq) ,1)           &
                      +ws(mg(g)%fine_e(iq),1)+ws(mg(g)%fine_ne(iq),1))
    end do

    ! ice residual
    call mgunpack_nsew(g,v(:,2,g),v_n,v_s,v_e,v_w)
    ws(1:ng,2) = -zzin(1:ng,g)*v_n(1:ng)-zzis(1:ng,g)*v_s(1:ng)     &
                 -zzie(1:ng,g)*v_e(1:ng)-zziw(1:ng,g)*v_w(1:ng)     &
                 -zzi(1:ng,g)*v(1:ng,2,g)+rhsi(1:ng,g)
    do iq = 1,ng4
      w(iq,13)=0.25*(ws(mg(g)%fine(iq)  ,2)+ws(mg(g)%fine_n(iq) ,2)  &
                        +ws(mg(g)%fine_e(iq),2)+ws(mg(g)%fine_ne(iq),2))
      ! special treatment of ice (neglect dfac)
      w(iq,14) = 0.25*dfaci*(zzi(mg(g)%fine(iq)  ,g)+zzi(mg(g)%fine_n(iq) ,g)    &
                            +zzi(mg(g)%fine_e(iq),g)+zzi(mg(g)%fine_ne(iq),g))
      w(iq,15) = 0.25*dfaci*(zzin(mg(g)%fine(iq)  ,g)+zzin(mg(g)%fine_n(iq) ,g) &
                            +zzin(mg(g)%fine_e(iq),g)+zzin(mg(g)%fine_ne(iq),g))
      w(iq,16) = 0.25*dfaci*(zzis(mg(g)%fine(iq)  ,g)+zzis(mg(g)%fine_n(iq) ,g) &
                            +zzis(mg(g)%fine_e(iq),g)+zzis(mg(g)%fine_ne(iq),g))
      w(iq,17) = 0.25*dfaci*(zzie(mg(g)%fine(iq)  ,g)+zzie(mg(g)%fine_n(iq) ,g) &
                            +zzie(mg(g)%fine_e(iq),g)+zzie(mg(g)%fine_ne(iq),g))
      w(iq,18) = 0.25*dfaci*(zziw(mg(g)%fine(iq)  ,g)+zziw(mg(g)%fine_n(iq) ,g) &
                            +zziw(mg(g)%fine_e(iq),g)+zziw(mg(g)%fine_ne(iq),g))
    end do
    
    ! merge grids if insufficent points on this processor
    if ( mg(g+1)%merge_len>1 ) then
      call mgcollect(g+1,w(:,1:18)) 
    end if
    
    if ( g+1<=mg_maxlevel_local ) then
      ng = mg(g+1)%ifull
      rhs(1:ng,g+1)   = w(1:ng,1)
      yyz(1:ng,g+1)   = w(1:ng,2)
      yyn(1:ng,g+1)   = w(1:ng,3)
      yys(1:ng,g+1)   = w(1:ng,4)
      yye(1:ng,g+1)   = w(1:ng,5)
      yyw(1:ng,g+1)   = w(1:ng,6)
      zz(1:ng,g+1)    = w(1:ng,7)
      zzn(1:ng,g+1)   = w(1:ng,8)
      zzs(1:ng,g+1)   = w(1:ng,9)
      zze(1:ng,g+1)   = w(1:ng,10)
      zzw(1:ng,g+1)   = w(1:ng,11)
      hh(1:ng,g+1)    = w(1:ng,12)
      rhsi(1:ng,g+1)  = w(1:ng,13)
      zzi(1:ng,g+1)   = w(1:ng,14)
      zzin(1:ng,g+1)  = w(1:ng,15)
      zzis(1:ng,g+1)  = w(1:ng,16)
      zzie(1:ng,g+1)  = w(1:ng,17)
      zziw(1:ng,g+1)  = w(1:ng,18)
    end if 

  end do

  ! solve for coarse grid with coloured SOR
  if ( mg_maxlevel_local==mg_maxlevel ) then
  
    g = mg_maxlevel  
    ng = mg(g)%ifull
  
    ! pack rhsc_c, helmc_c, zznc_c, zzec_c, zzwc_c and zzsc_c by colour
    ! pack yy by colour
    ! pack zz,hh and rhs by colour
    do nc = 1,3
      do iq = 1,mg_ifull_maxcolour
        ! ocean  
        yyzcu(iq,nc)  = yyz(col_iq(iq,nc),g)
        yyncu(iq,nc)  = yyn(col_iq(iq,nc),g)
        yyscu(iq,nc)  = yys(col_iq(iq,nc),g)
        yyecu(iq,nc)  = yye(col_iq(iq,nc),g)
        yywcu(iq,nc)  = yyw(col_iq(iq,nc),g)
        zzhhcu(iq,nc) = zz(col_iq(iq,nc),g) + hh(col_iq(iq,nc),g)
        zzncu(iq,nc)  = zzn(col_iq(iq,nc),g)
        zzscu(iq,nc)  = zzs(col_iq(iq,nc),g)
        zzecu(iq,nc)  = zze(col_iq(iq,nc),g)
        zzwcu(iq,nc)  = zzw(col_iq(iq,nc),g)
        rhscu(iq,nc)  = rhs(col_iq(iq,nc),g)
        ! ice
        zzicu(iq,nc)   = zzi(col_iq(iq,nc),g)
        zzincu(iq,nc)  = zzin(col_iq(iq,nc),g)
        zziscu(iq,nc)  = zzis(col_iq(iq,nc),g)
        zziecu(iq,nc)  = zzie(col_iq(iq,nc),g)
        zziwcu(iq,nc)  = zziw(col_iq(iq,nc),g)
        rhsicu(iq,nc)  = rhsi(col_iq(iq,nc),g)
      end do  
    end do  
  
    ! solve non-linear water free surface and solve for ice with coloured SOR
    ! first guess
    bu(1:ng) = zz(1:ng,g) + hh(1:ng,g)
    cu(1:ng) = -rhs(1:ng,g)
    v(:,1:2,g) = 0.
    v(1:ng,1,g) = -2.*cu(1:ng)/(bu(1:ng)+sqrt(max(bu(1:ng)**2-4.*yyz(1:ng,g)*cu(1:ng),0.1))) ! ocean
    v(1:ng,2,g) = rhsi(1:ng,g) / zzi(1:ng,g)                                        ! ice
    do itrc = 1,itr_mgice
      ! store previous guess for convegence test
      ws(1:ng,1:2) = v(1:ng,1:2,g)
      do nc = 1,3

        ! ocean
        bu(1:mg_ifull_maxcolour) = yyncu(:,nc)*v(col_iqn(:,nc),1,g) + yyscu(:,nc)*v(col_iqs(:,nc),1,g)     &
                                 + yyecu(:,nc)*v(col_iqe(:,nc),1,g) + yywcu(:,nc)*v(col_iqw(:,nc),1,g)     &
                                 + zzhhcu(:,nc)
        cu(1:mg_ifull_maxcolour) = zzncu(:,nc)*v(col_iqn(:,nc),1,g) + zzscu(:,nc)*v(col_iqs(:,nc),1,g)     &
                                 + zzecu(:,nc)*v(col_iqe(:,nc),1,g) + zzwcu(:,nc)*v(col_iqw(:,nc),1,g)     &
                                 - rhscu(:,nc)
        do iq = 1,mg_ifull_maxcolour
          v(col_iq(iq,nc),1,g) = -2.*cu(iq)/(bu(iq)+sqrt(max(bu(iq)**2-4.*yyzcu(iq,nc)*cu(iq),0.1)))
        end do  
        
        ! ice
        do iq = 1,mg_ifull_maxcolour
          v(col_iq(iq,nc),2,g) = ( - zzincu(iq,nc)*v(col_iqn(iq,nc),2,g) - zziscu(iq,nc)*v(col_iqs(iq,nc),2,g)     &
                                   - zziecu(iq,nc)*v(col_iqe(iq,nc),2,g) - zziwcu(iq,nc)*v(col_iqw(iq,nc),2,g)     &
                                  + rhsicu(iq,nc) ) / zzicu(iq,nc)
        end do

      end do
      ! test for convergence
      dsolmax(1:2) = maxval( abs( v(1:ng,1:2,g) - ws(1:ng,1:2) ) )
      if ( dsolmax(1)<tol .and. dsolmax(2)<itol ) exit
    end do
    
  end if
  
  ! downscale grid
  do g = gmax,2,-1

    call mgbcasta(g+1,v(:,1:2,g+1))

    ! interpolation
    ng = mg(g)%ifull
    do k = 1,2
      do iq = 1,ng
        ws(iq,k) =  0.5*v(mg(g+1)%coarse_a(iq),k,g+1)   &
                  + 0.25*v(mg(g+1)%coarse_b(iq),k,g+1)  &
                  + 0.25*v(mg(g+1)%coarse_c(iq),k,g+1)
      end do  
    end do  

    ! extension
    ! No mgbounds as the v halo has already been updated and
    ! the coarse interpolation also updates the w halo
    v(1:ng,1:2,g) = v(1:ng,1:2,g) + ws(1:ng,1:2)
   
    
    do i = 1,itrend
      call mgbounds(g,v(:,1:2,g))
      ! ocean - post smoothing
      call mgunpack_nsew(g,v(:,1,g),v_n,v_s,v_e,v_w)
      bu(1:ng) = yyn(1:ng,g)*v_n(1:ng)+yys(1:ng,g)*v_s(1:ng)     &
               + yye(1:ng,g)*v_e(1:ng)+yyw(1:ng,g)*v_w(1:ng)     &
               + zz(1:ng,g) + hh(1:ng,g)
      cu(1:ng) = zzn(1:ng,g)*v_n(1:ng)+zzs(1:ng,g)*v_s(1:ng)     &
               + zze(1:ng,g)*v_e(1:ng)+zzw(1:ng,g)*v_w(1:ng)     &
               - rhs(1:ng,g)
      v(1:ng,1,g) = -2.*cu(1:ng)/(bu(1:ng)+sqrt(max(bu(1:ng)**2-4.*yyz(1:ng,g)*cu(1:ng),0.1)))
      
      ! ice - post smoothing
      call mgunpack_nsew(g,v(:,2,g),v_n,v_s,v_e,v_w)  
      v(1:ng,2,g) = ( - zzin(1:ng,g)*v_n(1:ng) - zzis(1:ng,g)*v_s(1:ng)     &
                      - zzie(1:ng,g)*v_e(1:ng) - zziw(1:ng,g)*v_w(1:ng)     &
                      + rhsi(1:ng,g) ) / zzi(1:ng,g)
      
    end do
    
    call mgbounds(g,v(:,1:2,g)) ! for next mgbcast

  end do

  
  ! fine grid
  call mgbcasta(2,v(:,1:2,2))

  ! interpolation
  ng = mg(1)%ifull
  do k = 1,2
    do iq = 1,ng
      ws(iq,k) = 0.5*v(mg(2)%coarse_a(iq),k,2)  &
               + 0.25*v(mg(2)%coarse_b(iq),k,2) &
               + 0.25*v(mg(2)%coarse_c(iq),k,2)
    end do  
  end do  

end if

vduma = 0.
if ( mg(1)%merge_len>1 ) then
  call mgbcasta(1,ws(:,1:2),nobounds=.true.)
  ir = mod(mg(1)%merge_pos-1,mg(1)%merge_row)+1   ! index for proc row
  ic = (mg(1)%merge_pos-1)/mg(1)%merge_row+1      ! index for proc col
  do n = 1,npan
    do j = 1,jpan
      iq_a = 1 + (j-1)*ipan + (n-1)*ipan*jpan
      iq_c = 1 + (ir-1)*ipan + (j-1+(ic-1)*jpan)*ipan*mg(1)%merge_row + (n-1)*ipan*jpan*mg(1)%merge_len
      vduma(iq_a:iq_a+ipan-1,1:2) = ws(iq_c:iq_c+ipan-1,1:2)
    end do
  end do
else
  ! remap mg halo to normal halo 
  vduma(1:ifull,1:2) = ws(1:ifull,1:2)
end if

! extension
vduma(1:ifull,1) = max( -10., min( 10., vduma(1:ifull,1) ) )
dumc(1:ifull,1) = max( neta(1:ifull)+vduma(1:ifull,1), -dd(1:ifull) )*ee(1:ifull)
dumc(1:ifull,2) = max( min( ipice(1:ifull)+vduma(1:ifull,2), ipmax(1:ifull) ), 0. )*ee(1:ifull)
 
neta(1:ifull)  = dumc(1:ifull,1)
ipice(1:ifull) = dumc(1:ifull,2)

call bounds(dumc(:,1:2))

do i = 1,itrend
  
  ! post smoothing
  do nc = 1,maxcolour

    do k = 1,2
      do iq = 1,ifull_colour(nc)  
        dumc_n(iq,k) = dumc(iqn(iq,nc),k)
        dumc_s(iq,k) = dumc(iqs(iq,nc),k)
        dumc_e(iq,k) = dumc(iqe(iq,nc),k)
        dumc_w(iq,k) = dumc(iqw(iq,nc),k)
      end do
    end do  
      
    ! update halo
    isc = 1
    iec = ifull_colour_border(nc)
  
    ! ocean
    bu(isc:iec)=yync(isc:iec,nc)*dumc_n(isc:iec,1)+yysc(isc:iec,nc)*dumc_s(isc:iec,1)      &
               +yyec(isc:iec,nc)*dumc_e(isc:iec,1)+yywc(isc:iec,nc)*dumc_w(isc:iec,1)      &
               +zzhhc(isc:iec,nc)
    cu(isc:iec)=zznc(isc:iec,nc)*dumc_n(isc:iec,1)+zzsc(isc:iec,nc)*dumc_s(isc:iec,1)      &
               +zzec(isc:iec,nc)*dumc_e(isc:iec,1)+zzwc(isc:iec,nc)*dumc_w(isc:iec,1)      &
               -rhsc(isc:iec,nc)     
    do iq = isc,iec
      dumc(iqx(iq,nc),1) = eec(iq,nc)*max( -ddc(iq,nc),    &
         -2.*cu(iq)/(bu(iq)+sqrt(max(bu(iq)**2-4.*yyc(iq,nc)*cu(iq),0.1))) )
    end do  
    
    ! sea-ice (cavitating fluid)
    do iq = isc,iec
      dumc(iqx(iq,nc),2) = max(0.,min(ipmaxc(iq,nc),                       &
         ( -zzncice(iq,nc)*dumc_n(iq,2) - zzscice(iq,nc)*dumc_s(iq,2)      &
           -zzecice(iq,nc)*dumc_e(iq,2) - zzwcice(iq,nc)*dumc_w(iq,2)      &
          + rhscice(iq,nc) ) / zzcice(iq,nc) ))
    end do  

    call bounds_colour_send(dumc(:,1:2),nc)
    
    ! update interior
    isc = ifull_colour_border(nc) + 1
    iec = ifull_colour(nc)

    ! ocean
    bu(isc:iec)=yync(isc:iec,nc)*dumc_n(isc:iec,1)+yysc(isc:iec,nc)*dumc_s(isc:iec,1)      &
               +yyec(isc:iec,nc)*dumc_e(isc:iec,1)+yywc(isc:iec,nc)*dumc_w(isc:iec,1)      &
               +zzhhc(isc:iec,nc)
    cu(isc:iec)=zznc(isc:iec,nc)*dumc_n(isc:iec,1)+zzsc(isc:iec,nc)*dumc_s(isc:iec,1)      &
               +zzec(isc:iec,nc)*dumc_e(isc:iec,1)+zzwc(isc:iec,nc)*dumc_w(isc:iec,1)      &
               -rhsc(isc:iec,nc)   
    do iq = isc,iec
      dumc(iqx(iq,nc),1) = eec(iq,nc)*max( -ddc(iq,nc),              &
         -2.*cu(iq)/(bu(iq)+sqrt(max(bu(iq)**2-4.*yyc(iq,nc)*cu(iq),0.1))) )
    end do  
    
    ! sea-ice (cavitating fluid)
    do iq = isc,iec
      dumc(iqx(iq,nc),2) = max(0.,min(ipmaxc(iq,nc),                       &
         ( -zzncice(iq,nc)*dumc_n(iq,2) - zzscice(iq,nc)*dumc_s(iq,2)      &
           -zzecice(iq,nc)*dumc_e(iq,2) - zzwcice(iq,nc)*dumc_w(iq,2)      &
          + rhscice(iq,nc) ) / zzcice(iq,nc) ))
    end do  

    call bounds_colour_recv(dumc(:,1:2),nc)
    
  end do
end do

call END_LOG(mgsetup_end)

! Main loop
do itr = 2,itr_mgice
   
  call START_LOG(mgfine_begin)
  
  do i = 1,itrbgn
      
    ! store previous solution to test for convergence  
    neta(1:ifull)  = dumc(1:ifull,1)
    ipice(1:ifull) = dumc(1:ifull,2)
  
    do nc = 1,maxcolour

      do k = 1,2
        do iq = 1,ifull_colour(nc)  
          dumc_n(iq,k) = dumc(iqn(iq,nc),k)
          dumc_s(iq,k) = dumc(iqs(iq,nc),k)
          dumc_e(iq,k) = dumc(iqe(iq,nc),k)
          dumc_w(iq,k) = dumc(iqw(iq,nc),k)
        end do
      end do  
      
      ! update halo
      isc = 1
      iec = ifull_colour_border(nc)
  
      ! ocean
      bu(isc:iec) = yync(isc:iec,nc)*dumc_n(isc:iec,1) + yysc(isc:iec,nc)*dumc_s(isc:iec,1)      &
                  + yyec(isc:iec,nc)*dumc_e(isc:iec,1) + yywc(isc:iec,nc)*dumc_w(isc:iec,1)      &
                  + zzhhc(isc:iec,nc)
      cu(isc:iec) = zznc(isc:iec,nc)*dumc_n(isc:iec,1) + zzsc(isc:iec,nc)*dumc_s(isc:iec,1)      &
                  + zzec(isc:iec,nc)*dumc_e(isc:iec,1) + zzwc(isc:iec,nc)*dumc_w(isc:iec,1)      &
                  - rhsc(isc:iec,nc)   
      do iq = isc,iec
        dumc(iqx(iq,nc),1) = eec(iq,nc)*max( -ddc(iq,nc),           &
           -2.*cu(iq)/(bu(iq)+sqrt(max(bu(iq)**2-4.*yyc(iq,nc)*cu(iq),0.1))) )
      end do
        
      ! ice (cavitating fluid)
      do iq = isc,iec
        dumc(iqx(iq,nc),2) = max(0.,min(ipmaxc(iq,nc),                       &
           ( -zzncice(iq,nc)*dumc_n(iq,2) - zzscice(iq,nc)*dumc_s(iq,2)      &
             -zzecice(iq,nc)*dumc_e(iq,2) - zzwcice(iq,nc)*dumc_w(iq,2)      &
            + rhscice(iq,nc) ) / zzcice(iq,nc) ))
      end do  

      call bounds_colour_send(dumc(:,1:2),nc)
    
      ! update interior
      isc = ifull_colour_border(nc) + 1
      iec = ifull_colour(nc)

      ! ocean
      bu(isc:iec)=yync(isc:iec,nc)*dumc_n(isc:iec,1) + yysc(isc:iec,nc)*dumc_s(isc:iec,1)      &
                + yyec(isc:iec,nc)*dumc_e(isc:iec,1) + yywc(isc:iec,nc)*dumc_w(isc:iec,1)      &
                + zzhhc(isc:iec,nc)
      cu(isc:iec)=zznc(isc:iec,nc)*dumc_n(isc:iec,1)+zzsc(isc:iec,nc)*dumc_s(isc:iec,1)      &
                 +zzec(isc:iec,nc)*dumc_e(isc:iec,1)+zzwc(isc:iec,nc)*dumc_w(isc:iec,1)      &
                 -rhsc(isc:iec,nc)
      do iq = isc,iec
        dumc(iqx(iq,nc),1) = eec(iq,nc)*max( -ddc(iq,nc),           &
           -2.*cu(iq)/(bu(iq)+sqrt(max(bu(iq)**2-4.*yyc(iq,nc)*cu(iq),0.1))) )
      end do
        
      ! ice (cavitating fluid)
      do iq = isc,iec
        dumc(iqx(iq,nc),2) = max(0.,min(ipmaxc(iq,nc),                       &
           ( -zzncice(iq,nc)*dumc_n(iq,2) - zzscice(iq,nc)*dumc_s(iq,2)      &
             -zzecice(iq,nc)*dumc_e(iq,2) - zzwcice(iq,nc)*dumc_w(iq,2)      &
            + rhscice(iq,nc) ) / zzcice(iq,nc) ))
      end do  
      
      call bounds_colour_recv(dumc(:,1:2),nc)
    
    end do
  end do
  
  ! test for convergence
  dsolmax_g = 0.
  dsolmax_g(1)   = maxval( abs( dumc(1:ifull,1) - neta(1:ifull) ) )
  dsolmax_g(2)   = maxval( abs( dumc(1:ifull,2) - ipice(1:ifull) ) )
  neta(1:ifull)  = dumc(1:ifull,1)
  ipice(1:ifull) = dumc(1:ifull,2)  

  ! residual - ocean
  call unpack_nsew(dumc(:,1),dumc_n(:,1),dumc_s(:,1),dumc_e(:,1),dumc_w(:,1))
  w(1:ifull,1)=(-neta(1:ifull)*(iyy(1:ifull)*neta(1:ifull)+iyyn(1:ifull)*dumc_n(1:ifull,1)  &
                                                          +iyys(1:ifull)*dumc_s(1:ifull,1)  &
                                                          +iyye(1:ifull)*dumc_e(1:ifull,1)  &
                                                          +iyyw(1:ifull)*dumc_w(1:ifull,1)) &
               -(izz(1:ifull,1)*neta(1:ifull)                                               &
                +izzn(1:ifull,1)*dumc_n(1:ifull,1)+izzs(1:ifull,1)*dumc_s(1:ifull,1)        &
                +izze(1:ifull,1)*dumc_e(1:ifull,1)+izzw(1:ifull,1)*dumc_w(1:ifull,1))       &
               -ihh(1:ifull)*neta(1:ifull)+irhs(1:ifull,1))*ee(1:ifull)
  
  w(1:ifull,2)  =  izz(1:ifull,1) +  iyy(1:ifull)*neta(1:ifull)
  w(1:ifull,3)  = izzn(1:ifull,1) + iyyn(1:ifull)*neta(1:ifull)
  w(1:ifull,4)  = izzs(1:ifull,1) + iyys(1:ifull)*neta(1:ifull)
  w(1:ifull,5)  = izze(1:ifull,1) + iyye(1:ifull)*neta(1:ifull)
  w(1:ifull,6)  = izzw(1:ifull,1) + iyyw(1:ifull)*neta(1:ifull)
  w(1:ifull,7)  = iyyn*dumc_n(1:ifull,1) + iyys*dumc_s(1:ifull,1)     &
                + iyye*dumc_e(1:ifull,1) + iyyw*dumc_w(1:ifull,1)     &
                + iyy*neta(1:ifull) + ihh
  
  ! residual ice
  call unpack_nsew(dumc(:,2),dumc_n(:,2),dumc_s(:,2),dumc_e(:,2),dumc_w(:,2))
  w(1:ifull,8)  =(- izzn(1:ifull,2)*dumc_n(1:ifull,2) - izzs(1:ifull,2)*dumc_s(1:ifull,2)     &
                  - izze(1:ifull,2)*dumc_e(1:ifull,2) - izzw(1:ifull,2)*dumc_w(1:ifull,2)     &
                  - izz(1:ifull,2)*ipice(1:ifull) + irhs(1:ifull,2))*ee(1:ifull)
  where ( ipice(1:ifull)>=ipmax(1:ifull) )
    w(1:ifull,8) = 0. ! improves convergence
  end where
  
  ! For when the inital grid cannot be upscaled
  call mgcollect(1,w(:,1:8),dsolmax_g(1:8))
  
  call END_LOG(mgfine_end)
  
  
  if ( mg_maxlevel_local>0 ) then
  
    call START_LOG(mgfine_begin)  
  
    
    ! restriction
    ! (since this always operates within a panel, then ine = ien is always true)
    ng4 = mg(1)%ifull_fine
    do iq = 1,ng4
      rhs(iq,2)=0.25*(w(mg(1)%fine(iq)  ,1)+w(mg(1)%fine_n(iq) ,1)       &
                     +w(mg(1)%fine_e(iq),1)+w(mg(1)%fine_ne(iq),1))
      zz(iq,2) =0.25*dfac*(w(mg(1)%fine(iq)  ,2)+w(mg(1)%fine_n(iq) ,2)   &
                          +w(mg(1)%fine_e(iq),2)+w(mg(1)%fine_ne(iq),2))
      zzn(iq,2)=0.25*dfac*(w(mg(1)%fine(iq)  ,3)+w(mg(1)%fine_n(iq) ,3)   &
                          +w(mg(1)%fine_e(iq),3)+w(mg(1)%fine_ne(iq),3))
      zzs(iq,2)=0.25*dfac*(w(mg(1)%fine(iq)  ,4)+w(mg(1)%fine_n(iq) ,4)   &
                          +w(mg(1)%fine_e(iq),4)+w(mg(1)%fine_ne(iq),4))
      zze(iq,2)=0.25*dfac*(w(mg(1)%fine(iq)  ,5)+w(mg(1)%fine_n(iq) ,5)   &
                          +w(mg(1)%fine_e(iq),5)+w(mg(1)%fine_ne(iq),5))
      zzw(iq,2)=0.25*dfac*(w(mg(1)%fine(iq)  ,6)+w(mg(1)%fine_n(iq) ,6)   &
                           +w(mg(1)%fine_e(iq),6)+w(mg(1)%fine_ne(iq),6))
      hh(iq,2)    =0.25*(w(mg(1)%fine(iq)  ,7)+w(mg(1)%fine_n(iq) ,7)     &
                        +w(mg(1)%fine_e(iq),7)+w(mg(1)%fine_ne(iq),7))
      
      rhsi(iq,2)=0.25*(w(mg(1)%fine(iq)  ,8)+w(mg(1)%fine_n(iq) ,8)      &
                      +w(mg(1)%fine_e(iq),8)+w(mg(1)%fine_ne(iq),8))
    end do  

    ! merge grids if insufficent points on this processor
    if ( mg(2)%merge_len>1 ) then
      w(1:ng4,1)  = rhs(1:ng4,2)
      w(1:ng4,2)  = zz(1:ng4,2)
      w(1:ng4,3)  = zzn(1:ng4,2)
      w(1:ng4,4)  = zzs(1:ng4,2)
      w(1:ng4,5)  = zze(1:ng4,2)
      w(1:ng4,6)  = zzw(1:ng4,2)
      w(1:ng4,7)  = hh(1:ng4,2)
      w(1:ng4,8)  = rhsi(1:ng4,2)
      call mgcollect(2,w(:,1:8),dsolmax_g(1:8))
      if ( 2<=mg_maxlevel_local ) then
        ng = mg(2)%ifull
        rhs(1:ng,2)  = w(1:ng,1)
        zz(1:ng,2)   = w(1:ng,2)
        zzn(1:ng,2)  = w(1:ng,3)
        zzs(1:ng,2)  = w(1:ng,4)
        zze(1:ng,2)  = w(1:ng,5)
        zzw(1:ng,2)  = w(1:ng,6)
        hh(1:ng,2)   = w(1:ng,7)
        rhsi(1:ng,2) = w(1:ng,8)
      end if  
    end if
    
    call END_LOG(mgfine_end)

    call START_LOG(mgup_begin)

    ! upscale grid
    do g = 2,gmax
  
      ng = mg(g)%ifull

      ! update
      ! possibly use colours here, although v is reset to zero every iteration
      ! assume zero for first guess of residual (also avoids additional bounds call)
      bu(1:ng) = zz(1:ng,g) + hh(1:ng,g)
      cu(1:ng) = -rhs(1:ng,g)
      v(:,1:2,g) = 0.
      v(1:ng,1,g) = -2.*cu(1:ng)/(bu(1:ng)+sqrt(max(bu(1:ng)**2-4.*yyz(1:ng,g)*cu(1:ng),0.1)))
      v(1:ng,2,g) = rhsi(1:ng,g) / zzi(1:ng,g)
      call mgbounds(g,v(:,1:2,g))
      
      do i = 2,itrbgn
        ! ocean - post smoothing
        call mgunpack_nsew(g,v(:,1,g),v_n,v_s,v_e,v_w)
        bu(1:ng) = yyn(1:ng,g)*v_n(1:ng)+yys(1:ng,g)*v_s(1:ng)     &
                 + yye(1:ng,g)*v_e(1:ng)+yyw(1:ng,g)*v_w(1:ng)     &
                 + zz(1:ng,g) + hh(1:ng,g)
        cu(1:ng) = zzn(1:ng,g)*v_n(1:ng)+zzs(1:ng,g)*v_s(1:ng)     &
                 + zze(1:ng,g)*v_e(1:ng)+zzw(1:ng,g)*v_w(1:ng)     &
                 - rhs(1:ng,g)
        v(1:ng,1,g) = -2.*cu(1:ng)/(bu(1:ng)+sqrt(max(bu(1:ng)**2-4.*yyz(1:ng,g)*cu(1:ng),0.1)))
        
        ! ice - post smoothing
        call mgunpack_nsew(g,v(:,2,g),v_n,v_s,v_e,v_w)
        v(1:ng,2,g) = ( - zzin(1:ng,g)*v_n(1:ng) - zzis(1:ng,g)*v_s(1:ng)     &
                        - zzie(1:ng,g)*v_e(1:ng) - zziw(1:ng,g)*v_w(1:ng)     &
                        + rhsi(1:ng,g) ) / zzi(1:ng,g)
        call mgbounds(g,v(:,1:2,g))
      end do
    
      ! restriction
      ! (calculate finer grid before mgcollect as the messages sent/recv are shorter)

      ng4 = mg(g)%ifull_fine

      ! ocean residual
      call mgunpack_nsew(g,v(:,1,g),v_n,v_s,v_e,v_w)
      ws(1:ng,1)=-(zz(1:ng,g)*v(1:ng,1,g)+zzn(1:ng,g)*v_n(1:ng)+zzs(1:ng,g)*v_s(1:ng)      &
                                         +zze(1:ng,g)*v_e(1:ng)+zzw(1:ng,g)*v_w(1:ng))     &
                  -hh(1:ng,g)*v(1:ng,1,g)+rhs(1:ng,g)
      w(1:ng4,1)=0.25*(ws(mg(g)%fine  ,1)+ws(mg(g)%fine_n ,1) &
                      +ws(mg(g)%fine_e,1)+ws(mg(g)%fine_ne,1))
      
      ws(1:ng,1) = zz(1:ng,g) + yyz(1:ng,g)*v(1:ng,1,g)
      do iq = 1,ng4
        w(iq,2) = 0.25*dfac*(ws(mg(g)%fine(iq)  ,1)+ws(mg(g)%fine_n(iq) ,1) &
                            +ws(mg(g)%fine_e(iq),1)+ws(mg(g)%fine_ne(iq),1))
      end do
      
      ws(1:ng,1) = zzn(1:ng,g) + yyn(1:ng,g)*v(1:ng,1,g)
      do iq = 1,ng4
        w(iq,3) = 0.25*dfac*(ws(mg(g)%fine(iq)  ,1)+ws(mg(g)%fine_n(iq) ,1) &
                            +ws(mg(g)%fine_e(iq),1)+ws(mg(g)%fine_ne(iq),1))
      end do
      
      ws(1:ng,1) = zzs(1:ng,g) + yys(1:ng,g)*v(1:ng,1,g)
      do iq = 1,ng4
        w(iq,4) = 0.25*dfac*(ws(mg(g)%fine(iq)  ,1)+ws(mg(g)%fine_n(iq) ,1) &
                            +ws(mg(g)%fine_e(iq),1)+ws(mg(g)%fine_ne(iq),1))
      end do
      
      ws(1:ng,1) = zze(1:ng,g) + yye(1:ng,g)*v(1:ng,1,g)
      do iq = 1,ng4
        w(iq,5) = 0.25*dfac*(ws(mg(g)%fine(iq)  ,1)+ws(mg(g)%fine_n(iq) ,1) &
                            +ws(mg(g)%fine_e(iq),1)+ws(mg(g)%fine_ne(iq),1))
      end do
      
      ws(1:ng,1) = zzw(1:ng,g) + yyw(1:ng,g)*v(1:ng,1,g)
      do iq = 1,ng4
        w(iq,6) = 0.25*dfac*(ws(mg(g)%fine(iq)  ,1)+ws(mg(g)%fine_n(iq) ,1) &
                            +ws(mg(g)%fine_e(iq),1)+ws(mg(g)%fine_ne(iq),1))
      end do
      
      
      ws(1:ng,1) = yyn(1:ng,g)*v_n(1:ng)+yys(1:ng,g)*v_s(1:ng)   &
                  +yye(1:ng,g)*v_e(1:ng)+yyw(1:ng,g)*v_w(1:ng)   &
                  +yyz(1:ng,g)*v(1:ng,1,g) + hh(1:ng,g)
      do iq = 1,ng4
        w(iq,7) = 0.25*(ws(mg(g)%fine(iq)  ,1)+ws(mg(g)%fine_n(iq) ,1) &
                       +ws(mg(g)%fine_e(iq),1)+ws(mg(g)%fine_ne(iq),1))
      end do  
 
      
      ! ice residual
      call mgunpack_nsew(g,v(:,2,g),v_n,v_s,v_e,v_w)
      ws(1:ng,2) = -zzin(1:ng,g)*v_n(1:ng)-zzis(1:ng,g)*v_s(1:ng)     &
                   -zzie(1:ng,g)*v_e(1:ng)-zziw(1:ng,g)*v_w(1:ng)     &
                   -zzi(1:ng,g)*v(1:ng,2,g)+rhsi(1:ng,g)
      do iq = 1,ng4
        w(iq,8)=0.25*(ws(mg(g)%fine(iq)  ,2)+ws(mg(g)%fine_n(iq) ,2)  &
                     +ws(mg(g)%fine_e(iq),2)+ws(mg(g)%fine_ne(iq),2))
      end do

      ! merge grids if insufficent points on this processor
      if ( mg(g+1)%merge_len>1 ) then
        call mgcollect(g+1,w(:,1:8),dsolmax_g(1:8))
      end if

      if ( g+1<=mg_maxlevel_local ) then
        ng = mg(g+1)%ifull
        rhs(1:ng,g+1)    =w(1:ng,1)
        zz(1:ng,g+1)     =w(1:ng,2)
        zzn(1:ng,g+1)    =w(1:ng,3)
        zzs(1:ng,g+1)    =w(1:ng,4)
        zze(1:ng,g+1)    =w(1:ng,5)
        zzw(1:ng,g+1)    =w(1:ng,6)
        hh(1:ng,g+1)     =w(1:ng,7)
        rhsi(1:ng,g+1)   =w(1:ng,8)
      end if  
      
    end do

    call END_LOG(mgup_end)

    ! solve coarse grid    
    if ( mg_maxlevel==mg_maxlevel_local ) then

      call START_LOG(mgcoarse_begin)  
    
      g = mg_maxlevel
      ng = mg(g)%ifull

      ! pack rhsc_c by colour
      ! pack zz,hh and rhs by colour
      do nc = 1,3
        do iq = 1,mg_ifull_maxcolour  
          ! ocean  
          zzhhcu(iq,nc) = zz(col_iq(iq,nc),g) + hh(col_iq(iq,nc),g)
          zzncu(iq,nc)  = zzn(col_iq(iq,nc),g)
          zzscu(iq,nc)  = zzs(col_iq(iq,nc),g)
          zzecu(iq,nc)  = zze(col_iq(iq,nc),g)
          zzwcu(iq,nc)  = zzw(col_iq(iq,nc),g)
          rhscu(iq,nc)  = rhs(col_iq(iq,nc),g)
          ! ice
          rhsicu(iq,nc) = rhsi(col_iq(iq,nc),g)
        end do  
      end do  
  
      ! solve non-linear water free surface and solve for ice with coloured SOR
      ! first guess
      bu(1:ng) = zz(1:ng,g) + hh(1:ng,g)
      cu(1:ng) = -rhs(1:ng,g)
      v(:,1:2,g) = 0.
      v(1:ng,1,g) = -2.*cu(1:ng)/(bu(1:ng)+sqrt(max(bu(1:ng)**2-4.*yyz(1:ng,g)*cu(1:ng),0.1))) ! ocean
      v(1:ng,2,g) = rhsi(1:ng,g) / zzi(1:ng,g)                                        ! ice
      do itrc = 1,itr_mgice
        ! store previous guess for convegence test
        ws(1:ng,1:2) = v(1:ng,1:2,g)
        do nc = 1,3
        
          ! ocean
          bu(1:mg_ifull_maxcolour) = yyncu(:,nc)*v(col_iqn(:,nc),1,g) + yyscu(:,nc)*v(col_iqs(:,nc),1,g)     &
                                   + yyecu(:,nc)*v(col_iqe(:,nc),1,g) + yywcu(:,nc)*v(col_iqw(:,nc),1,g)     &
                                   + zzhhcu(:,nc)
          cu(1:mg_ifull_maxcolour) = zzncu(:,nc)*v(col_iqn(:,nc),1,g) + zzscu(:,nc)*v(col_iqs(:,nc),1,g)     &
                                   + zzecu(:,nc)*v(col_iqe(:,nc),1,g) + zzwcu(:,nc)*v(col_iqw(:,nc),1,g)     &
                                   - rhscu(:,nc)
          do iq = 1,mg_ifull_maxcolour
            v(col_iq(iq,nc),1,g) = -2.*cu(iq)/(bu(iq)+sqrt(max(bu(iq)**2-4.*yyzcu(iq,nc)*cu(iq),0.1)))
          end do  
      
          ! ice
          do iq = 1,mg_ifull_maxcolour
            v(col_iq(iq,nc),2,g) = ( - zzincu(iq,nc)*v(col_iqn(iq,nc),2,g) - zziscu(iq,nc)*v(col_iqs(iq,nc),2,g)     &
                                     - zziecu(iq,nc)*v(col_iqe(iq,nc),2,g) - zziwcu(iq,nc)*v(col_iqw(iq,nc),2,g)     &
                                    + rhsicu(iq,nc) ) / zzicu(iq,nc)
          end do
          
        end do
        ! test for convergence
        dsolmax(1:2) = maxval( abs( v(1:ng,1:2,g) - ws(1:ng,1:2) ) )
        if ( dsolmax(1)<tol .and. dsolmax(2)<itol ) exit
      end do
      
      call END_LOG(mgcoarse_end)
      
    end if
  
    
    call START_LOG(mgdown_begin)
    
    ! downscale grid
    do g = gmax,2,-1

      call mgbcast(g+1,v(:,1:2,g+1),dsolmax_g(1:2))

      ! interpolation
      ng = mg(g)%ifull
      do k = 1,2
        do iq = 1,ng
          ws(iq,k) =  0.5*v(mg(g+1)%coarse_a(iq),k,g+1)   &
                    + 0.25*v(mg(g+1)%coarse_b(iq),k,g+1)  &
                    + 0.25*v(mg(g+1)%coarse_c(iq),k,g+1)
        end do  
      end do  
      
      ! extension
      ! No mgbounds as the v halo has already been updated and
      ! the coarse interpolation also updates the w halo
      v(1:ng,1:2,g) = v(1:ng,1:2,g) + ws(1:ng,1:2)

      do i = 1,itrend
        call mgbounds(g,v(:,1:2,g))  
        ! ocean - post smoothing
        call mgunpack_nsew(g,v(:,1,g),v_n,v_s,v_e,v_w)
        bu(1:ng) = yyn(1:ng,g)*v_n(1:ng)+yys(1:ng,g)*v_s(1:ng)     &
                 + yye(1:ng,g)*v_e(1:ng)+yyw(1:ng,g)*v_w(1:ng)     &
                 + zz(1:ng,g) + hh(1:ng,g)
        cu(1:ng) = zzn(1:ng,g)*v_n(1:ng)+zzs(1:ng,g)*v_s(1:ng)     &
                 + zze(1:ng,g)*v_e(1:ng)+zzw(1:ng,g)*v_w(1:ng)     &
                 - rhs(1:ng,g)
        v(1:ng,1,g) = -2.*cu(1:ng)/(bu(1:ng)+sqrt(max(bu(1:ng)**2-4.*yyz(1:ng,g)*cu(1:ng),0.1)))
        ! ice - post smoothing
        call mgunpack_nsew(g,v(:,2,g),v_n,v_s,v_e,v_w)  
        v(1:ng,2,g) = ( - zzin(1:ng,g)*v_n(1:ng) - zzis(1:ng,g)*v_s(1:ng)     &
                        - zzie(1:ng,g)*v_e(1:ng) - zziw(1:ng,g)*v_w(1:ng)     &
                        + rhsi(1:ng,g) ) / zzi(1:ng,g)
      end do
      
      call mgbounds(g,v(:,1:2,g)) ! for next mgbcast

    end do

    
    ! fine grid
    call mgbcast(2,v(:,1:2,2),dsolmax_g(1:2))

    ! interpolation
    ng = mg(1)%ifull
    do k = 1,2
      do iq = 1,ng
        ws(iq,k) = 0.5*v(mg(2)%coarse_a(iq),k,2)  &
                 + 0.25*v(mg(2)%coarse_b(iq),k,2) &
                 + 0.25*v(mg(2)%coarse_c(iq),k,2)
      end do  
    end do  
    
    call END_LOG(mgdown_end)
    
  end if

  
  call START_LOG(mgfine_begin)

  vduma = 0.
  if ( mg(1)%merge_len>1 ) then
    call mgbcast(1,ws(:,1:2),dsolmax_g(1:2),nobounds=.true.)
    ir=mod(mg(1)%merge_pos-1,mg(1)%merge_row)+1   ! index for proc row
    ic=(mg(1)%merge_pos-1)/mg(1)%merge_row+1      ! index for proc col
    do n=1,npan
      do j=1,jpan
        iq_a=1+(j-1)*ipan+(n-1)*ipan*jpan
        iq_c=1+(ir-1)*ipan+(j-1+(ic-1)*jpan)*ipan*mg(1)%merge_row+(n-1)*ipan*jpan*mg(1)%merge_len
        vduma(iq_a:iq_a+ipan-1,1:2)=ws(iq_c:iq_c+ipan-1,1:2)
      end do
    end do
  else
    ! remap mg halo to normal halo 
    vduma(1:ifull,1:2) = ws(1:ifull,1:2)
  end if
 
  ! extension
  vduma(1:ifull,1) = max( -10., min( 10., vduma(1:ifull,1) ) )
  neta(1:ifull) = max( neta(1:ifull)+vduma(1:ifull,1), -dd(1:ifull) )*ee(1:ifull)
  ipice(1:ifull) = max( min( ipice(1:ifull)+vduma(1:ifull,2), ipmax(1:ifull) ), 0. )*ee(1:ifull)
  
  ! update fine spatial scales
  dumc(1:ifull,1) = neta(1:ifull)
  dumc(1:ifull,2) = ipice(1:ifull)
  call bounds(dumc(:,1:2))
  
  do i = 1,itrend
    do nc = 1,maxcolour

      do k = 1,2
        do iq = 1,ifull_colour(nc)  
          dumc_n(iq,k) = dumc(iqn(iq,nc),k)
          dumc_s(iq,k) = dumc(iqs(iq,nc),k)
          dumc_e(iq,k) = dumc(iqe(iq,nc),k)
          dumc_w(iq,k) = dumc(iqw(iq,nc),k)
        end do
      end do  
      
      ! update halo
      isc = 1
      iec = ifull_colour_border(nc)
  
      ! ocean
      bu(isc:iec)=yync(isc:iec,nc)*dumc_n(isc:iec,1)+yysc(isc:iec,nc)*dumc_s(isc:iec,1)      &
                 +yyec(isc:iec,nc)*dumc_e(isc:iec,1)+yywc(isc:iec,nc)*dumc_w(isc:iec,1)      &
                 +zzhhc(isc:iec,nc)
      cu(isc:iec)=zznc(isc:iec,nc)*dumc_n(isc:iec,1)+zzsc(isc:iec,nc)*dumc_s(isc:iec,1)      &
                 +zzec(isc:iec,nc)*dumc_e(isc:iec,1)+zzwc(isc:iec,nc)*dumc_w(isc:iec,1)      &
                 -rhsc(isc:iec,nc)    
      do iq = isc,iec
        dumc(iqx(iq,nc),1) = eec(iq,nc)*max( -ddc(iq,nc),              &
           -2.*cu(iq)/(bu(iq)+sqrt(max(bu(iq)**2-4.*yyc(iq,nc)*cu(iq),0.1))) )
      end do  
    
      ! ice (cavitating fluid)
      do iq = isc,iec
        dumc(iqx(iq,nc),2) = max(0.,min(ipmaxc(iq,nc),                       &
           ( -zzncice(iq,nc)*dumc_n(iq,2) - zzscice(iq,nc)*dumc_s(iq,2)      &
             -zzecice(iq,nc)*dumc_e(iq,2) - zzwcice(iq,nc)*dumc_w(iq,2)      &
            + rhscice(iq,nc) ) / zzcice(iq,nc) ))
      end do  

      call bounds_colour_send(dumc,nc)
    
      ! update interior
      isc = ifull_colour_border(nc) + 1
      iec = ifull_colour(nc)

      ! ocean
      bu(isc:iec)=yync(isc:iec,nc)*dumc_n(isc:iec,1)+yysc(isc:iec,nc)*dumc_s(isc:iec,1)      &
                 +yyec(isc:iec,nc)*dumc_e(isc:iec,1)+yywc(isc:iec,nc)*dumc_w(isc:iec,1)      &
                 +zzhhc(isc:iec,nc)
      cu(isc:iec)=zznc(isc:iec,nc)*dumc_n(isc:iec,1)+zzsc(isc:iec,nc)*dumc_s(isc:iec,1)      &
                 +zzec(isc:iec,nc)*dumc_e(isc:iec,1)+zzwc(isc:iec,nc)*dumc_w(isc:iec,1)      &
                 -rhsc(isc:iec,nc)   
      do iq = isc,iec
        dumc(iqx(iq,nc),1) = eec(iq,nc)*max( -ddc(iq,nc),              &
           -2.*cu(iq)/(bu(iq)+sqrt(max(bu(iq)**2-4.*yyc(iq,nc)*cu(iq),0.1))) )
      end do
      
      ! ice (cavitating fluid)
      do iq = isc,iec
        dumc(iqx(iq,nc),2) = max(0.,min(ipmaxc(iq,nc),                       &
           ( -zzncice(iq,nc)*dumc_n(iq,2) - zzscice(iq,nc)*dumc_s(iq,2)      &
             -zzecice(iq,nc)*dumc_e(iq,2) - zzwcice(iq,nc)*dumc_w(iq,2)      &
            + rhscice(iq,nc) ) / zzcice(iq,nc) ))
      end do

      call bounds_colour_recv(dumc,nc)
    
    end do
  end do
  
  call END_LOG(mgfine_end)
 
  ! test for convergence
  if ( dsolmax_g(1)<tol .and. dsolmax_g(2)<itol ) exit
  
end do

neta(1:ifull+iextra)  = dumc(1:ifull+iextra,1)
ipice(1:ifull+iextra) = dumc(1:ifull+iextra,2)

totits      = itr
maxglobseta = dsolmax_g(1)
maxglobip   = dsolmax_g(2)

call END_LOG(waterhelm_end)

return
end subroutine mgmlo


subroutine optmx(il,schmidt,dt,phibar,accel)
   
! Calculate the optimum acceleration factor for SOR on the cubic-conformal
! grid. Linearly interpolate tabulated values as a function of lambda,
! resolution and stretching. Tabulated values were calculated using
! adaptive SOR.

implicit none

integer, intent(in) :: il    ! Number of grid points on a panel
real, intent(in) :: schmidt  ! Stretch factor
real, intent(in) :: dt       ! Time step
real, intent(in) :: phibar   ! Equivalent depth
real, intent(out) :: accel
real, parameter :: radsq=6.371e6**2

integer :: i, j, k
real :: loglam, logstr, wl, wr, ws, schmidt_min
integer, parameter :: nlam=38, nres=14, nsch=9
real, dimension(nlam,nres,nsch), save :: omega
real, parameter, dimension(nres) :: resol = (/  20, 30, 40, 50, 60, 70, &
     80, 90, 100, 120, 140, 160, 180, 200 /)

data (omega(i, 1, 1),i=1,nlam) /  &
      1.8525,1.8302,1.7941,1.7407,1.6672,1.5785,1.4818,1.3834,    &
      1.2855,1.1957,1.1156,1.0611,1.0279,1.0000,1.0000,1.0000,    &
      1.0000,1.0000,1.0000,1.0000,1.0000,1.0000,1.0000,1.0000,    &
      1.0000,1.0000,1.0000,1.0000,1.0000,1.0000,1.0000,1.0000,    &
      1.0000,1.0000,1.0000,1.0000,1.0000,1.0000 /
data (omega(i, 2, 1),i=1,nlam) /  &
      1.9035,1.8877,1.8612,1.8222,1.7676,1.6982,1.6177,1.5296,    &
      1.4358,1.3382,1.2442,1.1550,1.0888,1.0440,1.0189,1.0000,    &
      1.0000,1.0000,1.0000,1.0000,1.0000,1.0000,1.0000,1.0000,    &
      1.0000,1.0000,1.0000,1.0000,1.0000,1.0000,1.0000,1.0000,    &
      1.0000,1.0000,1.0000,1.0000,1.0000,1.0000 /
data (omega(i, 3, 1),i=1,nlam) /  &
      1.9301,1.9179,1.8974,1.8662,1.8226,1.7658,1.6995,1.6241,    &
      1.5378,1.4453,1.3481,1.2537,1.1682,1.0948,1.0478,1.0208,    &
      1.0000,1.0000,1.0000,1.0000,1.0000,1.0000,1.0000,1.0000,    &
      1.0000,1.0000,1.0000,1.0000,1.0000,1.0000,1.0000,1.0000,    &
      1.0000,1.0000,1.0000,1.0000,1.0000,1.0000 /
data (omega(i, 4, 1),i=1,nlam) /  &
      1.9463,1.9362,1.9198,1.8942,1.8574,1.8095,1.7521,1.6873,    &
      1.6121,1.5246,1.4322,1.3345,1.2408,1.1523,1.0868,1.0428,    &
      1.0183,1.0000,1.0000,1.0000,1.0000,1.0000,1.0000,1.0000,    &
      1.0000,1.0000,1.0000,1.0000,1.0000,1.0000,1.0000,1.0000,    &
      1.0000,1.0000,1.0000,1.0000,1.0000,1.0000 /
data (omega(i, 5, 1),i=1,nlam) /  &
      1.9569,1.9487,1.9347,1.9129,1.8814,1.8400,1.7899,1.7322,    &
      1.6655,1.5852,1.4972,1.4017,1.3059,1.2144,1.1306,1.0714,    &
      1.0337,1.0139,1.0000,1.0000,1.0000,1.0000,1.0000,1.0000,    &
      1.0000,1.0000,1.0000,1.0000,1.0000,1.0000,1.0000,1.0000,    &
      1.0000,1.0000,1.0000,1.0000,1.0000,1.0000 /
data (omega(i, 6, 1),i=1,nlam) /  &
      1.9645,1.9575,1.9457,1.9269,1.8995,1.8624,1.8180,1.7663,    &
      1.7063,1.6342,1.5502,1.4599,1.3627,1.2676,1.1803,1.1039,    &
      1.0535,1.0238,1.0000,1.0000,1.0000,1.0000,1.0000,1.0000,    &
      1.0000,1.0000,1.0000,1.0000,1.0000,1.0000,1.0000,1.0000,    &
      1.0000,1.0000,1.0000,1.0000,1.0000,1.0000 /
data (omega(i, 7, 1),i=1,nlam) /  &
      1.9701,1.9640,1.9538,1.9373,1.9129,1.8797,1.8396,1.7932,    &
      1.7382,1.6723,1.5952,1.5066,1.4117,1.3160,1.2237,1.1381,    &
      1.0767,1.0367,1.0000,1.0000,1.0000,1.0000,1.0000,1.0000,    &
      1.0000,1.0000,1.0000,1.0000,1.0000,1.0000,1.0000,1.0000,    &
      1.0000,1.0000,1.0000,1.0000,1.0000,1.0000 /
data (omega(i, 8, 1),i=1,nlam) /  &
      1.9741,1.9691,1.9600,1.9456,1.9238,1.8938,1.8569,1.8145,    &
      1.7635,1.7032,1.6313,1.5470,1.4567,1.3594,1.2645,1.1776,    &
      1.1019,1.0522,1.0231,1.0000,1.0000,1.0000,1.0000,1.0000,    &
      1.0000,1.0000,1.0000,1.0000,1.0000,1.0000,1.0000,1.0000,    &
      1.0000,1.0000,1.0000,1.0000,1.0000,1.0000 /
data (omega(i, 9, 1),i=1,nlam) /  &
      1.9775,1.9727,1.9650,1.9521,1.9323,1.9049,1.8715,1.8321,    &
      1.7861,1.7292,1.6617,1.5834,1.4946,1.3983,1.3025,1.2114,    &
      1.1282,1.0697,1.0327,1.0135,1.0000,1.0000,1.0000,1.0000,    &
      1.0000,1.0000,1.0000,1.0000,1.0000,1.0000,1.0000,1.0000,    &
      1.0000,1.0000,1.0000,1.0000,1.0000,1.0000 /
data (omega(i,10, 1),i=1,nlam) /  &
      1.9825,1.9786,1.9721,1.9617,1.9455,1.9219,1.8931,1.8597,    &
      1.8194,1.7696,1.7100,1.6395,1.5579,1.4668,1.3697,1.2744,    &
      1.1863,1.1085,1.0565,1.0254,1.0000,1.0000,1.0000,1.0000,    &
      1.0000,1.0000,1.0000,1.0000,1.0000,1.0000,1.0000,1.0000,    &
      1.0000,1.0000,1.0000,1.0000,1.0000,1.0000 /
data (omega(i,11, 1),i=1,nlam) /  &
      1.9859,1.9828,1.9773,1.9683,1.9547,1.9345,1.9089,1.8795,    &
      1.8440,1.8003,1.7465,1.6824,1.6074,1.5212,1.4282,1.3306,    &
      1.2373,1.1493,1.0847,1.0415,1.0000,1.0000,1.0000,1.0000,    &
      1.0000,1.0000,1.0000,1.0000,1.0000,1.0000,1.0000,1.0000,    &
      1.0000,1.0000,1.0000,1.0000,1.0000,1.0000 /
data (omega(i,12, 1),i=1,nlam) /  &
      1.9880,1.9857,1.9811,1.9734,1.9617,1.9441,1.9209,1.8947,    &
      1.8629,1.8237,1.7751,1.7166,1.6472,1.5669,1.4767,1.3799,    &
      1.2844,1.1951,1.1153,1.0610,1.0278,1.0000,1.0000,1.0000,    &
      1.0000,1.0000,1.0000,1.0000,1.0000,1.0000,1.0000,1.0000,    &
      1.0000,1.0000,1.0000,1.0000,1.0000,1.0000 /
data (omega(i,13, 1),i=1,nlam) /  &
      1.9893,1.9876,1.9840,1.9773,1.9667,1.9514,1.9307,1.9066,    &
      1.8780,1.8424,1.7983,1.7443,1.6798,1.6047,1.5182,1.4251,    &
      1.3274,1.2343,1.1469,1.0829,1.0404,1.0000,1.0000,1.0000,    &
      1.0000,1.0000,1.0000,1.0000,1.0000,1.0000,1.0000,1.0000,    &
      1.0000,1.0000,1.0000,1.0000,1.0000,1.0000 /
data (omega(i,14, 1),i=1,nlam) /  &
      1.9908,1.9894,1.9862,1.9803,1.9710,1.9573,1.9383,1.9163,    &
      1.8901,1.8576,1.8171,1.7671,1.7072,1.6364,1.5539,1.4634,    &
      1.3663,1.2712,1.1835,1.1063,1.0551,1.0246,1.0000,1.0000,    &
      1.0000,1.0000,1.0000,1.0000,1.0000,1.0000,1.0000,1.0000,    &
      1.0000,1.0000,1.0000,1.0000,1.0000,1.0000 /
data (omega(i, 1, 2),i=1,nlam) /  &
      1.9333,1.9063,1.8777,1.8402,1.7914,1.7271,1.6457,1.5476,    &
      1.4419,1.3413,1.2440,1.1527,1.0873,1.0433,1.0188,1.0000,    &
      1.0000,1.0000,1.0000,1.0000,1.0000,1.0000,1.0000,1.0000,    &
      1.0000,1.0000,1.0000,1.0000,1.0000,1.0000,1.0000,1.0000,    &
      1.0000,1.0000,1.0000,1.0000,1.0000,1.0000 /
data (omega(i, 2, 2),i=1,nlam) /  &
      1.9563,1.9415,1.9223,1.8967,1.8615,1.8148,1.7515,1.6735,    &
      1.5849,1.4908,1.3945,1.2956,1.2047,1.1212,1.0652,1.0304,    &
      1.0000,1.0000,1.0000,1.0000,1.0000,1.0000,1.0000,1.0000,    &
      1.0000,1.0000,1.0000,1.0000,1.0000,1.0000,1.0000,1.0000,    &
      1.0000,1.0000,1.0000,1.0000,1.0000,1.0000 /
data (omega(i, 3, 2),i=1,nlam) /  &
      1.9695,1.9588,1.9448,1.9256,1.8989,1.8617,1.8122,1.7480,    &
      1.6733,1.5907,1.4997,1.4047,1.3055,1.2138,1.1284,1.0701,    &
      1.0332,1.0000,1.0000,1.0000,1.0000,1.0000,1.0000,1.0000,    &
      1.0000,1.0000,1.0000,1.0000,1.0000,1.0000,1.0000,1.0000,    &
      1.0000,1.0000,1.0000,1.0000,1.0000,1.0000 /
data (omega(i, 4, 2),i=1,nlam) /  &
      1.9773,1.9691,1.9582,1.9433,1.9219,1.8916,1.8499,1.7956,    &
      1.7315,1.6600,1.5782,1.4860,1.3909,1.2920,1.2016,1.1188,    &
      1.0636,1.0296,1.0000,1.0000,1.0000,1.0000,1.0000,1.0000,    &
      1.0000,1.0000,1.0000,1.0000,1.0000,1.0000,1.0000,1.0000,    &
      1.0000,1.0000,1.0000,1.0000,1.0000,1.0000 /
data (omega(i, 5, 2),i=1,nlam) /  &
      1.9824,1.9758,1.9670,1.9549,1.9373,1.9120,1.8767,1.8295,    &
      1.7729,1.7095,1.6355,1.5517,1.4583,1.3617,1.2640,1.1769,    &
      1.1001,1.0513,1.0229,1.0000,1.0000,1.0000,1.0000,1.0000,    &
      1.0000,1.0000,1.0000,1.0000,1.0000,1.0000,1.0000,1.0000,    &
      1.0000,1.0000,1.0000,1.0000,1.0000,1.0000 /
data (omega(i, 6, 2),i=1,nlam) /  &
      1.9854,1.9805,1.9732,1.9627,1.9481,1.9265,1.8957,1.8542,    &
      1.8046,1.7471,1.6807,1.6025,1.5148,1.4193,1.3200,1.2271,    &
      1.1390,1.0776,1.0376,1.0000,1.0000,1.0000,1.0000,1.0000,    &
      1.0000,1.0000,1.0000,1.0000,1.0000,1.0000,1.0000,1.0000,    &
      1.0000,1.0000,1.0000,1.0000,1.0000,1.0000 /
data (omega(i, 7, 2),i=1,nlam) /  &
      1.9880,1.9839,1.9777,1.9688,1.9562,1.9373,1.9107,1.8739,    &
      1.8283,1.7765,1.7159,1.6435,1.5609,1.4685,1.3721,1.2739,    &
      1.1856,1.1065,1.0555,1.0251,1.0000,1.0000,1.0000,1.0000,    &
      1.0000,1.0000,1.0000,1.0000,1.0000,1.0000,1.0000,1.0000,    &
      1.0000,1.0000,1.0000,1.0000,1.0000,1.0000 /
data (omega(i, 8, 2),i=1,nlam) /  &
      1.9893,1.9861,1.9812,1.9734,1.9625,1.9460,1.9221,1.8889,    &
      1.8481,1.8009,1.7450,1.6768,1.5995,1.5104,1.4161,1.3168,    &
      1.2241,1.1367,1.0759,1.0366,1.0000,1.0000,1.0000,1.0000,    &
      1.0000,1.0000,1.0000,1.0000,1.0000,1.0000,1.0000,1.0000,    &
      1.0000,1.0000,1.0000,1.0000,1.0000,1.0000 /
data (omega(i, 9, 2),i=1,nlam) /  &
      1.9909,1.9876,1.9838,1.9770,1.9674,1.9527,1.9316,1.9016,    &
      1.8638,1.8203,1.7687,1.7057,1.6322,1.5483,1.4550,1.3583,    &
      1.2607,1.1741,1.0980,1.0500,1.0222,1.0000,1.0000,1.0000,    &
      1.0000,1.0000,1.0000,1.0000,1.0000,1.0000,1.0000,1.0000,    &
      1.0000,1.0000,1.0000,1.0000,1.0000,1.0000 /
data (omega(i,10, 2),i=1,nlam) /  &
      1.9922,1.9905,1.9873,1.9823,1.9742,1.9627,1.9451,1.9201,    &
      1.8878,1.8506,1.8055,1.7511,1.6856,1.6081,1.5202,1.4263,    &
      1.3270,1.2336,1.1444,1.0814,1.0398,1.0000,1.0000,1.0000,    &
      1.0000,1.0000,1.0000,1.0000,1.0000,1.0000,1.0000,1.0000,    &
      1.0000,1.0000,1.0000,1.0000,1.0000,1.0000 /
data (omega(i,11, 2),i=1,nlam) /  &
      1.9938,1.9918,1.9894,1.9855,1.9793,1.9697,1.9548,1.9338,    &
      1.9058,1.8732,1.8327,1.7845,1.7254,1.6545,1.5739,1.4828,    &
      1.3869,1.2882,1.1983,1.1162,1.0619,1.0286,1.0000,1.0000,    &
      1.0000,1.0000,1.0000,1.0000,1.0000,1.0000,1.0000,1.0000,    &
      1.0000,1.0000,1.0000,1.0000,1.0000,1.0000 /
data (omega(i,12, 2),i=1,nlam) /  &
      1.9935,1.9933,1.9906,1.9881,1.9829,1.9748,1.9621,1.9437,    &
      1.9194,1.8901,1.8542,1.8104,1.7573,1.6931,1.6166,1.5318,    &
      1.4366,1.3396,1.2432,1.1593,1.0871,1.0433,1.0189,1.0000,    &
      1.0000,1.0000,1.0000,1.0000,1.0000,1.0000,1.0000,1.0000,    &
      1.0000,1.0000,1.0000,1.0000,1.0000,1.0000 /
data (omega(i,13, 2),i=1,nlam) /  &
      1.9945,1.9933,1.9921,1.9895,1.9857,1.9787,1.9676,1.9513,    &
      1.9298,1.9034,1.8715,1.8317,1.7823,1.7233,1.6519,1.5714,    &
      1.4796,1.3837,1.2851,1.1955,1.1141,1.0605,1.0278,1.0000,    &
      1.0000,1.0000,1.0000,1.0000,1.0000,1.0000,1.0000,1.0000,    &
      1.0000,1.0000,1.0000,1.0000,1.0000,1.0000 /
data (omega(i,14, 2),i=1,nlam) /  &
      1.9953,1.9943,1.9932,1.9910,1.9874,1.9814,1.9720,1.9575,    &
      1.9381,1.9146,1.8852,1.8487,1.8033,1.7488,1.6830,1.6051,    &
      1.5189,1.4229,1.3237,1.2305,1.1418,1.0796,1.0387,1.0000,    &
      1.0000,1.0000,1.0000,1.0000,1.0000,1.0000,1.0000,1.0000,    &
      1.0000,1.0000,1.0000,1.0000,1.0000,1.0000 /
data (omega(i, 1, 3),i=1,nlam) /  &
      1.9503,1.9368,1.9213,1.9003,1.8732,1.8363,1.7868,1.7224,    &
      1.6395,1.5399,1.4319,1.3303,1.2333,1.1432,1.0805,1.0393,    &
      1.0169,1.0000,1.0000,1.0000,1.0000,1.0000,1.0000,1.0000,    &
      1.0000,1.0000,1.0000,1.0000,1.0000,1.0000,1.0000,1.0000,    &
      1.0000,1.0000,1.0000,1.0000,1.0000,1.0000 /
data (omega(i, 2, 3),i=1,nlam) /  &
      1.9704,1.9622,1.9518,1.9385,1.9204,1.8947,1.8595,1.8117,    &
      1.7482,1.6675,1.5768,1.4808,1.3835,1.2843,1.1945,1.1128,    &
      1.0597,1.0275,1.0000,1.0000,1.0000,1.0000,1.0000,1.0000,    &
      1.0000,1.0000,1.0000,1.0000,1.0000,1.0000,1.0000,1.0000,    &
      1.0000,1.0000,1.0000,1.0000,1.0000,1.0000 /
data (omega(i, 3, 3),i=1,nlam) /  &
      1.9801,1.9746,1.9669,1.9573,1.9436,1.9248,1.8977,1.8602,    &
      1.8095,1.7435,1.6667,1.5823,1.4914,1.3937,1.2942,1.2034,    &
      1.1196,1.0643,1.0300,1.0000,1.0000,1.0000,1.0000,1.0000,    &
      1.0000,1.0000,1.0000,1.0000,1.0000,1.0000,1.0000,1.0000,    &
      1.0000,1.0000,1.0000,1.0000,1.0000,1.0000 /
data (omega(i, 4, 3),i=1,nlam) /  &
      1.9857,1.9816,1.9758,1.9684,1.9576,1.9429,1.9213,1.8906,    &
      1.8486,1.7931,1.7260,1.6529,1.5693,1.4776,1.3798,1.2807,    &
      1.1915,1.1105,1.0582,1.0267,1.0000,1.0000,1.0000,1.0000,    &
      1.0000,1.0000,1.0000,1.0000,1.0000,1.0000,1.0000,1.0000,    &
      1.0000,1.0000,1.0000,1.0000,1.0000,1.0000 /
data (omega(i, 5, 3),i=1,nlam) /  &
      1.9881,1.9856,1.9814,1.9755,1.9668,1.9547,1.9368,1.9115,    &
      1.8756,1.8272,1.7694,1.7035,1.6276,1.5424,1.4476,1.3506,    &
      1.2529,1.1673,1.0926,1.0467,1.0206,1.0000,1.0000,1.0000,    &
      1.0000,1.0000,1.0000,1.0000,1.0000,1.0000,1.0000,1.0000,    &
      1.0000,1.0000,1.0000,1.0000,1.0000,1.0000 /
data (omega(i, 6, 3),i=1,nlam) /  &
      1.9905,1.9881,1.9850,1.9805,1.9731,1.9627,1.9480,1.9262,    &
      1.8950,1.8530,1.8012,1.7419,1.6741,1.5940,1.5049,1.4084,    &
      1.3086,1.2165,1.1299,1.0713,1.0340,1.0000,1.0000,1.0000,    &
      1.0000,1.0000,1.0000,1.0000,1.0000,1.0000,1.0000,1.0000,    &
      1.0000,1.0000,1.0000,1.0000,1.0000,1.0000 /
data (omega(i, 7, 3),i=1,nlam) /  &
      1.9913,1.9895,1.9872,1.9836,1.9778,1.9688,1.9562,1.9372,    &
      1.9102,1.8724,1.8257,1.7730,1.7098,1.6357,1.5516,1.4579,    &
      1.3609,1.2627,1.1758,1.0987,1.0506,1.0226,1.0000,1.0000,    &
      1.0000,1.0000,1.0000,1.0000,1.0000,1.0000,1.0000,1.0000,    &
      1.0000,1.0000,1.0000,1.0000,1.0000,1.0000 /
data (omega(i, 8, 3),i=1,nlam) /  &
      1.9927,1.9912,1.9893,1.9862,1.9813,1.9735,1.9626,1.9458,    &
      1.9218,1.8882,1.8461,1.7972,1.7399,1.6716,1.5910,1.5003,    &
      1.4051,1.3054,1.2136,1.1276,1.0697,1.0331,1.0000,1.0000,    &
      1.0000,1.0000,1.0000,1.0000,1.0000,1.0000,1.0000,1.0000,    &
      1.0000,1.0000,1.0000,1.0000,1.0000,1.0000 /
data (omega(i, 9, 3),i=1,nlam) /  &
      1.9938,1.9925,1.9901,1.9877,1.9840,1.9772,1.9676,1.9527,    &
      1.9311,1.9011,1.8621,1.8179,1.7642,1.7001,1.6242,1.5388,    &
      1.4443,1.3471,1.2497,1.1646,1.0906,1.0455,1.0200,1.0000,    &
      1.0000,1.0000,1.0000,1.0000,1.0000,1.0000,1.0000,1.0000,    &
      1.0000,1.0000,1.0000,1.0000,1.0000,1.0000 /
data (omega(i,10, 3),i=1,nlam) /  &
      1.9953,1.9933,1.9925,1.9900,1.9870,1.9825,1.9749,1.9628,    &
      1.9450,1.9199,1.8871,1.8487,1.8019,1.7460,1.6791,1.5998,    &
      1.5102,1.4154,1.3156,1.2229,1.1350,1.0749,1.0361,1.0000,    &
      1.0000,1.0000,1.0000,1.0000,1.0000,1.0000,1.0000,1.0000,    &
      1.0000,1.0000,1.0000,1.0000,1.0000,1.0000 /
data (omega(i,11, 3),i=1,nlam) /  &
      1.9949,1.9947,1.9928,1.9920,1.9897,1.9856,1.9800,1.9699,    &
      1.9550,1.9334,1.9055,1.8719,1.8308,1.7803,1.7198,1.6486,    &
      1.5649,1.4724,1.3758,1.2770,1.1882,1.1080,1.0566,1.0258,    &
      1.0000,1.0000,1.0000,1.0000,1.0000,1.0000,1.0000,1.0000,    &
      1.0000,1.0000,1.0000,1.0000,1.0000,1.0000 /
data (omega(i,12, 3),i=1,nlam) /  &
      1.9958,1.9956,1.9941,1.9921,1.9907,1.9877,1.9832,1.9751,    &
      1.9622,1.9435,1.9194,1.8897,1.8524,1.8077,1.7523,1.6866,    &
      1.6084,1.5217,1.4257,1.3258,1.2324,1.1427,1.0803,1.0393,    &
      1.0000,1.0000,1.0000,1.0000,1.0000,1.0000,1.0000,1.0000,    &
      1.0000,1.0000,1.0000,1.0000,1.0000,1.0000 /
data (omega(i,13, 3),i=1,nlam) /  &
      1.9964,1.9963,1.9950,1.9933,1.9922,1.9897,1.9859,1.9790,    &
      1.9678,1.9513,1.9300,1.9032,1.8702,1.8290,1.7782,1.7177,    &
      1.6459,1.5623,1.4692,1.3725,1.2738,1.1854,1.1060,1.0552,    &
      1.0251,1.0000,1.0000,1.0000,1.0000,1.0000,1.0000,1.0000,    &
      1.0000,1.0000,1.0000,1.0000,1.0000,1.0000 /
data (omega(i,14, 3),i=1,nlam) /  &
      1.9969,1.9969,1.9957,1.9943,1.9920,1.9903,1.9874,1.9820,    &
      1.9722,1.9577,1.9385,1.9147,1.8846,1.8463,1.8005,1.7433,    &
      1.6763,1.5965,1.5070,1.4120,1.3123,1.2198,1.1326,1.0732,    &
      1.0351,1.0000,1.0000,1.0000,1.0000,1.0000,1.0000,1.0000,    &
      1.0000,1.0000,1.0000,1.0000,1.0000,1.0000 /
data (omega(i, 1, 4),i=1,nlam) /  &
      1.9630,1.9551,1.9462,1.9346,1.9197,1.8993,1.8721,1.8352,    &
      1.7857,1.7212,1.6380,1.5380,1.4294,1.3276,1.2305,1.1408,    &
      1.0788,1.0384,1.0000,1.0000,1.0000,1.0000,1.0000,1.0000,    &
      1.0000,1.0000,1.0000,1.0000,1.0000,1.0000,1.0000,1.0000,    &
      1.0000,1.0000,1.0000,1.0000,1.0000,1.0000 /
data (omega(i, 2, 4),i=1,nlam) /  &
      1.9789,1.9741,1.9686,1.9611,1.9512,1.9379,1.9197,1.8941,    &
      1.8589,1.8109,1.7472,1.6662,1.5748,1.4782,1.3807,1.2814,    &
      1.1920,1.1107,1.0583,1.0268,1.0000,1.0000,1.0000,1.0000,    &
      1.0000,1.0000,1.0000,1.0000,1.0000,1.0000,1.0000,1.0000,    &
      1.0000,1.0000,1.0000,1.0000,1.0000,1.0000 /
data (omega(i, 3, 4),i=1,nlam) /  &
      1.9856,1.9831,1.9791,1.9738,1.9666,1.9573,1.9433,1.9244,    &
      1.8972,1.8597,1.8088,1.7425,1.6650,1.5801,1.4889,1.3909,    &
      1.2913,1.2008,1.1175,1.0628,1.0293,1.0000,1.0000,1.0000,    &
      1.0000,1.0000,1.0000,1.0000,1.0000,1.0000,1.0000,1.0000,    &
      1.0000,1.0000,1.0000,1.0000,1.0000,1.0000 /
data (omega(i, 4, 4),i=1,nlam) /  &
      1.9892,1.9872,1.9847,1.9810,1.9759,1.9684,1.9578,1.9428,    &
      1.9209,1.8903,1.8481,1.7923,1.7247,1.6511,1.5670,1.4750,    &
      1.3770,1.2778,1.1889,1.1084,1.0569,1.0260,1.0000,1.0000,    &
      1.0000,1.0000,1.0000,1.0000,1.0000,1.0000,1.0000,1.0000,    &
      1.0000,1.0000,1.0000,1.0000,1.0000,1.0000 /
data (omega(i, 5, 4),i=1,nlam) /  &
      1.9917,1.9902,1.9879,1.9851,1.9815,1.9755,1.9669,1.9546,    &
      1.9369,1.9112,1.8753,1.8266,1.7683,1.7020,1.6256,1.5400,    &
      1.4449,1.3478,1.2501,1.1650,1.0907,1.0456,1.0201,1.0000,    &
      1.0000,1.0000,1.0000,1.0000,1.0000,1.0000,1.0000,1.0000,    &
      1.0000,1.0000,1.0000,1.0000,1.0000,1.0000 /
data (omega(i, 6, 4),i=1,nlam) /  &
      1.9925,1.9918,1.9897,1.9876,1.9850,1.9804,1.9733,1.9630,    &
      1.9480,1.9263,1.8947,1.8527,1.8003,1.7417,1.6724,1.5919,    &
      1.5024,1.4056,1.3057,1.2138,1.1276,1.0698,1.0332,1.0000,    &
      1.0000,1.0000,1.0000,1.0000,1.0000,1.0000,1.0000,1.0000,    &
      1.0000,1.0000,1.0000,1.0000,1.0000,1.0000 /
data (omega(i, 7, 4),i=1,nlam) /  &
      1.9939,1.9925,1.9915,1.9898,1.9873,1.9839,1.9780,1.9691,    &
      1.9562,1.9374,1.9100,1.8720,1.8258,1.7719,1.7083,1.6338,    &
      1.5493,1.4552,1.3581,1.2599,1.1733,1.0968,1.0494,1.0220,    &
      1.0000,1.0000,1.0000,1.0000,1.0000,1.0000,1.0000,1.0000,    &
      1.0000,1.0000,1.0000,1.0000,1.0000,1.0000 /
data (omega(i, 8, 4),i=1,nlam) /  &
      1.9949,1.9938,1.9920,1.9905,1.9886,1.9861,1.9814,1.9738,    &
      1.9625,1.9460,1.9220,1.8880,1.8456,1.7963,1.7387,1.6699,    &
      1.5888,1.4977,1.4023,1.3024,1.2109,1.1254,1.0682,1.0323,    &
      1.0000,1.0000,1.0000,1.0000,1.0000,1.0000,1.0000,1.0000,    &
      1.0000,1.0000,1.0000,1.0000,1.0000,1.0000 /
data (omega(i, 9, 4),i=1,nlam) /  &
      1.9941,1.9937,1.9931,1.9918,1.9902,1.9877,1.9840,1.9775,    &
      1.9675,1.9528,1.9313,1.9009,1.8617,1.8172,1.7631,1.6986,    &
      1.6223,1.5364,1.4416,1.3443,1.2469,1.1623,1.0888,1.0444,    &
      1.0000,1.0000,1.0000,1.0000,1.0000,1.0000,1.0000,1.0000,    &
      1.0000,1.0000,1.0000,1.0000,1.0000,1.0000 /
data (omega(i,10, 4),i=1,nlam) /  &
      1.9954,1.9952,1.9948,1.9937,1.9915,1.9898,1.9874,1.9827,    &
      1.9748,1.9628,1.9452,1.9201,1.8869,1.8482,1.8019,1.7448,    &
      1.6775,1.5976,1.5078,1.4126,1.3127,1.2202,1.1327,1.0733,    &
      1.0352,1.0000,1.0000,1.0000,1.0000,1.0000,1.0000,1.0000,    &
      1.0000,1.0000,1.0000,1.0000,1.0000,1.0000 /
data (omega(i,11, 4),i=1,nlam) /  &
      1.9963,1.9965,1.9951,1.9940,1.9931,1.9919,1.9896,1.9857,    &
      1.9798,1.9698,1.9550,1.9337,1.9053,1.8716,1.8302,1.7793,    &
      1.7185,1.6467,1.5625,1.4697,1.3729,1.2741,1.1856,1.1060,    &
      1.0553,1.0251,1.0000,1.0000,1.0000,1.0000,1.0000,1.0000,    &
      1.0000,1.0000,1.0000,1.0000,1.0000,1.0000 /
data (omega(i,12, 4),i=1,nlam) /  &
      1.9969,1.9958,1.9957,1.9949,1.9943,1.9919,1.9906,1.9883,    &
      1.9834,1.9749,1.9622,1.9438,1.9193,1.8895,1.8526,1.8069,    &
      1.7511,1.6850,1.6078,1.5193,1.4230,1.3258,1.2296,1.1403,    &
      1.0786,1.0383,1.0000,1.0000,1.0000,1.0000,1.0000,1.0000,    &
      1.0000,1.0000,1.0000,1.0000,1.0000,1.0000 /
data (omega(i,13, 4),i=1,nlam) /  &
      1.9974,1.9966,1.9962,1.9954,1.9951,1.9931,1.9921,1.9894,    &
      1.9854,1.9788,1.9678,1.9516,1.9299,1.9031,1.8699,1.8284,    &
      1.7772,1.7163,1.6442,1.5600,1.4665,1.3697,1.2710,1.1829,    &
      1.1039,1.0540,1.0245,1.0000,1.0000,1.0000,1.0000,1.0000,    &
      1.0000,1.0000,1.0000,1.0000,1.0000,1.0000 /
data (omega(i,14, 4),i=1,nlam) /  &
      1.9977,1.9971,1.9967,1.9961,1.9944,1.9940,1.9932,1.9901,    &
      1.9875,1.9818,1.9721,1.9578,1.9384,1.9146,1.8845,1.8466,    &
      1.7996,1.7420,1.6748,1.5943,1.5045,1.4093,1.3093,1.2172,    &
      1.1303,1.0716,1.0342,1.0000,1.0000,1.0000,1.0000,1.0000,    &
      1.0000,1.0000,1.0000,1.0000,1.0000,1.0000 /
data (omega(i, 1, 5),i=1,nlam) /  &
      1.9710,1.9660,1.9605,1.9539,1.9453,1.9341,1.9193,1.8991,    &
      1.8718,1.8350,1.7854,1.7209,1.6376,1.5376,1.4288,1.3269,    &
      1.2299,1.1402,1.0784,1.0381,1.0000,1.0000,1.0000,1.0000,    &
      1.0000,1.0000,1.0000,1.0000,1.0000,1.0000,1.0000,1.0000,    &
      1.0000,1.0000,1.0000,1.0000,1.0000,1.0000 /
data (omega(i, 2, 5),i=1,nlam) /  &
      1.9836,1.9809,1.9777,1.9735,1.9682,1.9611,1.9510,1.9378,    &
      1.9196,1.8939,1.8587,1.8107,1.7469,1.6658,1.5743,1.4776,    &
      1.3800,1.2807,1.1913,1.1102,1.0580,1.0266,1.0000,1.0000,    &
      1.0000,1.0000,1.0000,1.0000,1.0000,1.0000,1.0000,1.0000,    &
      1.0000,1.0000,1.0000,1.0000,1.0000,1.0000 /
data (omega(i, 3, 5),i=1,nlam) /  &
      1.9888,1.9873,1.9851,1.9826,1.9788,1.9737,1.9668,1.9572,    &
      1.9437,1.9243,1.8971,1.8596,1.8086,1.7422,1.6646,1.5795,    &
      1.4883,1.3902,1.2906,1.2001,1.1169,1.0625,1.0291,1.0000,    &
      1.0000,1.0000,1.0000,1.0000,1.0000,1.0000,1.0000,1.0000,    &
      1.0000,1.0000,1.0000,1.0000,1.0000,1.0000 /
data (omega(i, 4, 5),i=1,nlam) /  &
      1.9914,1.9901,1.9889,1.9871,1.9847,1.9809,1.9757,1.9683,    &
      1.9577,1.9428,1.9209,1.8902,1.8480,1.7921,1.7244,1.6507,    &
      1.5664,1.4743,1.3763,1.2771,1.1883,1.1079,1.0565,1.0258,    &
      1.0000,1.0000,1.0000,1.0000,1.0000,1.0000,1.0000,1.0000,    &
      1.0000,1.0000,1.0000,1.0000,1.0000,1.0000 /
data (omega(i, 5, 5),i=1,nlam) /  &
      1.9926,1.9916,1.9908,1.9896,1.9878,1.9854,1.9812,1.9754,    &
      1.9669,1.9546,1.9369,1.9111,1.8753,1.8265,1.7681,1.7016,    &
      1.6251,1.5394,1.4442,1.3470,1.2494,1.1644,1.0903,1.0453,    &
      1.0200,1.0000,1.0000,1.0000,1.0000,1.0000,1.0000,1.0000,    &
      1.0000,1.0000,1.0000,1.0000,1.0000,1.0000 /
data (omega(i, 6, 5),i=1,nlam) /  &
      1.9938,1.9936,1.9926,1.9917,1.9895,1.9878,1.9850,1.9802,    &
      1.9732,1.9630,1.9480,1.9263,1.8947,1.8525,1.8001,1.7414,    &
      1.6719,1.5913,1.5018,1.4049,1.3050,1.2131,1.1271,1.0694,    &
      1.0329,1.0000,1.0000,1.0000,1.0000,1.0000,1.0000,1.0000,    &
      1.0000,1.0000,1.0000,1.0000,1.0000,1.0000 /
data (omega(i, 7, 5),i=1,nlam) /  &
      1.9950,1.9938,1.9926,1.9919,1.9914,1.9891,1.9871,1.9837,    &
      1.9778,1.9691,1.9562,1.9374,1.9100,1.8725,1.8257,1.7717,    &
      1.7080,1.6333,1.5487,1.4545,1.3574,1.2592,1.1727,1.0963,    &
      1.0491,1.0220,1.0000,1.0000,1.0000,1.0000,1.0000,1.0000,    &
      1.0000,1.0000,1.0000,1.0000,1.0000,1.0000 /
data (omega(i, 8, 5),i=1,nlam) /  &
      1.9959,1.9950,1.9938,1.9931,1.9928,1.9908,1.9891,1.9863,    &
      1.9813,1.9737,1.9625,1.9459,1.9220,1.8879,1.8455,1.7961,    &
      1.7384,1.6694,1.5882,1.4972,1.4016,1.3017,1.2103,1.1248,    &
      1.0678,1.0321,1.0000,1.0000,1.0000,1.0000,1.0000,1.0000,    &
      1.0000,1.0000,1.0000,1.0000,1.0000,1.0000 /
data (omega(i, 9, 5),i=1,nlam) /  &
      1.9966,1.9958,1.9947,1.9943,1.9925,1.9920,1.9897,1.9877,    &
      1.9839,1.9774,1.9675,1.9528,1.9313,1.9008,1.8616,1.8171,    &
      1.7628,1.6982,1.6217,1.5358,1.4409,1.3436,1.2462,1.1619,    &
      1.0883,1.0442,1.0000,1.0000,1.0000,1.0000,1.0000,1.0000,    &
      1.0000,1.0000,1.0000,1.0000,1.0000,1.0000 /
data (omega(i,10, 5),i=1,nlam) /  &
      1.9957,1.9964,1.9960,1.9957,1.9937,1.9938,1.9919,1.9905,    &
      1.9871,1.9830,1.9747,1.9628,1.9452,1.9201,1.8868,1.8481,    &
      1.8017,1.7445,1.6770,1.5971,1.5071,1.4119,1.3120,1.2195,    &
      1.1321,1.0729,1.0350,1.0000,1.0000,1.0000,1.0000,1.0000,    &
      1.0000,1.0000,1.0000,1.0000,1.0000,1.0000 /
data (omega(i,11, 5),i=1,nlam) /  &
      1.9966,1.9958,1.9967,1.9950,1.9953,1.9948,1.9934,1.9913,    &
      1.9897,1.9861,1.9797,1.9698,1.9550,1.9337,1.9053,1.8715,    &
      1.8300,1.7790,1.7181,1.6462,1.5620,1.4690,1.3722,1.2733,    &
      1.1850,1.1055,1.0549,1.0250,1.0000,1.0000,1.0000,1.0000,    &
      1.0000,1.0000,1.0000,1.0000,1.0000,1.0000 /
data (omega(i,12, 5),i=1,nlam) /  &
      1.9972,1.9968,1.9972,1.9958,1.9961,1.9940,1.9945,1.9927,    &
      1.9903,1.9880,1.9832,1.9749,1.9622,1.9438,1.9193,1.8894,    &
      1.8525,1.8067,1.7508,1.6845,1.6073,1.5187,1.4223,1.3254,    &
      1.2289,1.1397,1.0782,1.0380,1.0000,1.0000,1.0000,1.0000,    &
      1.0000,1.0000,1.0000,1.0000,1.0000,1.0000 /
data (omega(i,13, 5),i=1,nlam) /  &
      1.9976,1.9974,1.9962,1.9963,1.9966,1.9946,1.9951,1.9938,    &
      1.9916,1.9898,1.9859,1.9788,1.9678,1.9516,1.9299,1.9031,    &
      1.8698,1.8282,1.7769,1.7159,1.6437,1.5594,1.4658,1.3690,    &
      1.2702,1.1823,1.1034,1.0536,1.0243,1.0000,1.0000,1.0000,    &
      1.0000,1.0000,1.0000,1.0000,1.0000,1.0000 /
data (omega(i,14, 5),i=1,nlam) /  &
      1.9979,1.9980,1.9967,1.9967,1.9971,1.9958,1.9957,1.9945,    &
      1.9927,1.9901,1.9873,1.9822,1.9721,1.9578,1.9384,1.9146,    &
      1.8844,1.8464,1.7994,1.7417,1.6743,1.5954,1.5039,1.4086,    &
      1.3086,1.2165,1.1297,1.0712,1.0340,1.0000,1.0000,1.0000,    &
      1.0000,1.0000,1.0000,1.0000,1.0000,1.0000 /
data (omega(i, 1, 6),i=1,nlam) /  &
      1.9763,1.9729,1.9694,1.9652,1.9601,1.9536,1.9452,1.9340,    &
      1.9192,1.8990,1.8717,1.8349,1.7853,1.7209,1.6375,1.5374,    &
      1.4286,1.3267,1.2297,1.1400,1.0783,1.0381,1.0000,1.0000,    &
      1.0000,1.0000,1.0000,1.0000,1.0000,1.0000,1.0000,1.0000,    &
      1.0000,1.0000,1.0000,1.0000,1.0000,1.0000 /
data (omega(i, 2, 6),i=1,nlam) /  &
      1.9867,1.9847,1.9831,1.9806,1.9774,1.9734,1.9681,1.9610,    &
      1.9510,1.9378,1.9196,1.8939,1.8587,1.8107,1.7469,1.6657,    &
      1.5741,1.4774,1.3798,1.2805,1.1912,1.1100,1.0579,1.0265,    &
      1.0000,1.0000,1.0000,1.0000,1.0000,1.0000,1.0000,1.0000,    &
      1.0000,1.0000,1.0000,1.0000,1.0000,1.0000 /
data (omega(i, 3, 6),i=1,nlam) /  &
      1.9900,1.9893,1.9884,1.9869,1.9847,1.9825,1.9788,1.9737,    &
      1.9668,1.9572,1.9437,1.9243,1.8976,1.8596,1.8086,1.7421,    &
      1.6645,1.5794,1.4881,1.3900,1.2904,1.2000,1.1168,1.0624,    &
      1.0290,1.0000,1.0000,1.0000,1.0000,1.0000,1.0000,1.0000,    &
      1.0000,1.0000,1.0000,1.0000,1.0000,1.0000 /
data (omega(i, 4, 6),i=1,nlam) /  &
      1.9930,1.9912,1.9909,1.9897,1.9889,1.9868,1.9846,1.9809,    &
      1.9756,1.9683,1.9577,1.9428,1.9209,1.8902,1.8480,1.7921,    &
      1.7243,1.6506,1.5663,1.4742,1.3761,1.2769,1.1881,1.1078,    &
      1.0564,1.0258,1.0000,1.0000,1.0000,1.0000,1.0000,1.0000,    &
      1.0000,1.0000,1.0000,1.0000,1.0000,1.0000 /
data (omega(i, 5, 6),i=1,nlam) /  &
      1.9936,1.9931,1.9930,1.9919,1.9907,1.9899,1.9877,1.9853,    &
      1.9812,1.9753,1.9669,1.9546,1.9369,1.9111,1.8752,1.8265,    &
      1.7680,1.7015,1.6250,1.5392,1.4441,1.3469,1.2492,1.1645,    &
      1.0901,1.0453,1.0000,1.0000,1.0000,1.0000,1.0000,1.0000,    &
      1.0000,1.0000,1.0000,1.0000,1.0000,1.0000 /
data (omega(i, 6, 6),i=1,nlam) /  &
      1.9950,1.9947,1.9931,1.9933,1.9924,1.9911,1.9893,1.9878,    &
      1.9850,1.9802,1.9732,1.9630,1.9480,1.9263,1.8947,1.8525,    &
      1.8001,1.7413,1.6718,1.5912,1.5016,1.4047,1.3048,1.2130,    &
      1.1269,1.0693,1.0329,1.0000,1.0000,1.0000,1.0000,1.0000,    &
      1.0000,1.0000,1.0000,1.0000,1.0000,1.0000 /
data (omega(i, 7, 6),i=1,nlam) /  &
      1.9956,1.9943,1.9942,1.9943,1.9935,1.9924,1.9910,1.9899,    &
      1.9871,1.9837,1.9778,1.9691,1.9562,1.9374,1.9100,1.8725,    &
      1.8256,1.7716,1.7079,1.6331,1.5485,1.4543,1.3572,1.2590,    &
      1.1726,1.0962,1.0490,1.0220,1.0000,1.0000,1.0000,1.0000,    &
      1.0000,1.0000,1.0000,1.0000,1.0000,1.0000 /
data (omega(i, 8, 6),i=1,nlam) /  &
      1.9964,1.9954,1.9950,1.9952,1.9943,1.9934,1.9922,1.9905,    &
      1.9891,1.9863,1.9812,1.9738,1.9625,1.9459,1.9220,1.8879,    &
      1.8454,1.7960,1.7383,1.6693,1.5881,1.4970,1.4014,1.3015,    &
      1.2101,1.1247,1.0677,1.0320,1.0000,1.0000,1.0000,1.0000,    &
      1.0000,1.0000,1.0000,1.0000,1.0000,1.0000 /
data (omega(i, 9, 6),i=1,nlam) /  &
      1.9970,1.9967,1.9957,1.9958,1.9950,1.9942,1.9931,1.9916,    &
      1.9906,1.9877,1.9839,1.9774,1.9675,1.9527,1.9313,1.9008,    &
      1.8615,1.8170,1.7628,1.6980,1.6216,1.5356,1.4407,1.3434,    &
      1.2460,1.1536,1.0882,1.0441,1.0000,1.0000,1.0000,1.0000,    &
      1.0000,1.0000,1.0000,1.0000,1.0000,1.0000 /
data (omega(i,10, 6),i=1,nlam) /  &
      1.9963,1.9962,1.9966,1.9951,1.9960,1.9952,1.9944,1.9933,    &
      1.9917,1.9905,1.9871,1.9831,1.9747,1.9628,1.9452,1.9201,    &
      1.8868,1.8480,1.8017,1.7444,1.6769,1.5969,1.5070,1.4118,    &
      1.3118,1.2194,1.1320,1.0728,1.0349,1.0000,1.0000,1.0000,    &
      1.0000,1.0000,1.0000,1.0000,1.0000,1.0000 /
data (omega(i,11, 6),i=1,nlam) /  &
      1.9965,1.9972,1.9971,1.9958,1.9966,1.9958,1.9952,1.9944,    &
      1.9931,1.9912,1.9896,1.9861,1.9797,1.9698,1.9550,1.9337,    &
      1.9053,1.8715,1.8300,1.7790,1.7180,1.6461,1.5618,1.4707,    &
      1.3720,1.2732,1.1849,1.1053,1.0549,1.0250,1.0000,1.0000,    &
      1.0000,1.0000,1.0000,1.0000,1.0000,1.0000 /
data (omega(i,12, 6),i=1,nlam) /  &
      1.9972,1.9980,1.9976,1.9958,1.9971,1.9963,1.9956,1.9950,    &
      1.9941,1.9925,1.9903,1.9879,1.9832,1.9749,1.9622,1.9438,    &
      1.9193,1.8894,1.8525,1.8067,1.7508,1.6844,1.6072,1.5185,    &
      1.4221,1.3220,1.2288,1.1395,1.0781,1.0380,1.0000,1.0000,    &
      1.0000,1.0000,1.0000,1.0000,1.0000,1.0000 /
data (omega(i,13, 6),i=1,nlam) /  &
      1.9977,1.9969,1.9980,1.9960,1.9974,1.9967,1.9960,1.9954,    &
      1.9947,1.9935,1.9916,1.9898,1.9859,1.9788,1.9678,1.9516,    &
      1.9299,1.9036,1.8698,1.8282,1.7768,1.7158,1.6436,1.5593,    &
      1.4677,1.3688,1.2701,1.1822,1.1033,1.0536,1.0000,1.0000,    &
      1.0000,1.0000,1.0000,1.0000,1.0000,1.0000 /
data (omega(i,14, 6),i=1,nlam) /  &
      1.9979,1.9973,1.9958,1.9970,1.9977,1.9970,1.9963,1.9958,    &
      1.9951,1.9941,1.9925,1.9912,1.9873,1.9822,1.9721,1.9578,    &
      1.9384,1.9146,1.8844,1.8464,1.7994,1.7416,1.6742,1.5953,    &
      1.5057,1.4084,1.3084,1.2163,1.1296,1.0711,1.0339,1.0000,    &
      1.0000,1.0000,1.0000,1.0000,1.0000,1.0000 /
data (omega(i, 1, 7),i=1,nlam) /  &
      1.9812,1.9780,1.9751,1.9723,1.9690,1.9650,1.9600,1.9535,    &
      1.9451,1.9340,1.9192,1.8990,1.8718,1.8349,1.7854,1.7209,    &
      1.6376,1.5375,1.4287,1.3268,1.2297,1.1401,1.0783,1.0381,    &
      1.0000,1.0000,1.0000,1.0000,1.0000,1.0000,1.0000,1.0000,    &
      1.0000,1.0000,1.0000,1.0000,1.0000,1.0000 /
data (omega(i, 2, 7),i=1,nlam) /  &
      1.9884,1.9875,1.9860,1.9843,1.9829,1.9805,1.9774,1.9734,    &
      1.9681,1.9610,1.9510,1.9378,1.9196,1.8944,1.8587,1.8107,    &
      1.7469,1.6658,1.5742,1.4774,1.3799,1.2806,1.1912,1.1100,    &
      1.0579,1.0266,1.0000,1.0000,1.0000,1.0000,1.0000,1.0000,    &
      1.0000,1.0000,1.0000,1.0000,1.0000,1.0000 /
data (omega(i, 3, 7),i=1,nlam) /  &
      1.9920,1.9911,1.9901,1.9890,1.9884,1.9868,1.9852,1.9824,    &
      1.9787,1.9737,1.9668,1.9572,1.9438,1.9243,1.8976,1.8596,    &
      1.8086,1.7422,1.6646,1.5794,1.4882,1.3900,1.2904,1.2000,    &
      1.1168,1.0624,1.0290,1.0000,1.0000,1.0000,1.0000,1.0000,    &
      1.0000,1.0000,1.0000,1.0000,1.0000,1.0000 /
data (omega(i, 4, 7),i=1,nlam) /  &
      1.9941,1.9934,1.9927,1.9918,1.9908,1.9896,1.9888,1.9868,    &
      1.9846,1.9809,1.9756,1.9683,1.9578,1.9428,1.9209,1.8902,    &
      1.8480,1.7921,1.7244,1.6506,1.5663,1.4742,1.3761,1.2770,    &
      1.1882,1.1078,1.0565,1.0258,1.0000,1.0000,1.0000,1.0000,    &
      1.0000,1.0000,1.0000,1.0000,1.0000,1.0000 /
data (omega(i, 5, 7),i=1,nlam) /  &
      1.9942,1.9947,1.9943,1.9937,1.9929,1.9919,1.9907,1.9899,    &
      1.9877,1.9853,1.9812,1.9753,1.9669,1.9546,1.9370,1.9111,    &
      1.8753,1.8265,1.7681,1.7015,1.6250,1.5393,1.4441,1.3469,    &
      1.2493,1.1646,1.0902,1.0453,1.0000,1.0000,1.0000,1.0000,    &
      1.0000,1.0000,1.0000,1.0000,1.0000,1.0000 /
data (omega(i, 6, 7),i=1,nlam) /  &
      1.9952,1.9956,1.9952,1.9947,1.9941,1.9933,1.9923,1.9910,    &
      1.9893,1.9878,1.9850,1.9802,1.9732,1.9630,1.9480,1.9264,    &
      1.8952,1.8525,1.8001,1.7413,1.6719,1.5912,1.5016,1.4048,    &
      1.3048,1.2130,1.1270,1.0693,1.0329,1.0000,1.0000,1.0000,    &
      1.0000,1.0000,1.0000,1.0000,1.0000,1.0000 /
data (omega(i, 7, 7),i=1,nlam) /  &
      1.9960,1.9963,1.9957,1.9954,1.9949,1.9942,1.9934,1.9924,    &
      1.9910,1.9899,1.9871,1.9837,1.9778,1.9691,1.9562,1.9374,    &
      1.9100,1.8726,1.8256,1.7717,1.7080,1.6332,1.5486,1.4563,    &
      1.3572,1.2590,1.1727,1.0962,1.0491,1.0000,1.0000,1.0000,    &
      1.0000,1.0000,1.0000,1.0000,1.0000,1.0000 /
data (omega(i, 8, 7),i=1,nlam) /  &
      1.9967,1.9969,1.9963,1.9960,1.9956,1.9950,1.9943,1.9934,    &
      1.9921,1.9905,1.9891,1.9863,1.9812,1.9738,1.9625,1.9460,    &
      1.9220,1.8879,1.8455,1.7961,1.7383,1.6694,1.5881,1.4970,    &
      1.4015,1.3016,1.2101,1.1247,1.0678,1.0320,1.0000,1.0000,    &
      1.0000,1.0000,1.0000,1.0000,1.0000,1.0000 /
data (omega(i, 9, 7),i=1,nlam) /  &
      1.9973,1.9973,1.9968,1.9964,1.9961,1.9956,1.9949,1.9941,    &
      1.9931,1.9916,1.9906,1.9877,1.9839,1.9774,1.9675,1.9528,    &
      1.9313,1.9008,1.8616,1.8171,1.7628,1.6980,1.6217,1.5357,    &
      1.4408,1.3435,1.2461,1.1537,1.0882,1.0441,1.0000,1.0000,    &
      1.0000,1.0000,1.0000,1.0000,1.0000,1.0000 /
data (omega(i,10, 7),i=1,nlam) /  &
      1.9980,1.9965,1.9972,1.9966,1.9965,1.9962,1.9957,1.9952,    &
      1.9944,1.9933,1.9917,1.9905,1.9871,1.9831,1.9747,1.9628,    &
      1.9452,1.9201,1.8868,1.8481,1.8017,1.7445,1.6770,1.5970,    &
      1.5088,1.4118,1.3118,1.2194,1.1320,1.0728,1.0349,1.0000,    &
      1.0000,1.0000,1.0000,1.0000,1.0000,1.0000 /
data (omega(i,11, 7),i=1,nlam) /  &
      1.9974,1.9973,1.9976,1.9970,1.9968,1.9965,1.9961,1.9957,    &
      1.9952,1.9943,1.9931,1.9913,1.9896,1.9861,1.9806,1.9698,    &
      1.9550,1.9337,1.9053,1.8715,1.8300,1.7790,1.7181,1.6461,    &
      1.5619,1.4689,1.3721,1.2732,1.1849,1.1054,1.0549,1.0251,    &
      1.0000,1.0000,1.0000,1.0000,1.0000,1.0000 /
data (omega(i,12, 7),i=1,nlam) /  &
      1.9963,1.9979,1.9980,1.9973,1.9971,1.9968,1.9965,1.9961,    &
      1.9956,1.9950,1.9940,1.9925,1.9903,1.9880,1.9833,1.9749,    &
      1.9622,1.9438,1.9193,1.8894,1.8525,1.8067,1.7508,1.6845,    &
      1.6072,1.5185,1.4221,1.3221,1.2288,1.1396,1.0781,1.0380,    &
      1.0000,1.0000,1.0000,1.0000,1.0000,1.0000 /
data (omega(i,13, 7),i=1,nlam) /  &
      1.9970,1.9983,1.9954,1.9975,1.9970,1.9969,1.9967,1.9964,    &
      1.9959,1.9954,1.9946,1.9935,1.9916,1.9898,1.9859,1.9793,    &
      1.9678,1.9516,1.9300,1.9036,1.8698,1.8282,1.7768,1.7158,    &
      1.6436,1.5593,1.4657,1.3689,1.2701,1.1823,1.1033,1.0536,    &
      1.0000,1.0000,1.0000,1.0000,1.0000,1.0000 /
data (omega(i,14, 7),i=1,nlam) /  &
      1.9976,1.9978,1.9957,1.9976,1.9969,1.9968,1.9967,1.9965,    &
      1.9962,1.9957,1.9951,1.9941,1.9925,1.9912,1.9873,1.9822,    &
      1.9721,1.9578,1.9384,1.9146,1.8844,1.8464,1.7994,1.7416,    &
      1.6743,1.5953,1.5057,1.4084,1.3084,1.2164,1.1296,1.0711,    &
      1.0339,1.0000,1.0000,1.0000,1.0000,1.0000 /
data (omega(i, 1, 8),i=1,nlam) /  &
      1.9875,1.9840,1.9804,1.9775,1.9749,1.9721,1.9689,1.9650,    &
      1.9600,1.9535,1.9451,1.9340,1.9192,1.8994,1.8718,1.8349,    &
      1.7854,1.7209,1.6376,1.5375,1.4286,1.3268,1.2297,1.1400,    &
      1.0783,1.0381,1.0000,1.0000,1.0000,1.0000,1.0000,1.0000,    &
      1.0000,1.0000,1.0000,1.0000,1.0000,1.0000 /
data (omega(i, 2, 8),i=1,nlam) /  &
      1.9912,1.9893,1.9878,1.9872,1.9858,1.9847,1.9828,1.9804,    &
      1.9774,1.9734,1.9681,1.9610,1.9517,1.9378,1.9196,1.8944,    &
      1.8587,1.8107,1.7469,1.6658,1.5742,1.4774,1.3798,1.2806,    &
      1.1912,1.1100,1.0579,1.0266,1.0000,1.0000,1.0000,1.0000,    &
      1.0000,1.0000,1.0000,1.0000,1.0000,1.0000 /
data (omega(i, 3, 8),i=1,nlam) /  &
      1.9923,1.9925,1.9917,1.9909,1.9900,1.9889,1.9884,1.9868,    &
      1.9852,1.9824,1.9787,1.9737,1.9668,1.9572,1.9439,1.9243,    &
      1.8976,1.8596,1.8086,1.7422,1.6646,1.5794,1.4882,1.3900,    &
      1.2904,1.2000,1.1168,1.0624,1.0290,1.0000,1.0000,1.0000,    &
      1.0000,1.0000,1.0000,1.0000,1.0000,1.0000 /
data (omega(i, 4, 8),i=1,nlam) /  &
      1.9938,1.9931,1.9925,1.9933,1.9926,1.9918,1.9908,1.9896,    &
      1.9888,1.9868,1.9846,1.9809,1.9756,1.9683,1.9578,1.9429,    &
      1.9209,1.8902,1.8480,1.7921,1.7244,1.6506,1.5663,1.4742,    &
      1.3761,1.2770,1.1881,1.1078,1.0565,1.0258,1.0000,1.0000,    &
      1.0000,1.0000,1.0000,1.0000,1.0000,1.0000 /
data (omega(i, 5, 8),i=1,nlam) /  &
      1.9948,1.9943,1.9938,1.9947,1.9942,1.9936,1.9928,1.9919,    &
      1.9906,1.9899,1.9877,1.9853,1.9812,1.9753,1.9669,1.9546,    &
      1.9370,1.9111,1.8753,1.8265,1.7681,1.7015,1.6250,1.5393,    &
      1.4441,1.3469,1.2493,1.1645,1.0902,1.0453,1.0000,1.0000,    &
      1.0000,1.0000,1.0000,1.0000,1.0000,1.0000 /
data (omega(i, 6, 8),i=1,nlam) /  &
      1.9966,1.9963,1.9959,1.9955,1.9951,1.9946,1.9940,1.9933,    &
      1.9923,1.9910,1.9893,1.9878,1.9850,1.9802,1.9732,1.9630,    &
      1.9480,1.9263,1.8952,1.8525,1.8001,1.7413,1.6719,1.5912,    &
      1.5016,1.4048,1.3048,1.2130,1.1270,1.0693,1.0329,1.0000,    &
      1.0000,1.0000,1.0000,1.0000,1.0000,1.0000 /
data (omega(i, 7, 8),i=1,nlam) /  &
      1.9970,1.9968,1.9965,1.9961,1.9958,1.9954,1.9949,1.9942,    &
      1.9934,1.9924,1.9910,1.9899,1.9871,1.9837,1.9778,1.9691,    &
      1.9562,1.9374,1.9100,1.8726,1.8256,1.7716,1.7080,1.6332,    &
      1.5486,1.4544,1.3572,1.2590,1.1727,1.0962,1.0491,1.0000,    &
      1.0000,1.0000,1.0000,1.0000,1.0000,1.0000 /
data (omega(i, 8, 8),i=1,nlam) /  &
      1.9974,1.9972,1.9969,1.9966,1.9963,1.9960,1.9955,1.9950,    &
      1.9943,1.9933,1.9921,1.9905,1.9891,1.9863,1.9831,1.9738,    &
      1.9625,1.9460,1.9220,1.8879,1.8455,1.7961,1.7383,1.6693,    &
      1.5881,1.4991,1.4015,1.3016,1.2101,1.1247,1.0677,1.0320,    &
      1.0000,1.0000,1.0000,1.0000,1.0000,1.0000 /
data (omega(i, 9, 8),i=1,nlam) /  &
      1.9976,1.9975,1.9972,1.9970,1.9967,1.9964,1.9960,1.9955,    &
      1.9949,1.9941,1.9931,1.9916,1.9906,1.9877,1.9839,1.9774,    &
      1.9675,1.9528,1.9313,1.9008,1.8622,1.8171,1.7628,1.6980,    &
      1.6217,1.5357,1.4408,1.3434,1.2461,1.1537,1.0882,1.0441,    &
      1.0000,1.0000,1.0000,1.0000,1.0000,1.0000 /
data (omega(i,10, 8),i=1,nlam) /  &
      1.9979,1.9974,1.9974,1.9972,1.9970,1.9968,1.9965,1.9961,    &
      1.9957,1.9952,1.9944,1.9933,1.9917,1.9905,1.9871,1.9833,    &
      1.9747,1.9628,1.9452,1.9201,1.8868,1.8481,1.8017,1.7445,    &
      1.6770,1.5970,1.5088,1.4118,1.3118,1.2194,1.1320,1.0728,    &
      1.0349,1.0000,1.0000,1.0000,1.0000,1.0000 /
data (omega(i,11, 8),i=1,nlam) /  &
      1.9983,1.9977,1.9976,1.9974,1.9972,1.9970,1.9968,1.9965,    &
      1.9961,1.9957,1.9951,1.9943,1.9931,1.9912,1.9896,1.9861,    &
      1.9809,1.9698,1.9550,1.9337,1.9053,1.8721,1.8300,1.7790,    &
      1.7181,1.6461,1.5618,1.4689,1.3721,1.2732,1.1850,1.1053,    &
      1.0549,1.0000,1.0000,1.0000,1.0000,1.0000 /
data (omega(i,12, 8),i=1,nlam) /  &
      1.9958,1.9979,1.9973,1.9976,1.9974,1.9972,1.9970,1.9967,    &
      1.9964,1.9961,1.9956,1.9949,1.9940,1.9925,1.9903,1.9879,    &
      1.9832,1.9749,1.9622,1.9438,1.9193,1.8894,1.8525,1.8067,    &
      1.7508,1.6845,1.6072,1.5185,1.4221,1.3220,1.2288,1.1395,    &
      1.0781,1.0380,1.0000,1.0000,1.0000,1.0000 /
data (omega(i,13, 8),i=1,nlam) /  &
      1.9965,1.9980,1.9972,1.9974,1.9973,1.9971,1.9970,1.9968,    &
      1.9966,1.9963,1.9959,1.9954,1.9946,1.9934,1.9916,1.9898,    &
      1.9859,1.9794,1.9678,1.9516,1.9300,1.9036,1.8698,1.8282,    &
      1.7768,1.7172,1.6436,1.5593,1.4657,1.3689,1.2701,1.1825,    &
      1.1033,1.0536,1.0000,1.0000,1.0000,1.0000 /
data (omega(i,14, 8),i=1,nlam) /  &
      1.9974,1.9981,1.9972,1.9971,1.9970,1.9968,1.9967,1.9967,    &
      1.9966,1.9965,1.9962,1.9957,1.9950,1.9941,1.9925,1.9912,    &
      1.9873,1.9822,1.9721,1.9578,1.9384,1.9146,1.8844,1.8464,    &
      1.7994,1.7430,1.6743,1.5953,1.5059,1.4084,1.3084,1.2163,    &
      1.1296,1.0711,1.0339,1.0000,1.0000,1.0000 /
data (omega(i, 1, 9),i=1,nlam) /  &
      1.9931,1.9905,1.9872,1.9837,1.9802,1.9773,1.9748,1.9721,    &
      1.9689,1.9650,1.9599,1.9535,1.9451,1.9340,1.9192,1.8994,    &
      1.8723,1.8348,1.7851,1.7206,1.6372,1.5371,1.4283,1.3264,    &
      1.2294,1.1398,1.0781,1.0380,1.0000,1.0000,1.0000,1.0000,    &
      1.0000,1.0000,1.0000,1.0000,1.0000,1.0000 /
data (omega(i, 2, 9),i=1,nlam) /  &
      1.9943,1.9920,1.9908,1.9891,1.9877,1.9871,1.9858,1.9847,    &
      1.9828,1.9804,1.9774,1.9734,1.9681,1.9610,1.9518,1.9377,    &
      1.9195,1.8943,1.8586,1.8105,1.7466,1.6655,1.5738,1.4771,    &
      1.3795,1.2802,1.1909,1.1098,1.0578,1.0265,1.0000,1.0000,    &
      1.0000,1.0000,1.0000,1.0000,1.0000,1.0000 /
data (omega(i, 3, 9),i=1,nlam) /  &
      1.9948,1.9933,1.9921,1.9923,1.9916,1.9908,1.9899,1.9889,    &
      1.9884,1.9867,1.9852,1.9824,1.9787,1.9737,1.9668,1.9571,    &
      1.9440,1.9242,1.8975,1.8595,1.8084,1.7419,1.6643,1.5791,    &
      1.4878,1.3897,1.2901,1.1997,1.1166,1.0622,1.0289,1.0000,    &
      1.0000,1.0000,1.0000,1.0000,1.0000,1.0000 /
data (omega(i, 4, 9),i=1,nlam) /  &
      1.9953,1.9943,1.9936,1.9930,1.9938,1.9932,1.9926,1.9918,    &
      1.9908,1.9896,1.9888,1.9868,1.9846,1.9808,1.9756,1.9682,    &
      1.9577,1.9428,1.9208,1.8900,1.8478,1.7919,1.7241,1.6503,    &
      1.5660,1.4738,1.3758,1.2766,1.1878,1.1075,1.0563,1.0257,    &
      1.0000,1.0000,1.0000,1.0000,1.0000,1.0000 /
data (omega(i, 5, 9),i=1,nlam) /  &
      1.9959,1.9951,1.9946,1.9941,1.9937,1.9947,1.9942,1.9936,    &
      1.9928,1.9919,1.9906,1.9899,1.9877,1.9853,1.9812,1.9753,    &
      1.9669,1.9546,1.9369,1.9115,1.8751,1.8263,1.7678,1.7013,    &
      1.6247,1.5389,1.4458,1.3436,1.2489,1.1560,1.0899,1.0452,    &
      1.0000,1.0000,1.0000,1.0000,1.0000,1.0000 /
data (omega(i, 6, 9),i=1,nlam) /  &
      1.9959,1.9968,1.9965,1.9962,1.9958,1.9955,1.9951,1.9946,    &
      1.9940,1.9933,1.9923,1.9910,1.9893,1.9878,1.9850,1.9802,    &
      1.9732,1.9629,1.9480,1.9263,1.8951,1.8524,1.7999,1.7411,    &
      1.6716,1.5909,1.5013,1.4044,1.3045,1.2127,1.1267,1.0691,    &
      1.0328,1.0000,1.0000,1.0000,1.0000,1.0000 /
data (omega(i, 7, 9),i=1,nlam) /  &
      1.9975,1.9972,1.9969,1.9967,1.9964,1.9961,1.9958,1.9954,    &
      1.9948,1.9942,1.9934,1.9924,1.9909,1.9899,1.9871,1.9837,    &
      1.9778,1.9691,1.9562,1.9373,1.9099,1.8724,1.8255,1.7714,    &
      1.7077,1.6329,1.5482,1.4560,1.3568,1.2587,1.1642,1.0960,    &
      1.0489,1.0000,1.0000,1.0000,1.0000,1.0000 /
data (omega(i, 8, 9),i=1,nlam) /  &
      1.9978,1.9975,1.9973,1.9971,1.9969,1.9966,1.9963,1.9960,    &
      1.9955,1.9950,1.9943,1.9933,1.9921,1.9905,1.9891,1.9863,    &
      1.9837,1.9737,1.9625,1.9459,1.9219,1.8878,1.8453,1.7959,    &
      1.7381,1.6691,1.5878,1.4987,1.4011,1.3012,1.2098,1.1244,    &
      1.0676,1.0319,1.0000,1.0000,1.0000,1.0000 /
data (omega(i, 9, 9),i=1,nlam) /  &
      1.9980,1.9977,1.9976,1.9974,1.9972,1.9969,1.9967,1.9964,    &
      1.9960,1.9955,1.9949,1.9941,1.9931,1.9916,1.9906,1.9877,    &
      1.9839,1.9773,1.9675,1.9527,1.9312,1.9007,1.8614,1.8169,    &
      1.7626,1.6979,1.6214,1.5353,1.4404,1.3432,1.2457,1.1534,    &
      1.0880,1.0440,1.0000,1.0000,1.0000,1.0000 /
data (omega(i,10, 9),i=1,nlam) /  &
      1.9980,1.9979,1.9977,1.9975,1.9974,1.9972,1.9970,1.9967,    &
      1.9965,1.9961,1.9957,1.9952,1.9944,1.9933,1.9917,1.9905,    &
      1.9871,1.9833,1.9747,1.9628,1.9451,1.9200,1.8867,1.8479,    &
      1.8015,1.7442,1.6767,1.5967,1.5086,1.4114,1.3115,1.2191,    &
      1.1318,1.0726,1.0348,1.0000,1.0000,1.0000 /
data (omega(i,11, 9),i=1,nlam) /  &
      1.9981,1.9980,1.9979,1.9977,1.9975,1.9974,1.9972,1.9970,    &
      1.9968,1.9965,1.9961,1.9957,1.9951,1.9943,1.9931,1.9912,    &
      1.9896,1.9861,1.9807,1.9698,1.9549,1.9336,1.9052,1.8714,    &
      1.8299,1.7798,1.7178,1.6459,1.5615,1.4705,1.3717,1.2728,    &
      1.1846,1.1051,1.0547,1.0251,1.0000,1.0000 /
data (omega(i,12, 9),i=1,nlam) /  &
      1.9982,1.9981,1.9979,1.9978,1.9976,1.9975,1.9973,1.9972,    &
      1.9970,1.9967,1.9964,1.9961,1.9956,1.9949,1.9940,1.9925,    &
      1.9903,1.9879,1.9832,1.9749,1.9622,1.9437,1.9192,1.8898,    &
      1.8523,1.8065,1.7506,1.6842,1.6069,1.5182,1.4218,1.3217,    &
      1.2285,1.1393,1.0779,1.0379,1.0000,1.0000 /
data (omega(i,13, 9),i=1,nlam) /  &
      1.9980,1.9979,1.9978,1.9976,1.9975,1.9973,1.9972,1.9971,    &
      1.9970,1.9968,1.9966,1.9963,1.9959,1.9954,1.9946,1.9934,    &
      1.9916,1.9898,1.9859,1.9794,1.9677,1.9515,1.9299,1.9035,    &
      1.8703,1.8288,1.7778,1.7156,1.6433,1.5590,1.4654,1.3685,    &
      1.2697,1.1736,1.1031,1.0534,1.0000,1.0000 /
data (omega(i,14, 9),i=1,nlam) /  &
      1.9974,1.9976,1.9974,1.9972,1.9971,1.9970,1.9968,1.9968,    &
      1.9967,1.9967,1.9966,1.9964,1.9962,1.9957,1.9950,1.9941,    &
      1.9925,1.9912,1.9873,1.9822,1.9721,1.9577,1.9383,1.9146,    &
      1.8848,1.8471,1.8003,1.7427,1.6740,1.5950,1.5053,1.4081,    &
      1.3081,1.2160,1.1293,1.0709,1.0339,1.0000 /

loglam = log(4 * radsq / ( phibar * dt**2 ))/log(2.)

if ( il < resol(1) .or. il >= resol(nres) ) then
   write(6,*) "Error: resolution out of range in optmx"
   stop
end if
schmidt_min = 0.5**(nsch-1)
if ( schmidt > 1.0 .or. schmidt <= schmidt_min ) then
   write(6,*) "Warning: schmidt out of range in optmx", schmidt
end if
if ( loglam < 1.0 ) then
   write(6,*) "Error: lambda out of range in optmx"
   stop
end if

! Table size chosen so that this is true even for the highest resolution
! and stretching.
if ( loglam > nlam-1 ) then
   accel = 1.0
   return
end if

! Find i, such that lam(i) <= lambda < lam(i+1)
!      j, such that res(j) <= il < res(j+1)
!      k, such that sch(k) >= lambda > sch(k+1)
   
i = 1 + floor(loglam)
if ( il <= 100 ) then
   j = il/10 - 1
else
   j = (il-100)/20 + 9
end if
logstr = -log(max(schmidt,schmidt_min))/log(2.)
k = 1 + floor(logstr)

! Interpolation weights
wl = 1. - ( loglam - floor(loglam) )
ws = 1. - ( logstr - floor(logstr) )
wr = 1. - ( il - resol(j) )/ real(resol(j+1)-resol(j))

accel = wl*wr*ws*omega(i,j,k) +                    &
        (1-wl)*wr*ws*omega(i+1,j,k) +              &
        wl*(1-wr)*ws*omega(i,j+1,k) +              &
        (1-wl)*(1-wr)*ws*omega(i+1,j+1,k) +        &
        wl*wr*(1-ws)*omega(i,j,k+1) +              &
        (1-wl)*wr*(1-ws)*omega(i+1,j,k+1) +        & 
        wl*(1-wr)*(1-ws)*omega(i,j+1,k+1) +        &
        (1-wl)*(1-wr)*(1-ws)*omega(i+1,j+1,k+1)

end subroutine optmx
                 

! Initialise multi-grid arrays
subroutine mgsor_init

use cc_mpi
use newmpar_m
use parm_m
use parmdyn_m

implicit none

integer g, gp, np, iq, iqq, nn, ii, jj
integer mipan, mjpan, hipan, hjpan, mil_g, iia, jja
integer i, j, n, mg_npan, mxpr, mypr, sii, eii, sjj, ejj
integer cid, ix, jx, colour, rank, ncol, nrow
integer npanx, na, nx, ny, drow, dcol, mmx, mmy
logical lglob

if ( .not.sorfirst ) return

! Begin full initialisation
mg_maxsize = ifull + iextra ! first guess of maximum grid size
lglob = (nproc==1)          ! Global gather flag
mipan = ipan                ! local number of rows
mjpan = jpan                ! local number of columns
mil_g = il_g                ! global grid size
mxpr = il_g/ipan            ! number of processors over rows
mypr = il_g/jpan            ! number of processors over columns

! calculate number of levels
mg_maxlevel = 1
g = 1
gp = 2
do while ( mod( il_g, gp )==0 )
  g = g + 1
  gp = 2*gp
end do
mg_maxlevel = g
mg_maxlevel_local = mg_maxlevel

if ( myid==0 ) then
  write(6,*) "Initialise multi-grid arrays"
end if

allocate( mg(mg_maxlevel) )
mg(:)%comm_merge = 0

allocate( mg_bnds(0:nproc-1,mg_maxlevel) )

hipan = mipan
hjpan = mjpan

! calculate fine grid for finest grid level
g = 1
mg_ifull_maxcolour = 0
mg_npan = npan
mg(1)%npanx = npan
mg(1)%merge_len = 1
mg(1)%merge_row = 1
mg(1)%nmax = 1
mg(1)%ifull_coarse = 0

allocate( mg(1)%procmap(0:nproc-1) )
do n = 0,nproc-1
  mg(1)%procmap(n) = n ! default for now
end do


! check if coarse grid needs mgcollect
if ( mod( mipan, 2 )/=0 .or. mod( mjpan, 2 )/=0 .or. g==mg_maxlevel ) then

  if ( mod( mxpr, 2 )==0 .and. mod( mypr, 2 )==0 .and. g<mg_maxlevel ) then
   
    ! This case occurs when there are multiple processors on a panel.
    ! Consequently, npan should be 1 or 6.
    if ( npan>1 .and. npan<6 ) then
      write(6,*) "ERROR: Invalid gather4"
      call ccmpi_abort(-1)
    end if
    
    mg(1)%merge_len = 4
    mg(1)%merge_row = 2
    mg(1)%nmax = 4
    mxpr = mxpr/2
    mypr = mypr/2
    mipan = 2*mipan
    mjpan = 2*mjpan

    if ( myid==0 ) then
      write(6,*) "--> Multi-grid gather4 at level           ",g,mipan,mjpan
    end if

    allocate( mg(1)%merge_list(4) )

    ! find my processor and surrounding members of the gather
    nn = 1 - noff

    ix = -1
    jx = -1
    do jj = 1,mil_g,hjpan
      do ii = 1,mil_g,hipan
        if ( mg_fproc(1,ii,jj,nn)==myid ) then
          ix = ii
          jx = jj
          exit
        end if
      end do
      if ( ix>0 ) exit
    end do
    if ( ix<1 ) then
      write(6,*) "ERROR: Cannot locate processor in gather4"
      call ccmpi_abort(-1)
    end if
    ii = ix
    jj = jx
    if ( mod( (ii-1)/hipan, 2 )/=0 ) ii = ii - hipan
    if ( mod( (jj-1)/hjpan, 2 )/=0 ) jj = jj - hjpan
    ix = ix - ii ! offset for myid from ii
    jx = jx - jj ! offset for myid from jj
   
    mg(1)%merge_list(1) = mg_fproc(1,ii,      jj      ,nn)
    mg(1)%merge_list(2) = mg_fproc(1,ii+hipan,jj      ,nn)
    mg(1)%merge_list(3) = mg_fproc(1,ii,      jj+hjpan,nn)
    mg(1)%merge_list(4) = mg_fproc(1,ii+hipan,jj+hjpan,nn)
       
    do j = 1,mil_g,mjpan
      do i = 1,mil_g,mipan
        cid = mg_fproc(1,i+ix,j+jx,nn) ! processor in same merge position as myid
                                       ! we will maintain communications with this processor
        do jja = 1,mjpan,hjpan
          do iia = 1,mipan,hipan
            ! update fproc map with processor that owns this data
            mg(1)%procmap(mg_fproc_1(1,i+iia-1,j+jja-1,nn)) = cid  
            !mg(1)%fproc(i+iia-1,j+jja-1,nn) = cid
          end do
        end do
      end do
    end do

    mg(1)%merge_pos = -1
    do j = 1,4
      if ( mg(1)%merge_list(j)==myid ) then
        mg(1)%merge_pos = j
        exit
      end if
    end do
    if ( mg(g)%merge_pos<1 ) then
      write(6,*) "ERROR: Invalid merge_pos g,pos ",g,mg(g)%merge_pos
      call ccmpi_abort(-1)
    end if
      
    ! fix any remaining panels
    if ( npan==1 ) then
      ! must loop over all panels
      do nn = 0,npanels
        do j = 1,mil_g,mjpan
          do i = 1,mil_g,mipan
            cid = mg_fproc(1,i+ix,j+jx,nn) ! processor in same merge position as myid
                                           ! we will maintain communications with this processor
            do jja = 1,mjpan
              do iia = 1,mipan
                ! update fproc map with processor that owns this data
                mg(1)%procmap(mg_fproc_1(1,i+iia-1,j+jja-1,nn)) = cid
                !mg(1)%fproc(i+iia-1,j+jja-1,nn) = cid
              end do
            end do
          end do
        end do
      end do
    end if

  else if ( .not.lglob ) then ! collect all data to one processor
    lglob = .true.
    if ( uniform_decomp ) then
      mg(1)%merge_len = mxpr*mypr
    else
      mg(1)%merge_len = min( 6*mxpr*mypr, nproc )
      mg(1)%npanx = 1
      mg_npan = 6
    end if

    mg(1)%merge_row = mxpr
    mg(1)%nmax = mg(1)%merge_len
    mipan = mipan*mxpr
    mjpan = mjpan*mypr
    mxpr = 1
    mypr = 1

    if ( myid==0 ) then
      write(6,*) "--> Multi-grid gatherall at level         ",1,mipan,mjpan
    end if
      
    ! find gather members
    if ( uniform_decomp ) then
      allocate( mg(1)%merge_list(mg(1)%merge_len) )
      iqq = 0
      do jj = 1,mil_g,hjpan
        do ii = 1,mil_g,hipan
          iqq = iqq + 1
          mg(1)%merge_list(iqq) = mg_fproc(1,ii,jj,0)
        end do
      end do
      if ( iqq/=mg(1)%merge_len ) then
        write(6,*) "ERROR: merge_len mismatch ",iqq,mg(1)%merge_len,1
        call ccmpi_abort(-1)
      end if
      mg(1)%merge_pos = -1
      do i = 1,mg(1)%merge_len
        if ( mg(1)%merge_list(i)==myid ) then
          mg(1)%merge_pos = i
          exit
        end if
      end do
      if ( mg(1)%merge_pos<1 ) then
        write(6,*) "ERROR: Invalid merge_pos g,pos ",1,mg(1)%merge_pos
        call ccmpi_abort(-1)
      end if
    else
      allocate( mg(1)%merge_list(mg(1)%merge_len) )      
      iqq = 0
      do n = 1,6/npan
        nn = (n-1)*npan
        do jj = 1,mil_g,hjpan
          do ii = 1,mil_g,hipan
            iqq = iqq + 1
            mg(1)%merge_list(iqq) = mg_fproc(1,ii,jj,nn)
          end do
        end do
      end do
      if ( iqq/=mg(1)%merge_len ) then
        write(6,*) "ERROR: merge_len mismatch ",iqq,mg(1)%merge_len,1
        stop
      end if
      mg(1)%merge_pos = -1
      do i = 1,mg(1)%merge_len
        if ( mg(1)%merge_list(i)==myid ) then
          mg(1)%merge_pos=i
          exit
        end if
      end do
      if ( mg(1)%merge_pos<1 ) then
        write(6,*) "ERROR: Invalid merge_pos g,pos ",1,mg(1)%merge_pos
        call ccmpi_abort(-1)
      end if
    end if

    ! modify mg_fproc for remaining processor
    mg(1)%procmap(:) = myid
      
  else ! all data is already on one processor
    if ( 1/=mg_maxlevel ) then
      write(6,*) "ERROR: g/=mg_maxlevel ",1,mg_maxlevel
      call ccmpi_abort(-1)
    end if
    if ( myid==0 ) then
      write(6,*) "--> Multi-grid toplevel                   ",1,mipan,mjpan
    end if
    if ( .not.uniform_decomp ) then
      mg(1)%npanx = 1  
    end if
    mg(1)%merge_pos = 1
  end if
    
  ! define split and local comm
  if ( mg(1)%merge_len>1 ) then
    colour = mg(1)%merge_list(1)
    rank = mg(1)%merge_pos-1
    call ccmpi_commsplit(mg(1)%comm_merge,comm_world,colour,rank)
    if ( rank/=0 ) then
      mg_maxlevel_local = 0
      mg(1)%nmax = 0
    end if
    deallocate( mg(1)%merge_list )
  end if  
    
else
  if ( myid==0 ) then
    write(6,*) "--> Multi-grid fine level                 ",1,mipan,mjpan
  end if
  mg(1)%merge_pos = 1
end if

mg(1)%ipan = mipan
mg(1)%ifull = mipan*mjpan*mg_npan
mg(1)%ifull_fine = mg(1)%ifull/4

np = mg(1)%ifull
allocate( mg(1)%in(np), mg(1)%ie(np), mg(1)%is(np), mg(1)%iw(np) )
allocate( mg(1)%ine(np), mg(1)%ien(np), mg(1)%inw(np), mg(1)%iwn(np) )
allocate( mg(1)%ise(np), mg(1)%ies(np), mg(1)%isw(np), mg(1)%iws(np) )

call mg_index(1,mil_g,mipan,mjpan)

if ( mg_maxlevel>1 ) then
  allocate( mg(1)%fine(mg(1)%ifull_fine), mg(1)%fine_n(mg(1)%ifull_fine) )
  allocate( mg(1)%fine_e(mg(1)%ifull_fine), mg(1)%fine_ne(mg(1)%ifull_fine) )
  iqq = 0
  do n = 1,mg_npan
    do jj = 1,mjpan,2
      do ii = 1,mipan,2
        iq = indx(ii,jj,n-1,mipan,mjpan)
        iqq = iqq + 1
        mg(1)%fine(iqq) = iq
      end do
    end do
  end do
  mg(1)%fine_n = mg(1)%in(mg(1)%fine)
  mg(1)%fine_e = mg(1)%ie(mg(1)%fine)
  mg(1)%fine_ne = mg(1)%in(mg(g)%ie(mg(1)%fine))
end if

mg_maxsize = max( mg_maxsize, mg(1)%ifull+mg(1)%iextra )

! loop over upscaled grids
do g = 2,mg_maxlevel

  mipan = mipan/2
  mjpan = mjpan/2
  mil_g = mil_g/2
  hipan = mipan
  hjpan = mjpan

  ! Calculate size of grid at this level
  allocate( mg(g)%procmap(0:nproc-1) )
  mg(g)%procmap(:) = mg(g-1)%procmap(:)
  deallocate( mg(g-1)%procmap )

  ! default if no gather for upscaled grid
  mg(g)%npanx = npan
  mg(g)%merge_len = 1
  mg(g)%merge_row = 1
  mg(g)%nmax = 1
  
  ! check for multi-grid gather
  if ( mod(mipan,2)/=0 .or. mod(mjpan,2)/=0 .or. g==mg_maxlevel ) then ! grid cannot be subdivided on current processor
 
    if ( mod(mxpr,2)==0 .and. mod(mypr,2)==0 .and. g<mg_maxlevel ) then ! collect data over adjacent processors (proc=4)

      ! This case occurs when there are multiple processors on a panel.
      ! Consequently, npan should be 1 or 6.
      if ( npan>1 .and. npan<6 ) then
        write(6,*) "ERROR: Invalid gather4"
        call ccmpi_abort(-1)
      end if

      mg(g)%merge_len = 4
      mg(g)%merge_row = 2
      mg(g)%nmax = 4
      mxpr = mxpr/2
      mypr = mypr/2
      mipan = 2*mipan
      mjpan = 2*mjpan

      if ( myid==0 ) then
        write(6,*) "--> Multi-grid gather4 at level           ",g,mipan,mjpan
      end if

      allocate( mg(g)%merge_list(4) )
      
      nn = 1 - noff

      ! find my processor and surrounding members of the gather
      ix = -1
      jx = -1
      do jj = 1,mil_g,hjpan
        do ii = 1,mil_g,hipan
          if ( mg_fproc(g,ii,jj,nn)==myid ) then
            ix = ii
            jx = jj
            exit
          end if
        end do
        if ( ix>0 ) exit
      end do
      if ( ix<1 ) then
        write(6,*) "ERROR: Cannot locate processor in gather4 ",myid
        call ccmpi_abort(-1)
      end if
      ii = ix
      jj = jx
      mmx = mod( (ii-1)/hipan, 2 )
      mmy = mod( (jj-1)/hjpan, 2 )
      ii = ii - mmx*hipan ! left corner of merge
      jj = jj - mmy*hjpan ! bottom corner of merge
      ix = ix - ii ! offset for myid from ii
      jx = jx - jj ! offset for myid from jj
       
      mg(g)%merge_list(1) = mg_fproc(g,ii,      jj      ,nn)
      mg(g)%merge_list(2) = mg_fproc(g,ii+hipan,jj      ,nn)
      mg(g)%merge_list(3) = mg_fproc(g,ii,      jj+hjpan,nn)
      mg(g)%merge_list(4) = mg_fproc(g,ii+hipan,jj+hjpan,nn)
 
      do j = 1,mil_g,mjpan
        do i = 1,mil_g,mipan
          cid = mg_fproc(g,i+ix,j+jx,nn) ! processor in same merge position as myid
                                         ! we will maintain communications with this processor
          do jja = 1,mjpan
            do iia = 1,mipan
              ! update fproc map with processor that owns this data
              mg(g)%procmap(mg_fproc_1(g,i+iia-1,j+jja-1,nn)) = cid  
            end do
          end do
        end do
      end do

      mg(g)%merge_pos = -1
      do j = 1,4
        if ( mg(g)%merge_list(j)==myid ) then
          mg(g)%merge_pos = j
          exit
        end if
      end do
       
      if ( mg(g)%merge_pos<1 ) then
        write(6,*) "ERROR: Invalid merge_pos g,pos ",g,mg(g)%merge_pos
        call ccmpi_abort(-1)
      end if

      ! fix any remaining panels
      if ( npan==1 ) then
        ! must loop over all panels
        do nn = 0,npanels
          do j = 1,mil_g,mjpan
            do i = 1,mil_g,mipan
              cid = mg_fproc(g,i+ix,j+jx,nn) ! processor in same merge position as myid
                                             ! we will maintain communications with this processor
              do jja = 1,mjpan
                do iia = 1,mipan
                  ! update fproc map with processor that owns this data
                  mg(g)%procmap(mg_fproc_1(g,i+iia-1,j+jja-1,nn)) = cid  
                end do
              end do
            end do
          end do
        end do
      end if
    
    else if ( .not.lglob ) then ! collect all data to one processor
      lglob = .true.
      if ( uniform_decomp ) then
        mg(g)%merge_len = mxpr*mypr
      else
        mg(g)%merge_len = min( 6*mxpr*mypr, nproc )
        mg(g)%npanx = 1
        mg_npan = 6
      end if

      mg(g)%merge_row = mxpr
      mg(g)%nmax = mg(g)%merge_len
      mipan = mipan*mxpr
      mjpan = mjpan*mypr
      mxpr = 1
      mypr = 1

      if ( myid==0 ) then
        write(6,*) "--> Multi-grid gatherall at level         ",g,mipan,mjpan
      end if
      
      ! find gather members
      if ( uniform_decomp ) then
        allocate( mg(g)%merge_list(mg(g)%merge_len) )
        iqq = 0
        do jj = 1,mil_g,hjpan
          do ii = 1,mil_g,hipan
            iqq = iqq + 1
            mg(g)%merge_list(iqq) = mg_fproc(g,ii,jj,0)
          end do
        end do
        if ( iqq/=mg(g)%merge_len ) then
          write(6,*) "ERROR: merge_len mismatch ",iqq,mg(g)%merge_len,g
          call ccmpi_abort(-1)
        end if
        mg(g)%merge_pos = -1
        do i = 1,mg(g)%merge_len
          if ( mg(g)%merge_list(i)==myid ) then
            mg(g)%merge_pos = i
            exit
          end if
        end do
        if ( mg(g)%merge_pos<1 ) then
          write(6,*) "ERROR: Invalid merge_pos g,pos ",g,mg(g)%merge_pos
          call ccmpi_abort(-1)
        end if
      else
        allocate( mg(g)%merge_list(mg(g)%merge_len) )      
        iqq = 0
        do n = 1,6/npan
          nn = (n-1)*npan
          do jj = 1,mil_g,hjpan
            do ii = 1,mil_g,hipan
              iqq = iqq + 1
              mg(g)%merge_list(iqq) = mg_fproc(g,ii,jj,nn)
            end do
          end do
        end do
        if ( iqq/=mg(g)%merge_len ) then
          write(6,*) "ERROR: merge_len mismatch ",iqq,mg(g)%merge_len,g
          stop
        end if
        mg(g)%merge_pos = -1
        do i = 1,mg(g)%merge_len
          if ( mg(g)%merge_list(i)==myid ) then
            mg(g)%merge_pos=i
            exit
          end if
        end do
        if ( mg(g)%merge_pos<1 ) then
          write(6,*) "ERROR: Invalid merge_pos g,pos ",g,mg(g)%merge_pos
          call ccmpi_abort(-1)
        end if
      end if

      ! modify mg_fproc for remaining processor
      mg(g)%procmap(:) = myid
      
    else ! all data is already on one processor
      if ( g/=mg_maxlevel ) then
        write(6,*) "ERROR: g/=mg_maxlevel ",g,mg_maxlevel
        call ccmpi_abort(-1)
      end if
      if ( myid==0 ) then
        write(6,*) "--> Multi-grid toplevel                   ",g,mipan,mjpan
      end if
      if ( .not.uniform_decomp ) then
        mg(g)%npanx = 1  
      end if
      mg(g)%merge_pos = 1
    end if

    ! define split and local comm
    if ( mg(g)%merge_len>1 ) then
      colour = mg(g)%merge_list(1)
      rank = mg(g)%merge_pos-1
      call ccmpi_commsplit(mg(g)%comm_merge,comm_world,colour,rank)
      if ( rank/=0 ) then
        mg_maxlevel_local = min( g-1, mg_maxlevel_local )
        mg(g)%nmax = 0
      end if
      deallocate( mg(g)%merge_list )
    end if
  
  else
    if ( myid==0 ) then
      write(6,*) "--> Multi-grid local subdivision at level ",g,mipan,mjpan
    end if
    ! no messages sent, but we allocate this array for the
    ! coarse calculation below
    mg(g)%merge_pos = 1
  end if
  

  ! total number of points over all panels on this processor
  mg(g)%ipan = mipan
  mg(g)%ifull = mipan*mjpan*mg_npan
  mg(g)%ifull_fine = mg(g)%ifull/4  
  nrow = mg(g)%merge_row
  ncol = mg(g)%merge_len/nrow
  drow = mipan/nrow
  dcol = mjpan/ncol
  npanx = mg_npan
  
  if ( .not.uniform_decomp .and. lglob ) then
    npanx = 1
    dcol = 6*mjpan/ncol
  end if
  
  np = mg(g)%ifull
  allocate( mg(g)%in(np), mg(g)%ie(np), mg(g)%is(np), mg(g)%iw(np) )
  allocate( mg(g)%ine(np), mg(g)%ien(np), mg(g)%inw(np), mg(g)%iwn(np) )
  allocate( mg(g)%ise(np), mg(g)%ies(np), mg(g)%isw(np), mg(g)%iws(np) )
  
  call mg_index(g,mil_g,mipan,mjpan)

  gmax = min( mg_maxlevel-1, mg_maxlevel_local )
  
  ! ifine is the index of the SW point on the next finer grid.
  if ( g<=gmax ) then
    np = mg(g)%ifull_fine
    allocate( mg(g)%fine(np), mg(g)%fine_n(np) )
    allocate( mg(g)%fine_e(np), mg(g)%fine_ne(np) )
    ! mipan and mjpan should always be an even number here
    iqq = 0
    do n = 1,mg_npan
      do jj = 1,mjpan,2
        do ii = 1,mipan,2
          iq = indx(ii,jj,n-1,mipan,mjpan)
          iqq = iqq + 1
          mg(g)%fine(iqq) = iq
        end do
      end do
    end do
    mg(g)%fine_n = mg(g)%in(mg(g)%fine)
    mg(g)%fine_e = mg(g)%ie(mg(g)%fine)
    mg(g)%fine_ne = mg(g)%in(mg(g)%ie(mg(g)%fine))
  end if
  
  mg_maxsize = max( mg_maxsize, mg(g)%ifull+mg(g)%iextra )

  ! Set up pointers for neighbours on the next coarser grid.
  ! ic1 is the nearest neigbour, ic2 and ic3 are equally distant and
  ! ic4 is the most distant.
  ! These are best defined separately on each face. The code here
  ! requires an even number of points to ensure the coarser grid
  ! always sits inside the fine grid.

  if ( g<=gmax+1 ) then
    mg(g)%ifull_coarse = mg(g-1)%ifull
    np = mg(g)%ifull_coarse
    allocate( mg(g)%coarse_a(np), mg(g)%coarse_b(np), mg(g)%coarse_c(np) )
    mg(g)%coarse_a = 0 ! unassigned
 
    do n = 1,npanx
  
      na = mg(g)%merge_pos
      nx = mod( na-1, nrow ) + 1
      ny = (na-1)/nrow + 1
      sii = (nx-1)*drow + 1
      eii = nx*drow
      sjj = (ny-1)*dcol + 1
      ejj = ny*dcol
    
      do jj = sjj,ejj
        jja = 2*(jj-sjj) + 1
        do ii = sii,eii
          iia = 2*(ii-sii) + 1
       
          iqq = indx(ii,jj,n-1,mipan,mjpan)   ! coarse grid
        
          ! odd, odd          
          iq = indx(iia,jja,n-1,2*drow,2*dcol)     ! fine grid
          mg(g)%coarse_a(iq) =          iqq
          mg(g)%coarse_b(iq) = mg(g)%is(iqq)
          mg(g)%coarse_c(iq) = mg(g)%iw(iqq)
        
          ! odd, even
          iq = indx(iia,jja+1,n-1,2*drow,2*dcol)   ! fine grid
          mg(g)%coarse_a(iq) =          iqq
          mg(g)%coarse_b(iq) = mg(g)%in(iqq)
          mg(g)%coarse_c(iq) = mg(g)%iw(iqq)
          
          ! even, odd
          iq = indx(iia+1,jja,n-1,2*drow,2*dcol)   ! fine grid 
          mg(g)%coarse_a(iq) =          iqq
          mg(g)%coarse_b(iq) = mg(g)%is(iqq)
          mg(g)%coarse_c(iq) = mg(g)%ie(iqq)

          ! even, even
          iq = indx(iia+1,jja+1,n-1,2*drow,2*dcol) ! fine grid
          mg(g)%coarse_a(iq) =          iqq
          mg(g)%coarse_b(iq) = mg(g)%in(iqq)
          mg(g)%coarse_c(iq) = mg(g)%ie(iqq)
      
        end do
      end do
    
    end do
  
  end if

  ! free some memory
  if ( g>=gmax+2 ) then
    deallocate( mg(g)%in, mg(g)%ie, mg(g)%is, mg(g)%iw )
    deallocate( mg(g)%ine, mg(g)%ien, mg(g)%inw, mg(g)%iwn )
    deallocate( mg(g)%ise, mg(g)%ies, mg(g)%isw, mg(g)%iws )
  end if
  
end do


if ( myid==0 ) then
  mg_maxlevel_decomp = mg_maxlevel    
  mg_minsize = 6*mil_g*mil_g
else
  mg_maxlevel_decomp = mg_maxlevel_local
  mg_minsize = 0
end if

! free some memory
deallocate( mg(mg_maxlevel)%procmap )

sorfirst = .false.

return
end subroutine mgsor_init

subroutine mgzz_init(zz,zzn,zze,zzw,zzs)

use cc_mpi
use newmpar_m

implicit none

integer g, np
real, dimension(ifull), intent(in) :: zz, zzn, zze, zzw, zzs
real, dimension(mg_maxsize,5) :: dum

  ! Later equations use zz/delta**2. It's more efficient to put
  ! this scaling in here. delta = delta_0 * 2**(m-1) so increases
  ! by a factor of 2 each grid. Therefore using a factor of 1/4
  ! in the averaging calculation will compensate.

if ( myid==0 ) then
  write(6,*) "Initialise atmosphere multi-grid coupling arrays"
end if

if ( zzfirst ) then
  do g = 1,gmax+1
    np = mg(g)%ifull
    allocate( mg(g)%zzn(np),mg(g)%zze(np), mg(g)%zzs(np),mg(g)%zzw(np),mg(g)%zz(np) )
  end do
  zzfirst = .false.
end if

! no need for bounds call as all indices lie within this processor
g = 1

dum(1:ifull,1) = zzn(1:ifull)
dum(1:ifull,2) = zzs(1:ifull)
dum(1:ifull,3) = zze(1:ifull)
dum(1:ifull,4) = zzw(1:ifull)
dum(1:ifull,5) = zz(1:ifull)
call mgcollect(1,dum)
np = mg(1)%ifull
mg(1)%zzn(1:np) = dum(1:np,1)
mg(1)%zzs(1:np) = dum(1:np,2)
mg(1)%zze(1:np) = dum(1:np,3)
mg(1)%zzw(1:np) = dum(1:np,4)
mg(1)%zz(1:np)  = dum(1:np,5)

do g = 2,gmax+1
  np = mg(g-1)%ifull_fine
  dum(1:np,1) = 0.25*dfac*(mg(g-1)%zzn(mg(g-1)%fine(1:np)  )+mg(g-1)%zzn(mg(g-1)%fine_n(1:np)) &
                          +mg(g-1)%zzn(mg(g-1)%fine_e(1:np))+mg(g-1)%zzn(mg(g-1)%fine_ne(1:np)))
  dum(1:np,2) = 0.25*dfac*(mg(g-1)%zzs(mg(g-1)%fine(1:np)  )+mg(g-1)%zzs(mg(g-1)%fine_n(1:np)) &
                          +mg(g-1)%zzs(mg(g-1)%fine_e(1:np))+mg(g-1)%zzs(mg(g-1)%fine_ne(1:np)))
  dum(1:np,3) = 0.25*dfac*(mg(g-1)%zze(mg(g-1)%fine(1:np)  )+mg(g-1)%zze(mg(g-1)%fine_n(1:np)) &
                          +mg(g-1)%zze(mg(g-1)%fine_e(1:np))+mg(g-1)%zze(mg(g-1)%fine_ne(1:np)))
  dum(1:np,4) = 0.25*dfac*(mg(g-1)%zzw(mg(g-1)%fine(1:np)  )+mg(g-1)%zzw(mg(g-1)%fine_n(1:np)) &
                          +mg(g-1)%zzw(mg(g-1)%fine_e(1:np))+mg(g-1)%zzw(mg(g-1)%fine_ne(1:np)))
  dum(1:np,5) = 0.25*dfac*(mg(g-1)%zz(mg(g-1)%fine(1:np)  ) +mg(g-1)%zz(mg(g-1)%fine_n(1:np))  &
                          +mg(g-1)%zz(mg(g-1)%fine_e(1:np)) +mg(g-1)%zz(mg(g-1)%fine_ne(1:np)))
  call mgcollect(g,dum)
  np = mg(g)%ifull
  mg(g)%zzn(1:np) = dum(1:np,1)
  mg(g)%zzs(1:np) = dum(1:np,2)
  mg(g)%zze(1:np) = dum(1:np,3)
  mg(g)%zzw(1:np) = dum(1:np,4)
  mg(g)%zz(1:np)  = dum(1:np,5)
end do

return
end subroutine mgzz_init

subroutine mgunpack_nsew(g,data_in,data_n,data_s,data_e,data_w)

use cc_mpi
use indices_m

implicit none

integer, intent(in) :: g
integer ng, mg_npan, mg_ipan, mg_jpan
integer i, j, n, iq
real, dimension(:), intent(in) :: data_in
real, dimension(:), intent(out) :: data_n, data_s, data_e, data_w

ng = mg(g)%ifull

mg_npan = mg(g)%npanx
mg_ipan = mg(g)%ipan
mg_jpan = ng/(mg(g)%ipan*mg_npan)

data_e(1:ng-1)       = data_in(2:ng)
data_w(2:ng)         = data_in(1:ng-1)
data_n(1:ng-mg_ipan) = data_in(mg_ipan+1:ng)
data_s(mg_ipan+1:ng) = data_in(1:ng-mg_ipan)
do n = 1,mg_npan
  do j = 1,mg_jpan
    iq = 1 + (j-1)*mg_ipan + (n-1)*mg_ipan*mg_jpan
    data_w(iq) = data_in(mg(g)%iw(iq))
    iq = mg_ipan + (j-1)*mg_ipan + (n-1)*mg_ipan*mg_jpan
    data_e(iq) = data_in(mg(g)%ie(iq))
  end do
  do i = 1,mg_ipan
    iq = i + (n-1)*mg_ipan*mg_jpan
    data_s(iq) = data_in(mg(g)%is(iq))
    iq = i + (mg_jpan-1)*mg_ipan + (n-1)*mg_ipan*mg_jpan
    data_n(iq) = data_in(mg(g)%in(iq))
  end do
end do
    
return
end subroutine mgunpack_nsew

end module helmsolve
