! Conformal Cubic Atmospheric Model
    
! Copyright 2015-2024 Commonwealth Scientific Industrial Research Organisation (CSIRO)
    
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
    
subroutine upglobal

use aerointerface          ! Aerosol interface
use arrays_m               ! Atmosphere dyamics prognostic arrays
use cc_mpi                 ! CC MPI routines
use cfrac_m                ! Cloud fraction
use const_phys             ! Physical constants
use diag_m                 ! Diagnostic routines
use epst_m                 ! Off-centre terms
use indices_m              ! Grid index arrays
use kuocom_m               ! JLM convection
use liqwpar_m              ! Cloud water mixing ratios
use map_m                  ! Grid map arrays
use newmpar_m              ! Grid parameters
use nharrs_m               ! Non-hydrostatic atmosphere arrays
use nlin_m                 ! Atmosphere non-linear dynamics
use parm_m                 ! Model configuration
use parmdyn_m              ! Dynamics parameters
use parmhor_m              ! Horizontal advection parameters
use sbar_m                 ! Saved dynamic arrays
use sigs_m                 ! Atmosphere sigma levels
use staguvmod              ! Reversible grid staggering
use tkeeps, only : tke,eps ! TKE-EPS boundary layer
use tracers_m              ! Tracer data
use unn_m                  ! Saved dynamic arrays
use vadv                   ! Vertical advection
use vecsuv_m               ! Map to cartesian coordinates
use vvel_m                 ! Additional vertical velocity
use work3f_m               ! Grid work arrays
use xarrs_m                ! Saved dynamic arrays
use xyzinfo_m              ! Grid coordinate arrays

implicit none

integer, parameter :: ntest=0       ! ~8+ for diagnostic stability tests
integer ii, intsch, iq, jj, k, kk
integer idjdd, nstart
integer, save :: numunstab = 0
integer, dimension(ifull) :: nits
integer, dimension(max(naero,ntrac,6)) :: nfield
real, dimension(ifull) :: nvadh_inv_pass
real, dimension(ifull+iextra,kl,8) :: bb
real, dimension(ifull+iextra,kl,3) :: uvw
real, dimension(ifull+iextra,kl) :: dd
real, dimension(ifull+iextra) :: aa
real, dimension(ifull,kl) :: theta
real, dimension(ifull) :: sdmx, tempry
real denb, vecdot, vdot1, vdot2
real vec1x, vec1y, vec1z
real vec2x, vec2y, vec2z
real vec3x, vec3y, vec3z
real, dimension(kl) :: factr
real(kind=8), dimension(ifull,kl) :: x3d, y3d, z3d
#ifdef debug
integer, save :: num_hight = 0
real, dimension(kl) :: diag_temp
#endif

call START_LOG(upglobal_begin)


select case( intsch_mode )
  case(1)
    intsch = 1
  case(0)
    intsch = 0
  case default
    intsch = mod(ktau, 2)
end select    

do k = 1,kl
  ! finish off RHS terms; this coriolis term was once in nonlin
  ux(1:ifull,k) = ux(1:ifull,k) + 0.5*dt*(1.-epsf)*f(1:ifull)*v(1:ifull,k) ! end of Eq. 129
  vx(1:ifull,k) = vx(1:ifull,k) - 0.5*dt*(1.-epsf)*f(1:ifull)*u(1:ifull,k) ! end of Eq. 130
end do      ! k loop

call depts1(x3d,y3d,z3d,intsch)
      
!     calculate factr for choice of nt_adv, as usually used
select case(nt_adv)
  case(0)
    factr(:) = 0.
  case(3)  ! 1. up to sig=.3
    where ( sig(1:kl) >= 0.3 )
      factr(1:kl) = stdlapse*(rdry*nritch_t/grav)
    elsewhere
      factr(1:kl) = 0.      
    end where
  case(4) ! (1, .5, 0) for sig=(1, .75, .5)
    factr(1:kl) = max(2.*sig(1:kl)-1., 0.)*stdlapse*(rdry*300./grav)
  case(5) ! 1 to 0 for sig=1 to 0
    factr(1:kl) = sig(1:kl)*stdlapse*(rdry*nritch_t/grav)
  case(6) ! (1, .5625, 0) for sig=(1, .5, .2)
    factr(1:kl) = max(0.,1.25*(sig(1:kl)-.2)*(2.-sig(1:kl)))*stdlapse*(rdry*nritch_t/grav)
  case(7) ! 1 up to .4, then lin decr. to .2, then 0
    factr(1:kl) = max(0.,min(1.,(sig(1:kl)-.2)/(.4-.2)))*stdlapse*(rdry*nritch_t/grav)
  case(8) ! .8 up to .4, then lin decr. to .2, then 0
    factr(1:kl) = .8*max(0.,min(1.,(sig(1:kl)-.2)/(.4-.2)))*stdlapse*(rdry*nritch_t/grav)
  case(9) ! (1,1,.84375,.5,.15625,0) for sig=(1,.6,.5,.4,.3,.2)  
    where ( sig(1:kl) > 0.6 )
      factr(1:kl) = stdlapse*(rdry*nritch_t/grav)
    elsewhere ( sig(1:kl) >= 0.2 )
      factr(1:kl) = 3.*((sig(1:kl)-.2)/.4)**2 -2.*((sig(1:kl)-.2)/.4)**3*stdlapse*(rdry*nritch_t/grav)
    elsewhere
      factr(1:kl) = 0.
    end where
  case(10) ! (1,1,.741,.259,0) for sig=(1,.5,.4,.3,.2) 
    where ( sig(1:kl) > 0.5 )
      factr(1:kl) = stdlapse*(rdry*nritch_t/grav)
    elsewhere ( sig(1:kl) >= 0.2 )
      factr(1:kl) = 3.*((sig(1:kl)-.2)/.3)**2 -2.*((sig(1:kl)-.2)/.3)**3*stdlapse*(rdry*nritch_t/grav)  
    elsewhere
      factr(1:kl) = 0.
    end where
end select

#ifdef debug  
if ( mydiag ) then
  if ( tx(idjd,kl)>264. ) then  !cb
    write(6,*)
    write(6,*) 'upglobal a',ktau,id,jd,tx(idjd,kl)
  end if    ! (tx(idjd,kl)>264.)
end if

if ( num_hight<100 ) then
  do iq = 1,ifull
    if ( tx(iq,kl)>264. ) then  !cb
      write(6,*) 'upglobal ktau,myid,iq,large_tx  ',ktau,myid,iq,tx(iq,kl)
      diag_temp(1:kl) = sdot(iq,1:kl)
      write (6,"('sdot_iq',9f7.3/7x,9f7.3)") diag_temp(1:kl)
      num_hight=num_hight+1
    end if
  end do
end if 
#endif

aa(1:ifull) = zs(1:ifull)/(rdry*nritch_t)    ! save zs/(r*t) for nt_adv schemes 
do k = 1,kl   
  dd(1:ifull,k) = aa(1:ifull)
end do     ! k loop


!-------------------------moved up here May 06---------------------------
! N.B. this moved one is doing vadv on just extra pslx terms    
sdmx(:) = maxval(abs(sdot), 2)
nits(:) = int(1.+sdmx(:)/2.)
nvadh_inv_pass(:) = 0.5/real(nits(:)) ! use - for nvadu
call vadvtvd(tx,ux,vx,nvadh_inv_pass,nits)
if ( (diag.or.nmaxpr==1) .and. mydiag ) then
  write(6,*) 'in upglobal after vadv1'
  write (6,"('tx_a',9f8.2)")   tx(idjd,:)
  write (6,"('ux_a',9f8.3)")   ux(idjd,:)
  write (6,"('vx_a',9f8.3)")   vx(idjd,:)
  write (6,"('qg_a',9f8.3)")   qg(idjd,:)
end if
!------------------------------------------------------------------

do k = 1,kl   
  ! N.B. [D + dsigdot/dsig] saved in adjust5 (or updps) as pslx
  pslx(1:ifull,k) = psl(1:ifull) - pslx(1:ifull,k)*dt*.5*(1.-epst(:)) + aa(1:ifull)
  tx(1:ifull,k)   = tx(1:ifull,k) + aa(1:ifull)*factr(k)   !cy  
end do   ! k

if ( nmaxpr==1 .and. nproc==1 ) then
  write(6,*) 'pslx_3p before advection'
  write (6,"('pslx_b',9f8.4)") pslx(idjd,:)
  write (6,"(i6,8i8)") (ii,ii=id-4,id+4)
  write (6,"(9f8.4)") ((pslx(max(min(ii+jj*il,ifull),1),nlv),ii=idjd-4,idjd+4),jj=2,-2,-1)
end if


! call bounds before calling ints
call START_LOG(ints_begin)
if ( mup/=0 ) then
  call bounds(dd,corner=.true.)
  if ( nh/=0 ) then
    call bounds(pslx,nrows=2)      
    call bounds(h_nh,nrows=2)
  else
    call bounds(pslx,nrows=2)      
  end if ! nh/=0
  call bounds(tx,nrows=2)
end if    ! mup/=0
do k = 1,kl
  uvw(1:ifull,k,1) = ax(1:ifull)*ux(1:ifull,k) + bx(1:ifull)*vx(1:ifull,k)
  uvw(1:ifull,k,2) = ay(1:ifull)*ux(1:ifull,k) + by(1:ifull)*vx(1:ifull,k)
  uvw(1:ifull,k,3) = az(1:ifull)*ux(1:ifull,k) + bz(1:ifull)*vx(1:ifull,k)
end do
if ( mup/=0 ) then
  call bounds(uvw,nrows=2)
end if
if ( mspec==1 .and. mup/=0 ) then
  if ( ldr/=0 .and. ncloud>=100 .and. ncloud<200 ) then
    ! Lin microphysics  
    call bounds(qg,nrows=2)
    call bounds(qlg,nrows=2)
    call bounds(qfg,nrows=2)
    call bounds(stratcloud,nrows=2)
    call bounds(ni,nrows=2)   
    call bounds(qrg,nrows=2)
    call bounds(qsng,nrows=2)
    call bounds(qgrg,nrows=2)
  else if ( ldr/=0 ) then
    ! Leon microphysics  
    call bounds(qg,nrows=2)
    call bounds(qlg,nrows=2)
    call bounds(qfg,nrows=2)
    call bounds(stratcloud,nrows=2)
    call bounds(qrg,nrows=2)
    call bounds(qsng,nrows=2)
    call bounds(qgrg,nrows=2)
  else
    ! diagnostic (depreciated)  
    call bounds(qg,nrows=2)
  end if    ! ldr/=0
  if ( ntrac>0 ) then
    call bounds(tr,nrows=2)
  end if
  if ( nvmix==6 .or. nvmix==9 ) then
    call bounds(tke,nrows=2)
    call bounds(eps,nrows=2)
  endif                 ! nvmix==6 .or. nvmix==9
  if ( abs(iaero)>=2 .and. nhstest>=0 ) then
    call bounds(xtg,nrows=2)
  end if
end if     ! mspec==1
call END_LOG(ints_end)


#if GPU
!$acc data create(xg,yg,nface)
!$acc update device(xg,yg,nface)
#endif


if ( mup/=0 ) then
  call ints_bl(dd,intsch,nface,xg,yg)  ! advection on all levels
  if ( nh/=0 ) then
    ! non-hydrostatic version
    bb(:,:,1) = pslx(:,:)
    bb(:,:,2) = h_nh(:,:)
    call ints(bb(:,:,1:2),2,intsch,nface,xg,yg,1)      
    pslx(:,:) = bb(:,:,1)
    h_nh(:,:) = bb(:,:,2)
  else
    call ints(pslx,1,intsch,nface,xg,yg,1)    
  end if ! nh/=0
  call ints(tx,1,intsch,nface,xg,yg,3)
end if    ! mup/=0

do k = 1,kl
  pslx(1:ifull,k) = pslx(1:ifull,k) - dd(1:ifull,k)      
  tx(1:ifull,k)   = tx(1:ifull,k)   - dd(1:ifull,k)*factr(k)
end do

!------------------------------------------------------------------
if ( nmaxpr==1 .and. nproc==1 ) then
  write(6,*) 'pslx_3p & dd after advection'
  write (6,"('pslx_a',9f8.4)") pslx(idjd,:)
  write (6,"('aa#',9f8.4)") ((aa(ii+jj*il),ii=idjd-1,idjd+1),jj=-1,1)
  write (6,"('dd1#',9f8.4)") ((dd(ii+jj*il,1),ii=idjd-1,idjd+1),jj=-1,1)
  write (6,"('dd_a',9f8.4)") dd(idjd,:)
  write (6,"('nface',18i4)") nface(idjd,:)
  write (6,"('xg',9f8.4)") xg(idjd,:)
  write (6,"('yg',9f8.4)") yg(idjd,:)
  write (6,"(i6,8i8)") (ii,ii=id-4,id+4)
  idjdd=max(5+2*il,min(idjd,ifull-4-2*il))  ! for following prints
  write (6,"(39f8.4)") ((pslx(ii+jj*il,nlv),ii=idjdd-4,idjdd+4),jj=2,-2,-1)
! uc(1:ifull,1)=-pslx(1:ifull,1)*dsig(1) 
! do k=2,kl
!    uc(1:ifull,1)=uc(1:ifull,1)-pslx(1:ifull,k)*dsig(k)
! enddo
!  write(6,*) 'integ pslx after advection'
!  write (6,"(i6,8i8)") (ii,ii=id-4,id+4)
!  write (6,"(9f8.4)") ((uc(ii+jj*il,1),ii=idjdd-4,idjdd+4),jj=2,-2,-1)
!  write(6,*) 'corresp integ ps after advection'
!  write (6,"(i6,8i8)") (ii,ii=id-4,id+4)
!  write (6,"(9f8.2)") ((1.e5*exp(uc(ii+jj*il,1)),ii=idjdd-4,idjdd+4),jj=2,-2,-1)
end if

!     now comes ux & vx section
!if ( diag ) then
!  if ( mydiag ) then
!    write(6,*) 'unstaggered now as uavx and vavx: globpea uses ux & vx'
!    write(6,*) 'ux ',(ux(idjd,kk),kk=1,nlv)
!  end if
!  call printa('uavx',ux,ktau,nlv,ia,ib,ja,jb,0.,1.)
!  if ( mydiag ) write(6,*) 'vx ',(vx(idjd,kk),kk=1,nlv)
!  call printa('vavx',vx,ktau,nlv,ia,ib,ja,jb,0.,1.)
!  if ( mydiag ) write(6,*)'unstaggered u and v as uav and vav: globpea uses u & v'
!  call printa('uav ',u,ktau,nlv,ia,ib,ja,jb,0.,1.)
!  call printa('vav ',v,ktau,nlv,ia,ib,ja,jb,0.,1.)
!end if

if ( diag ) then
  if ( mydiag ) then
    write(6,*) 'uc,vc,wc before advection'
    write (6,'(a,18e20.10)') 'uc,vc,wc ',uvw(idjd,nlv,1),uvw(idjd,nlv,2),uvw(idjd,nlv,3)
  end if
  call printa('uc  ',uvw(:,:,1),ktau,nlv,ia,ib,ja,jb,0.,1.)
  call printa('vc  ',uvw(:,:,2),ktau,nlv,ia,ib,ja,jb,0.,1.)
  call printa('wc  ',uvw(:,:,3),ktau,nlv,ia,ib,ja,jb,0.,1.)
  call printa('xg  ',xg,ktau,nlv,ia,ib,ja,jb,0.,1.)
  call printa('yg  ',yg,ktau,nlv,ia,ib,ja,jb,0.,1.)
  if ( mydiag ) write(6,*) 'nface ',nface(idjd,:)
end if

if ( mup/=0 ) then
  call ints(uvw,3,intsch,nface,xg,yg,2)
end if

if ( diag ) then
  if ( mydiag ) write(6,*) 'uc,vc,wc after advection'
  call printa('uc  ',uvw(:,:,1),ktau,nlv,ia,ib,ja,jb,0.,1.)
  call printa('vc  ',uvw(:,:,2),ktau,nlv,ia,ib,ja,jb,0.,1.)
  call printa('wc  ',uvw(:,:,3),ktau,nlv,ia,ib,ja,jb,0.,1.)
end if

! rotate wind vector to arrival point
do k = 1,kl
  do iq = 1,ifull
    ! the following normalization may be done, but has ~zero effect
    ! dena=sqrt(x3d(iq)**2+y3d(iq)**2+z3d(iq)**2)
    ! x3d(iq)=x3d(iq)/dena
    ! y3d(iq)=y3d(iq)/dena
    ! z3d(iq)=z3d(iq)/dena
    ! cross product n1xn2 into vec1
    vec1x = real(y3d(iq,k)*z(iq) - y(iq)*z3d(iq,k))
    vec1y = real(z3d(iq,k)*x(iq) - z(iq)*x3d(iq,k))
    vec1z = real(x3d(iq,k)*y(iq) - x(iq)*y3d(iq,k))
    denb = vec1x**2 + vec1y**2 + vec1z**2
    ! N.B. rotation formula is singular for small denb,
    ! but the rotation is unnecessary in this case
    if ( denb>1.e-4 ) then
      vecdot = real(x3d(iq,k)*x(iq) + y3d(iq,k)*y(iq) + z3d(iq,k)*z(iq))
      vec2x = real(x3d(iq,k)*vecdot - x(iq))
      vec2y = real(y3d(iq,k)*vecdot - y(iq))
      vec2z = real(z3d(iq,k)*vecdot - z(iq))
      vec3x = real(x3d(iq,k) - vecdot*x(iq))
      vec3y = real(y3d(iq,k) - vecdot*y(iq))
      vec3z = real(z3d(iq,k) - vecdot*z(iq))
      vdot1 = (vec1x*uvw(iq,k,1) + vec1y*uvw(iq,k,2) + vec1z*uvw(iq,k,3))/denb
      vdot2 = (vec2x*uvw(iq,k,1) + vec2y*uvw(iq,k,2) + vec2z*uvw(iq,k,3))/denb
      uvw(iq,k,1) = vdot1*vec1x + vdot2*vec3x
      uvw(iq,k,2) = vdot1*vec1y + vdot2*vec3y
      uvw(iq,k,3) = vdot1*vec1z + vdot2*vec3z
    end if ! (denb>1.e-4)
    ! convert back to conformal-cubic velocity components (unstaggered)
    ! globpea: this can be sped up later
    ux(iq,k) = ax(iq)*uvw(iq,k,1) + ay(iq)*uvw(iq,k,2) + az(iq)*uvw(iq,k,3)
    vx(iq,k) = bx(iq)*uvw(iq,k,1) + by(iq)*uvw(iq,k,2) + bz(iq)*uvw(iq,k,3)
  end do   ! iq
end do     ! k


if ( diag .and. k==nlv ) then
  if ( mydiag ) write(6,*) 'after advection in upglobal; unstaggered ux and vx:'
  call printa('ux  ',ux,ktau,nlv,ia,ib,ja,jb,0.,1.)
  call printa('vx  ',vx,ktau,nlv,ia,ib,ja,jb,0.,1.)
end if

if ( mspec==1 .and. mup/=0 ) then   ! advect qg after preliminary step
  if ( ldr/=0 .and. ncloud>=100 .and. ncloud<200 ) then
    ! Lin microphysics  
    bb(:,:,1) = qg(:,:)
    bb(:,:,2) = qlg(:,:)
    bb(:,:,3) = qfg(:,:)
    bb(:,:,4) = stratcloud(:,:)
    bb(:,:,5) = qrg(:,:)
    bb(:,:,6) = qsng(:,:)
    bb(:,:,7) = qgrg(:,:)
    bb(:,:,8) = ni(:,:)
    call ints(bb(:,:,1:8),8,intsch,nface,xg,yg,4)
    qg(1:ifull,1:kl) = bb(1:ifull,1:kl,1)
    qlg(1:ifull,1:kl) = bb(1:ifull,1:kl,2)
    qfg(1:ifull,1:kl) = bb(1:ifull,1:kl,3)
    stratcloud(1:ifull,1:kl) = bb(1:ifull,1:kl,4)
    qrg(1:ifull,1:kl) = bb(1:ifull,1:kl,5)
    qsng(1:ifull,1:kl) = bb(1:ifull,1:kl,6)
    qgrg(1:ifull,1:kl) = bb(1:ifull,1:kl,7)
    ni(1:ifull,1:kl) = bb(1:ifull,1:kl,8)  
  else if ( ldr/=0 ) then
    ! Leon microphysics  
    bb(:,:,1) = qg(:,:)
    bb(:,:,2) = qlg(:,:)
    bb(:,:,3) = qfg(:,:)
    bb(:,:,4) = stratcloud(:,:)
    bb(:,:,5) = qrg(:,:)
    bb(:,:,6) = qsng(:,:)
    bb(:,:,7) = qgrg(:,:)
    call ints(bb(:,:,1:7),7,intsch,nface,xg,yg,4)
    qg(1:ifull,1:kl) = bb(1:ifull,1:kl,1)
    qlg(1:ifull,1:kl) = bb(1:ifull,1:kl,2)
    qfg(1:ifull,1:kl) = bb(1:ifull,1:kl,3)
    stratcloud(1:ifull,1:kl) = bb(1:ifull,1:kl,4)
    qrg(1:ifull,1:kl) = bb(1:ifull,1:kl,5)
    qsng(1:ifull,1:kl) = bb(1:ifull,1:kl,6)
    qgrg(1:ifull,1:kl) = bb(1:ifull,1:kl,7)
  else
    ! diagnostic microphysics 
    call ints(qg,1,intsch,nface,xg,yg,3)
  end if    ! ldr/=0
  if ( ntrac>0 ) then
    !if ( nmaxpr==1 .and. mydiag ) then
    !  write (6,"('xg#',9f8.2)") diagvals(xg(:,nlv))
    !  write (6,"('yg#',9f8.2)") diagvals(yg(:,nlv))
    !  write (6,"('nface#',9i8)") diagvals(nface(:,nlv))
    !  write (6,"('xlat#',9f8.2)") diagvals(tr(:,nlv,ngas+1))
    !  write (6,"('xlon#',9f8.2)") diagvals(tr(:,nlv,ngas+2))
    !  write (6,"('xpre#',9f8.2)") diagvals(tr(:,nlv,ngas+3))
    !end if
    call ints(tr(:,:,1:ntrac),ntrac,intsch,nface,xg,yg,5)
    !if ( nmaxpr==1 .and. mydiag ) then
    !  write (6,"('ylat#',9f8.2)") diagvals(tr(:,nlv,ngas+1))
    !  write (6,"('ylon#',9f8.2)") diagvals(tr(:,nlv,ngas+2))
    !  write (6,"('ypre#',9f8.2)") diagvals(tr(:,nlv,ngas+3))
    !endif
  endif  ! (ntrac>0)
  if ( nvmix==6 .or. nvmix==9 ) then
    bb(:,:,1) = tke(:,:)
    bb(:,:,2) = eps(:,:)
    call ints(bb(:,:,1:2),2,intsch,nface,xg,yg,4)
    tke(1:ifull,1:kl) = bb(1:ifull,1:kl,1)
    eps(1:ifull,1:kl) = bb(1:ifull,1:kl,2)
  endif                 ! nvmix==6 .or. nvmix==9
  if ( abs(iaero)>=2 .and. nhstest>=0 ) then
    call ints(xtg(:,:,1:naero),naero,intsch,nface,xg,yg,5)
  end if
end if     ! mspec==1


#ifdef GPU
!$acc end data
#endif


do k = 2,kl
  sdot(:,k) = sbar(:,k)
end do  

if ( mod(ktau, nmaxpr)==0 .and. mydiag ) then
  write(6,*) 'upglobal ktau,sdmx,nits,nvadh_pass ',ktau,sdmx(idjd),nits(idjd),nint(1./nvadh_inv_pass(idjd))
endif
if ( (diag.or.nmaxpr==1) .and. mydiag ) then
  write(6,*) 'in upglobal before vadv2'
  write (6,"('tx_b',9f8.2)")   tx(idjd,:)
  write (6,"('ux_b',9f8.3)")   ux(idjd,:)
  write (6,"('vx_b',9f8.3)")   vx(idjd,:)
  write (6,"('qg_b',9f8.3)")   qg(idjd,:)
endif

call vadvtvd(tx,ux,vx,nvadh_inv_pass,nits)


if ( (diag.or.nmaxpr==1) .and. mydiag ) then
  write(6,*) 'in upglobal after vadv2'
  write (6,"('tx_c',9f8.2)")   tx(idjd,:)
  write (6,"('ux_c',9f8.3)")   ux(idjd,:)
  write (6,"('vx_c',9f8.3)")   vx(idjd,:)
  write (6,"('qg_c',9f8.3)")   qg(idjd,:)
endif

do k = 1,kl
  ! adding later (after 2nd vadv) than for npex=1
  ux(1:ifull,k) = ux(1:ifull,k) + 0.5*dt*un(1:ifull,k) ! dyn contrib
  vx(1:ifull,k) = vx(1:ifull,k) + 0.5*dt*vn(1:ifull,k) ! dyn contrib
      
  ! second part of usual m=6 coriolis treatment (after 2nd vadv)
  ! incorporate coriolis terms (done here as for m=6 instead of in adjust5)
  tempry(1:ifull) = ux(1:ifull,k) + .5*dt*(1.+epsf)*f(1:ifull)*vx(1:ifull,k) ! Eq. 133
  vx(1:ifull,k)   = vx(1:ifull,k) - .5*dt*(1.+epsf)*f(1:ifull)*ux(1:ifull,k) ! Eq. 134
  ux(1:ifull,k)   = tempry(1:ifull)
  
  tx(1:ifull,k) = tx(1:ifull,k) + .5*dt*tn(1:ifull,k) 
end do

!     now interpolate ux,vx to the staggered grid
call staguv(ux,vx,ux,vx)


if ( diag ) then
  if ( mydiag ) then
    write(6,*) 'near end of upglobal staggered ux and vx:'
    write(6,*) 'un_u ',un(idjd,:)
    write(6,*) 'vn_u ',vn(idjd,:)
    write(6,*) 'tn_u ',tn(idjd,:)
    write (6,"('tx_u1',9f8.2/5x,9f8.2)") tx(idjd,:)
  end if
  call printa('ux  ',ux,ktau,nlv,ia,ib,ja,jb,0.,1.)
  call printa('vx  ',vx,ktau,nlv,ia,ib,ja,jb,0.,1.)
end if

if ( ntest>4 ) then
  ! diagnostic check for unstable layers
  do iq = 1,ifull
    theta(iq,1:kl) = tx(iq,1:kl)*sig(1:kl)**(-rdry/cp)
    do k = ntest,kl   ! e.g. 8,kl
      if ( theta(iq,k)<theta(iq,k-1) ) then  ! based on tx
        write(6,*)"unstable layer in upglobal for ktau,iq,k's,del ",ktau,iq,k-1,k,theta(iq,k-1)-theta(iq,k)
        write (6,"('theta',9f7.2/5x,9f7.2)") theta(iq,:)
        write (6,"('sdot',9f7.3/4x,9f7.3)")  sdot(iq,1:kl)
        numunstab = numunstab + 1
      end if
    end do
  end do
end if

if( (diag.or.nmaxpr==1) .and. mydiag ) then
  write(6,*) 'near end of upglobal for ktau= ',ktau
  write (6,"('tx_u2',9f8.2/5x,9f8.2)")  tx(idjd,:)
  write (6,"('qg_u',9f8.3/4x,9f8.3)")   qg(idjd,:)
  write (6,"('ql_u',9f8.3/4x,9f8.3)")   qlg(idjd,:)
  write (6,"('qf_u',9f8.3/4x,9f8.3)")   qfg(idjd,:)
end if 

call END_LOG(upglobal_end)
      
return
end subroutine upglobal
