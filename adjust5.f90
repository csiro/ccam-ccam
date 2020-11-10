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


subroutine adjust5

use aerosolldr             ! LDR prognostic aerosols
use arrays_m               ! Atmosphere dyamics prognostic arrays
use cc_mpi                 ! CC MPI routines
use cfrac_m                ! Cloud fraction
use const_phys             ! Physical constants
use diag_m                 ! Diagnostic routines
use dpsdt_m                ! Vertical velocity
use epst_m                 ! Off-centre terms
use helmsolve              ! Implicit solver for Helmholtz equation
use indices_m              ! Grid index arrays
use liqwpar_m              ! Cloud water mixing ratios
use map_m                  ! Grid map arrays
use morepbl_m              ! Additional boundary layer diagnostics
use newmpar_m              ! Grid parameters
use nharrs_m               ! Non-hydrostatic atmosphere arrays
use nlin_m                 ! Atmosphere non-linear dynamics
use parm_m                 ! Model configuration
use parmdyn_m              ! Dynamics parameters
use pbl_m                  ! Boundary layer arrays
use sigs_m                 ! Atmosphere sigma levels
use staguvmod              ! Reversible grid staggering
use tbar2d_m               ! Atmosphere dynamics reference temperature
use tracers_m              ! Tracer data
use vadv                   ! Vertical advection
use vecsuv_m               ! Map to cartesian coordinates
use vecs_m                 ! Eigenvectors for atmosphere dynamics
use vvel_m                 ! Additional vertical velocity
use work3sav_m             ! Water and tracer saved arrays
use xarrs_m                ! Saved dynamic arrays
use xyzinfo_m              ! Grid coordinate arrays

implicit none

include 'kuocom.h'         ! Convection parameters

integer k, l, iq
integer its
integer, save :: precon_in = -99999
real, dimension(:), allocatable, save :: zz, zzn, zze, zzw, zzs
real, dimension(:), allocatable, save :: pfact, alff, alf, alfe
real, dimension(:), allocatable, save :: alfn, alfu, alfv
real, dimension(ifull+iextra,kl) :: p, cc, dd, pe
real, dimension(ifull,kl) :: omgfnl, wrk1, wrk2, wrk3, wrk4
real, dimension(ifull,kl) :: helm, rhsl, omgf, e
real, dimension(ifull+iextra,kl) :: d
real, dimension(ifull) :: ps_sav, pslxint, pslsav
real, dimension(ifull) :: delps, bb
real, dimension(ifull) :: ptemp, ptempsav
real hdt, hdtds
real alph_p, alph_pm, delneg, delpos
real const_nh
real sumx
real, save :: dtsave = 0.

call START_LOG(adjust_begin)

hdt   = 0.5*dt
hdtds = hdt/ds
alph_p  = 0.
alph_pm = 0.
const_nh = 0.
cc(:,:) = 0.
dd(:,:) = 0.
p(:,:) = 0.
pe(:,:) = 0.

! time step can change during initialisation
if ( abs(dt-dtsave)>=1.e-20 ) then ! dt/=dtsave
  if ( .not.allocated(zz) ) then
    allocate( zz(ifull), zzn(ifull), zze(ifull), zzw(ifull) )
    allocate( zzs(ifull), pfact(ifull), alff(ifull+iextra) )
    allocate( alf(ifull+iextra), alfe(ifull+iextra) )
    allocate( alfn(ifull+iextra), alfu(ifull), alfv(ifull) )
    alff(:) = 0.
    alf(:) = 0.
    alfe(:) = 0.
    alfn(:) = 0.
  end if    
  call adjust_init(zz,zzn,zze,zzw,zzs,pfact,alff,alf,alfe,alfn,alfu,alfv)
  precon_in = precon
  precon = min(precon, 0)  ! 22/4/07
  if ( precon<-9999 ) then
    call mgsor_init
    call mgzz_init(zz,zzn,zze,zzw,zzs)
  end if
end if

if ( diag .or. nmaxpr==1 ) then
  if ( mydiag ) then
    write(6,*) 'entering adjust5'
    write (6,"('tx_a1',10f8.2)") tx(idjd,:)
    write (6,"('ux_stag',10f8.2)") ux(idjd,:)
    write (6,"('vx_stag',10f8.2)") vx(idjd,:)
    write (6,"('qg_a1',3p10f7.3)") qg(idjd,:)
    write (6,"('pslx_3p',3p9f8.4)") pslx(idjd,:)
  end if
  call printa('pslx',pslx,ktau,nlv,ia,ib,ja,jb,0.,1000.)
  call printa('tx  ',tx,ktau,nlv,ia,ib,ja,jb,200.,1.)
  call maxmin(alf,'a ',ktau,1.,1)
  call maxmin(alfe,'ae',ktau,1.,1)
  call maxmin(alfn,'an',ktau,1.,1)
end if

! recompute nonlinear sigma-dot contribution for updating tn, tx
! e contains intgrl{-pslx x dsigma} from 0 to sigma (dsig is -ve)
! vert. integ. nonlin part of mass weighted div into e
do k = 1,kl
  pslx(1:ifull,k) = pslx(1:ifull,k)*2./dt  ! i.e. [RHS of Eq. 115]*2/dt, i.e. M
end do    ! k  loop
e(1:ifull,kl) = dsig(kl)*pslx(1:ifull,kl)
wrk3(1:ifull,kl) = -sig(kl)*pslx(1:ifull,kl) ! integration following eig.f
do k = kl-1,1,-1
  e(1:ifull,k) = e(1:ifull,k+1) + dsig(k)*pslx(1:ifull,k)  ! i.e. -M_bar_(sig-.5)
  wrk3(1:ifull,k) = e(1:ifull,k+1) + (sigmh(k+1)-sig(k))*pslx(1:ifull,k) ! integration following eig.f
end do     ! k loop
pslxint(1:ifull) = -e(1:ifull,1)*dt/2. ! pslxint holds integrated pslx, Eq. 116

! full-level (1.+epsp)*omega/ps into omgfnl (nonlin. part only), Eq. 118
do k = 1,kl
  omgfnl(1:ifull,k) = -wrk3(1:ifull,k) - sig(k)*pslx(1:ifull,k)
end do     ! k loop
     
! redefine ux, vx
do k = 1,kl
  ! N.B. alfu is alf/(1+epsu)
  cc(1:ifull,k) = ux(1:ifull,k)/emu(1:ifull)*alfu(1:ifull)  ! Eq. 136
  dd(1:ifull,k) = vx(1:ifull,k)/emv(1:ifull)*alfv(1:ifull)  ! Eq. 137
end do


call boundsuv(cc,dd,stag=-9) ! only update isv and iwu


do k = 1,kl
  ! N.B. the omgfnl term on LHS of Eq. 121 not yet added in
  tx(1:ifull,k) = tx(1:ifull,k) + hdt*tbar2d(1:ifull)*omgfnl(1:ifull,k)*roncp/sig(k)  ! with correct epsp
end do ! k loop

! calculate hydrostatic heights from the tx array (add zs and pslxint terms after NHS terms)
p(1:ifull,1) = zs(1:ifull) + bet(1)*(tx(1:ifull,1)-280.) + rdry*tbar2d(1:ifull)*pslxint(1:ifull) ! Eq. 146
do k = 2,kl
  p(1:ifull,k) = p(1:ifull,k-1) + bet(k)*(tx(1:ifull,k)-280.) + betm(k)*(tx(1:ifull,k-1)-280.)
end do ! k loop

if ( nh/=0 .and. (ktau>=knh.or.lrestart) ) then
  ! add in departure values of p-related nh terms  & omgfnl terms    
  if ( abs(epsp)<=1. ) then
    const_nh = 2.*rdry/(dt*grav*grav*(1.+epsh))
  else
    const_nh=2.*rdry/(dt*grav*grav)
  end if  
  ! note that linear part of omega/ps for tau+1 is included in eig.f90
  ! wrk2 contains the tstar and non-linear part of omega/ps at tau+1
  select case(nh)
    case (8)
      if ( abs(epsp)<=1. ) then
        do k = 1,kl
          ! omgfnl already includes (1+epsp)  
          wrk2(:,k) = const_nh*(tbar(1)*omgfnl(:,k)/((1.-epsp)*sig(k))-h_nh(1:ifull,k))*t(1:ifull,k)**2/tbar(1)
        end do    
      else
        do k = 1,kl
          wrk2(:,k) = const_nh*(tbar(1)*omgfnl(:,k)/sig(k)-h_nh(1:ifull,k))*t(1:ifull,k)**2/tbar(1)
        end do          
      end if
    case default ! usually nh=5 or 6
      if ( abs(epsp)<=1. ) then
        do k = 1,kl
          ! omgfnl already includes (1+epsp)
          wrk2(:,k) = const_nh*(tbar(1)*omgfnl(1:ifull,k)/((1.-epsp)*sig(k))-h_nh(1:ifull,k))*tbar2d(1:ifull)
        end do
      else
        do k = 1,kl
          wrk2(:,k) = const_nh*(tbar(1)*omgfnl(1:ifull,k)/sig(k)-h_nh(1:ifull,k))*tbar2d(1:ifull)
        end do
      end if
  end select    
  wrk1(:,1) = bet(1)*wrk2(:,1)
  do k = 2,kl
    wrk1(:,k) = wrk1(:,k-1) + bet(k)*wrk2(:,k) + betm(k)*wrk2(:,k-1)
  end do   ! k loop
  if ( (diag.or.nmaxpr==1) .and. mydiag )then
    write(6,*) 'adjust5 omgfnl ',(omgfnl(idjd,k),k=1,kl)
    write(6,*) 'adjust5 h_nh ',(h_nh(idjd,k),k=1,kl)
    write(6,*) 'adjust5 pa ',(p(idjd,k),k=1,kl)
    write(6,*) 'adjust5 wrk1 ',(wrk1(idjd,k),k=1,kl)
  end if
  do k = 1,kl
    p(1:ifull,k) = p(1:ifull,k) + wrk1(1:ifull,k)  ! nh
  end do
end if     ! (nh/=0)

! form divergence of rhs (xu & xv) terms
do k = 1,kl
  ! d is xd in Eq. 157, divided by em**2/ds
  d(1:ifull,k) = cc(1:ifull,k) - cc(iwu,k) + dd(1:ifull,k) - dd(isv,k)
end do

! transform p & d to eigenvector space
do k = 1,kl
  pe(1:ifull,k) = einv(k,1)*p(1:ifull,1)
  rhsl(1:ifull,k) = einv(k,1)*d(1:ifull,1)
  do l = 2,kl
    pe(1:ifull,k) = pe(1:ifull,k) + einv(k,l)*p(1:ifull,l)     ! xp in eig space
    rhsl(1:ifull,k) = rhsl(1:ifull,k) + einv(k,l)*d(1:ifull,l) ! xd in eig space
  end do  ! l loop

  ! N.B.   pfact(iq)=4.*( ds/(dt*em(iq)) )**2    
  helm(1:ifull,k) = pfact(1:ifull)*tbar(1)/(bam(k)*(1.+epst(1:ifull))*tbar2d(1:ifull))
  rhsl(1:ifull,k) = rhsl(1:ifull,k)/hdtds - helm(1:ifull,k)*pe(1:ifull,k) ! Eq. 161 mult
end do    ! k loop


if ( precon<-9999 ) then
  ! Multi-grid
  call mghelm(zz,zzn,zze,zzw,zzs,helm,pe,rhsl)
else if ( precon<0 ) then
  ! SOR
  call helmsor(zz,zzn,zze,zzw,zzs,helm,pe,rhsl)
else
  ! Congujate gradient
  call bounds(pe)
  call helmsol(zz,zzn,zze,zzw,zzs,helm,pe,rhsl)
end if ! (precon<-9999) .. else ..


if ( diag .or. nmaxpr==1 ) then   !  only for last k of loop (i.e. 1)
  ! Some diagnostics requiring extra bounds calls have been removed
  if ( mydiag ) then
    write (6,"('tx_a2',10f8.2)") tx(idjd,:)
    write(6,*) 'omgfnl_a2 ',omgfnl(idjd,:)
    write(6,*) 'adjust5 cc ',cc(idjd,:)
    write(6,*) 'adjust5 dd ',dd(idjd,:)
    write(6,*) 'adjust5 d_in ',d(idjd,:)
    write(6,*) 'adjust5 p_in ',p(idjd,:)
    write(6,*) 'adjust5 pe ',pe(idjd,:)
    write(6,*) 'adjust5 pe_e ',pe(ie(idjd),:)
    write(6,*) 'adjust5 pe_w ',pe(iw(idjd),:)
    write(6,*) 'adjust5 rhsl ',rhsl(idjd,:)
    write(6,*) 'adjust5 helm ',helm(idjd,:)
  end if
  call printa('psnt',pslxint,ktau,0,ia,ib,ja,jb,0.,1000.)
  call printa('tx  ',tx,ktau,nlv,ia,ib,ja,jb,200.,1.)
  call printa('rhsl',rhsl,ktau,1,ia,ib,ja,jb,0.,0.)
  call printa('pe  ',pe,ktau,1,ia,ib,ja,jb,0.,0.)
end if

do k = 1,kl
  ! first p from pe
  p(1:ifull,k) = emat(k,1)*pe(1:ifull,1)
  do l = 2,kl              ! this is the remaining expensive loop
    p(1:ifull,k) = p(1:ifull,k) + emat(k,l)*pe(1:ifull,l)
  end do                 !  l loop
end do

call bounds(p,corner=.true.)

!      now u & v
if ((diag.or.nmaxpr==1).and.mydiag) then
  write(6,*) 'iq,k,fu,alfu,alfu*ux(iq,k) ',idjd,nlv,fu(idjd),alfu(idjd),alfu(idjd)*ux(idjd,nlv)
  write(6,*) 'alff & n e w s (in(iq)),alff(ine(iq)),alff(is(iq)),alff(ise(iq)),alfe(iq) ', &
      alff(in(idjd)),alff(ine(idjd)),alff(is(idjd)),alff(ise(idjd)),alfe(idjd)
  sumx = alf(ie(idjd))-alf(idjd)+.25*(alff(in(idjd))+alff(ine(idjd))-alff(is(idjd))-alff(ise(idjd)))-alfe(idjd)
  write(6,*) 'sum  ',sumx
  write(6,*) 'p & n e w s ne se ',p(idjd,nlv),p(in(idjd),nlv),p(ie(idjd),nlv),p(iw(idjd),nlv),p(is(idjd),nlv), &
                                  p(ine(idjd),nlv),p(ise(idjd),nlv)
  write(6,*) 'p & direct n s ',p(idjd,nlv),p(idjd+il,nlv),p(idjd-il,nlv)
end if

do k = 1,kl
  do iq = 1,ifull
    cc(iq,k) = alfu(iq)*ux(iq,k) - hdtds*emu(iq)*(                                        &
               alf(ie(iq))*p(ie(iq),k)-alf(iq)*p(iq,k)-.5*alfe(iq)*(p(iq,k)+p(ie(iq),k))  &
               +.25*(alff(in(iq))*p(in(iq),k) +alff(ine(iq))*p(ine(iq),k)                 &
               -alff(is(iq))*p(is(iq),k) -alff(ise(iq))*p(ise(iq),k)) ) ! Eq. 139
    dd(iq,k) = alfv(iq)*vx(iq,k) - hdtds*emv(iq)*(                                        &
               alf(in(iq))*p(in(iq),k)-alf(iq)*p(iq,k)-.5*alfn(iq)*(p(iq,k)+p(in(iq),k))  &
               -.25*(alff(ien(iq))*p(ien(iq),k) +alff(ie(iq))*p(ie(iq),k)                 &
               -alff(iwn(iq))*p(iwn(iq),k) -alff(iw(iq))*p(iw(iq),k)) ) ! Eq. 140
  end do  
end do     !  k loop 

call boundsuv(cc,dd,stag=-9) ! only update isv and iwu

! calculate linear part only of sigma-dot and omega/ps
do k = 1,kl
  d(1:ifull,k) = (cc(1:ifull,k)/emu(1:ifull)-cc(iwu,k)/emu(iwu)   &
                 +dd(1:ifull,k)/emv(1:ifull)-dd(isv,k)/emv(isv))  &
                 *em(1:ifull)**2/ds ! Eq. 101
end do     ! k  loop
if ( nmaxpr==1 ) then
  call maxmin(d,'dv',ktau,0.,kl)
  if ( nproc==1 ) write(6,*) 'cc,cc-,dd,dd-', cc(idjd,nlv)/emu(idjd),cc(iwu(idjd),nlv)/emu(iwu(idjd)), &
                                              dd(idjd,nlv)/emv(idjd),dd(isv(idjd),nlv)/emv(isv(idjd))      
end if   ! (nmaxpr==1)


! following is JLM sliding divergence damping (N.B. d dimensions with iextra)
if ( nh/=0 .and. ktau<=-knh .and. .not.lrestart ) then  ! e.g. knh=-10, divdamp=450.
  if ( nmaxpr==1 .and. mydiag ) then
    write(6,*) 'div damping'
    write(6,*) 'cc,cc-,dd,dd-', cc(idjd,nlv)/emu(idjd),cc(iwu(idjd),nlv)/emu(iwu(idjd)), &
                                            dd(idjd,nlv)/emv(idjd),dd(isv(idjd),nlv)/emv(isv(idjd))      
    write(6,*) 'd',d(idjd,:)
  end if
  call bounds(d)
  do k=1,kl
    cc(1:ifull,k)=cc(1:ifull,k)+(d(ie,k)-d(1:ifull,k))*divdamp*real(1-knh-ktau)*emu(1:ifull)*dt/(ds*real(abs(knh))) 
    dd(1:ifull,k)=dd(1:ifull,k)+(d(in,k)-d(1:ifull,k))*divdamp*real(1-knh-ktau)*emv(1:ifull)*dt/(ds*real(abs(knh)))
  end do
  call boundsuv(cc,dd,stag=-9) ! only update isv and iwu
  do k=1,kl
    d(1:ifull,k)=(cc(1:ifull,k)/emu(1:ifull)-cc(iwu,k)/emu(iwu)   &
                 +dd(1:ifull,k)/emv(1:ifull)-dd(isv,k)/emv(isv))  &
                 *em(1:ifull)**2/ds ! Eq. 101
  end do     ! k  loop
  if ( nmaxpr==1 .and. mydiag ) then
    write(6,*) 'after div damping'
    write(6,*) 'cc,cc-,dd,dd-', cc(idjd,nlv)/emu(idjd),cc(iwu(idjd),nlv)/emu(iwu(idjd)), &
                                            dd(idjd,nlv)/emv(idjd),dd(isv(idjd),nlv)/emv(isv(idjd))      
    write(6,*) 'd',d(idjd,:)
  end if
end if  ! (nh/=0.and.ktau<=-knh.and..not.lrestart)



! straightforward rev. cubic interp of u and v (i.e. nuv=10)
if ( nstag==0 ) then
  call staguv(u,v,wrk1,wrk2)
  its = ifull + iextra
  if ( nmaxpr==1 .and. nproc==1 ) then
    write(6,"('u_u0 ',10f8.2)") (u(iq,nlv),iq=idjd-3,idjd+3)
    write(6,"('v_u0 ',10f8.2)") (v(iq,nlv),iq=idjd-3*its,idjd+3*its,its)
    write(6,"('u_s0 ',10f8.2)") (wrk1(iq,nlv),iq=idjd-3,idjd+3)
    write(6,"('v_s0 ',10f8.2)") (wrk2(iq,nlv),iq=idjd-3*ifull,idjd+3*ifull,ifull)
    write(6,"('u_s1 ',10f8.2)") (cc(iq,nlv),iq=idjd-3,idjd+3)
    write(6,"('v_s1 ',10f8.2)") (dd(iq,nlv),iq=idjd-3*its,idjd+3*its,its)       
  endif
  wrk1(:,:) = cc(1:ifull,:) - wrk1(:,:) ! staggered increment
  wrk2(:,:) = dd(1:ifull,:) - wrk2(:,:) ! staggered increment
  call unstaguv(wrk1,wrk2,wrk3,wrk4) 
  u(1:ifull,:) = u(1:ifull,:) + wrk3(:,:)
  v(1:ifull,:) = v(1:ifull,:) + wrk4(:,:)
  if ( nmaxpr==1 .and. nproc==1 ) then
    write(6,"('u_u1 ',10f8.2)") (u(iq,nlv),iq=idjd-3,idjd+3)
    write(6,"('v_u1 ',10f8.2)") (v(iq,nlv),iq=idjd-3*its,idjd+3*its,its)
 endif
else
  call unstaguv(cc,dd,u,v) ! usual
endif

! vert. integ. div into e
wrk2(1:ifull,kl) = -dsig(kl)*d(1:ifull,kl)
wrk3(1:ifull,kl) = sig(kl)*d(1:ifull,kl)
do k = kl-1,1,-1
  wrk2(1:ifull,k) = wrk2(1:ifull,k+1) - dsig(k)*d(1:ifull,k)
  wrk3(1:ifull,k) = wrk2(1:ifull,k+1) - (sigmh(k+1)-sig(k))*d(1:ifull,k)
end do     ! k  loop

! full-level omega/ps into omgf (linear part only)
do k = 1,kl
  omgf(1:ifull,k) = -wrk3(1:ifull,k) ! in Eq. 110
end do     ! k  loop
pslsav(1:ifull) = psl(1:ifull) ! saved for gas fixers below, and diags
ps_sav(1:ifull) = ps(1:ifull)  ! saved for gas fixers below, and diags
psl(1:ifull)    = pslxint(1:ifull) - hdt*wrk2(1:ifull,1)*(1.+epst(1:ifull))  ! Eq. 116

if ( mod(ktau, nmaxpr)==0 ) vx(1:ifull,:) = sdot(1:ifull,1:kl) ! for accln

do k = 1,kl
  ! save [D + dsigdot/dsig] in pslx for next use in nonlin
  pslx(1:ifull,k) = (pslx(1:ifull,k)-psl(1:ifull)*2./dt)/(1.+epst(1:ifull)) ! from Eq. 115
end do   !  k loop
do k = kl,2,-1
  ! calculate latest sdot (at level k-.5)
  sdot(1:ifull,k) = sdot(1:ifull,k+1) - dsig(k)*(pslx(1:ifull,k)-d(1:ifull,k))
end do   !  k loop
      
do k = 2,kl
  ! and convert sdot (at level k-.5) to units of grid-steps/timestep
  ! dtin is used to be ready for next full timestep
  sdot(1:ifull,k) = sdot(1:ifull,k)*dtin/(sig(k)-sig(k-1))
end do    ! k  loop

if ( mod(ktau,nmaxpr)==0 ) then
  vx(1:ifull,:) = sdot(1:ifull,1:kl) - vx(1:ifull,:)
  ! convert to approx m/s/s
  do k = 2,kl
    vx(1:ifull,k) = vx(1:ifull,k)*rdry*(sig(k-1)-sig(k))*.5*(t(1:ifull,k)+t(1:ifull,k-1))/(sigmh(k)*grav*dtin*dt)
  end do
  call maxmin(vx,'ac',ktau,100.,kl)  ! max min of accln * 100
end if
if ( (diag.or.nmaxpr==1) .and. mydiag ) then
  write(6,"('diva5p ',5p10f8.2)") d(idjd,:)
  write(6,"('omgf_l*dt',10f8.4)") omgf(idjd,:)*dt
  write(6,"('omgfnl*dt',10f8.4)") omgfnl(idjd,:)*dt
  write(6,"('u_a2 ',10f8.2)") u(idjd,:)
  write(6,"('v_a2 ',10f8.2)") v(idjd,:)
  write(6,*) 'ps,psl ',ps(idjd),psl(idjd)
  write(6,"('pslx_3p',3p9f8.4)") pslx(idjd,:)
  write(6,"('sdot_a2',10f8.3)") sdot(idjd,1:kl)
end if

do k = 1,kl
  ! save full omega/ps in dpsldt for use in nonlin next time step (& outfile)
  ! N.B. omgfnl part already incorp. into tx above
  dpsldt(1:ifull,k) = omgfnl(1:ifull,k)/(1.+epst(1:ifull)) + omgf(1:ifull,k)
  t(1:ifull,k) = tx(1:ifull,k) + hdt*(1.+epst(1:ifull))*tbar2d(1:ifull)*omgf(1:ifull,k)*roncp/sig(k) ! Eq 121 F26
end do     ! k  loop

if ( (diag.or.nmaxpr==1) .and. mydiag ) then
  write(6,"('omgf_a2',10f8.3)") ps(idjd)*dpsldt(idjd,1:kl)
end if

if ( nh/=0 .and. (ktau>knh.or.lrestart) ) then
   
  ! original method for estimating phi_nh
  do k = 1,kl
    phi(:,k) = p(1:ifull,k) - rdry*tbar2d(:)*psl(1:ifull)
  end do
  ! extract NHS component
  bb(:) = zs(1:ifull) + bet(1)*(t(1:ifull,1)-280.)
  phi_nh(:,1) = phi(:,1) - bb(:)
  do k = 2,kl
    bb(:) = bb(:) + bet(k)*(t(1:ifull,k)-280.) + betm(k)*(t(1:ifull,k-1)-280.)
    phi_nh(:,k) = phi(:,k) - bb(:)
  end do
    
  !! alternative method for estimating phi_nh
  !do k = 1,kl
  !  wrk2(1:ifull,k) = const_nh*tbar2d(1:ifull)*(tbar(1)*dpsldt(1:ifull,k)/sig(k)-h_nh(1:ifull,k))
  !end do
  !phi_nh(1:ifull,1) = bet(1)*wrk2(1:ifull,1)
  !do k = 2,kl
  !  phi_nh(1:ifull,k) = phi_nh(1:ifull,k-1) + bet(k)*wrk2(1:ifull,k) + betm(k)*wrk2(1:ifull,k-1)
  !end do   ! k loop 
  
  ! update phi for use in next time step
  phi(1:ifull,1) = zs(1:ifull) + bet(1)*t(1:ifull,1)
  do k = 2,kl
    phi(1:ifull,k) = phi(1:ifull,k-1) + bet(k)*t(1:ifull,k) + betm(k)*t(1:ifull,k-1)
  end do
  phi(1:ifull,1:kl) = phi(1:ifull,1:kl) + phi_nh(1:ifull,1:kl)

  if ( nmaxpr==1 .and. mydiag ) then
    write(6,*) 'phi_nh ',(phi_nh(idjd,k),k=1,kl)
  end if

end if  ! (nh/=0.and.(ktau>knh.or.lrestart))


call START_LOG(mfix_begin)

if ( mfix==-1 ) then   ! perform conservation fix on psl
  ! delpos is the sum of all positive changes over globe
  ! delneg is the sum of all negative changes over globe
  ! alph_p is chosen to satisfy alph_p*delpos + delneg/alph_p = 0
  ! _l means local to this processor        
  delps(1:ifull) = psl(1:ifull) - pslsav(1:ifull)
  call ccglobal_posneg(delps,delpos,delneg)
  alph_p = sqrt( -delneg/max(1.e-30,delpos))
  alph_pm=1./max(1.e-30,alph_p)
  psl(1:ifull) = pslsav(1:ifull) + alph_p*max(0.,delps(1:ifull)) + alph_pm*min(0.,delps(1:ifull))
end if    !  (mfix==-1)

if ( mfix/=3 ) then
  ps(1:ifull) = 1.e5*exp(psl(1:ifull))     
end if
if ( mfix==1 .or. mfix==2 ) then   ! perform conservation fix on ps
  ! fix is on ps (not psl) from 24/1/06      
  ! delpos is the sum of all positive changes over globe
  ! delneg is the sum of all negative changes over globe
  ! alph_p is chosen to satisfy alph_p*delpos + delneg/alph_p = 0
  ! _l means local to this processor
  bb(1:ifull)    = ps(1:ifull)   
  delps(1:ifull) = ps(1:ifull) - ps_sav(1:ifull)
  call ccglobal_posneg(delps,delpos,delneg)
  if ( mfix==1 ) then
    alph_p  = sqrt(-delneg/max(1.e-30, delpos))
    alph_pm = 1./max(1.e-30, alph_p)
  else if ( mfix==2 ) then
    if ( delpos>-delneg ) then
      alph_p  = 1.
      alph_pm = -delpos/delneg
    else
      alph_p  = -delneg/delpos
      alph_pm = 1.
    end if
  end if                  ! (mfix==2)
  ps(1:ifull) = ps_sav(1:ifull) + alph_p*max(0., delps(1:ifull)) + alph_pm*min(0., delps(1:ifull))
  ! psl(1:ifull)=log(1.e-5*ps(1:ifull))
  ! following is cheaper and maintains full precision of psl
  psl(1:ifull) = psl(1:ifull) + (ps(1:ifull)/bb(1:ifull)-1.)     
end if !  (mfix==1.or.mfix==2)
      
! newer scheme - preferred
if ( mfix==3 ) then   ! perform conservation fix on ps (best for 32-bit)
  ! fix is on ps (not psl) from 24/1/06      
  ! delpos is the sum of all positive changes over globe
  ! delneg is the sum of all negative changes over globe
  ! alph_p is chosen to satisfy alph_p*delpos + delneg/alph_p = 0
  ! _l means local to this processor     
  delps(1:ifull) = psl(1:ifull) - pslsav(1:ifull)
  delps(1:ifull) = ps_sav(1:ifull)*delps(1:ifull)*(1.+0.5*delps(1:ifull))         
  call ccglobal_posneg(delps,delpos,delneg)
  alph_p  = sqrt( -delneg/max(1.e-30,delpos) )
  alph_pm = 1./max(1.e-30,alph_p)
  delps(1:ifull) = alph_p*max(0.,delps(1:ifull)) + alph_pm*min(0.,delps(1:ifull))
  delps(1:ifull) = delps(1:ifull)/ps_sav(1:ifull)
  psl(1:ifull)   = pslsav(1:ifull) + delps(1:ifull)*(1.-0.5*delps(1:ifull))
  ps(1:ifull)    = 1.e5*exp(psl(1:ifull))
end if    ! (mfix==3)
      
! following dpsdt diagnostic is in hPa/day
if ( epsp>1. .and. epsp<2. ) then
  dpsdtbb(1:ifull)     = dpsdtb(1:ifull)    
  dpsdtb(1:ifull)      = dpsdt(1:ifull)    
end if
dpsdt(1:ifull) = (ps(1:ifull)-ps_sav(1:ifull))*24.*3600./(100.*dt)

!------------------------------------------------------------------------
! Air temperature conservation
if ( mfix_t/=0 ) then
  ptemp(1:ifull) = ps(1:ifull)**.286
  ptempsav(1:ifull) = ps_sav(1:ifull)**.286
  do k = 1,kl
    t(1:ifull,k) = t(1:ifull,k)/ptemp(1:ifull)
    tsav(1:ifull,k) = tsav(1:ifull,k)/ptempsav(1:ifull)
  end do  
  call massfix(mfix_t,1,t,tsav,ps,ps_sav,.false.)
  do k = 1,kl
    t(1:ifull,k) = t(1:ifull,k)*ptemp(1:ifull)
  end do
end if

!------------------------------------------------------------------------
! Cloud water conservation
if ( mfix_qg/=0 .and. mspec==1 .and. ldr/=0 ) then
  do k = 1,kl
    stratcloud(1:ifull,k) = min( max( stratcloud(1:ifull,k), 0. ), 1. )
    qg(1:ifull,k) = max( qg(1:ifull,k), qgmin-qfg(1:ifull,k)-qlg(1:ifull,k), 0. )
    qfg(1:ifull,k) = max( qfg(1:ifull,k), 0. )
    qlg(1:ifull,k) = max( qlg(1:ifull,k), 0. )
  end do
  call massfix(mfix_qg,1,qg,qgsav,ps,ps_sav,.false.) 
  call massfix(mfix_qg,1,qfg,qfgsav,ps,ps_sav,.true.) 
  call massfix(mfix_qg,1,qlg,qlgsav,ps,ps_sav,.true.) 
  do k = 1,kl
    qlg(1:ifull,k) = max( qlg(1:ifull,k), 0. )
    qfg(1:ifull,k) = max( qfg(1:ifull,k), 0. )    
    qg(1:ifull,k)  = max( qg(1:ifull,k), qgmin-qfg(1:ifull,k)-qlg(1:ifull,k), 0. )
  end do
else if ( mfix_qg/=0 .and. mspec==1 ) then
  qg(1:ifull,1:kl) = max( qg(1:ifull,1:kl), qgmin )
  call massfix(mfix_qg,1,qg,qgsav,ps,ps_sav,.false.)
  qg(1:ifull,1:kl) = max( qg(1:ifull,1:kl), qgmin )
end if !  (mfix_qg/=0.and.mspec==1.and.ldr/=0) ..else..

!------------------------------------------------------------------------
! Tracer conservation
if ( mfix_tr/=0 .and. mspec==1 .and. ngas>0 ) then
  call massfix(mfix_tr,ngas,tr,trsav,ps,ps_sav,.true.)
end if !  (mfix_tr/=0.and.mspec==1.and.ngas>0)

!--------------------------------------------------------------
! Aerosol conservation
if ( mfix_aero/=0 .and. mspec==1 .and. abs(iaero)>=2 ) then
  xtg(1:ifull,1:kl,1:naero) = max( xtg(1:ifull,1:kl,1:naero), 0. )
  call massfix(mfix_aero,naero,xtg,xtgsav,ps,ps_sav,.true.)
  xtg(1:ifull,1:kl,1:naero) = max( xtg(1:ifull,1:kl,1:naero), 0. )
end if ! (mfix_aero/=0.and.mspec==1.and.abs(iaero)>=2)
!--------------------------------------------------------------
      
call END_LOG(mfix_end)


if ( (diag.or.nmaxpr==1) .and. mydiag ) then
  write(6,*) 'at end of adjust5 for ktau= ',ktau
  write(6,*) 'ps_sav,ps ',ps_sav(idjd),ps(idjd)
  write(6,"('dpsdt# mb/d ',9f8.1)") diagvals(dpsdt) 
  write(6,"('qg_a4 ',3p10f8.3)") qg(idjd,:)
  write(6,"('qgs',3p10f8.3)") (qgsav(idjd,k)/ps_sav(idjd),k=1,kl)
  write(6,"('qf_a4',3p10f8.3)") qfg(idjd,:)
  write(6,"('ql_a4',3p10f8.3)") qlg(idjd,:)
end if
if ( diag ) then
  call bounds(psl)
  if ( mydiag ) then
    iq=idjd
    write(6,*) 'adjust5 d ',d(iq,:)
    write(6,*) 'adjust5 wrk2 ',idjd, wrk2(idjd,:)
    write(6,"('adjust5 psl_3p & n e w s',3p9f8.4)") psl(idjd),psl(in(iq)),psl(ie(iq)),psl(iw(iq)),psl(is(iq))
  end if
  call printa('psl ',psl,ktau,0,ia,ib,ja,jb,0.,1000.)
  if ( mydiag ) write(6,*) 'adjust5 t ',t(idjd,:)
  call printa('t   ',t,ktau,nlv,ia,ib,ja,jb,200.,1.)
  if ( mydiag ) write(6,*) 'adjust5 u ',u(idjd,:)
  call printa('u   ',u,ktau,nlv,ia,ib,ja,jb,0.,1.)
  if ( mydiag ) write(6,*) 'adjust5 v ',v(idjd,:)
  call printa('v   ',v,ktau,nlv,ia,ib,ja,jb,0.,1.)
  if ( mydiag ) write(6,*) 'ps diagnostic from end of adjust:'
  bb(1:ifull)=ps(1:ifull)-ps_sav(1:ifull)
  call printa('dps ',bb,ktau,0,ia,ib,ja,jb,0.,.01)
  call printa('ps  ',ps,ktau,0,ia,ib,ja,jb,1.e5,.01)
  if ( sig(nlv)<.3 ) then
    call printa('qg  ',qg,ktau,nlv,ia,ib,ja,jb,0.,1.e6)
  else
    call printa('qg  ',qg,ktau,nlv,ia,ib,ja,jb,0.,1.e3)
  end if
end if

dtsave = dt
      
call END_LOG(adjust_end)

end subroutine adjust5

subroutine adjust_init(zz,zzn,zze,zzw,zzs,pfact,alff,alf,alfe,alfn,alfu,alfv)

use cc_mpi
use indices_m
use map_m
use newmpar_m
use parm_m
use parmdyn_m

implicit none

integer :: iq, n
real, dimension(ifull), intent(out) :: zz, zzn, zze, zzw, zzs, pfact
real, dimension(ifull), intent(out) :: alfu, alfv
real, dimension(ifull+iextra), intent(out) :: alff, alf, alfe, alfn
real :: hdt

hdt = 0.5*dt
alf(1:ifull) = (1.+epsu)/(1.+(hdt*(1.+epsf)*f(1:ifull))**2) ! Eq. 138
alff(1:ifull) = alf(1:ifull)*f(1:ifull)*hdt*(1.+epsf)       ! now includes hdt in alff
alfu(1:ifull) = 1./(1.+(hdt*(1.+epsf)*fu(1:ifull))**2)      ! i.e. alf/(1+epsu)
alfv(1:ifull) = 1./(1.+(hdt*(1.+epsf)*fv(1:ifull))**2)      ! i.e. alf/(1+epsu)
! These really only need to be recomputed when time step changes.
call bounds(alf)
call bounds(alff,corner=.true.)
pfact(1:ifull) = 4.*( ds/(dt*em(1:ifull)) )**2    
alfe(1:ifull) = alf(ie)-alf(1:ifull)+0.25*(alff(ine)+alff(in)-alff(ise)-alff(is)) ! alf3 Eq. 141 times ds
alfn(1:ifull) = alf(in)-alf(1:ifull)-0.25*(alff(ien)+alff(ie)-alff(iwn)-alff(iw)) ! alf4 Eq. 142 times ds
call boundsuv(alfe,alfn)

! see Eq 158/159
! need care with vector quantities on w (odd) & s (even) panel boundaries
zz(1:ifull) = 0.5*(alfe(iwu)-alfe(1:ifull)+alfn(isv)-alfn(1:ifull)) - 4.*alf(1:ifull) ! i,j   coeff
zzn(1:ifull) = alf(in) - 0.5*alfn(1:ifull)                                            ! i,j+1 coeff
zzw(1:ifull) = alf(iw) + 0.5*alfe(iwu)                                                ! i-1,j coeff
zze(1:ifull) = alf(ie) - 0.5*alfe(1:ifull)                                            ! i+1,j coeff
zzs(1:ifull) = alf(is) + 0.5*alfn(isv)                                                ! i,j-1 coeff
! N.B. there are some special z values at the 8 vertices
if ( edge_s .and. edge_w ) then
  do n = 1,npan            ! 1,6    
    iq = indp(1,1,n)
    zzs(iq) = zzs(iq) + 0.25*alff(is(iq)) ! i,j-1 coeff
    zzw(iq) = zzw(iq) - 0.25*alff(iw(iq)) ! i-1,j coeff
  end do   ! n loop
end if
if( edge_n .and. edge_e ) then
  do n = 1,npan            ! 1,6        
    iq = indp(ipan,jpan,n)
    zzn(iq) = zzn(iq) + 0.25*alff(in(iq)) ! i,j+1 coeff
    zze(iq) = zze(iq) - 0.25*alff(ie(iq)) ! i+1,j coeff
  end do   ! n loop
end if
if ( edge_s .and. edge_e ) then
  do n = 1,npan            ! 1,6          
    iq = indp(ipan,1,n)
    zzs(iq) = zzs(iq) - 0.25*alff(is(iq)) ! i,j-1 coeff
    zze(iq) = zze(iq) + 0.25*alff(ie(iq)) ! i+1,j coeff
  end do   ! n loop
end if
if ( edge_n .and. edge_w ) then
  do n = 1,npan            ! 1,6    
    iq = indp(1,jpan,n)
    zzn(iq) = zzn(iq) - 0.25*alff(in(iq)) ! i,j+1 coeff
    zzw(iq) = zzw(iq) + 0.25*alff(iw(iq)) ! i-1,j coeff
  end do   ! n loop
end if

return
end subroutine adjust_init

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Mass fixer subroutine
    
subroutine massfix(mfix,ntr,s,ssav_in,ps,pssav,llim)
      
use cc_mpi
use newmpar_m
use sigs_m
      
implicit none
      
integer, intent(in) :: mfix, ntr
integer k, i
real, dimension(ifull+iextra,kl,ntr), intent(inout) :: s
real, dimension(ifull,kl,ntr), intent(in) :: ssav_in
real, dimension(ifull), intent(in) :: ps, pssav
real, dimension(ntr) :: delpos, delneg, ratio, alph_g
logical, intent(in) :: llim

ratio = 0. ! for cray compiler

select case(mfix)
  !case(4)  
  !  do i = 1,ntr
  !    do k = 1,kl
  !      wrk1(1:ifull,k,i) = s(1:ifull,k,i)*ps(1:ifull) - ssav_in(1:ifull,k,i)*pssav(1:ifull)
  !    end do
  !  end do
  !  call ccglobal_posneg(wrk1(:,:,1:ntr),delpos3,delneg3,dsig)
  !  do i = 1,ntr
  !    do k = 1,kl
  !      wrk1(1:ifull,k,i) = (s(1:ifull,k,i)-ssav_in(1:ifull,k,i))*ps(1:ifull)
  !    end do
  !  end do
  !  call ccglobal_posneg(wrk1(:,:,1:ntr),delpos,delneg,dsig)
  !  alph_g(1:ntr) = -(delpos3(1:ntr)+delneg3(1:ntr))/(max(delpos(1:ntr),1.e-30)-delneg(1:ntr))
  !  do i = 1,ntr
  !    do k = 1,kl
  !      s(1:ifull,k,i) = s(1:ifull,k,i)                                                           &
  !          + alph_g(i)*(max(0.,wrk1(1:ifull,k,i))-min(0.,wrk1(1:ifull,k,i)))/ps(1:ifull)
  !    end do
  !  end do
  
  case(3) ! newer scheme - preferred
    do i = 1,ntr
      do k = 1,kl
        s(1:ifull,k,i) = s(1:ifull,k,i) - ssav_in(1:ifull,k,i)*pssav(1:ifull)/ps(1:ifull) 
        s(1:ifull,k,i) = pssav(1:ifull)*s(1:ifull,k,i)*(1.+0.5*s(1:ifull,k,i))
      end do   ! k loop
    end do
    call ccglobal_posneg(s,delpos,delneg,dsig)
    if ( llim ) then
      ratio(1:ntr) = -delneg(1:ntr)/max(delpos(1:ntr), 1.e-30)
    else
      ratio(1:ntr) = -delneg(1:ntr)/delpos(1:ntr)
    end if
    alph_g(1:ntr) = min(ratio(1:ntr), sqrt(ratio(1:ntr)))
    do i = 1,ntr
      do k = 1,kl
        s(1:ifull,k,i) = alph_g(i)*max(0., s(1:ifull,k,i))+min(0., s(1:ifull,k,i))/max(1., alph_g(i))
        s(1:ifull,k,i) = s(1:ifull,k,i)/pssav(1:ifull)
        s(1:ifull,k,i) = ssav_in(1:ifull,k,i)*pssav(1:ifull)/ps(1:ifull) + s(1:ifull,k,i)*(1.-0.5*s(1:ifull,k,i))
      end do    ! k  loop
    end do

  case(2)
    do i = 1,ntr
      do k = 1,kl
        s(1:ifull,k,i) = s(1:ifull,k,i) - ssav_in(1:ifull,k,i)*pssav(1:ifull)/ps(1:ifull) 
      end do   ! k loop
    end do
    call ccglobal_posneg(s,delpos,delneg,dsig)
    if ( llim ) then
      ratio(1:ntr) = -delneg(1:ntr)/max(delpos(1:ntr), 1.e-30)
      alph_g(1:ntr) = max(sqrt(ratio(1:ntr)), 1.e-30)
    else
      ratio(1:ntr) = -delneg(1:ntr)/delpos(1:ntr)
      alph_g(1:ntr) = sqrt(ratio(1:ntr))
    end if
    do i = 1,ntr
      do k = 1,kl
        s(1:ifull,k,i) = ssav_in(1:ifull,k,i)*pssav(1:ifull)/ps(1:ifull) &
            + alph_g(i)*max(0., s(1:ifull,k,i))+min(0., s(1:ifull,k,i))/max(1., alph_g(i))
      end do    ! k  loop
    end do

  case(1) ! original scheme
    do i = 1,ntr
      do k = 1,kl
        s(1:ifull,k,i) = s(1:ifull,k,i) - ssav_in(1:ifull,k,i)*pssav(1:ifull)/ps(1:ifull) 
      end do   ! k loop
    end do
    call ccglobal_posneg(s,delpos,delneg,dsig)
    if ( llim ) then
      ratio(1:ntr) = -delneg(1:ntr)/max(delpos(1:ntr), 1.e-30)
    else
      ratio(1:ntr) = -delneg(1:ntr)/delpos(1:ntr)
    end if
    alph_g(1:ntr) = min(ratio(1:ntr), sqrt(ratio(1:ntr)))
    do i = 1,ntr
      do k = 1,kl
        s(1:ifull,k,i) = ssav_in(1:ifull,k,i)*pssav(1:ifull)/ps(1:ifull) &
            + alph_g(i)*max(0., s(1:ifull,k,i))+min(0., s(1:ifull,k,i))/max(1., alph_g(i))
      end do    ! k  loop
    end do
   
  case default    
    write(6,*) "ERROR: Unknow mfix option ",mfix
    call ccmpi_abort(-1)
  
end select
  
return
end subroutine massfix
    
