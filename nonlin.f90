! Conformal Cubic Atmospheric Model
    
! Copyright 2015-2019 Commonwealth Scientific Industrial Research Organisation (CSIRO)
    
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

subroutine nonlin

use aerosolldr          ! LDR prognostic aerosols
use arrays_m            ! Atmosphere dyamics prognostic arrays
use cc_mpi              ! CC MPI routines
use const_phys          ! Physical constants
use diag_m              ! Diagnostic routines
use epst_m              ! Off-centre terms
use hs_phys_m           ! Held & Suarez
use indices_m           ! Grid index arrays
use latlong_m           ! Lat/lon coordinates
use liqwpar_m           ! Cloud water mixing ratios
use map_m               ! Grid map arrays
use morepbl_m           ! Additional boundary layer diagnostics
use newmpar_m           ! Grid parameters
use nharrs_m            ! Non-hydrostatic atmosphere arrays
use nlin_m              ! Atmosphere non-linear dynamics
use parm_m              ! Model configuration
use parmdyn_m           ! Dynamics parameters
use savuvt_m            ! Saved dynamic arrays
use sigs_m              ! Atmosphere sigma levels
use staguvmod           ! Reversible grid staggering
use tbar2d_m            ! Atmosphere dynamics reference temperature
use tracers_m           ! Tracer data
use unn_m               ! Saved dynamic arrays
use vadv                ! Vertical advection
use vecsuv_m            ! Map to cartesian coordinates
use vvel_m              ! Additional vertical velocity
use work3sav_m          ! Water and tracer saved arrays
use xarrs_m             ! Saved dynamic arrays
use xyzinfo_m           ! Grid coordinate arrays

implicit none

include 'kuocom.h'      ! Convection parameters

integer iq, k, n
integer, save :: num = 0
real invconst_nh, contv
real, dimension(ifull,kl) :: aa, bb
real, dimension(ifull+iextra,kl) :: p, tv, phiv
real, dimension(ifull+iextra,kl+1) :: duma
real, dimension(ifull) :: t_nh, spmax2, termlin
real, allocatable, save, dimension(:) :: epstsav
      
call START_LOG(nonlin_begin)

if ( epsp<-2. ) then
  if ( .not.allocated(epstsav) ) then
    allocate(epstsav(ifull))
    epstsav(:) = epst(:)
  end if
  spmax2(1:ifull) = max(u(1:ifull,3*kl/4)**2+v(1:ifull,3*kl/4)**2, u(1:ifull,kl)**2+v(1:ifull,kl)**2)
  where ( spmax2>(.8*ds/(em(1:ifull)*dt))**2 )
    ! setting epst for Courant number > .8        
    epst(1:ifull) = epstsav(1:ifull)
  elsewhere
    epst(1:ifull) = 0.
  end where
end if  ! (epsp<-2.)

! *** save t to avoid overshoot
do k = 1,kl
  tsav(1:ifull,k) = t(1:ifull,k)  
end do

! *** following qgsav should be before first vadv call
do k = 1,kl
  qgsav(1:ifull,k) = qg(1:ifull,k)      ! for qg  conservation in adjust5
end do  
if ( ldr/=0 ) then
  do k = 1,kl  
    qfgsav(1:ifull,k)  = qfg(1:ifull,k)
    qlgsav(1:ifull,k)  = qlg(1:ifull,k)
  end do  
end if   ! (ldr/=0)
      
if ( abs(iaero)>=2 ) then
  do n = 1,naero 
    do k = 1,kl 
      xtgsav(1:ifull,k,n) = xtg(1:ifull,k,n)
    end do  
  end do  
end if   ! abs(iaero)==2

if ( ngas>0 ) then
  do n = 1,naero
    do k = 1,kl  
      trsav(1:ifull,k,n) = tr(1:ifull,k,n) ! for tr conservation in adjust5
    end do
  end do  
end if   ! (ngas>=1)
 
if ( diag .or. nmaxpr==1 ) then
  call bounds(ps)
  if ( mydiag ) then
    write(6,*) "qgsav ",qgsav(idjd,nlv)
    write(6,*) 'at beginning of nonlin'
    write(6,*) 'roncp ',roncp
    write (6,"('tn0*dt',9f8.3/6x,9f8.3)") tn(idjd,:)*dt
    write (6,"('un0*dt',9f8.3/6x,9f8.3)") un(idjd,:)*dt
    write (6,"('vn0*dt',9f8.3/6x,9f8.3)") vn(idjd,:)*dt
    write (6,"('tbar',9f8.3/4x,9f8.3)") tbar(:)
    write (6,"('sig ',9f8.5/4x,9f8.5)") sig(:)
    write (6,"('rata',9f8.5/4x,9f8.5)") rata(:)
    write (6,"('ratb',9f8.5/4x,9f8.5)") ratb(:)
    write (6,"('em & news',5f10.4)") em(idjd),em(in(idjd)),em(ie(idjd)),em(iw(idjd)),em(is(idjd))
    write (6,"('emu,emu_w,emv,emv_s',4f10.4)") emu(idjd),emu(iwu(idjd)),emv(idjd),emv(isv(idjd))
    write (6,"('psl & news ',5f9.5)") psl(idjd),psl(in(idjd)),psl(ie(idjd)),psl(iw(idjd)),psl(is(idjd))
    write (6,"('ps  & news ',5f9.3)") ps(idjd), ps(in(idjd)),ps(ie(idjd)),ps(iw(idjd)),ps(is(idjd))
  endif
  call printa('u   ',u,ktau,nlv,ia,ib,ja,jb,0.,1.)
  call printa('v   ',v,ktau,nlv,ia,ib,ja,jb,0.,1.)
  call printa('t   ',t,ktau,nlv,ia,ib,ja,jb,200.,1.)

  if ( mydiag ) then
    write(6,*) 'in nonlin before possible vertical advection',ktau
    write (6,"('epst#  ',9f8.2)") diagvals(epst) 
    write (6,"('sdot#  ',9f8.3)") diagvals(sdot(:,nlv)) 
    write (6,"('sdotn  ',9f8.3/7x,9f8.3)") sdot(idjd,1:kl)
    write (6,"('omgf#  ',9f8.3)") ps(idjd)*dpsldt(idjd,nlv)
    write (6,"('omgfn  ',9f8.3/7x,9f8.3)") ps(idjd)*dpsldt(idjd,:)
    write (6,"('t   ',9f8.3/4x,9f8.3)")     t(idjd,:)
    write (6,"('u   ',9f8.3/4x,9f8.3)")     u(idjd,:)
    write (6,"('v   ',9f8.3/4x,9f8.3)")     v(idjd,:)
    write (6,"('qg  ',9f8.3/4x,9f8.3)")     qg(idjd,:)
  end if
end if

if ( nhstest==1 ) then ! Held and Suarez test case
!$omp parallel
  call hs_phys
!$omp end parallel
endif

if ( diag ) then
  call printa('sdot',sdot,ktau,nlv+1,ia,ib,ja,jb,0.,10.)
  call printa('omgf',dpsldt,ktau,nlv,ia,ib,ja,jb,0.,1.e5)
  aa(1:ifull,1) = rata(nlv)*sdot(1:ifull,nlv+1) + ratb(nlv)*sdot(1:ifull,nlv)
  if ( mydiag ) write(6,*) 'k,aa,emu,emv',nlv,aa(idjd,1),emu(idjd),emv(idjd)
  call printa('sgdf',aa(:,1),ktau,nlv,ia,ib,ja,jb,0.,10.)
end if   ! (diag)

! extra qfg & qlg terms included in tv from April 04
do k = 1,kl
  tv(1:ifull,k) = (.61*qg(1:ifull,k)-qfg(1:ifull,k)-qlg(1:ifull,k))*t(1:ifull,k)  ! just add-on at this stage 
end do
contv = (1.61-cpv/cp)/.61      ! about -.26/.61
if ( ktau==1 .and. myid==0 ) then
  write(6,*)'in nonlin ntbar =',ntbar 
end if
! Note that ntbar=-1 is altered in globpe to a +ve value
if ( ntbar==-2 .and. num==0 ) then
  tbar2d(1:ifull) = t(1:ifull,1) + contv*tv(1:ifull,1)
else if ( ntbar==0 ) then
  tbar2d(1:ifull) = tbar(1)
else if ( ntbar>0 ) then
  tbar2d(1:ifull) = t(1:ifull,ntbar)
else if ( ntbar==-3 ) then
  tbar2d(1:ifull) = max( t(1:ifull,1), t(1:ifull,2), t(1:ifull,3), t(1:ifull,kl) )
else if ( ntbar==-4 ) then
  tbar2d(1:ifull) = max( t(1:ifull,1), t(1:ifull,2), t(1:ifull,4), t(1:ifull,kl) )
end if     ! (ntbar==-4)
      
! update (linearized) hydrostatic geopotential phi
phi(:,1) = zs(1:ifull) + bet(1)*t(1:ifull,1) 
do k = 2,kl
  phi(:,k) = phi(:,k-1) + bet(k)*t(1:ifull,k) + betm(k)*t(1:ifull,k-1)
end do    ! k  loop
      
! update non-hydrostatic terms from Miller-White height equation
if ( nh/=0 ) then
  invconst_nh = 0.5*dt*grav*grav*(1.-epsh)/rdry  
  do k = 1,kl  
    phi(:,k) = phi(:,k) + phi_nh(:,k)
    h_nh(1:ifull,k) = tbar(1)*dpsldt(:,k)/sig(k)
  enddo
  if ( nmaxpr==1 ) then
    if ( mydiag ) write(6,*) 'h_nh.a ',(h_nh(idjd,k),k=1,kl)
  end if
  select case(nh)
    case(2) ! was -2 add in other term explicitly, more consistently
      ! N.B. nh=2 needs lapsbot=3        
      do k = 2,kl
        h_nh(1:ifull,k) = h_nh(1:ifull,k) - ((phi(:,k)-phi(:,k-1))/bet(k)+t(1:ifull,k))*invconst_nh/tbar2d(:)
      enddo
      k = 1
      h_nh(1:ifull,k) = h_nh(1:ifull,k) - ((phi(:,k)-zs(1:ifull))/bet(k)+t(1:ifull,k))*invconst_nh/tbar2d(:)
    case(3)
      do k = 2,kl-1
        ! now includes epst
        h_nh(1:ifull,k) = h_nh(1:ifull,k) - (sig(k)*(phi(:,k+1)-phi(:,k-1))/(rdry*(sig(k+1)-sig(k-1))) &
                        + t(1:ifull,k))*invconst_nh/tbar2d(:)
      enddo
      k = 1
      h_nh(1:ifull,k) = h_nh(1:ifull,k) - (sig(k)*(phi(:,k+1)-zs(1:ifull))/(rdry*(sig(k+1)-1.))        &
                      + t(1:ifull,k))*invconst_nh/tbar2d(:)
      k = kl
      h_nh(1:ifull,k) = h_nh(1:ifull,k) - (sig(k)*(phi(:,k)-phi(:,k-1))/(rdry*(sig(k)-sig(k-1)))       &
                      + t(1:ifull,k))*invconst_nh/tbar2d(:)
    case(4) ! was -3 add in other term explicitly, more accurately?
      do k = 2,kl-1
        h_nh(1:ifull,k) = h_nh(1:ifull,k) - (((sig(k)-sig(k-1))*(phi(:,k+1)-phi(:,k))/(sig(k+1)-sig(k))+                  &
                           ((sig(k+1)-sig(k))*(phi(:,k)-phi(:,k-1))/(sig(k)-sig(k-1))))*sig(k)/(rdry*(sig(k+1)-sig(k-1))) &
                         + t(1:ifull,k))*invconst_nh/tbar2d(:)
      end do
      k = 1
      h_nh(1:ifull,k) = h_nh(1:ifull,k) - (((sig(k)-1.)*(phi(:,k+1)-phi(:,k))/(sig(k+1)-sig(k))+               &
                          ((sig(k+1)-sig(k))*(phi(:,k)-zs(1:ifull))/(sig(k)-1.)))*sig(k)/(rdry*(sig(k+1)-1.))  &
                      + t(1:ifull,k))*invconst_nh/tbar2d(:)
    case(5,6)
      ! MJT - This method is compatible with bet(k) and betm(k)
      ! This is the similar to nh==2, but works for all lapsbot
      ! and only involves phi_nh as the hydrostatic component
      ! is eliminated.
      ! t_nh = -(sig/rdry)*d(phi_nh)/d(sig)
      ! so phi_nh(k) - phi_nh(k-1) = bet(k)*t_nh(k) + betm(k)*t_nh(k-1)
      t_nh(:) = phi_nh(:,1)/bet(1)
      h_nh(1:ifull,1) = h_nh(1:ifull,1) + t_nh(:)*invconst_nh/tbar2d(:)
      do k = 2,kl
        t_nh(:) = (phi_nh(:,k)-phi_nh(:,k-1)-betm(k)*t_nh(:))/bet(k)
        h_nh(1:ifull,k) = h_nh(1:ifull,k) + t_nh(:)*invconst_nh/tbar2d(:)
      end do
    case(8)
      t_nh(:) = phi_nh(:,1)/bet(1)
      h_nh(1:ifull,1) = h_nh(1:ifull,1) + t_nh(:)*invconst_nh*tbar(1)/t(1:ifull,1)**2
      do k = 2,kl
        t_nh(:) = (phi_nh(:,k)-phi_nh(:,k-1)-betm(k)*t_nh(:))/bet(k)
        h_nh(1:ifull,k) = h_nh(1:ifull,k) + t_nh(:)*invconst_nh*tbar(1)/t(1:ifull,k)**2
      end do
  end select
  if ( nmaxpr==1 ) then
    if ( mydiag ) then
      write(6,*) 'h_nh.b ',(h_nh(idjd,k),k=1,kl)
      write(6,*) 'phi    ',(phi(idjd,k),k=1,kl)
      write(6,*) 'phi_nh ',(phi_nh(idjd,k),k=1,kl)
    end if
    call maxmin(h_nh,'h_',ktau,1.,kl)
  end if
else
  do k = 1,kl  
    phi_nh(:,k) = 0. ! set to hydrostatic approximation
    h_nh(:,k) = 0.
  end do  
end if      ! (nh/=0) ..else..

do k = 1,kl
  tn(1:ifull,k) = tn(1:ifull,k) + (t(1:ifull,k)+contv*tv(1:ifull,k)-tbar2d(1:ifull))*dpsldt(1:ifull,k)*roncp/sig(k) 
  termlin(1:ifull) = tbar2d(1:ifull)*dpsldt(1:ifull,k)*roncp/sig(k) ! full dpsldt used here  
  ! add in  cnon*dt*tn(iq,k)  term at bottom
  tx(1:ifull,k) = t(1:ifull,k) + 0.5*dt*(1.-epst(1:ifull))*termlin(1:ifull)  
end do      ! k  loop

if ( (diag.or.nmaxpr==1) .and. mydiag ) then
  iq = idjd
  k = nlv
  write(6,*) 'dpsldt,roncp,sig ',dpsldt(iq,k),roncp,sig(k)
  write(6,*) 'contv,tbar2d,termlin_nlv ',contv,tbar2d(iq),tbar2d(iq)*dpsldt(iq,k)*roncp/sig(k)
  write(6,*) 'tv,tn ',tv(iq,k),tn(iq,k)
end if
               
! calculate (linearized) augmented geopotential height terms and save in p
do k = 1,kl
  p(1:ifull,k) = phi(1:ifull,k) + rdry*tbar2d(1:ifull)*psl(1:ifull)
end do      ! k  loop

! calculate virtual temp extra terms (mainly -ve): phi -phi_v +r*tbar*psl
phiv(1:ifull,1) = rdry*tbar2d(1:ifull)*psl(1:ifull) - bet(1)*tv(1:ifull,1)
do k = 2,kl
  phiv(1:ifull,k) = phiv(1:ifull,k-1) - bet(k)*tv(1:ifull,k) - betm(k)*tv(1:ifull,k-1)
end do    ! k  loop

! also need full Tv
do k = 1,kl
  tv(1:ifull,k) = t(1:ifull,k) + tv(1:ifull,k)  
end do


! MJT notes - This is the first bounds call after the physics
! routines, so load balance is a significant issue
duma(ifull+1:ifull+iextra,:) = 0. ! avoids float invalid errors
call bounds(p,nehalf=.true.)
call bounds(tv,nehalf=.true.)
duma(1:ifull,1:kl) = phiv(1:ifull,:)
duma(1:ifull,kl+1) = psl(1:ifull)
call bounds(duma(:,1:kl+1))
phiv(ifull+1:ifull+iextra,1:kl) = duma(ifull+1:ifull+iextra,1:kl)
psl(ifull+1:ifull+iextra)       = duma(ifull+1:ifull+iextra,kl+1)


do k = 1,kl
  ! calculate staggered ux,vx first
  aa(1:ifull,k) = -0.5*dt*emu(1:ifull)*(p(ie,k)-p(1:ifull,k))*(1.-epsu)/ds
  bb(1:ifull,k) = -0.5*dt*emv(1:ifull)*(p(in,k)-p(1:ifull,k))*(1.-epsu)/ds

  ! calculate staggered dyn residual contributions first
  un(1:ifull,k) = emu(1:ifull)*(phiv(ie,k)-phiv(1:ifull,k)-0.5*rdry*(tv(ie,k)+tv(1:ifull,k))*(psl(ie)-psl(1:ifull)))/ds
  vn(1:ifull,k) = emv(1:ifull)*(phiv(in,k)-phiv(1:ifull,k)-0.5*rdry*(tv(in,k)+tv(1:ifull,k))*(psl(in)-psl(1:ifull)))/ds
  aa(1:ifull,k) = aa(1:ifull,k) + 0.5*dt*un(1:ifull,k) ! still staggered
  bb(1:ifull,k) = bb(1:ifull,k) + 0.5*dt*vn(1:ifull,k) ! still staggered
end do    ! k loop

if ( diag ) then
  if ( mydiag ) then
    write(6,*) 'tv ',tv(idjd,:)
    write (6,"('tn1*dt',9f8.3/6x,9f8.3)") tn(idjd,:)*dt
    write (6,"('un1*dt',9f8.3/6x,9f8.3)") un(idjd,:)*dt
    write (6,"('vn1*dt',9f8.3/6x,9f8.3)") vn(idjd,:)*dt
  end if
end if                     ! (diag)


! Bounds call for unstaggering winds
! Bounds call for unstaggering winds
call unstaguv(aa,bb,ux,vx) ! convert to unstaggered positions
call unstaguv(un,vn,un,vn) ! convert to unstaggered positions


if ( diag ) then
  call printa('aa  ',aa,ktau,nlv,ia,ib,ja,jb,0.,1.)
  call printa('bb  ',bb,ktau,nlv,ia,ib,ja,jb,0.,1.)
end if                     ! (diag)

do k = 1,kl
  ux(1:ifull,k) = u(1:ifull,k) + ux(1:ifull,k)
  vx(1:ifull,k) = v(1:ifull,k) + vx(1:ifull,k)
     
  tx(1:ifull,k) = tx(1:ifull,k) + 0.5*dt*tn(1:ifull,k)
end do  

if ( diag ) then
  if ( mydiag ) then
    write(6,*) 'at end of nonlin; idjd,dt', idjd,dt
    write(6,*) 'p1 . & e ',p(idjd,nlv),p(ie(idjd),nlv)
    write(6,*) 'p1 . & n ',p(idjd,nlv),p(in(idjd),nlv)
    write(6,*) 'tx ',tx(idjd,:)
    write (6,"('tn2*dt',9f8.3/6x,9f8.3)")   tn(idjd,:)*dt
    write (6,"('un2*dt',9f8.3/6x,9f8.3)")   un(idjd,:)*dt
    write (6,"('vn2*dt',9f8.3/6x,9f8.3)")   vn(idjd,:)*dt
    write (6,"('ux  ',9f8.2/4x,9f8.2)")     ux(idjd,:)
    write (6,"('vx  ',9f8.2/4x,9f8.2)")     vx(idjd,:)
  end if
  call printa('psl ',psl,ktau,0,ia,ib,ja,jb,0.,100.)
  call printa('pslx',pslx,ktau,nlv,ia,ib,ja,jb,0.,100.)
  call printa('tn  ',tn,ktau,nlv,ia,ib,ja,jb,0.,100.*dt)
  call printa('un  ',un,ktau,nlv,ia,ib,ja,jb,0.,100.*dt)
  call printa('vn  ',vn,ktau,nlv,ia,ib,ja,jb,0.,100.*dt)
  call printa('tx  ',tx,ktau,nlv,ia,ib,ja,jb,200.,1.)
  call printa('ux  ',ux,ktau,nlv,ia,ib,ja,jb,0.,1.)
  call printa('vx  ',vx,ktau,nlv,ia,ib,ja,jb,0.,1.)
end if

num = 1

call END_LOG(nonlin_end)

return
end subroutine nonlin
