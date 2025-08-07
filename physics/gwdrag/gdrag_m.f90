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
    
module gdrag_m

implicit none

private
public he,helo
public gdrag_init,gdrag_sbl,gdrag_end,gwdrag

real, dimension(:), allocatable, save :: helo
real, dimension(:), allocatable, save :: he
integer, save :: kbot_gwd

contains

subroutine gdrag_init(ifull)

implicit none

integer, intent(in) :: ifull

allocate(he(ifull),helo(ifull))

return
end subroutine gdrag_init

subroutine gdrag_sbl

use parm_m, only : sigbot_gwd
use sigs_m, only : sig
use cc_mpi, only : mydiag

implicit none

integer, dimension(1) :: kpos

! set bottom level, above which permit wave breaking
! set sigbot_gw<.5 to give kbot_gwd=1 and use bvng, very similar to older scheme
kbot_gwd = 1
if ( sigbot_gwd>=.5 ) then
  kpos = minloc(abs(sig-sigbot_gwd)) ! finds k value closest to sigbot_gwd    
  kbot_gwd = kpos(1) ! JLM
end if
if ( mydiag ) write(6,*) 'in gwdrag sigbot_gwd,kbot:',sigbot_gwd,kbot_gwd
! MJT notes - he defined in indata.f90 before calling gdrag_sbl

return
end subroutine gdrag_sbl

subroutine gdrag_end

implicit none

deallocate(he,helo)

return
end subroutine gdrag_end

subroutine gwdrag

use cc_mpi, only : mydiag
use arrays_m
use newmpar_m
use parm_m, only : idjd
use pbl_m

implicit none

integer tile, js, je
integer idjd_t
real, dimension(imax,kl) :: lt, lu, lv
logical mydiag_t

do tile = 1,ntiles
  js = (tile-1)*imax + 1
  je = tile*imax
  
  idjd_t = mod(idjd-1,imax) + 1
  mydiag_t = ((idjd-1)/imax==tile-1).and.mydiag
  
  lt = t(js:je,:)
  lu = u(js:je,:)
  lv = v(js:je,:)
  
  call gwdrag_work(lt,lu,lv,tss(js:je),he(js:je),idjd_t,mydiag_t)

  u(js:je,:) = lu
  v(js:je,:) = lv
 
end do

return
end subroutine gwdrag

!  this is vectorized jlm version with kbot_gwd generalization July 2015
!  Parameters and suggested values (jlm July 2015):
!  ngwd  -20  (similar to Chouinard; previously we used weaker gwdrag with -5)
!  helim 1600 limit of launching height (similar to Chouinard; previously we used weaker 800 m)
!  fc2  0.5 limit of Froude number**2 (previously we had tighter limit of 1.)
!       -ve value forces wave breaking at top level, even if fc2 condn not satisfied
!  sigbot_gwd 0.8 breaking may only occur from this sigma level up (previously 1.)
    
subroutine gwdrag_work(t,u,v,tss,he,idjd,mydiag)

use const_phys, only : grav, rdry, cp
use parm_m, only : vmodmin, sigbot_gwd, fc2, dt, alphaj, ngwd
use sigs_m, only : sig, dsig

implicit none

integer, parameter :: ntest = 0 ! ntest= 0 for diags off; ntest= 1 for diags on
integer, intent(in) :: idjd
integer iq, k, imax, kl
real, dimension(:,:), intent(in)    :: t
real, dimension(size(t,1),size(t,2)), intent(inout) :: u, v
real, dimension(size(t,1),size(t,2)) :: uu,fni,bvnf
real, dimension(size(t,1),size(t,2)) :: theta_full
real, dimension(size(t,1),size(t,2)) :: dtheta_dz_kmh
real, dimension(size(t,1)), intent(in) :: tss, he
real, dimension(size(t,1)) :: temp,fnii
real, dimension(size(t,1)) :: bvng ! to be depreciated
real, dimension(size(t,1)) :: apuw,apvw,alambda,wmag
real, dimension(size(t,2)) :: dsk,sigk
real dzx, dzi, froude2_inv, uux, xxx
logical, intent(in) :: mydiag

imax = size(t,1)
kl = size(t,2)

! older values:  
!   ngwd=-5  helim=800.  fc2=1.  sigbot_gw=0. alphaj=1.E-6 (almost equiv to 0.0075)
! new JLM suggestion to resemble Chouinard et al values (apart from alphaj):
!   ngwd=-20 helim=1600. fc2=-.5 sigbot_gw=1. alphaj=0.05
! If desire to tune, only need to vary alphaj (increase for stronger GWD)

do k = 1,kl
  dsk(k) = -dsig(k)
  sigk(k) = sig(k)**(rdry/cp)
end do

! put theta in theta_full()
do k = 1,kl
  do iq = 1,imax
    theta_full(iq,k) = t(iq,k)/sigk(k)                ! gwdrag
  end do  
end do

!  calc d(theta)/dz  at half-levels , using 1/dz at level k-.5
dzx = .5*grav*(1.+sig(1))/((1.-sig(1))*rdry)
do iq = 1,imax    
  dzi = dzx/t(iq,1)
  dtheta_dz_kmh(iq,1) = (theta_full(iq,1)-tss(iq))*dzi    
end do
do k = 2,kl
  do iq = 1,imax  
    dzx = grav*(sig(k-1)+sig(k))/((sig(k-1)-sig(k))*rdry)    
    dzi = dzx/(t(iq,k-1)+t(iq,k)) 
    dtheta_dz_kmh(iq,k) = (theta_full(iq,k)-theta_full(iq,k-1))*dzi
  end do
end do    ! k loop          

!     form wmag at surface
wmag(:) = sqrt(max(u(:,1)**2+v(:,1)**2, vmodmin**2)) ! MJT suggestion


!**** calculate Brunt-Vaisala frequency at full levels (bvnf)
!      if unstable reference level then no gwd 
!            - happens automatically via effect on bvnf & alambda
do k = 1,kl-1
  do iq = 1,imax
    bvnf(iq,k) = sqrt(max(1.e-20, grav*0.5*(dtheta_dz_kmh(iq,k)+dtheta_dz_kmh(iq,k+1))/theta_full(iq,k)) ) ! MJT fixup
  end do    
end do    ! k loop
bvnf(:,kl) = sqrt(max(1.e-20, grav*dtheta_dz_kmh(:,kl)/theta_full(:,kl)))    ! jlm fixup

!**    calc (t*/n*/wmag)/he**2
if ( sigbot_gwd<.5 ) then !  to be depreciated
  bvng(:) = sqrt(max(1.e-20, grav*dtheta_dz_kmh(:,1)/tss(:))) ! tries to use a sfce value rather than level 1 
  temp(:) = tss(:)/max(bvng(:)*wmag(:)*he(:)**2, 1.e-10) 
else
  temp(:) = theta_full(:,1)/max(bvnf(:,1)*wmag(:)*he(:)**2, 1.e-10)      
end if

do k = 1,2 ! uu is +ve wind compt in dirn of (u_1,v_1)
  do iq = 1,imax
    uu(iq,k) = max(0., u(iq,k)*u(iq,1)+v(iq,k)*v(iq,1))/wmag(iq)
  end do
end do    ! k loop

!**** set uu() to zero above if uu() zero below
!**** uu>0 at k=1, uu>=0 at k=1+1 - only set for k=1+2 to kl  
do k = 3,kl
  do iq = 1,imax
    if ( uu(iq,k-1)<1.e-20 ) then
      uu(iq,k) = 0.
    else
      uu(iq,k) = max(0., u(iq,k)*u(iq,1)+v(iq,k)*v(iq,1))/wmag(iq)
    end if
  end do  
end do    ! k loop

!       calc max(1-Fc**2/F**2,0) : put in fni()
do k = kbot_gwd,kl
  do iq = 1,imax
    froude2_inv = sig(k)*temp(iq)*uu(iq,k)**3/(sigk(k)*bvnf(iq,k)*theta_full(iq,k))
    fni(iq,k) = max(0., 1.-abs(fc2)*froude2_inv)
  end do
end do    ! k loop
if ( fc2<0. ) then
  fni(:,kl) = 1.
end if

! form integral of above*uu**2 from sig=sigbot_gwd to sig=0
fnii(:) = -fni(:,kbot_gwd)*dsig(kbot_gwd)*uu(:,kbot_gwd)*uu(:,kbot_gwd)
do k = kbot_gwd+1,kl
  do iq = 1,imax
    fnii(iq) = fnii(iq)-fni(iq,k)*dsig(k)*uu(iq,k)*uu(iq,k)
  end do  
end do    ! k loop

!     Chouinard et al. use alpha=.01
!alphaj=0.01*1.e-4  ! jlm   .01 *rhos*g/ps
!      if integral=0., reset to some +ve value
!      form alambda=(g/p*).alpha.rhos.he.N*.wmag/integral(above)
if ( alphaj<1.e-5 ) then  ! for backward compatibility - will depreciate
  alambda(:) = alphaj*he(:)*bvnf(:,kbot_gwd)*wmag(:)/max(fnii(:), 1.e-9)
else  ! newer usage with alphaj around 0.0075 (similar to resemble Hal's value)
  alambda(:) = alphaj*he(:)*bvnf(:,kbot_gwd)*wmag(:)*grav/(rdry*tss(:)*max(fnii(:), 1.e-9))
end if  
!      define apuw=alambda.u1/wmag , apvw=alambda.v1/wmag
apuw(:) = alambda(:)*u(:,1)/wmag(:)
apvw(:) = alambda(:)*v(:,1)/wmag(:)

do k = kbot_gwd,kl
  do iq = 1,imax
    !**** form fni=alambda*max(--,0) and
    !**** solve for uu at t+1 (implicit solution)
    uux = 2.*uu(iq,k)/(1.+sqrt(1. + 4.*dt*alambda(iq)*fni(iq,k)*uu(iq,k)))
    !       N.B. 4.*dt*alambda(iq)*fni(iq,k)*uu(iq,k)) can be ~300
    !**** form dv/dt due to gw-drag at each level
    !**** = -alambda.v*/wmag.uu(t+1)**2.max(--,0)
    xxx = uux**2*fni(iq,k)
    u(iq,k) = u(iq,k) - apuw(iq)*xxx*dt
    v(iq,k) = v(iq,k) - apvw(iq)*xxx*dt
  end do
end do     ! k loop

if ( ntest==1 .and. mydiag ) then ! JLM
  do iq = idjd-1,idjd+1
    write(6,*) 'from gwdrag, iq,ngwd,alambda,fnii,apuw,apvw,wmag',  &
    iq,ngwd,alambda(iq),fnii(iq),apuw(iq),apvw(iq),wmag(iq)
    write(6,*) 'temp,bvnf_1,he,tss',   &
      temp(iq),bvnf(iq,1),he(iq),tss(iq)
    if ( sigbot_gwd<.5 ) write(6,*) 'bvng',bvng(iq) !  to be depreciated
    write(6,*) 't',t(iq,kbot_gwd:kl)
    write(6,*) 'theta_full',theta_full(iq,1:kl)
    write(6,*) 'bvnf',bvnf(iq,1:kl)
    write(6,*) 'dtheta_dz_kmh',dtheta_dz_kmh(iq,1:kl)
    write(6,*) 'uu',uu(iq,kbot_gwd:kl)
    !write(6,*) 'froude2_inv',froude2_inv(iq,kbot_gwd:kl)
    write(6,*) 'fni',fni(iq,kbot_gwd:kl)
    !write(6,*) 'uux',uux(iq,kbot_gwd:kl)
!   following in reverse k order to assist viewing by grep	 
    !write(6,*) 'uincr',(-apuw(iq)*xxx(iq,k)*dt,k=kl,kbot_gwd,-1)
    !write(6,*) 'vincr',(-apvw(iq)*xxx(iq,k)*dt,k=kl,kbot_gwd,-1)
  end do
end if

return
end subroutine gwdrag_work

end module gdrag_m
