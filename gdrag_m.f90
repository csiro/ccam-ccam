! Conformal Cubic Atmospheric Model
    
! Copyright 2015 Commonwealth Scientific Industrial Research Organisation (CSIRO)
    
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

real, dimension(:), allocatable, save :: he,helo
integer, save :: nb, imax
integer, save :: kbot

type blocked_data_1d
  real, dimension(:), allocatable :: data
end type blocked_data_1d
type blocked_data_2d
  real, dimension(:,:), allocatable :: data
end type blocked_data_2d
type(blocked_data_1d), dimension(:), allocatable, save :: b_he
type(blocked_data_1d), dimension(:), allocatable, save :: b_tss
type(blocked_data_2d), dimension(:), allocatable, save :: b_phi_nh
type(blocked_data_2d), dimension(:), allocatable, save :: b_t
type(blocked_data_2d), dimension(:), allocatable, save :: b_u
type(blocked_data_2d), dimension(:), allocatable, save :: b_v

contains

subroutine gdrag_init(ifull,nbin)
use newmpar_m, only : kl
use arrays_m, only : t,u,v
use nharrs_m, only : phi_nh
use pbl_m, only : tss

implicit none

integer, intent(in) :: ifull,nbin
integer :: i

allocate(he(ifull),helo(ifull))
nb=nbin
imax=ifull/nb

allocate(b_he(nb))
allocate(b_tss(nb))
allocate(b_phi_nh(nb))
allocate(b_t(nb))
allocate(b_u(nb))
allocate(b_v(nb))
do i=1,nb
  allocate(b_he(i)%data(imax))
  allocate(b_tss(i)%data(imax))
  allocate(b_phi_nh(i)%data(imax,kl))
  allocate(b_t(i)%data(imax,kl))
  allocate(b_u(i)%data(imax,kl))
  allocate(b_v(i)%data(imax,kl))
end do

return
end subroutine gdrag_init

subroutine gdrag_sbl
use parm_m, only : sigbot_gwd
use sigs_m, only : sig
use cc_mpi, only : mydiag

implicit none

integer, dimension(1) :: kpos

! set bottom level, above which permit wave breaking
! set sigbot_gw<.5 to give kbot=1 and use bvng, very similar to older scheme
kbot = 1
if ( sigbot_gwd>=.5 ) then
  kpos = minloc(abs(sig-sigbot_gwd)) ! finds k value closest to sigbot_gwd    
  kbot = kpos(1) ! JLM
end if
if ( mydiag ) write(6,*) 'in gwdrag sigbot_gwd,kbot:',sigbot_gwd,kbot

return
end subroutine gdrag_sbl

subroutine gdrag_end
use arrays_m, only : t,u,v
use nharrs_m, only : phi_nh
use pbl_m, only : tss

implicit none
integer :: i

deallocate(he,helo)
do i=1,nb
  deallocate(b_he(i)%data)
  deallocate(b_tss(i)%data)
  deallocate(b_phi_nh(i)%data)
  deallocate(b_t(i)%data)
  deallocate(b_u(i)%data)
  deallocate(b_v(i)%data)
end do
deallocate(b_he)
deallocate(b_tss)
deallocate(b_phi_nh)
deallocate(b_t)
deallocate(b_u)
deallocate(b_v)

return
end subroutine gdrag_end

subroutine gwdrag
use cc_mpi
use arrays_m, only : t,u,v
use nharrs_m, only : phi_nh
use pbl_m, only : tss

implicit none
integer :: i,is,ie

do i=1,nb
  is=(i-1)*imax+1
  ie=i*imax
  b_he(i)%data=he(is:ie)
  b_tss(i)%data=tss(is:ie)
  b_phi_nh(i)%data=phi_nh(is:ie,:)
  b_t(i)%data=t(is:ie,:)
  b_u(i)%data=u(is:ie,:)
  b_v(i)%data=v(is:ie,:)
end do
call start_log(gdrag_begin)
!$omp parallel do
do i=1,nb
  call gwdrag_work(i)
end do
call end_log(gdrag_end)
do i=1,nb
  is=(i-1)*imax+1
  ie=i*imax
  u(is:ie,:)=b_u(i)%data
  v(is:ie,:)=b_v(i)%data
end do

end subroutine gwdrag

subroutine gwdrag_work(tile)   ! globpea/darlam (but not staggered)
!  this is vectorized jlm version with kbot generalization July 2015
!  Parameters and suggested values (jlm July 2015):
!  ngwd  -20  (similar to Chouinard; previously we used weaker gwdrag with -5)
!  helim 1600 limit of launching height (similar to Chouinard; previously we used weaker 800 m)
!  fc2  0.5 limit of Froude number**2 (previously we had tighter limit of 1.)
!       -ve value forces wave breaking at top level, even if fc2 condn not satisfied
!  sigbot_gwd 0.8 breaking may only occur from this sigma level up (previously 1.)
use arrays_m
use cc_mpi, only : mydiag
use const_phys
use liqwpar_m
use newmpar_m
use nharrs_m
use parm_m
use pbl_m
use sigs_m
#ifdef _OPENMP
use omp_lib, only : omp_in_parallel
#endif
implicit none
integer, intent(in) :: tile
integer, parameter :: ntest = 0 ! ntest= 0 for diags off; ntest= 1 for diags on
integer iq,k
real dzx
real, dimension(1:imax,kl) :: uu,fni,bvnf
real, dimension(1:imax,kl) :: theta_full
real, dimension(1:imax) :: dzi, uux, xxx, froude2_inv
real, dimension(1:imax,kl) :: tnhs
real, dimension(1:imax,kl) :: dtheta_dz_kmh
real, dimension(1:imax) :: temp,fnii
real, dimension(1:imax) :: bvng ! to be depreciated
real, dimension(1:imax) :: apuw,apvw,alambda,wmag
real, dimension(kl) :: dsk,sigk
logical :: serial

#ifdef _OPENMP
serial=.not.omp_in_parallel()
#else
serial=.true.
#endif

! older values:  
!   ngwd=-5  helim=800.  fc2=1.  sigbot_gw=0. alphaj=1.E-6 (almost equiv to 0.0075)
! new JLM suggestion to resemble Chouinard et al values (apart from alphaj):
!   ngwd=-20 helim=1600. fc2=-.5 sigbot_gw=1. alphaj=0.05
! If desire to tune, only need to vary alphaj (increase for stronger GWD)

! Non-hydrostatic terms
tnhs(:,1) = b_phi_nh(tile)%data(:,1)/bet(1)
do k = 2,kl
  ! representing non-hydrostatic term as a correction to air temperature
  tnhs(:,k) = (b_phi_nh(tile)%data(:,k)-b_phi_nh(tile)%data(:,k-1)-betm(k)*tnhs(:,k-1))/bet(k)
end do      

do k = 1,kl
  dsk(k) = -dsig(k)
  sigk(k) = sig(k)**(rdry/cp)
  ! put theta in theta_full()
  theta_full(:,k) = b_t(tile)%data(1:imax,k)/sigk(k)                ! gwdrag
end do

!  calc d(theta)/dz  at half-levels , using 1/dz at level k-.5
dzx = .5*grav*(1.+sig(1))/((1.-sig(1))*rdry)    
dzi(:) = dzx/(b_t(tile)%data(1:imax,1)+tnhs(:,1))
dtheta_dz_kmh(:,1) = (theta_full(:,1)-b_tss(tile)%data(:))*dzi(:)    
do k = 2,kl
 dzx = grav*(sig(k-1)+sig(k))/((sig(k-1)-sig(k))*rdry)  
 dzi(:) = dzx/(b_t(tile)%data(1:imax,k-1)+b_t(tile)%data(1:imax,k)+tnhs(:,k-1)+tnhs(:,k)) 
 dtheta_dz_kmh(:,k) = (theta_full(:,k)-theta_full(:,k-1))*dzi(:)
end do    ! k loop          

!     form wmag at surface
wmag(1:imax) = sqrt(max(b_u(tile)%data(1:imax,1)**2+b_v(tile)%data(1:imax,1)**2, vmodmin**2)) ! MJT suggestion


!**** calculate Brunt-Vaisala frequency at full levels (bvnf)
!      if unstable reference level then no gwd 
!            - happens automatically via effect on bvnf & alambda
do k = 1,kl-1
  bvnf(:,k) = sqrt(max(1.e-20, grav*0.5*(dtheta_dz_kmh(:,k)+dtheta_dz_kmh(:,k+1))/theta_full(:,k)) ) ! MJT fixup
end do    ! k loop
bvnf(:,kl) = sqrt(max(1.e-20, grav*dtheta_dz_kmh(:,kl)/theta_full(:,kl)))    ! jlm fixup

!**    calc (t*/n*/wmag)/he**2
if ( sigbot_gwd<.5 ) then !  to be depreciated
  bvng(1:imax) = sqrt(max(1.e-20, grav*dtheta_dz_kmh(1:imax,1)/b_tss(tile)%data(1:imax))) ! tries to use a sfce value rather than level 1 
  temp(1:imax) = b_tss(tile)%data(1:imax)/max(bvng(1:imax)*wmag(1:imax)*b_he(tile)%data(1:imax)**2, 1.e-10) 
else
  temp(1:imax) = theta_full(1:imax,1)/max(bvnf(1:imax,1)*wmag(1:imax)*b_he(tile)%data(1:imax)**2, 1.e-10)      
end if

do k = 1,2 ! uu is +ve wind compt in dirn of (u_1,v_1)
  uu(1:imax,k) = max(0., b_u(tile)%data(1:imax,k)*b_u(tile)%data(1:imax,1)+b_v(tile)%data(1:imax,k)*b_v(tile)%data(1:imax,1))/wmag(1:imax)
end do    ! k loop

!**** set uu() to zero above if uu() zero below
!**** uu>0 at k=1, uu>=0 at k=1+1 - only set for k=1+2 to kl  
do k = 3,kl
  where ( uu(1:imax,k-1)<1.e-20 )
    uu(1:imax,k) = 0.
  elsewhere
    uu(1:imax,k) = max(0., b_u(tile)%data(1:imax,k)*b_u(tile)%data(1:imax,1)+b_v(tile)%data(1:imax,k)*b_v(tile)%data(1:imax,1))/wmag(1:imax)
  end where
end do    ! k loop

do k = kbot,kl
  !       calc max(1-Fc**2/F**2,0) : put in fni()
  froude2_inv(:) = sig(k)*temp(:)*uu(:,k)**3/(sigk(k)*bvnf(:,k)*theta_full(:,k))
  fni(:,k) = max(0., 1.-abs(fc2)*froude2_inv(:))
end do    ! k loop
if ( fc2<0. ) then
  fni(:,kl) = 1.
end if

! form integral of above*uu**2 from sig=sigbot_gwd to sig=0
fnii(:) = -fni(:,kbot)*dsig(kbot)*uu(:,kbot)*uu(:,kbot)
do k = kbot+1,kl
  fnii(:) = fnii(:)-fni(:,k)*dsig(k)*uu(:,k)*uu(:,k)
end do    ! k loop

!     Chouinard et al. use alpha=.01
!alphaj=0.01*1.e-4  ! jlm   .01 *rhos*g/ps
!      if integral=0., reset to some +ve value
!      form alambda=(g/p*).alpha.rhos.he.N*.wmag/integral(above)
if ( alphaj<1.e-5 ) then  ! for backward compatibility - will depreciate
  alambda(1:imax) = alphaj*b_he(tile)%data(1:imax)*bvnf(1:imax,kbot)*wmag(1:imax)/max(fnii(1:imax), 1.e-9)
else  ! newer usage with alphaj around 0.0075 (similar to resemble Hal's value)
  alambda(1:imax) = alphaj*b_he(tile)%data(1:imax)*bvnf(1:imax,kbot)*wmag(1:imax)*grav/(rdry*b_tss(tile)%data(1:imax)*max(fnii(1:imax), 1.e-9))
end if  
!      define apuw=alambda.u1/wmag , apvw=alambda.v1/wmag
apuw(1:imax) = alambda(1:imax)*b_u(tile)%data(1:imax,1)/wmag(1:imax)
apvw(1:imax) = alambda(1:imax)*b_v(tile)%data(1:imax,1)/wmag(1:imax)

do k = kbot,kl
  !**** form fni=alambda*max(--,0) and
  !**** solve for uu at t+1 (implicit solution)
  uux(:) = 2.*uu(:,k)/(1.+sqrt(1. + 4.*dt*alambda(:)*fni(:,k)*uu(:,k)))
  !       N.B. 4.*dt*alambda(iq)*fni(iq,k)*uu(iq,k)) can be ~300
  !**** form dv/dt due to gw-drag at each level
  !**** = -alambda.v*/wmag.uu(t+1)**2.max(--,0)
  xxx(:) = uux(:)*uux(:)*fni(:,k)
  b_u(tile)%data(1:imax,k) = b_u(tile)%data(1:imax,k) - apuw(:)*xxx(:)*dt
  b_v(tile)%data(1:imax,k) = b_v(tile)%data(1:imax,k) - apvw(:)*xxx(:)*dt
end do     ! k loop


if ( ntest==1 .and. mydiag .and. serial .and. nb==1 ) then ! JLM
  do iq = idjd-1,idjd+1
    write(6,*) 'from gwdrag, iq,ngwd,alambda,fnii,apuw,apvw,wmag',  &
    iq,ngwd,alambda(iq),fnii(iq),apuw(iq),apvw(iq),wmag(iq)
    write(6,*) 'temp,bvnf_1,he,tss',   &
      temp(iq),bvnf(iq,1),b_he(tile)%data(iq),b_tss(tile)%data(iq)
    if ( sigbot_gwd<.5 ) write(6,*) 'bvng',bvng(iq) !  to be depreciated
    write(6,*) 't',b_t(tile)%data(iq,kbot:kl)
    write(6,*) 'theta_full',theta_full(iq,1:kl)
    write(6,*) 'bvnf',bvnf(iq,1:kl)
    write(6,*) 'dtheta_dz_kmh',dtheta_dz_kmh(iq,1:kl)
    write(6,*) 'uu',uu(iq,kbot:kl)
    !write(6,*) 'froude2_inv',froude2_inv(iq,kbot:kl)
    write(6,*) 'fni',fni(iq,kbot:kl)
    !write(6,*) 'uux',uux(iq,kbot:kl)
!   following in reverse k order to assist viewing by grep	 
    !write(6,*) 'uincr',(-apuw(iq)*xxx(iq,k)*dt,k=kl,kbot,-1)
    !write(6,*) 'vincr',(-apvw(iq)*xxx(iq,k)*dt,k=kl,kbot,-1)
  end do
end if

return
end subroutine gwdrag_work

end module gdrag_m! Conformal Cubic Atmospheric Model
