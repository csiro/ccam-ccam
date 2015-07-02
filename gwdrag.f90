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
    
subroutine gwdrag   ! globpea/darlam (but not staggered)
!     this is vectorized jlm version
use arrays_m
use gdrag_m
use liqwpar_m
use morepbl_m
use nharrs_m
use nlin_m
use pbl_m
use sigs_m
use soil_m
implicit none
!integer, parameter :: ndzx = 0 !  as per Hal      gwd1a
integer, parameter :: ndzx = 1  !  as per jlm      gwd1b
integer, parameter :: ntest = 0 ! ntest= 0 for diags off; ntest= 1 for diags on
include 'newmpar.h'
include 'const_phys.h'
include 'parm.h'
integer iq,i,j,k
real dzx
real uu(ifull,kl),fni(ifull,kl),bvnf(ifull,kl)
real, dimension(ifull,kl) :: thf
real, dimension(ifull) :: dzi, xxx
real tnhs(ifull,kl),tv(ifull,kl)
real dthdz(ifull,kl)
real bvng(ifull),temp(ifull),fnii(ifull)
real apuw(ifull),apvw(ifull),alam(ifull),wmag(ifull),frsav(ifull)
real dsk(kl),sigk(kl)

! interface mods for darlam & globpe
do k = 1,kl
 dsk(k)=-dsig(k)
 sigk(k)=sig(k)**(rdry/cp)
end do
      
! Non-hydrostatic terms
tv=t(1:ifull,:)*(1.+0.61*qg(1:ifull,:)-qlg(1:ifull,:)-qfg(1:ifull,:))
tnhs(:,1)=phi_nh(:,1)/bet(1)
do k = 2,kl
  ! representing non-hydrostatic term as a correction to air temperature
  tnhs(:,k)=(phi_nh(:,k)-phi_nh(:,k-1)-betm(k)*tnhs(:,k-1))/bet(k)
end do      

do k = 1,kl
!       put theta in thf()
  thf(:,k)=t(1:ifull,k)/sigk(k)                ! gwdrag
end do    ! k loop

!     calc d(theta)/dz  at half-levels , using 1/dz at level k-.5
if (ndzx==0) then
  dzx=-grav*sig(1)/(dsig(1)*rdry)                ! Hal's
  dzi(:)=dzx/(tv(:,1)+tnhs(:,1))
  dthdz(:,1)=max(thf(:,1)-tss(:),0.)*dzi(:)      ! was wrong dzi - jlm
  do k=2,kl
    dzx=-2.*grav*sig(k)/(dsig(k)*rdry)                      ! Hal's
    dzi(:)=dzx/(tv(:,k-1)+tv(:,k)+tnhs(:,k-1)+tnhs(:,k))    ! gwdrag
    dthdz(:,k)=(thf(:,k)-thf(:,k-1))*dzi(:)
  end do    ! k loop          
else if (ndzx==1) then
  dzx=.5*grav*(1.+sig(1))/((1.-sig(1))*rdry)     ! fixup by jlm
  dzi(:)=dzx/(tv(:,1)+tnhs(:,1))
  dthdz(:,1)=max(thf(:,1)-tss(:),0.)*dzi(:)      ! was wrong dzi - jlm
  do k=2,kl
    dzx=grav*(sig(k-1)+sig(k))/((sig(k-1)-sig(k))*rdry)   ! fixup by jlm
    dzi(:)=dzx/(tv(:,k-1)+tv(:,k)+tnhs(:,k-1)+tnhs(:,k))  ! gwdrag
    dthdz(:,k)=(thf(:,k)-thf(:,k-1))*dzi(:)
  end do    ! k loop          
end if

!     form new wmag at surface
wmag(:)=sqrt(max(u(:,1)**2+v(:,1)**2,vmodmin**2)) ! MJT suggestion


!**** calculate Brunt-Vaisala frequency
!**** surface bvng() and full levels bvnf(,)

!     calculate bvng,  (Brunt-Vaisala-N-ground)
!     limit value of bvng to about 50 degrees per 166 m
!     if unstable (surface to bl) then no gwd (set bvng=0)
!           - happens automatically via next line & effect on bvng & alam
bvng(:)=min(.1,sqrt(grav*dthdz(:,1)/tss(:)))
!****  froude number calcs
!****  calculate (fc/f)**2 where fc**2=fc2=0.5  (for Hal & Ch., 1. for jlm)
!**    calc fc2*(t*/n*/wmag)/he**2
temp(:)=fc2*tss(:)/max( bvng(:)*wmag(:)*he(:)**2,1.e-10)  !jlm

!     calculate bvnf at other levels,  (Brunt-Vaisala-N-full)
do k = 1,kl-1
  !bvnf(:,k)=sqrt( max(1.e-20,grav*(dthdz(:,k)+dthdz(:,k+1))/(thf(:,k)+thf(:,k+1))) ) ! jlm fixup
  bvnf(:,k)=sqrt( max(1.e-20,grav*0.5*(dthdz(:,k)+dthdz(:,k+1))/thf(:,k)) ) ! MJT suggestion
end do    ! k loop
bvnf(:,kl)=sqrt(max(1.e-20,grav*dthdz(:,kl)/thf(:,kl)))    ! jlm fixup

do k = 1,2
  uu(:,k)=max(0. , u(:,k)*u(:,1)+v(:,k)*v(:,1))/wmag(:)
end do    ! k loop

!**** set uu() to zero above if uu() zero below
!**** uu>0 at k=1, uu>=0 at k=2 - only set for k=3 to kl  OK
do k = 3,kl
  where ( uu(:,k-1)==0. )
    uu(:,k) = 0.
  elsewhere
    uu(:,k) = max(0. , u(:,k)*u(:,1)+v(:,k)*v(:,1))/wmag(:)
  end where
end do    ! k loop

do k=1,kl
!       calc max(1-fc**2/f**2,0) : put in fni()
  fni(:,k) = max(0. , 1.-sig(k)*temp(:)*uu(:,k)**3/(sigk(k)*bvnf(:,k)*thf(:,k)))
end do    ! k loop

!     form integral of above*uu**2 from sig=1 to sig=0
fnii(:) = -fni(:,1)*dsig(1)*uu(:,1)*uu(:,1)
do k=2,kl
  fnii(:) = fnii(:)-fni(:,k)*dsig(k)*uu(:,k)*uu(:,k)
end do    ! k loop

!     Chouinard et al. use alpha=.01
!alphaj=0.01*1.e-4  ! jlm   .01 *rhos*g/ps
!hal  alphah=0.0075*g/r  ! actually alpha*(g/ps)*rhos  *tss  ! Ch et al
!      if integral=0., reset to some +ve value
!      form alam=(g/p*).alpha.rhos.he.N*.wmag/integral(above)
alam(:) = alphaj*he(:)*bvng(:)*wmag(:)/max(fnii(:),1.e-20)
!      define apuw=alam.u1/wmag , apvw=alam.v1/wmag
apuw(:) = alam(:)*u(:,1)/wmag(:)
apvw(:) = alam(:)*v(:,1)/wmag(:)

!**** form fni=alam*max(--,0) and
!**** solve for uu at t+1 (implicit solution)
do k = 1,kl
  uu(:,k) = 2.*uu(:,k)/(1.+sqrt(1. + 4.*dt*alam(:)*fni(:,k)*uu(:,k)))
!       N.B. 4.*dt*alam(iq)*fni(iq,k)*uu(iq,k)) can be ~300
end do    ! k loop

!**** form dv/dt due to gw-drag at each level
!**** = -alam.v*/wmag.uu(t+1)**2.max(--,0)
do k = 1,kl
  xxx(:) = uu(:,k)*uu(:,k)*fni(:,k)
  u(1:ifull,k) = u(1:ifull,k)-apuw(:)*xxx(:)*dt
  v(1:ifull,k) = v(1:ifull,k)-apvw(:)*xxx(:)*dt
end do     ! k loop

#ifdef debug      
if(ntest==1)then
  write(6,*) 'from gwdrag, ngwd,alam,fnii,apuw,apvw,wmag',ngwd,alam(idjd),fnii(idjd),apuw(idjd),apvw(idjd),wmag(idjd)
  write(6,*) 'temp,bvng,he,tss',temp(idjd),bvng(idjd),he(idjd),tss(idjd)
  write(6,*) 't',t(idjd,:)
  write(6,*) 'thf',thf(idjd,:)
  write(6,*) 'bvnf',bvnf(idjd,:)
  write(6,*) 'dthdz',dthdz(idjd,:)
  write(6,*) 'fni',fni(idjd,:)
  write(6,*) 'uu',uu(idjd,:)
  write(6,*) 'un',vn(idjd,:)
  write(6,*) 'vn',vn(idjd,:)
endif
#endif

return
end subroutine gwdrag
