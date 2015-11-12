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
    
subroutine pbldif(theta,rkh,rkm,uav,vav,cgmap)
! vectorized version      

use arrays_m   !t
use cc_mpi, only : mydiag, myid
use cfrac_m
use diag_m
use extraout_m !ustar
use map_m
use morepbl_m  !fg,eg
use nharrs_m
use sigs_m     !sig,sigmh
use soil_m     !land

implicit none

include 'newmpar.h'
include 'const_phys.h'
include 'kuocom.h'
include 'parm.h'         !dtin

integer, parameter :: ntest=0
integer, parameter :: nrkmin=1   ! 1 original (& from 0510); 2 new; 3 newer
integer, parameter :: npblmin=4  ! 1 original (best for Oz); 2 new ; 3,4 newer
integer kmax,iq
integer k                 ! level index
integer, dimension(ifull) :: iflag

!------------------------------------------------------------------------
! 
! Atmospheric boundary layer computation.
!
! Nonlocal scheme that determines eddy diffusivities based on a
! diagnosed boundary layer height and a turbulent velocity scale;
! also, countergradient effects for heat and moisture, and constituents
! are included, along with temperature and humidity perturbations which 
! measure the strength of convective thermals in the lower part of the 
! atmospheric boundary layer.
!
! For more information, see Holtslag, A.A.M., and B.A. Boville, 1993:
! Local versus Nonlocal Boundary-Layer Diffusion in a Global Climate
! Model. J. Clim., vol. 6., p. 1825--1842.
!
! Updated by Holtslag and Hack to exclude the surface layer from the
! definition of the boundary layer Richardson number. Ri is now defined
! across the outer layer of the pbl (between the top of the surface
! layer and the pbl top) instead of the full pbl (between the surface and
! the pbl top). For simplicity, the surface layer is assumed to be the
! region below the first model level (otherwise the boundary layer depth 
! determination would require iteration).
!
! NOTE that all calculation in this module is at temperature points (DARLAM)
!------------------------------Code history--------------------------------
!
! Original version:  B. Boville
! Standardized:      J. Rosinski, June 1992
! Reviewed:          B. Boville, P. Rasch, August 1992
! Reviewed:          B. Boville, P. Rasch, April 1996
!
! Modified for boundary layer height diagnosis: Bert Holtslag, june 1994
! >>>>>>>>>  (Use ricr = 0.3 in this formulation)
!
!------------------------------Arguments--------------------------------
!
! Input arguments:u,v,fg,eg,theta,ustar,uav,vav
      
! Input & Output arguments
real, dimension(ifull,kl) :: rkm           ! eddy diffusivity for momentum [m2/s]
real, dimension(ifull,kl) :: rkh           ! eddy diffusivity for heat [m2/s]
real, dimension(ifull,kl) :: theta         ! potential temperature [K]
!     also qg                              ! mixing ratio [kg/kg}

real, dimension(ifull,kl) :: cgh           ! counter-gradient term for heat [K/m]
real, dimension(ifull,kl) :: cgq           ! counter-gradient term for constituents
real, dimension(ifull,kl) :: zg
real ztodtgor,delsig,tmp1,sigotbk,sigotbkm1
real cgs                     ! counter-gradient star (cg/flux)
!
!---------------------------Local parameters----------------------------
!
real, parameter :: tiny=1.e-36             ! lower bound for wind magnitude
!
!---------------------------Local workspace-----------------------------
!

real, dimension(ifull) :: heatv         ! surface virtual heat flux
real, dimension(ifull) :: thvref        ! reference level virtual temperature
real, dimension(ifull) :: phiminv       ! inverse phi function for momentum
real, dimension(ifull) :: phihinv       ! inverse phi function for heat 
real, dimension(ifull) :: wm            ! turbulent velocity scale for momentum
real, dimension(ifull) :: rkhfs         ! surface kinematic heat flux [mK/s]
real, dimension(ifull) :: rkqfs         ! sfc kinematic constituent flux [m/s]
real, dimension(ifull,kl) :: rino       ! bulk Richardson no. from level to ref lev
real, dimension(ifull) :: tlv           ! ref. level pot tmp + tmp excess
real, dimension(ifull) :: wstr          ! w*, convective velocity scale
real, dimension(ifull) :: obklen        ! Obukhov length
real, dimension(ifull,kl) :: tnhs       ! Non-hydrostatic term represented as temperature adjustment
real tkv                                ! model level potential temperature
real therm                              ! thermal virtual temperature excess
real pmid                               ! midpoint pressures
real vvk                                ! velocity magnitude squared
real zmzp                               ! level height halfway between zm and zp
real fak1                               ! k*ustar*pblh
real fak2                               ! k*wm*pblh
real fak3                               ! fakn*wstr/wm 
real pblk                               ! level eddy diffusivity for momentum
real pr                                 ! Prandtl number for eddy diffusivities
real zm                                 ! current level height
real zp                                 ! current level height + one level up
real zl                                 ! zmzp / Obukhov length
real zh                                 ! zmzp / pblh      at half levels
real zzh                                ! (1-(zmzp/pblh))**2      at half levels
real rrho                               ! 1./bottom level density (temporary)
real term                               ! intermediate calculation
real fac                                ! interpolation factor

!------------------------------Commons----------------------------------
real, dimension(ifull,kl) :: uav,vav
real, dimension(ifull) :: cgmap

real, parameter :: c1     = 0.61
real, parameter :: betam  = 15.0  ! Constant in wind gradient expression
real, parameter :: betas  = 5.0   ! Constant in surface layer gradient expression
real, parameter :: betah  = 15.0  ! Constant in temperature gradient expression
real, parameter :: fak    = 8.5   ! Constant in surface temperature excess
real, parameter :: fakn   = 7.2   ! Constant in turbulent prandtl number
real, parameter :: ricr   = 0.25  ! Critical richardson number
real, parameter :: sffrac = 0.1   ! Surface layer fraction of boundary layer
real, parameter :: vk     = 0.4   ! Von Karman's constant
real, parameter :: zkmin = 0.01   ! Minimum kneutral*f(ri)
real ccon    ! fak * sffrac * vk
real binm    ! betam * sffrac
real binh    ! betah * sffrac
real rkmin ! minimum eddy coeffs based on Hourdin et. al. (2001)

      
kmax=kl-1
if (nlocal==6) then
  fac=10.
else
  fac=100.
end if
binh   = betah*sffrac
binm   = betam*sffrac
ccon   = fak*sffrac*vk
!****************************************************************
zg(1:ifull,1) = bet(1)*t(1:ifull,1)/grav
do k = 2,kl
  zg(1:ifull,k) = zg(1:ifull,k-1)+(bet(k)*t(1:ifull,k)+betm(k)*t(1:ifull,k-1))/grav
enddo         ! k  loop
! Non-hydrostatic terms
zg(:,:) = zg(:,:)+phi_nh(:,:)/grav
tnhs(:,1) = phi_nh(:,1)/bet(1)
do k = 2,kl
  ! representing non-hydrostatic term as a correction to air temperature
  tnhs(:,k) = (phi_nh(:,k)-phi_nh(:,k-1)-betm(k)*tnhs(:,k-1))/bet(k)
end do
cgh(:,:) = 0.   ! 3D
cgq(:,:) = 0.   ! 3D
if ( ktau==1 .and. myid==0 ) then
  write(6,*) 'in pbldif nrkmin,npblmin: ',nrkmin,npblmin 
end if
      
! Compute kinematic surface fluxes
do iq=1,ifull
  pmid=ps(iq)*sigmh(1) 
  rrho = rdry*t(iq,1)/pmid
  ustar(iq) = max(ustar(iq),0.01)
  rkhfs(iq) = fg(iq)*rrho/cp           !khfs=w'theta'
  rkqfs(iq) = eg(iq)*rrho/hl           !kqfs=w'q'

  ! Compute various arrays for use later:

  thvref(iq) = theta(iq,1)*(1.0 + 0.61*qg(iq,1))
  heatv(iq)  = rkhfs(iq) + 0.61*theta(iq,1)*rkqfs(iq)
  wm(iq)     = 0.
  ! obklen at t point
  obklen(iq) = -thvref(iq)*ustar(iq)**3/(grav*vk*(heatv(iq) + sign(1.e-10,heatv(iq))))
         
  ! >>>> Define first a new factor fac=100 for use in Richarson number
  !      Calculate virtual potential temperature first level
  !      and initialize pbl height to z1 i.e  1st full level

  pblh(iq) = zg(iq,1)    
  rino(iq,1) = 0.
enddo

! PBL height calculation:
! Search for level of pbl. Scan upward until the Richardson number between
! the first level and the current level exceeds the "critical" value.
! Richardson no. is computed using eq. (4.d.18) NCAR technical report, CCM3)

iflag(:)=0
do k=2,kmax
  do iq=1,ifull
    vvk = (uav(iq,k) - uav(iq,1))**2 + (vav(iq,k) - vav(iq,1))**2 + fac*ustar(iq)**2
    tkv = theta(iq,k)*(1. + 0.61*qg(iq,k))
    rino(iq,k) = grav*(tkv - thvref(iq))*(zg(iq,k)-zg(iq,1))/max(thvref(iq)*vvk,tiny)
    if(rino(iq,k)>=ricr.and.iflag(iq)==0)then
      pblh(iq) = zg(iq,k-1) + (ricr - rino(iq,k-1))/(rino(iq,k) - rino(iq,k-1))*(zg(iq,k) - zg(iq,k-1))
      iflag(iq)=1
    endif  ! (rino(iq,k)>=ricr.and.iflag(iq)==0)
  enddo  ! iq loop
enddo   ! k loop
if(nmaxpr==1.and.mydiag)then
  write (6,"('zg',9f8.1)") zg(idjd,1:kmax)
  write (6,"('rino_pa',9f8.3)") rino(idjd,1:kmax)
endif

! Set pbl height to maximum value where computation exceeds number of
! layers allowed
 
! Improve estimate of pbl height for the unstable points.
! Find unstable points (virtual heat flux is positive):
!
phiminv=0. ! MJT bug fix for IBM compiler
tlv=thvref

do iq=1,ifull
  if(heatv(iq)>0.)then  ! unstable case
    phiminv(iq) = (1. - binm*pblh(iq)/obklen(iq))**(1./3.)
    wm(iq)= ustar(iq)*phiminv(iq)
    ! therm: 2nd term in eq. (4.d.19):
    ! temperature excess due to convective thermal
    therm = heatv(iq)*fak/wm(iq)
    ! eq. (4.d.19) : tlv then used in eq. (4.d.18) to improve pblh
    tlv(iq) = thvref(iq) + therm
  end if  
end do

! Improve pblh estimate for unstable conditions using the
! convective temperature excess:

do k=2,kmax
  do iq=1,ifull
    vvk = (uav(iq,k) - uav(iq,1))**2 + (vav(iq,k) - vav(iq,1))**2 + fac*ustar(iq)**2
    vvk = max(vvk,tiny)
    tkv = theta(iq,k)*(1. + 0.61*qg(iq,k))
    rino(iq,k) = grav*(tkv - tlv(iq))*(zg(iq,k)-zg(iq,1))/max(thvref(iq)*vvk,tiny)     ! (see (4.d.18)
  enddo  !  i loop
enddo   !  k loop

iflag(:)=0
do k=2,kmax
  do iq=1,ifull
    if(heatv(iq)>0..and.iflag(iq)==0)then  ! unstable case
      pblh(iq) = zg(iq,kl)    ! large default for unstable case
      if(rino(iq,k)>=ricr)then
        pblh(iq) = zg(iq,k-1) + (ricr - rino(iq,k-1))/(rino(iq,k) - rino(iq,k-1))*(zg(iq,k) - zg(iq,k-1))
        iflag(iq)=1  ! i.e. found it
      endif  ! (rino(iq,k)>=ricr)
    endif    ! (heatv(iq)>0..and.iflag(iq)==0)
  enddo     ! i loop
enddo      ! k loop

! Points for which pblh exceeds number of pbl layers allowed;
! set to maximum
 
! PBL height must be greater than some minimum mechanical mixing depth
! Several investigators have proposed minimum mechanical mixing depth
! relationships as a function of the local friction velocity, u*.  We 
! make use of a linear relationship of the form h = c u* where c=700.
! The scaling arguments that give rise to this relationship most often 
! represent the coefficient c as some constant over the local coriolis
! parameter.  Here we make use of the experimental results of Koracin 
! and Berkowicz (1988) [BLM, Vol 43] for which they recommend 0.07/f
! where f was evaluated at 39.5 N and 52 N.  Thus we use a typical mid
! latitude value for f so that c = 0.07/f = 700.
 
if(npblmin==1)pblh(:) = max(pblh(:),min(200.,700.*ustar(:)))
if(npblmin==2)pblh(:) = max(pblh(:),.07*ustar(:)/max(.5e-4,abs(f(1:ifull))))
if(npblmin==3)pblh(:) = max(pblh(:),.07*ustar(:)/max(1.e-4,abs(f(1:ifull)))) ! to ~agree 39.5N
if(npblmin==4)pblh(1:ifull) = max(pblh(1:ifull),50.)

! pblh is now available; do preparation for diffusivity calculation:

! Do additional preparation for unstable cases only, set temperature
! and moisture perturbations depending on stability.

do iq=1,ifull
  if(heatv(iq)>0.)then  ! unstable case
    phiminv(iq) =     (1. - binm*pblh(iq)/obklen(iq))**(1./3.)
    phihinv(iq) = sqrt(1. - binh*pblh(iq)/obklen(iq))
    wm(iq)      = ustar(iq)*phiminv(iq)
    wstr(iq)    = (heatv(iq)*grav*pblh(iq)/thvref(iq))**(1./3.)
  end if
end do

! Main level loop to compute the diffusivities and 
! counter-gradient terms:

if(nlocal==3)then
  ! suppress nonlocal scheme over the sea   jlm
  do iq=1,ifull
    if(.not.land(iq))pblh(iq)=0.
  enddo
endif  !  (nlocal==3)

if(nlocal==4)then
  ! suppress nonlocal scheme for column if cloudy layer in pbl   jlm
  do k=1,kl/2
    do iq=1,ifull
      if(zg(iq,k)<pblh(iq).and.cfrac(iq,k)>0.)pblh(iq)=0.
    enddo
  enddo
endif  !  (nlocal==4)

if(nlocal==5)then
  ! suppress nonlocal scheme for column if cloudy layer in pbl   jlm
  ! restores pblh at the bottom to have it available in convjlm/vertmix
  do k=1,kl/2
    do iq=1,ifull
      if(zg(iq,k)<pblh(iq).and.cfrac(iq,k)>0.) pblh(iq)=-pblh(iq)  
    enddo
  enddo
endif  !  (nlocal==5)

if(nlocal==2)then
  do k=1,kmax-1        
    ! suppress nonlocal scheme if cloudy layers in pbl   jlm
    ! note this allows layers below to be done as nonlocal
    do iq=1,ifull
      if(zg(iq,k)<pblh(iq).and.cfrac(iq,k)>0.)pblh(iq)=0.
    enddo
  end do
endif  !  (nlocal==2)

! Find levels within boundary layer:
! This is where Kh is at half model levels 
! zmzp = 0.5*(zm + zp)

do k=1,kmax-1
  do iq=1,ifull
    fak1 = ustar(iq)*pblh(iq)*vk
    zm = zg(iq,k)
    zp = zg(iq,k+1)
    if (zm < pblh(iq)) then
      zmzp = 0.5*(zm + zp)
      zh = zmzp/pblh(iq)
      zl = zmzp/obklen(iq)
      zzh= 0.
      if (zh<=1.0) zzh = (1. - zh)**2

! stblev for points zm < plbh and stable and neutral
! unslev for points zm < plbh and unstable
      if(heatv(iq)>0.)then  ! unstable case
        fak2   = wm(iq)*pblh(iq)*vk
! unssrf, unstable within surface layer of pbl
! Unstable for surface layer; counter-gradient terms zero
        if (zh<sffrac) then
          term =     (1. - betam*zl)**(1./3.)
          pblk = fak1*zh*zzh*term
          pr = term/sqrt(1. - betah*zl)
        else
! unsout, unstable within outer   layer of pbl
! Unstable for outer layer; counter-gradient terms non-zero:
          pblk = fak2*zh*zzh
          fak3 = fakn*wstr(iq)/wm(iq)
          cgs     = fak3/(pblh(iq)*wm(iq))
          cgh(iq,k) = rkhfs(iq)*cgs                 !eq. (4.d.17)
          cgq(iq,k) = rkqfs(iq)*cgs                 !eq. (4.d.17)
          pr = phiminv(iq)/phihinv(iq) + ccon*fak3/fak
        end if
        rkm(iq,k) = max(pblk,rkm(iq,k))
        rkh(iq,k) = max(pblk/pr,rkh(iq,k))
      elseif(nlocal>0)then    ! following are stable or neutral
! Stable and neutral points; set diffusivities; counter-gradient
! terms zero for stable case:
! term: pblk is Kc in eq. (4.d.16)
! but reverts to Louis stable treatment for nlocal=-1
        if (zl<=1.) then   ! 0 < z/L < 1.
          pblk = fak1*zh*zzh/(1. + betas*zl)
        else
          pblk = fak1*zh*zzh/(betas + zl)
        endif
        if(nrkmin==2)rkmin=vk*ustar(iq)*zmzp*zzh
        if(nrkmin==3)rkmin=max(rkh(iq,k),vk*ustar(iq)*zmzp*zzh)
        if(nrkmin==1.or.nlocal==6)rkmin=rkh(iq,k)
        if(ntest==1.and.mydiag)then
          if(iq==idjd)then
            write(6,*) 'in pbldif k,ustar,zmzp,zh,zl,zzh ',k,ustar(iq),zmzp,zh,zl,zzh
            write(6,*) 'rkh_L,rkmin,pblk,fak1,pblh ',rkh(iq,k),rkmin,pblk,fak1,pblh(iq)
          endif  ! (iq==idjd)
        endif    ! (ntest==1)
        rkm(iq,k) = max(pblk,rkmin)        
        rkh(iq,k) = rkm(iq,k)
      endif      ! (heatv(iq)>0.)    unstbl(i)
    endif        ! zm < pblh(iq)
  enddo         ! iq=1,ifull
enddo             !end of k loop
if(diag.and.mydiag)then
  if(heatv(idjd)>0.) write (6,"('rino_pb',9f8.3)") rino(idjd,1:kmax) ! not meaningful or used otherwise
  write (6,*) 'ricr,obklen,heatv,pblh ',ricr,obklen(idjd),heatv(idjd),pblh(idjd)
  write (6,"('rkh_p',9f9.3/5x,9f9.3)") rkh(idjd,1:kl-2)
endif


! turn off CG term for small grid spacing
do k=1,kl
  cgh(:,k)=cgh(:,k)*cgmap(:)
  cgq(:,k)=cgq(:,k)*cgmap(:)
end do


ztodtgor = dtin*grav/rdry
!     update theta and qtg due to counter gradient
do k=2,kmax-1
  do iq=1,ifull
    delsig = sigmh(k+1)-sigmh(k)
    tmp1 = ztodtgor/delsig
    sigotbk=sigmh(k+1)/(0.5*(t(iq,k+1) + t(iq,k) + tnhs(iq,k+1) + tnhs(iq,k)))
    sigotbkm1=sigmh(k)/(0.5*(t(iq,k-1) + t(iq,k) + tnhs(iq,k-1) + tnhs(iq,k)))
    theta(iq,k) = theta(iq,k) + tmp1*(sigotbk*rkh(iq,k)*cgh(iq,k) - sigotbkm1*rkh(iq,k-1)*cgh(iq,k-1))
    qg(iq,k) = qg(iq,k) + tmp1*(sigotbk*rkh(iq,k)*cgq(iq,k) - sigotbkm1*rkh(iq,k-1)*cgq(iq,k-1))
  end do
end do
k=1
do iq=1,ifull
  delsig = sigmh(k+1)-sigmh(k)
  tmp1 = ztodtgor/delsig
  sigotbk=sigmh(k+1)/(0.5*(t(iq,k+1) + t(iq,k) + tnhs(iq,k+1) + tnhs(iq,k)))
  theta(iq,k) = theta(iq,k) + tmp1*sigotbk*rkh(iq,k)*cgh(iq,k)
  qg(iq,k) = qg(iq,k) + tmp1*sigotbk*rkh(iq,k)*cgq(iq,k)
end do

if (ntest>0.and.mydiag) then
  write(6,*) 'pbldif'
  write(6,*) 'rkh= ',(rkh(idjd,k),k=1,kl)
  write(6,*) 'theta= ',(theta(idjd,k),k=1,kl)
  write(6,*) 'qg= ',(qg(idjd,k),k=1,kl)
  write(6,*) 'cgh= ',(cgh(idjd,k),k=1,kl)
  write(6,*) 'cgq= ',(cgq(idjd,k),k=1,kl)
endif
!
!    Check for neg qtg's and put the original vertical
!    profile back if a neg value is found. A neg value implies that the
!    quasi-equilibrium conditions assumed for the countergradient term are
!    strongly violated.
!    Original code rewritten by Rosinski 7/8/91 to vectorize in longitude.

if(nlocal==5)then
  ! restoring pblh to have it available in convjlm/vertmix jlm
  pblh(:)=abs(pblh(:))
endif  !  (nlocal==5)

return
end subroutine pbldif
