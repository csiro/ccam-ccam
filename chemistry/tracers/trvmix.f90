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
    
module trvmix

private
public tracervmix

contains
    
subroutine tracervmix

use arrays_m                        ! Atmosphere dyamics prognostic arrays
use carbpools_m                     ! Carbon pools
use const_phys                      ! Physical constants
use morepbl_m                       ! Additional boundary layer diagnostics
use newmpar_m                       ! Grid parameters
use parm_m, only : idjd, dt         ! Model configuration
use pbl_m                           ! Boundary layer arrays
use sigs_m                          ! Atmosphere sigma levels
use tracermodule                    ! Tracer routines
use tracers_m                       ! Tracer data

implicit none

integer tile, is, ie
integer iq, k
real, dimension(imax,numtracer) :: lco2em
real, dimension(imax,kl,ntrac) :: ltr
real, dimension(imax,kl) :: loh, lstrloss, ljmcf
real, dimension(imax,kl) :: lt
real, dimension(imax,kl) :: lat, lct
real, dimension(imax,kl) :: lrkhsave
real, dimension(imax) :: lfnee, lfpn, lfrp, lfrs, lmcfdep
real tmnht, dz, gt, rlogs1, rlogs2, rlogh1, rlog12, rong

!$omp do schedule(static) private(is,ie,iq,k),      &
!$omp private(lt,lat,lct),                          &
!$omp private(ltr,lco2em,loh,lstrloss,ljmcf),       &
!$omp private(lrkhsave,rong,rlogs1,rlogs2),         &
!$omp private(rlogh1,rlog12,tmnht,dz,gt,lfnee),     &
!$omp private(lfpn,lfrp,lfrs,lmcfdep)
do tile = 1,ntiles
  is = (tile-1)*imax + 1
  ie = tile*imax

  lt       = t(is:ie,:)
  lrkhsave = rkhsave(is:ie,:)
  
  rong = rdry/grav
  lat(:,1) = 0.
  lct(:,kl) = 0.
  rlogs1=log(sig(1))
  rlogs2=log(sig(2))
  rlogh1=log(sigmh(2))
  rlog12=1./(rlogs1-rlogs2)
  do iq = 1,imax
    tmnht = (lt(iq,2)*rlogs1-lt(iq,1)*rlogs2+(lt(iq,1)-lt(iq,2))*rlogh1)*rlog12
    dz = -tmnht*rong*((sig(2)-sig(1))/sigmh(2))  ! this is z(k+1)-z(k)
    gt = lrkhsave(iq,1)*dt*(sig(2)-sig(1))/(dz**2)
    lat(iq,2) = -gt/dsig(2)  
    lct(iq,1) = -gt/dsig(1)
  end do
  do k = 2,kl-1
    do iq = 1,imax
      ! Calculate half level heights and temperatures
      ! n.b. an approximate zh (in m) is quite adequate for this routine
      tmnht = ratha(k)*lt(iq,k+1) + rathb(k)*lt(iq,k)
      dz = -tmnht*rong*((sig(k+1)-sig(k))/sigmh(k+1))  ! this is z(k+1)-z(k)
      gt = lrkhsave(iq,k)*dt*(sig(k+1)-sig(k))/(dz**2)
      lat(iq,k+1) = -gt/dsig(k+1)  
      lct(iq,k) = -gt/dsig(k)
    end do
  end do
  
  ltr = tr(is:ie,:,:)
  lco2em = co2em(is:ie,:)
  if ( allocated(oh) ) loh = oh(is:ie,:)
  if ( allocated(strloss) ) lstrloss = strloss(is:ie,:)
  if ( allocated(jmcf) ) ljmcf = jmcf(is:ie,:)
  if ( allocated(fnee) ) lfnee = fnee(is:ie)
  if ( allocated(fpn) ) lfpn = fpn(is:ie)
  if ( allocated(frp) ) lfrp = frp(is:ie)
  if ( allocated(frs) ) lfrp = frs(is:ie)
  if ( allocated(mcfdep) ) lmcfdep = mcfdep(is:ie)
  lt = t(is:ie,:)
  ! Tracers
  call tracervmix_work(lat,lct,lt,ps(is:ie),cdtq(is:ie),ltr,lfnee,lfpn,       &
                       lfrp,lfrs,lco2em,loh,lstrloss,ljmcf,lmcfdep,tile,imax)
  tr(is:ie,:,:) = ltr

end do ! tile = 1,ntiles
!$omp end do nowait

return
end subroutine tracervmix
    
! ***************************************************************************
! Tracer emission, deposition, settling and turbulent mixing routines
subroutine tracervmix_work(at,ct,t,ps,cdtq,tr,fnee,fpn,frp,frs,co2em,oh,strloss,jmcf,mcfdep,tile,imax)

use const_phys
use diag_m
use newmpar_m, only : kl
use parm_m
use sigs_m
use tracermodule, only : tracunit,numtracer
use tracers_m, only : tracname,ntrac,ngas

implicit none

integer, intent(in) :: tile,imax
integer igas, k
real, dimension(imax,kl) :: updtr
real, intent(in), dimension(imax,kl) :: at, ct
real, dimension(imax,kl), intent(in) :: t
real, dimension(imax), intent(in) :: ps
real, dimension(imax), intent(in) :: cdtq
real, dimension(imax,kl,ntrac), intent(inout) :: tr
real, dimension(imax), intent(in) :: fnee
real, dimension(imax), intent(in) :: fpn
real, dimension(imax), intent(in) :: frp
real, dimension(imax), intent(in) :: frs
real, dimension(imax,numtracer), intent(in) :: co2em
real, dimension(imax,kl), intent(in) :: oh
real, dimension(imax,kl), intent(in) :: strloss
real, dimension(imax,kl), intent(in) :: jmcf
real, dimension(imax), intent(in) :: mcfdep
real, dimension(imax,kl) :: prf, dz, rhoa
real, dimension(imax,kl) :: trsrc
real molfact, radfact, co2fact, gasfact
logical decay, methloss, mcfloss

! Setup
molfact = 1000.*fair_molm          ! factor for units in mol/m2/s
co2fact = 1000.*fair_molm/fc_molm
radfact = 1.293                    ! test factor for radon units in Bq/m2/s, conc in Bq/m3

do k=1,kl
  dz(:,k) = -rdry*dsig(k)*t(1:imax,k)/(grav*sig(k))
  rhoa(:,k) = ps(1:imax)*sig(k)/(rdry*t(1:imax,k)) ! density of air (kg/m**3)
  prf(:,k) = ps(1:imax)*sig(k)
end do

! Tracer settling
call trsettling(rhoa,t,dz,prf,tr,imax)

do igas=1,ngas                  

  ! Tracer emission
  call trgassflux(igas,trsrc,fnee,fpn,frp,frs,co2em,tile,imax)
  
  ! change gasfact to be depend on tracer flux units
  if (trim(tracunit(igas))=='gC/m2/s') then
    gasfact = co2fact
    decay = .false.
  elseif (trim(tracunit(igas))=='mol/m2/s') then
    gasfact = molfact
    decay = .false.
  elseif (trim(tracunit(igas))=='Bq/m2/s') then
    gasfact = radfact
    decay = .true.
  else
!   assume no surface flux so gasfact could be anything but we'll 
!   set it to zero
    gasfact = 0.
    decay = .false.
  endif
  
  ! also set decay for tracer name 'radon' in case not in Bq/m2/s
  if ( trim(tracname(igas))=='radon' .or. tracname(igas)(1:2)=='Rn' ) then
    decay = .true.
  end if
  
  methloss = tracname(igas)(1:7)=='methane'  ! check for methane tracers to set flag to do loss
  mcfloss  = tracname(igas)(1:3)=='mcf'      ! check for mcf tracers to set flag to do loss

  ! deposition and decay terms
  call gasvmix(updtr,gasfact,igas,decay,trsrc,methloss,mcfloss,cdtq,dz(:,1), &
               tr,t,oh,strloss,jmcf,mcfdep,ps,imax)
  
  call trimt(at,ct,updtr,imax)
  tr(1:imax,:,igas) = updtr  
  
end do

return
end subroutine tracervmix_work

! ***************************************************************************
!     this routine put the correct tracer surface flux into trsrc
subroutine trgassflux(igas,trsrc,fnee,fpn,frp,frs,co2em,tile,imax)

use cable_ccam, only : cbmemiss
use dates_m
use newmpar_m, only : kl
use tracermodule, only : tracdaytime,traclevel,numtracer
use tracers_m, only : tracname,tractype

implicit none

integer, intent(in) :: tile,imax
integer igas, ierr, k
real, dimension(imax,kl), intent(out) :: trsrc
real, dimension(imax), intent(in) :: fnee
real, dimension(imax), intent(in) :: fpn
real, dimension(imax), intent(in) :: frp
real, dimension(imax), intent(in) :: frs
real, dimension(imax,numtracer), intent(in) :: co2em
integer nchar, mveg

!     initialise (to allow for ocean gridpoints for cbm fluxes)      
!     and non surface layers
trsrc = 0.

select case(trim(tractype(igas)))
    
  case('online')
    if (trim(tracname(igas)(1:3))=='cbm') then
      select case (trim(tracname(igas)))
        case('cbmnep'); trsrc(:,1) = fnee
        case('cbmpn');  trsrc(:,1) = fpn
        case('cbmrp');  trsrc(:,1) = frp
        case('cbmrs');  trsrc(:,1) = frs
        case default;   stop 'unknown online tracer name'
      end select
    else
      nchar = len_trim(tracname(igas))
      read(tracname(igas)(nchar-1:nchar),'(i2)',iostat=ierr) mveg
      if (ierr/=0) then
        write(6,*) 'unknown online tracer name or veg type number'
        write(6,*) trim(tracname(igas)),ierr
        stop
      end if
      !if (mveg<1.or.mveg>maxval(ivegt)) stop 'tracer selection: veg type out of range'
      select case (tracname(igas)(1:nchar-2))
        case('gpp');    call cbmemiss(trsrc(:,1),mveg,1,tile,imax)
        case('plresp'); call cbmemiss(trsrc(:,1),mveg,2,tile,imax)
        case('slresp'); call cbmemiss(trsrc(:,1),mveg,3,tile,imax)
        case default;   stop 'unknown online tracer name'
      end select
    endif
    
  case ('daypulseon')
    ! only add flux during day time
    if (tracdaytime(igas,1)<tracdaytime(igas,2) .and. tracdaytime(igas,1)<=timeg .and. tracdaytime(igas,2)>=timeg) then
      trsrc(:,1) = co2em(:,igas)
    elseif (tracdaytime(igas,1)>tracdaytime(igas,2) .and. (tracdaytime(igas,1)<=timeg .or. tracdaytime(igas,2)>=timeg)) then
      trsrc(:,1) = co2em(:,igas)
    else
      trsrc(:,1) = 0.
    endif
    
  case default
    ! emissions from file over levels
    do k=1,traclevel(igas)
      trsrc(:,k) = co2em(:,igas)/real(traclevel(igas))
    end do
    
end select

return
end subroutine trgassflux

! *****************************************************************
subroutine gasvmix(temptr,fluxfact,igas,decay,trsrc,methloss,mcfloss,vt,dz1, &
                   tr,t,oh,strloss,jmcf,mcfdep,ps,imax)

use const_phys
use newmpar_m, only : kl
use parm_m
use sigs_m 
use tracermodule, only : trdep
use tracers_m, only : ntrac

implicit none

integer, intent(in) :: imax
integer, intent(in) :: igas
integer k
real, dimension(imax,kl), intent(out) :: temptr
real, dimension(imax,kl) :: loss
real, dimension(imax,kl), intent(in) :: trsrc
real, dimension(imax), intent(in) :: vt, dz1
real, dimension(imax) :: dep
real, intent(in) :: fluxfact
real drate
logical, intent(in) :: decay, methloss, mcfloss
real, dimension(imax,kl,ntrac), intent(in) :: tr
real, dimension(imax,kl), intent(in) :: t
real, dimension(imax,kl), intent(in) :: oh
real, dimension(imax,kl), intent(in) :: strloss
real, dimension(imax,kl), intent(in) :: jmcf
real, dimension(imax), intent(in) :: mcfdep
real, dimension(imax), intent(in) :: ps

real, parameter :: koh    = 2.45e-12
real, parameter :: kohmcf = 1.64e-12

! decay rate for radon (using units of source, Bq/m2/s, to
! indicate that radon and need decay
if (decay) then
  drate = exp(-dt*2.11e-6)
else
  drate = 1.
endif

! rml 16/2/10 methane loss by OH and in stratosphere
if (methloss) then
  loss = tr(1:imax,:,igas)*dt*(koh*exp(-1775./t(1:imax,:))*oh(:,:) + strloss(:,:))
  dep = 0.
elseif (mcfloss) then
  loss = tr(1:imax,:,igas)*dt*(kohmcf*exp(-1520./t(1:imax,:))*oh(:,:) + jmcf(:,:))
  ! deposition
  dep  = exp(-mcfdep*dt/dz1)
elseif (trdep(igas)>0.) then
  loss = 0.
  dep  = exp(-vt*dt/dz1)
else
  loss = 0.
  dep  = 1.
endif

! implicit version due to potentially high transfer velocity relative to dz1
temptr(:,1) = tr(1:imax,1,igas)*drate*dep - fluxfact*grav*dt*trsrc(:,1)/(dsig(1)*ps(1:imax)) - loss(:,1)
do k = 2,kl
  temptr(:,k) = tr(1:imax,k,igas)*drate - fluxfact*grav*dt*trsrc(:,k)/(dsig(k)*ps(1:imax)) - loss(:,k)
end do

return
end subroutine gasvmix

! *********************************************************************
!     This is a copy of trim.f but trying to do all tracers at once.  
!     u initially now contains rhs; leaves with answer u (jlm)
!     n.b. we now always assume b = 1-a-c

subroutine trimt(a,c,rhs,imax)

use newmpar_m, only : kl

implicit none

!     N.B.  e, g, temp are just work arrays (not passed through at all)     

integer k
integer, intent(in) :: imax
real, dimension(imax,kl), intent(inout) :: rhs
real, dimension(imax,kl) :: g
real, dimension(imax,kl), intent(in) :: a, c
real, dimension(imax,kl) :: e, temp
real, dimension(imax) :: b

!     this routine solves the system
!       a(k)*u(k-1)+b(k)*u(k)+c(k)*u(k+1)=rhs(k)    for k=2,kl-1
!       with  b(k)*u(k)+c(k)*u(k+1)=rhs(k)          for k=1
!       and   a(k)*u(k-1)+b(k)*u(k)=rhs(k)          for k=kl

!     the Thomas algorithm is used
!     save - only needed if common/work removed

b(:)=1.-a(:,1)-c(:,1)
e(:,1)=c(:,1)/b(:)
do k=2,kl-1
  b(:)=1.-a(:,k)-c(:,k)
  temp(:,k)= 1./(b(:)-a(:,k)*e(:,k-1))
  e(:,k)=c(:,k)*temp(:,k)
enddo

!     use precomputed values of e array when available
b(:)=1.-a(:,1)-c(:,1)
g(:,1)=rhs(:,1)/b(:)
do k=2,kl-1
  g(:,k)=(rhs(:,k)-a(:,k)*g(:,k-1))*temp(:,k)
end do

!     do back substitution to give answer now
b(:)=1.-a(:,kl)-c(:,kl)
rhs(:,kl)=(rhs(:,kl)-a(:,kl)*g(:,kl-1))/(b(:)-a(:,kl)*e(:,kl-1))
do k=kl-1,1,-1
  rhs(:,k)=g(:,k)-e(:,k)*rhs(:,k+1)
end do
      
return
end subroutine trimt

! *********************************************************************
! Calculate settling
! Based on dust settling in aerosolldr.f90
subroutine trsettling(rhoa,tmp,delz,prf,tr,imax)

use const_phys
use newmpar_m, only : kl
use parm_m
use tracermodule, only : trden, trreff
use tracers_m, only : ngas,ntrac

implicit none

integer, intent(in) :: imax
real, dimension(imax,kl), intent(in) :: rhoa    !air density (kg/m3)
real, dimension(:,:), intent(in) :: tmp         !temperature (K)
real, dimension(imax,kl), intent(in) :: delz    !Layer thickness (m)
real, dimension(imax,kl), intent(in) :: prf     !Pressure (hPa)
real, dimension(imax,kl,ntrac), intent(inout) :: tr
real, dimension(imax) :: c_stokes, corr, c_cun
real, dimension(imax) :: newtr, b, dfall
real, dimension(imax,kl) :: vd_cor
integer nt,k

do nt = 1, ngas
    
  if ( trden(nt)>0. .and. trreff(nt)>0. ) then
    
    ! Settling velocity (m/s)  (Stokes Law)
    ! TRDEN       soil class density             (kg/m3)
    ! TRREFF      effective radius               (m)
    ! grav        gravity                        (m/s2)

    ! Solve at the model top
    ! Dynamic viscosity
    C_Stokes = 1.458E-6 * TMP(1:imax,kl)**1.5/(TMP(1:imax,kl)+110.4) 
    ! Cuningham correction
    Corr = 6.6E-8*prf(:,kl)/1013.*TMP(1:imax,kl)/293.15
    C_Cun = 1. + 1.249*corr/trreff(nt)
    ! Settling velocity
    Vd_cor(:,kl) =2./9.*grav*trden(nt)*trreff(nt)**2/C_Stokes*C_Cun
    ! Solve each vertical layer successively (layer l)
    do k = kl-1,1,-1
      ! Dynamic viscosity
      C_Stokes = 1.458E-6*TMP(1:imax,k)**1.5/(TMP(1:imax,k)+110.4) 
      ! Cuningham correction
      Corr = 6.6E-8*prf(:,k)/1013.*TMP(1:imax,k)/293.15
      C_Cun = 1. + 1.249*corr/trreff(nt)
      ! Settling velocity
      Vd_cor(:,k) = 2./9.*grav*trden(nt)*trreff(nt)*trreff(nt)/C_Stokes*C_Cun
    end do
  
    ! Update mixing ratio
    b = dt*VD_cor(:,kl)/DELZ(:,kl)
    newtr = tr(1:imax,kl,nt)*exp(-b)
    newtr = max( newtr, 0. )
    dfall = max( tr(1:imax,kl,nt) - newtr, 0. )
    tr(1:imax,kl,nt) = newtr
    ! Solve each vertical layer successively (layer l)
    do k = kl-1,1,-1
      ! Update mixing ratio
      b = dt*Vd_cor(:,k)/DELZ(:,k)
      dfall = dfall * delz(:,k+1)*rhoa(:,k+1)/(delz(:,k)*rhoa(:,k))
      ! Fout  = 1.-exp(-b)
      ! Fthru = 1.-Fout/b
      newtr = tr(1:imax,k,nt)*exp(-b) + dfall*(1.-exp(-b))/b
      newtr = max( newtr, 0. )
      dfall = max( tr(1:imax,k,nt) + dfall - newtr, 0. )
      tr(1:imax,k,nt) = newtr
    end do
  end if
end do

return
end subroutine trsettling

! ********************************************************************
end module trvmix
