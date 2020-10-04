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
    
! This module is the Rotstayn 1997 cloud microphysics parameterisation

! The scheme has been modifed by MJT for max/rnd cloud overlap and to include prognostic rainfall, snow and
! graupel.  The snow and graupel components are based on Lin et al 1983, with modifications to the slope
! and intercept to be consistent with Rotstayn 97. There is also an optional prognostic cloud fraction option
! based on Tiedtke (see cloudmod.f90).

! ldr    = 0    Diagnosed cloud scheme (depreciated)
! ldr   /= 0    Prognostic cloud condensate (different ice fall speed options)
    
! ncloud = 0    Standard LDR cloud microphysics with water vapour, liquid cloud and ice cloud
! ncloud = 2    Same as ncloud=0, but with prognostic rain and modified cfrac
! ncloud = 3    Same as ncloud=2, but with prognostic graupel and snow, as well as modified cfrac
! ncloud = 4    Use prognostic cloud fraction based on Tiedtke from GFDL-CM3, but autoconversion from ncloud=0
   
!                            Water vapour (qg)
!
!   Cloud water (qlg,cfrac)                      Cloud ice (qfg,cfrac)
!
!   Rain (qrg,rfrac)                             Snow (qsg,sfrac)         Graupel (qgrg,gfrac)

! qg, qlg, qfg, qrg, qsg and qgrg are mixing ratios (g/g) and cfrac, rfrac, sfrac, gfrac are area cover
! fractions
    
module leoncld_mod

private
public leoncld
public rhow, rhoice, um, Dva, rKa
public ti, tice, aa, bb
public Ec, Aurate
public rhosno, Eac
public Nor, rk2, rho0, Ecol
public wlc, wls, ticon
public aice, bice

real, parameter :: maxlintime = 120. ! time-step for Lin et al 83 cloud microphysics

! Physical constants
real, parameter :: rhow=1000.  !Density of water
real, parameter :: rhoice=917. !Density of ice
real, parameter :: um=1.717e-5 !Dynamic viscosity of air (at 0 deg. C)
real, parameter :: Dva=2.21    !Diffusivity of qv in air (0 deg. and 1 Pa)
real, parameter :: rKa=2.4e-2  !Thermal conductivity of air (0 deg)

! Tunable parameters for qcloud scheme
real, parameter :: ti = -40.               ! Min T for liquid water clouds in Celsius
real, parameter :: tice=273.15+ti          !Convert ti to Kelvin
real, parameter :: aa=-2/ti**3, bb=3/ti**2 ! Coeffs for cubic interp of fracice

! The following are used in the Manton-Cotton rain parameterization
real, parameter :: Ec=0.55                 !Mean collection efficiency for cloud drops
real, parameter :: Aurate=0.104*9.80616*Ec/um !Part of rate constant

! Parameters related to snow
real, parameter :: rhosno=100. !Assumed density of snow in kg/m^3
real, parameter :: Eac=0.7     !Mean collection efficiency of ql by snow

! Parameters related to rain
real, parameter :: Nor=8.0e6 !Marshall-Palmer intercept parameter
real, parameter :: rk2=142.  !Fall speed of rain: V(D)=rk2*sqrt(D)*sqrt(rho0/rhoa)
real, parameter :: rho0=1.2  !Standard air density
real, parameter :: Ecol=0.7  !Mean collection efficiency of ql by rain

! Parameters related to diagnosed convective cloud
real, parameter :: wlc=0.2e-3   !LWC of deep conv cloud (kg/m**3)
real, parameter :: wls=0.35e-3  !LWC of shallow (non-preciptating) conv cloud
real, parameter :: ticon=238.15 !Temp at which conv cloud becomes ice

! Parameters related to cloud radiative properties
real, parameter :: aice=1.016 !Constant in Platt optical depth for ice (SI units)
real, parameter :: bice=0.68  !Constant in Platt optical depth for ice (SI units)

interface pow75
  module procedure pow75_s, pow75_v
end interface

contains
    
subroutine leoncld

use aerointerface                 ! Aerosol interface
use arrays_m                      ! Atmosphere dyamics prognostic arrays
use cc_mpi, only : mydiag         ! CC MPI routines
use cc_omp                        ! CC OpenMP routines
use cfrac_m                       ! Cloud fraction
use cloudmod                      ! Prognostic cloud fraction
use const_phys                    ! Physical constants
use kuocomb_m                     ! JLM convection
use liqwpar_m                     ! Cloud water mixing ratios
use map_m                         ! Grid map arrays
use morepbl_m                     ! Additional boundary layer diagnostics
use newmpar_m                     ! Grid parameters
use nharrs_m                      ! Non-hydrostatic atmosphere arrays 
use parm_m, only : idjd, iaero
                                  ! Model configuration
use prec_m                        ! Precipitation
use sigs_m                        ! Atmosphere sigma levels
use soil_m                        ! Soil and surface data
use work3f_m                      ! Grid work arrays
use vvel_m                        ! Additional vertical velocity

implicit none

include 'kuocom.h'                ! Convection parameters

integer tile, is, ie, k
integer idjd_t
real, dimension(imax,kl) :: lcfrac, lgfrac, lppfevap, lppfmelt, lppfprec, lppfsnow
real, dimension(imax,kl) :: lppfstayice, lppfstayliq, lppfsubl, lpplambs, lppmaccr, lppmrate
real, dimension(imax,kl) :: lppqfsedice, lpprfreeze, lpprscav, lqccon, lqfg, lqfrad
real, dimension(imax,kl) :: lqg, lqgrg, lqlg, lqlrad, lqrg, lqsng, lrfrac, lsfrac, lt
real, dimension(imax,kl) :: ldpsldt, lnettend, lstratcloud, lclcon, lcdrop, lrhoa
real, dimension(ifull,kl) :: clcon, cdrop
logical mydiag_t

!$omp do schedule(static) private(is,ie),                                             &
!$omp private(k,lrhoa,lcdrop,lclcon)
do tile = 1,ntiles
  is = (tile-1)*imax + 1
  ie = tile*imax

  ! Calculate droplet concentration from aerosols (for non-convective faction of grid-box)
  do concurrent (k = 1:kl)
    lrhoa(1:imax,k) = ps(is:ie)*sig(k)/(rdry*t(is:ie,k))  
  end do
  call aerodrop(is,lcdrop,lrhoa,outconv=.true.)
  cdrop(is:ie,1:kl) = lcdrop(1:imax,1:kl)

  ! Calculate convective cloud fraction
  call convectivecloudfrac(lclcon,kbsav(is:ie),ktsav(is:ie),condc(is:ie))
  clcon(is:ie,1:kl) = lclcon(1:imax,1:kl)
end do
!$omp end do nowait

!$omp do schedule(static) private(is,ie),                                             &
!$omp private(lcfrac,lgfrac,lrfrac,lsfrac),                                           &
!$omp private(lppfevap,lppfmelt,lppfprec,lppfsnow,lppfstayice,lppfstayliq,lppfsubl),  &
!$omp private(lpplambs,lppmaccr,lppmrate,lppqfsedice,lpprfreeze,lpprscav),            &
!$omp private(lqccon,lqfg,lqfrad,lqg,lqgrg,lqlg,lqlrad,lqrg,lqsng,lt),                &
!$omp private(ldpsldt,lnettend,lstratcloud,lclcon,lcdrop,idjd_t,mydiag_t)
!$acc parallel loop copy(stratcloud,gfrac,rfrac,sfrac,t,qg,qgrg,qlg,qfg,qrg,qsng,      &
!$acc   nettend,condg,conds,condx,precip)                                              &
!$acc copyin(dpsldt,clcon,cdrop,kbsav,ktsav,land,ps,em,clcon,cdrop)                    &
!$acc copyout(cfrac,qlrad,qfrad,qccon,ppfevap,ppfmelt,ppfprec,ppfsnow,ppfstayice,      &
!$acc   ppfstayliq,ppfsubl,pplambs,ppmaccr,ppmrate,ppqfsedice,pprfreeze,pprscav)       &
!$acc present(sig)                                                                     &
!$acc private(lcfrac,lgfrac,lppfevap,lppfmelt,lppfprec,lppfsnow,lppfstayice,           &
!$acc   lppfstayliq,lppfsubl,lpplambs,lppmaccr,lppmrate,lppqfsedice,lpprfreeze,        &
!$acc   lpprscav,lqccon,lqfg,lqfrad,lqg,lqgrg,lqlg,lqlrad,lqrg,lqsng,lrfrac,lsfrac,lt, &
!$acc   ldpsldt,lnettend,lstratcloud,lclcon,lcdrop)
do tile = 1,ntiles
  is = (tile-1)*imax + 1
  ie = tile*imax
  
  idjd_t = mod(idjd-1,imax) + 1
  mydiag_t = ((idjd-1)/imax==tile-1).and.mydiag
  
  lcfrac   = cfrac(is:ie,:)
  lgfrac   = gfrac(is:ie,:)
  lrfrac   = rfrac(is:ie,:)
  lsfrac   = sfrac(is:ie,:)
  lqg      = qg(is:ie,:)
  lqgrg    = qgrg(is:ie,:)
  lqlg     = qlg(is:ie,:)
  lqfg     = qfg(is:ie,:)
  lqrg     = qrg(is:ie,:)
  lqsng    = qsng(is:ie,:)
  lqlrad   = qlrad(is:ie,:)
  lqfrad   = qfrad(is:ie,:)  
  lt       = t(is:ie,:)
  ldpsldt  = dpsldt(is:ie,:)
  lclcon   = clcon(is:ie,:)
  lcdrop   = cdrop(is:ie,:)
  lstratcloud = stratcloud(is:ie,:)
  if ( ncloud>=4 ) then
    lnettend    = nettend(is:ie,:)
  end if

  call leoncld_work(lcfrac,condg(is:ie),conds(is:ie),condx(is:ie),lgfrac,                           &
                    kbsav(is:ie),ktsav(is:ie),land(is:ie),                                          &
                    lppfevap,lppfmelt,lppfprec,lppfsnow,lppfstayice,lppfstayliq,lppfsubl,           &
                    lpplambs,lppmaccr,lppmrate,lppqfsedice,lpprfreeze,lpprscav,precip(is:ie),       &
                    ps(is:ie),lqccon,lqfg,lqfrad,lqg,lqgrg,lqlg,lqlrad,lqrg,lqsng,lrfrac,lsfrac,lt, &
                    ldpsldt,lnettend,lstratcloud,lclcon,lcdrop,em(is:ie),idjd_t,mydiag_t,           &
                    ncloud,nclddia,nevapls,ldr,rcrit_l,rcrit_s,rcm,imax,kl)

  cfrac(is:ie,:) = lcfrac
  gfrac(is:ie,:) = lgfrac
  rfrac(is:ie,:) = lrfrac
  sfrac(is:ie,:) = lsfrac
  qccon(is:ie,:) = lqccon
  qg(is:ie,:)    = lqg
  qlg(is:ie,:)   = lqlg
  qfg(is:ie,:)   = lqfg
  qrg(is:ie,:)   = lqrg
  qsng(is:ie,:)  = lqsng
  qgrg(is:ie,:)  = lqgrg
  qlrad(is:ie,:) = lqlrad
  qfrad(is:ie,:) = lqfrad
  t(is:ie,:)     = lt
  stratcloud(is:ie,:) = lstratcloud
  if ( abs(iaero)>=2 ) then
    ppfevap(is:ie,:)    = lppfevap
    ppfmelt(is:ie,:)    = lppfmelt
    ppfprec(is:ie,:)    = lppfprec
    ppfsnow(is:ie,:)    = lppfsnow
    ppfstayice(is:ie,:) = lppfstayice
    ppfstayliq(is:ie,:) = lppfstayliq
    ppfsubl(is:ie,:)    = lppfsubl
    pplambs(is:ie,:)    = lpplambs
    ppmaccr(is:ie,:)    = lppmaccr
    ppmrate(is:ie,:)    = lppmrate
    ppqfsedice(is:ie,:) = lppqfsedice
    pprfreeze(is:ie,:)  = lpprfreeze
    pprscav(is:ie,:)    = lpprscav
  end if
  if ( ncloud>=4 ) then
    nettend(is:ie,:)    = lnettend
  end if
  
end do
!$acc end parallel
!$omp end do nowait

return
end subroutine leoncld

! This subroutine is the interface for the LDR cloud microphysics
subroutine leoncld_work(cfrac,condg,conds,condx,gfrac,kbsav,ktsav,land,                 &
                        ppfevap,ppfmelt,ppfprec,ppfsnow,ppfstayice,ppfstayliq,ppfsubl,  &
                        pplambs,ppmaccr,ppmrate,ppqfsedice,pprfreeze,pprscav,precip,    &
                        ps,qccon,qfg,qfrad,qg,qgrg,qlg,qlrad,qrg,qsng,rfrac,sfrac,t,    &
                        dpsldt,nettend,stratcloud,clcon,cdrop,em,idjd,mydiag,           &
                        ncloud,nclddia,nevapls,ldr,rcrit_l,rcrit_s,rcm,imax,kl)
!$acc routine vector

use const_phys                    ! Physical constants
use estab                         ! Liquid saturation function
use parm_m, only : iaero, nmaxpr, dt
                                  ! Model configuration
use sigs_m                        ! Atmosphere sigma levels

implicit none

integer, intent(in) :: idjd, ncloud, nclddia, nevapls, ldr
integer, intent(in) :: imax, kl
integer, dimension(imax), intent(in) :: kbsav
integer, dimension(imax), intent(in) :: ktsav
real, dimension(imax,kl), intent(inout) :: cfrac, gfrac, rfrac, sfrac
real, dimension(imax,kl), intent(inout) :: qg, qlg, qfg, qrg, qsng, qgrg
real, dimension(imax,kl), intent(inout) :: qlrad, qfrad
real, dimension(imax,kl), intent(inout) :: t
real, dimension(imax,kl), intent(inout) :: nettend
real, dimension(imax,kl), intent(inout) :: stratcloud, clcon, cdrop
real, dimension(imax,kl), intent(out) :: qccon
real, dimension(imax,kl), intent(out) :: ppfevap
real, dimension(imax,kl), intent(out) :: ppfmelt
real, dimension(imax,kl), intent(out) :: ppfprec
real, dimension(imax,kl), intent(out) :: ppfsnow
real, dimension(imax,kl), intent(out) :: ppfstayice
real, dimension(imax,kl), intent(out) :: ppfstayliq
real, dimension(imax,kl), intent(out) :: ppfsubl
real, dimension(imax,kl), intent(out) :: pplambs
real, dimension(imax,kl), intent(out) :: ppmaccr
real, dimension(imax,kl), intent(out) :: ppmrate
real, dimension(imax,kl), intent(out) :: ppqfsedice
real, dimension(imax,kl), intent(out) :: pprfreeze
real, dimension(imax,kl), intent(out) :: pprscav
real, dimension(imax,kl), intent(in) :: dpsldt
real, dimension(imax), intent(inout) :: condg
real, dimension(imax), intent(inout) :: conds
real, dimension(imax), intent(inout) :: condx
real, dimension(imax), intent(inout) :: precip
real, dimension(imax), intent(in) :: ps
real, dimension(imax), intent(in) :: em
real, intent(in) :: rcrit_l, rcrit_s, rcm
logical, intent(in) :: mydiag
logical, dimension(imax), intent(in) :: land

integer, dimension(imax) :: kbase,ktop          !Bottom and top of convective cloud 
real, dimension(imax,kl) :: prf      !Pressure on full levels (hPa)
real, dimension(imax,kl) :: dprf     !Pressure thickness (hPa)
real, dimension(imax,kl) :: rhoa     !Air density (kg/m3)
real, dimension(imax,kl) :: dz       !Layer thickness (m)
real, dimension(imax,kl) :: ccov     !Cloud cover (may differ from cloud frac if vertically subgrid)
real, dimension(imax,kl) :: qsatg    !Saturation mixing ratio
real, dimension(imax,kl) :: qcl      !Vapour mixing ratio inside convective cloud
real, dimension(imax,kl) :: qenv     !Vapour mixing ratio outside convective cloud
real, dimension(imax,kl) :: tenv     !Temperature outside convective cloud
real, dimension(imax) :: precs                  !Amount of stratiform precipitation in timestep (mm)
real, dimension(imax) :: preci                  !Amount of stratiform snowfall in timestep (mm)
real, dimension(imax) :: precg                  !Amount of stratiform graupel in timestep (mm)
real, dimension(imax) :: wcon                   !Convective cloud water content (in-cloud, prescribed)

integer k, iq
real, dimension(imax,kl) :: qevap, qsubl, qauto, qcoll, qaccr, qaccf
real, dimension(imax,kl) :: fluxr, fluxi, fluxs, fluxg, fluxm, fluxf
real, dimension(imax,kl) :: pqfsedice, pfstayice, pfstayliq, pslopes, prscav
real, dimension(imax) :: prf_temp
real, dimension(imax) :: diag_temp
real invdt, fl

! meterological fields
do k = 1,kl
  prf_temp(:) = ps*sig(k)
  prf(:,k)    = 0.01*prf_temp    !ps is SI units
  dprf(:,k)   = -0.01*ps*dsig(k) !dsig is -ve
  rhoa(:,k)   = prf_temp/(rdry*t(:,k))             ! air density
  qsatg(:,k)  = qsat(prf_temp,t(:,k),imax)         ! saturated mixing ratio
  dz(:,k)     = -rdry*dsig(k)*t(:,k)/(grav*sig(k)) ! level thickness in metres 
  dz(:,k)     = min( max(dz(:,k), 1.), 2.e4 )
end do
 
! default values
kbase(:) = 0  ! default
ktop(:)  = 0  ! default
precs(:) = 0. ! rain
preci(:) = 0. ! snow
precg(:) = 0. ! graupel

!     Set up convective cloud column
where ( ktsav(:)<kl-1 )
  ktop(:)  = ktsav(:)
  kbase(:) = kbsav(:) + 1
  wcon(:)  = wlc
elsewhere
  wcon(:)  = 0.
end where


#ifndef GPU
if ( nmaxpr==1 .and. mydiag ) then
  !if ( ktau==1 ) then
  !  write(6,*)'in leoncloud Rcm ',Rcm
  !end if
  write(6,*) 'entering leoncld'
  diag_temp(:) = qg(idjd,:)
  write(6,"('qv  ',9f8.3/4x,9f8.3)") diag_temp(:)
  diag_temp(:) = qfg(idjd,:)
  write(6,"('qf  ',9f8.3/4x,9f8.3)") diag_temp(:)
  diag_temp(:) = qlg(idjd,:)
  write(6,"('ql  ',9f8.3/4x,9f8.3)") diag_temp(:)
  diag_temp(:) = qrg(idjd,:)
  write(6,"('qr  ',9f8.3/4x,9f8.3)") diag_temp(:)
  diag_temp(:) = qsng(idjd,:)
  write(6,"('qs  ',9f8.3/4x,9f8.3)") diag_temp(:)
  diag_temp(:) = qgrg(idjd,:) 
  write(6,"('qg  ',9f8.3/4x,9f8.3)") diag_temp(:)
endif
#endif


! Calculate convective cloud fraction and adjust moisture variables before calling newcloud
do concurrent (k = 1:kl)
  where ( clcon(:,k)>0. )  
    !ccw=wcon(iq)/rhoa(iq,k)  !In-cloud l.w. mixing ratio
    qccon(:,k)  = clcon(:,k)*wcon(:)/rhoa(:,k)  
    qenv(:,k)   = max( 1.e-8, (qg(:,k)-clcon(:,k)*max(qsatg(:,k),qg(:,k)))/(1.-clcon(:,k)) )
    qcl(:,k)    = (qg(:,k)-(1.-clcon(:,k))*qenv(:,k))/clcon(:,k)
    qlg(:,k)    = qlg(:,k)/(1.-clcon(:,k))  
    qfg(:,k)    = qfg(:,k)/(1.-clcon(:,k))  
    stratcloud(:,k) = stratcloud(:,k)/(1.-clcon(:,k)) 
  elsewhere
    qccon(:,k)  = 0.  
    qcl(:,k)    = qg(:,k)
    qenv(:,k)   = qg(:,k)
  end where    
  tenv(:,k)   = t(:,k) ! Assume T is the same in and out of convective cloud
end do


#ifndef GPU
if ( nmaxpr==1 .and. mydiag ) then
  write(6,*) 'before newcloud'
  diag_temp(:) = t(idjd,:)
  write(6,"('t   ',9f8.2/4x,9f8.2)") diag_temp
  diag_temp(:) = qg(idjd,:)
  write(6,"('qv  ',9f8.3/4x,9f8.3)") diag_temp
  diag_temp(:) = qfg(idjd,:)
  write(6,"('qf  ',9f8.3/4x,9f8.3)") diag_temp
  diag_temp(:) = qlg(idjd,:)
  write(6,"('ql  ',9f8.3/4x,9f8.3)") diag_temp
  diag_temp(:) = qrg(idjd,:)
  write(6,"('qr  ',9f8.3/4x,9f8.3)") diag_temp
  diag_temp(:) = qsng(idjd,:)
  write(6,"('qs  ',9f8.3/4x,9f8.3)") diag_temp
  diag_temp(:) = qgrg(idjd,:)
  write(6,"('qg  ',9f8.3/4x,9f8.3)") diag_temp
  diag_temp(:) = qenv(idjd,:)
  write(6,"('qnv ',9f8.3/4x,9f8.3)") diag_temp
  diag_temp(:) = qsatg(idjd,:)
  write(6,"('qsat',9f8.3/4x,9f8.3)") diag_temp
  diag_temp(:) = qcl(idjd,:)
  write(6,"('qcl ',9f8.3/4x,9f8.3)") diag_temp
  diag_temp(:) = clcon(idjd,:)
  write(6,"('clc ',9f8.3/4x,9f8.3)") diag_temp
  write(6,*) 'kbase,ktop ',kbase(idjd),ktop(idjd)
endif
#endif


!     Calculate cloud fraction and cloud water mixing ratios
call newcloud(dt,land,prf,rhoa,tenv,qenv,qlg,qfg,       &
              dpsldt,nettend,stratcloud,em,idjd,mydiag, &
              ncloud,nclddia,rcrit_l,rcrit_s,imax,kl)


! Vertically sub-grid cloud
ccov(1:imax,1:kl) = stratcloud(1:imax,1:kl)
do concurrent (k = 2:kl-1)
  where ( stratcloud(:,k-1)<1.e-10 .and. stratcloud(:,k)>1.e-2 .and. stratcloud(:,k+1)<1.e-10 )
    ccov(:,k) = sqrt(stratcloud(:,k))
  end where
end do
     

#ifndef GPU
if ( nmaxpr==1 .and. mydiag ) then
  write(6,*) 'after newcloud'
  diag_temp(:) = tenv(idjd,:)
  write (6,"('tnv ',9f8.2/4x,9f8.2)") diag_temp
  diag_temp(:) = qg(idjd,:) 
  write (6,"('qv0 ',9f8.3/4x,9f8.3)") diag_temp
  diag_temp(:) = qfg(idjd,:)
  write (6,"('qf  ',9f8.3/4x,9f8.3)") diag_temp
  diag_temp(:) = qlg(idjd,:)
  write (6,"('ql  ',9f8.3/4x,9f8.3)") diag_temp
  diag_temp(:) = qrg(idjd,:)
  write (6,"('qr  ',9f8.3/4x,9f8.3)") diag_temp
  diag_temp(:) = qsng(idjd,:)
  write (6,"('qs  ',9f8.3/4x,9f8.3)") diag_temp
  diag_temp(:) = qgrg(idjd,:)
  write (6,"('qg  ',9f8.3/4x,9f8.3)") diag_temp
  diag_temp(:) = qenv(idjd,:) ! really new qg
  write (6,"('qnv ',9f8.3/4x,9f8.3)") diag_temp
endif
#endif


!     Weight output variables according to non-convective fraction of grid-box            
do concurrent (k = 1:kl)
  t(:,k)  = clcon(:,k)*t(:,k) + (1.-clcon(:,k))*tenv(:,k)
  qg(:,k) = clcon(:,k)*qcl(:,k) + (1.-clcon(:,k))*qenv(:,k)
  where ( k>=kbase(:) .and. k<=ktop(:) )
    stratcloud(:,k) = stratcloud(:,k)*(1.-clcon(:,k))
    ccov(:,k) = ccov(:,k)*(1.-clcon(:,k))              
    qlg(:,k)  = qlg(:,k)*(1.-clcon(:,k))
    qfg(:,k)  = qfg(:,k)*(1.-clcon(:,k))
  end where  
end do


#ifndef GPU
if ( nmaxpr==1 .and. mydiag ) then
  write(6,*) 'before newsnowrain'
  diag_temp(:) = t(idjd,:)
  write (6,"('t   ',9f8.2/4x,9f8.2)") diag_temp
  diag_temp(:) = qg(idjd,:)
  write (6,"('qv  ',9f8.3/4x,9f8.3)") diag_temp
  diag_temp(:) = qfg(idjd,:)
  write (6,"('qf  ',9f8.3/4x,9f8.3)") diag_temp
  diag_temp(:) = qlg(idjd,:)
  write (6,"('ql  ',9f8.3/4x,9f8.3)") diag_temp
endif
!if ( diag .and. ntiles==1 ) then
!  call maxmin(t,' t',ktau,1.,kl)
!  call maxmin(qg,'qv',ktau,1.e3,kl)
!  call maxmin(qfg,'qf',ktau,1.e3,kl)
!  call maxmin(qlg,'ql',ktau,1.e3,kl)
!  call maxmin(qrg,'qr',ktau,1.e3,kl)
!  call maxmin(qsng,'qs',ktau,1.e3,kl)
!  call maxmin(qgrg,'qg',ktau,1.e3,kl)
!endif
#endif


! Add convective cloud water into fields for radiation
! done because sometimes newrain drops out all qlg, ending up with 
! zero cloud (although it will be rediagnosed as 1 next timestep)
do concurrent (k = 1:kl)
  do concurrent (iq = 1:imax)
    fl      = max(0., min(1., (t(iq,k)-ticon)/(273.15-ticon)))
    qlrad(iq,k) = qlg(iq,k) + fl*qccon(iq,k)
    qfrad(iq,k) = qfg(iq,k) + (1.-fl)*qccon(iq,k)
    cfrac(iq,k) = min( 1., ccov(iq,k)+clcon(iq,k) ) ! original
  end do
end do


!     Calculate precipitation and related processes
call newsnowrain(dt,rhoa,dz,prf,cdrop,t,qlg,qfg,qrg,qsng,qgrg,                    &
                 precs,qg,stratcloud,rfrac,sfrac,gfrac,preci,precg,qevap,qsubl,   &
                 qauto,qcoll,qaccr,qaccf,fluxr,fluxi,fluxs,fluxg,fluxm,           &
                 fluxf,pfstayice,pfstayliq,pqfsedice,pslopes,prscav,              &
                 condx,ktsav,idjd,mydiag,ncloud,nevapls,ldr,rcm,imax,kl)


#ifndef GPU
if ( nmaxpr==1 .and. mydiag ) then
  write(6,*) 'after newsnowrain'
  diag_temp(:) = t(idjd,:)
  write (6,"('t   ',9f8.2/4x,9f8.2)") diag_temp
  diag_temp(:) = qg(idjd,:)
  write (6,"('qv  ',9f8.3/4x,9f8.3)") diag_temp
  diag_temp(:) = qfg(idjd,:)
  write (6,"('qf  ',9f8.3/4x,9f8.3)") diag_temp
  diag_temp(:) = qlg(idjd,:)
  write (6,"('ql  ',9f8.3/4x,9f8.3)") diag_temp
  diag_temp(:) = qrg(idjd,:)
  write (6,"('qr  ',9f8.3/4x,9f8.3)") diag_temp
  diag_temp(:) = qsng(idjd,:)
  write (6,"('qs  ',9f8.3/4x,9f8.3)") diag_temp
  diag_temp(:) = qgrg(idjd,:)
  write (6,"('qg  ',9f8.3/4x,9f8.3)") diag_temp
end if
!if ( diag .and. ntiles==1 ) then
!  call maxmin(t,' t',ktau,1.,kl)
!  call maxmin(qg,'qv',ktau,1.e3,kl)
!  call maxmin(qfg,'qf',ktau,1.e3,kl)
!  call maxmin(qlg,'ql',ktau,1.e3,kl)
!  call maxmin(qrg,'qr',ktau,1.e3,kl)
!  call maxmin(qsng,'qs',ktau,1.e3,kl)
!  call maxmin(qgrg,'qg',ktau,1.e3,kl)
!endif
#endif


!--------------------------------------------------------------
! Store data needed by prognostic aerosol scheme
! MJT notes - invert levels for aerosol code
if ( abs(iaero)>=2 ) then
  invdt = 1./dt
  ppfprec(:,1) = 0.   !At TOA
  ppfmelt(:,1) = 0.   !At TOA
  ppfsnow(:,1) = 0.   !At TOA
  pprfreeze(:,1) = 0. !At TOA
  do concurrent (k = 1:kl-1)
    ppfprec(:,kl+1-k) = (fluxr(:,k+1)+fluxm(:,k)-fluxf(:,k))*invdt     !flux *entering* layer k
    ppfmelt(:,kl+1-k) = fluxm(:,k)*invdt                               !flux melting in layer k
    ppfsnow(:,kl+1-k) = (fluxi(:,k+1)+fluxs(:,k+1)+fluxg(:,k+1) &
                        -fluxm(:,k)+fluxf(:,k))*invdt                  !flux *entering* layer k
    pprfreeze(:,kl+1-k) = fluxf(:,k)*invdt                             !flux freezing in layer k
  end do
  do concurrent (k = 1:kl)
    ppfevap(:,kl+1-k)    = qevap(:,k)*rhoa(:,k)*dz(:,k)*invdt
    ppfsubl(:,kl+1-k)    = qsubl(:,k)*rhoa(:,k)*dz(:,k)*invdt !flux sublimating or staying in k
    pplambs(:,kl+1-k)    = pslopes(:,k)
    ppmrate(:,kl+1-k)    = (qauto(:,k)+qcoll(:,k))*invdt
    ppmaccr(:,kl+1-k)    = qaccr(:,k)*invdt
    ppfstayice(:,kl+1-k) = pfstayice(:,k)
    ppfstayliq(:,kl+1-k) = pfstayliq(:,k)
    ppqfsedice(:,kl+1-k) = pqfsedice(:,k)
    pprscav(:,kl+1-k)    = prscav(:,k)
  end do
end if
!--------------------------------------------------------------


!========================= Jack's diag stuff =========================
!if ( ncfrp==1 ) then  ! from here to near end; Jack's diag stuff
!  do iq = 1,icfrp
!    tautot(iq)  = 0.
!    cldmax(iq)  = 0.
!    ctoptmp(iq) = 0.
!    ctoppre(iq) = 0.
!    do k = 1,kl
!      fice(iq,k) = 0.
!    enddo
!    kcldfmax(iq) = 0.
!  enddo
!!      cfrp data
!  do k = 1,kl-1
!    do iq = 1,icfrp
!      taul(iq,k) = 0.
!      taui(iq,k) = 0.
!      Reffl = 0.
!      if ( stratcloud(iq,k)>0. ) then
!        tau_sfac = 1.
!        fice(iq,k) = qfrad(iq,k)/(qfrad(iq,k)+qlrad(iq,k)) ! 16/1/06
!!            Liquid water clouds
!        if ( qlg(iq,k)>1.0e-8 ) then
!          Wliq = rhoa(iq,k)*qlg(iq,k)/(stratcloud(iq,k)*(1-fice(iq,k))) !kg/m^3
!          if ( .not.land(iq) ) then !sea
!            rk = 0.8
!          else            !land
!            rk = 0.67
!          endif
!! Reffl is the effective radius at the top of the cloud (calculated following
!! Martin etal 1994, JAS 51, 1823-1842) due to the extra factor of 2 in the
!! formula for reffl. Use mid cloud value of Reff for emissivity.
!          Reffl = (3*2*Wliq/(4*rhow*pi*rk*cdrop(iq,k)))**(1./3)
!          qlpath = Wliq*dz(iq,k)
!          taul(iq,k) = tau_sfac*1.5*qlpath/(rhow*Reffl)
!        endif ! qlg
!! Ice clouds
!        if ( qfg(iq,k)>1.0e-8 ) then
!          Wice = rhoa(iq,k)*qfg(iq,k)/(stratcloud(iq,k)*fice(iq,k)) !kg/m**3
!          sigmai = aice*Wice**bice !visible ext. coeff. for ice
!          taui(iq,k) = sigmai*dz(iq,k) !visible opt. depth for ice
!          taui(iq,k) = tau_sfac*taui(iq,k)
!        endif ! qfg
!      endif !stratcloud
!    enddo ! iq
!  enddo ! k
!! Code to get vertically integrated value...
!! top down to get highest level with stratcloud=cldmax (kcldfmax)
!  do k = kl-1,1,-1
!    do iq = 1,icfrp
!      tautot(iq) = tautot(iq)+stratcloud(iq,k)*(fice(iq,k)*taui(iq,k)+(1.-fice(iq,k))*taul(iq,k))
!      if ( stratcloud(iq,k)>cldmax(iq) ) kcldfmax(iq) = k
!      cldmax(iq) = max(cldmax(iq),stratcloud(iq,k))
!    enddo ! iq
!  enddo ! k
!
!  do iq = 1,icfrp
!    if ( cldmax(iq)>1.e-10 ) then
!      tautot(iq) = tautot(iq)/cldmax(iq)
!
!      cfd = 0.
!      do k = kl,kcldfmax(iq),-1
!        fcf = max(0.,stratcloud(iq,k)-cfd) ! cld frac. from above
!        ctoptmp(iq) = ctoptmp(iq)+fcf*t(iq,k)/cldmax(iq)
!        ctoppre(iq) = ctoppre(iq)+fcf*prf(iq,k)/cldmax(iq)
!        cfd = max(stratcloud(iq,k),cfd)
!      enddo ! k=kl,kcldfmax(iq),-1
!
!    endif ! (cldmax(iq).gt.1.e-10) then
!  enddo   ! iq
!endif    ! ncfrp.eq.1
!========================= end of Jack's diag stuff ======================

condx(:)  = condx(:) + precs(:)
conds(:)  = conds(:) + preci(:)
condg(:)  = condg(:) + precg(:)
precip(:) = precip(:) + precs(:)

return
end subroutine leoncld_work


! from arguments
!      ttg - temperature (K)
!      qtg - water vapour mixing ratio (kg/kg) - called qenv in leoncld
!      qlg - cloud liquid water mixing ratio (kg/kg)
!      qfg - cloud ice mixing ratio (kg/kg)
!
! Output:
!
! from arguments
!      stratcloud - cloudy fraction of grid box
! 
!******************************************************************************

 subroutine newcloud(tdt,land,prf,rhoa,ttg,qtg,qlg,qfg,        &
                     dpsldt,nettend,stratcloud,em,idjd,mydiag, &
                     ncloud,nclddia,rcrit_l,rcrit_s,imax,kl)
!$acc routine vector
 
! This routine is part of the prognostic cloud water scheme

use const_phys                    ! Physical constants
use estab                         ! Liquid saturation function
use parm_m, only : diag, ds       ! Model configuration
use sigs_m                        ! Atmosphere sigma levels
 
implicit none

! Argument list
integer, intent(in) :: idjd, ncloud, nclddia
integer, intent(in) :: imax, kl
real, dimension(imax,kl), intent(in) :: prf
real, dimension(imax,kl), intent(in) :: rhoa
real, dimension(imax,kl), intent(inout) :: ttg
real, dimension(imax,kl), intent(inout) :: qtg
real, dimension(imax,kl), intent(inout) :: qlg
real, dimension(imax,kl), intent(inout) :: qfg
real, dimension(imax,kl), intent(in) :: dpsldt
real, dimension(imax,kl), intent(inout) :: nettend
real, dimension(imax,kl), intent(inout) :: stratcloud
real, dimension(imax), intent(in) :: em
real, intent(in) :: tdt
real, intent(in) :: rcrit_l, rcrit_s
logical, intent(in) :: mydiag
logical, dimension(imax), intent(in) :: land

! Local work arrays and variables
real, dimension(imax,kl) :: qsw
real, dimension(imax,kl) :: qcg, qtot, tliq
real, dimension(imax,kl) :: fice, qcold, rcrit
real, dimension(imax) :: tk
real es, Aprpr, Bprpr, Cice
real qi0, fd, Crate, Qfdep
real qsl, fl, qfnew
real, dimension(imax) :: pk, deles
real, dimension(imax) :: diag_temp
real, dimension(imax) :: qsi
real hlrvap, qs, dqsdt, al, qc, delq

integer k, iq

real decayfac
real, parameter :: rhoic = 700.
real, parameter :: cm0 = 1.e-12 !Initial crystal mass

! Start code : ----------------------------------------------------------


#ifndef GPU
if ( diag.and.mydiag ) then
  write(6,*) 'entering newcloud'
  diag_temp(:) = prf(idjd,:)
  write(6,'(a,30f10.3)') 'prf ',diag_temp
  diag_temp(:) = ttg(idjd,:)
  write(6,'(a,30f10.3)') 'ttg ',diag_temp
  diag_temp(:) = qtg(idjd,:)
  write(6,*) 'qtg ',diag_temp
  diag_temp(:) = qlg(idjd,:)
  write(6,*) 'qlg ',diag_temp
  diag_temp(:) = qfg(idjd,:)
  write(6,*) 'qfg ',diag_temp
end if
#endif

! First melt cloud ice or freeze cloud water to give correct ice fraction fice.
! Then calculate the cloud conserved variables qtot and tliq.
! Note that qcg is the total cloud water (liquid+frozen)

do concurrent (k = 1:kl)
  do concurrent (iq = 1:imax)
    if ( ttg(iq,k)>=tfrz ) then
      fice(iq,k) = 0.
    else if ( ttg(iq,k)>=tice .and. qfg(iq,k)>1.e-12 ) then
      fice(iq,k) = min(qfg(iq,k)/(qfg(iq,k)+qlg(iq,k)), 1.)
    else if ( ttg(iq,k)>=tice ) then
      fice(iq,k) = 0.
    else
      fice(iq,k) = 1.
    end if
    qcg(iq,k)   = qlg(iq,k) + qfg(iq,k)
    qcold(iq,k) = qcg(iq,k)
    qfnew       = fice(iq,k)*qcg(iq,k)
    ttg(iq,k)   = ttg(iq,k) + hlfcp*(qfnew-qfg(iq,k)) !Release L.H. of fusion
    qfg(iq,k)   = qfnew
    qlg(iq,k)   = max(0., qcg(iq,k)-qfg(iq,k))
    qtot(iq,k) = qtg(iq,k) + qcg(iq,k)
    tliq(iq,k) = ttg(iq,k) - hlcp*qcg(iq,k) - hlfcp*qfg(iq,k) 
  end do
end do

#ifndef GPU
if ( diag .and. mydiag ) then
  write(6,*) 'within newcloud'
  diag_temp = ttg(idjd,:)
  write(6,*) 'ttg ',diag_temp
  diag_temp = qcold(idjd,:)
  write(6,*) 'qcold ',diag_temp
  diag_temp = qcg(idjd,:)
  write(6,*) 'qcg ',diag_temp
  diag_temp = qlg(idjd,:)
  write(6,*) 'qlg ',diag_temp
  diag_temp = qfg(idjd,:)
  write(6,*) 'qfg ',diag_temp
  diag_temp = fice(idjd,:)
  write(6,*) 'fice ',diag_temp
end if
#endif


! Precompute the array of critical relative humidities 
if ( nclddia==-3 ) then
  do concurrent (k = 1:kl)
    where ( land(:) )
      rcrit(:,k)=max( rcrit_l, (1.-16.*(1.-sig(k))**3) )
    elsewhere
      rcrit(:,k)=max( rcrit_s, (1.-16.*(1.-sig(k))**3) )
    end where
  enddo
else if ( nclddia<0 ) then
  do concurrent (k = 1:kl)
    where ( land(:) )
      rcrit(:,k)=max( rcrit_l, (1.-4.*(1.-sig(k))**2) )
    elsewhere
      rcrit(:,k)=max( rcrit_s, (1.-4.*(1.-sig(k))**2) )
    end where
  enddo
else if ( nclddia==1 ) then
  do concurrent (k = 1:kl)
    where ( land(:) )
      rcrit(:,k)=max( rcrit_l, sig(k)**3 )
    elsewhere
      rcrit(:,k)=max( rcrit_s, sig(k)**3 )
    end where
  enddo
else if ( nclddia==2 ) then
  do concurrent (k = 1:kl)
    where ( land(:) )
      rcrit(:,k)=rcrit_l
    elsewhere
      rcrit(:,k)=rcrit_s
    end where
  enddo
else if ( nclddia==3 ) then
  do concurrent (k = 1:kl)
    where ( land(:) )
      rcrit(:,k)=max( rcrit_l, sig(k)**2 )          ! .75 for R21 Mk2
    elsewhere
      rcrit(:,k)=max( rcrit_s, sig(k)**2 )          ! .85 for R21 Mk2
    end where
  enddo
else if ( nclddia==4 ) then
  do concurrent (k = 1:kl)
    where ( land(:) )
      rcrit(:,k)=max( rcrit_l, sig(k) )             ! .75 for test Mk2/3
    elsewhere
      rcrit(:,k)=max( rcrit_s, sig(k) )             ! .9  for test Mk2/3
    end where
  enddo
else if ( nclddia==5 ) then  ! default till May 08
  do concurrent (k = 1:kl)
    where ( land(:) )
      rcrit(:,k)=max( rcrit_l, min(.99,sig(k)) )    ! .75 for same as T63
    elsewhere
      rcrit(:,k)=max( rcrit_s, min(.99,sig(k)) )    ! .85 for same as T63
    end where
  enddo
else if ( nclddia==6 ) then
  do concurrent (k = 1:kl)
    where ( land(:) )
      rcrit(:,k)=max(rcrit_l*(1.-.15*sig(k)),sig(k)**4)
    elsewhere
      rcrit(:,k)=max(rcrit_s*(1.-.15*sig(k)),sig(k)**4)
    end where
  enddo
else if ( nclddia==7 ) then
  do concurrent (k = 1:kl)
    where ( land(:) )
      rcrit(:,k)=max(rcrit_l*(1.-.2*sig(k)),sig(k)**4)
    elsewhere
      rcrit(:,k)=max(rcrit_s*(1.-.2*sig(k)),sig(k)**4)
    end where
  enddo
else if ( nclddia>7 ) then  ! e.g. 12    JLM
  ! MJT notes - Lopez (2002) "Implementation and validation of a new pronostic large-scale cloud
  ! and precipitation scheme for climate and data-assimilation purposes" Q J R Met Soc 128, 229-257,
  ! has a useful discussion of the dependence of RHcrit on grid spacing
  do concurrent (k = 1:kl) ! typically set rcrit_l=.75,  rcrit_s=.85
    do concurrent (iq = 1:imax)
      tk(iq) = ds/(em(iq)*208498.) ! MJT suggestion
      fl = (1.+real(nclddia))*tk(iq)/(1.+real(nclddia)*tk(iq))
      ! for rcit_l=.75 & nclddia=12 get rcrit=(0.751, 0.769, .799, .901, .940, .972, .985) for (200, 100, 50, 10, 5, 2, 1) km
      if ( land(iq) ) then
        rcrit(iq,k) = max( 1.-fl*(1.-rcrit_l), sig(k)**3 )
      else
        rcrit(iq,k) = max( 1.-fl*(1.-rcrit_s), sig(k)**3 )
      end if
    end do
  end do
end if  ! (nclddia<0)  .. else ..


if ( ncloud<=3 ) then
  ! usual diagnostic cloud fraction
      
  ! Calculate cloudy fraction of grid box (stratcloud) and gridbox-mean cloud water
  ! using the triangular PDF of Smith (1990)

  do concurrent (k = 1:kl)
    ! Calculate qs and gam=(L/cp)*dqsdt,  at temperature tliq
    pk(:) = 100.0*prf(:,k)
    qsi(:) = qsati(pk,tliq(:,k),imax)   !Ice value
    deles(:) = esdiffx(tliq(:,k),imax)  ! MJT suggestion
    do concurrent (iq = 1:imax)
      hlrvap = (hl+fice(iq,k)*hlf)/rvap   
      qsl = qsi(iq) + epsil*deles(iq)/pk(iq)             !qs over liquid
      qsw(iq,k) = fice(iq,k)*qsi(iq) + (1.-fice(iq,k))*qsl !Weighted qs at temperature Tliq
      qs = qsw(iq,k)
      dqsdt = qs*hlrvap/tliq(iq,k)**2
      al = 1./(1.+(hlcp+fice(iq,k)*hlfcp)*dqsdt)  !Smith's notation
      qc = qtot(iq,k) - qs
      delq = (1.-rcrit(iq,k))*qs     !UKMO style (equivalent to above)
      if ( qc<=-delq ) then
        stratcloud(iq,k) = 0.
        qcg(iq,k) = 0.
      else if ( qc<=0. ) then
        stratcloud(iq,k) = max( 1.e-6, 0.5*((qc+delq)/delq)**2 )  ! for roundoff
        qcg(iq,k) = max( 1.e-8, al*(qc+delq)**3/(6.*delq**2) )    ! for roundoff
      else if ( qc<delq ) then
        stratcloud(iq,k) = max( 1.e-6, 1.-0.5*((qc-delq)/delq)**2 ) ! for roundoff
        qcg(iq,k) = max( 1.e-8, al*(qc-(qc-delq)**3/(6.*delq**2)) ) ! for roundoff
      else
        stratcloud(iq,k) = 1.
        qcg(iq,k) = al*qc
      end if
    end do ! iq loop
  end do   ! k loop

#ifndef GPU
  if ( diag .and. mydiag ) then
    diag_temp(:) = rcrit(idjd,:)
    write(6,*) 'rcrit ',diag_temp
    diag_temp(:) = qtot(idjd,:)
    write(6,*) 'qtot ',diag_temp
    !diag_temp(:) = qsi(idjd,:)
    !write(6,*) 'qsi',diag_temp
    diag_temp(:) = tliq(idjd,:)
    write(6,*) 'tliq',diag_temp
    !diag_temp(:) = qsl(idjd,:)
    !write(6,*) 'qsl ',diag_temp
    diag_temp(:) = qsw(idjd,:)
    write(6,*) 'qsw ',diag_temp
    diag_temp(:) = stratcloud(idjd,:)
    write(6,*) 'stratcloud',diag_temp
    diag_temp(:) = qtot(idjd,:)-qsw(idjd,:)
    write(6,*) 'qc  ',diag_temp  
    diag_temp(:) = qcg(idjd,:)
    write(6,*) 'qcg ',diag_temp
    diag_temp(:) = (1.-rcrit(idjd,:))*qsw(idjd,:)
    write(6,*) 'delq ',diag_temp 
  endif
#endif

  ! Assume condensation or evaporation retains ice fraction fice.
  ! Introduce a time-decay factor for cirrus (as suggested by results of Khvorostyanov & Sassen,
  ! JAS, 55, 1822-1845, 1998). Their suggested range for the time constant is 0.5 to 2 hours.
  ! The grid-box-mean values of qtg and ttg are adjusted later on (below).
  decayfac = exp ( -tdt/7200. )      ! Try 2 hrs
  !decayfac = 0.                     ! Instant adjustment (old scheme)
  do k = 1,kl
    where( ttg(:,k)>=Tice )
      qfg(:,k) = fice(:,k)*qcg(:,k)
      qlg(:,k) = qcg(:,k) - qfg(:,k)
    elsewhere                                 ! Cirrus T range
      qfg(:,k) = qcold(:,k)*decayfac + qcg(:,k)*(1.-decayfac)
      qlg(:,k) = 0.
      qcg(:,k) = qfg(:,k)
    end where
  end do
  
else
  
  ! Tiedtke prognostic cloud fraction model
  ! MJT notes - we use ttg instead of tliq
  do concurrent (k = 1:kl)
    pk(:) = 100.*prf(:,k)
    qsi(:) = qsati(pk,ttg(:,k),imax) ! Ice value
    deles(:) = esdiffx(ttg(:,k),imax)
    do concurrent (iq = 1:imax)
      qsl = qsi(iq) + epsil*deles(iq)/pk(iq)               ! Liquid value
      qsw(iq,k) = fice(iq,k)*qsi(iq) + (1.-fice(iq,k))*qsl ! Weighted qs at temperature Tliq
    end do  
  end do
  
  call progcloud(tdt,qcg,qtot,prf,rhoa,fice,qsw,ttg,rcrit,  &
                 dpsldt,nettend,stratcloud,imax,kl)

  decayfac = exp ( -tdt/7200. )      ! Try 2 hrs
  !decayfac = 0.                     ! Instant adjustment (old scheme)
  do concurrent (k = 1:kl)
    where( ttg(:,k)>=Tice )
      qfg(:,k) = fice(:,k)*qcg(:,k)
      qlg(:,k) = qcg(:,k) - qfg(:,k)
    elsewhere                                 ! Cirrus T range
      qfg(:,k) = qcold(:,k)*decayfac + qcg(:,k)*(1.-decayfac)
      qlg(:,k) = 0.
      qcg(:,k) = qfg(:,k)
    end where
  end do  
  
end if ! ncloud<=3 ..else..


! Do the vapour deposition calculation in mixed-phase clouds:
! Calculate deposition on cloud ice, assuming es(T) is the weighted value of the 
! liquid and ice values.
do concurrent (k = 1:kl)  
  do concurrent (iq = 1:imax)
    if ( stratcloud(iq,k)>0.) then  
      Tk(iq) = tliq(iq,k) + hlcp*(qlg(iq,k)+qfg(iq,k))/stratcloud(iq,k) !T in liq cloud  
      if ( Tk(iq)<tfrz .and. qlg(iq,k)>1.e-8 ) then
        !fl = qlg(iq,k)/max(qfg(iq,k)+qlg(iq,k),1.e-30)
        pk(iq)    = 100.*prf(iq,k)
        qs        = qsati(pk(iq),Tk(iq))
        es        = qs*pk(iq)/0.622 !ice value
        Aprpr     = hl/(rKa*Tk(iq))*(hls/(rvap*Tk(iq))-1.)
        Bprpr     = rvap*Tk(iq)/((Dva/pk(iq))*es)
        deles(iq) = (1.-fice(iq,k))*esdiffx(Tk(iq))
        Cice      = 1.e3*exp(12.96*deles(iq)/es - 0.639) !Meyers et al 1992
        qi0       = cm0*Cice/rhoa(iq,k) !Initial ice mixing ratio
        ! Next 2 lines are for assumption of fully mixed ql and qf (also a line further down).
        qi0       = max(qi0, qfg(iq,k)/stratcloud(iq,k)) !Assume all qf and ql are mixed
        fd        = 1.       !Fraction of cloud in which deposition occurs
        !fd        = fl      !Or, use option of adjacent ql,qf
        Crate     = 7.8*((Cice/rhoa(iq,k))**2/rhoic)**(1./3.)*deles(iq)/((Aprpr+Bprpr)*es)
        qfdep     = fd*stratcloud(iq,k)*sqrt(((2./3.)*Crate*tdt+qi0**(2./3.))**3)
        ! Also need this line for fully-mixed option...
        qfdep     = qfdep - qfg(iq,k)
        qfdep      = min(qfdep, qlg(iq,k))
        qlg(iq,k) = qlg(iq,k) - qfdep
        qfg(iq,k) = qfg(iq,k) + qfdep
        fice(iq,k) = qfg(iq,k)/max(qfg(iq,k)+qlg(iq,k),1.e-30)
      end if
    end if  
  end do
end do    

! Calculate new values of vapour mixing ratio and temperature
do concurrent (k = 1:kl)
  qtg(:,k) = qtot(:,k) - qcg(:,k)
  ttg(:,k) = tliq(:,k) + hlcp*qcg(:,k) + hlfcp*qfg(:,k)
end do

#ifndef GPU
if ( diag .and. mydiag ) then
   write(6,*) 'at end of newcloud'
   diag_temp(:) = ttg(idjd,:)
   write(6,*) 'ttg ',diag_temp
   diag_temp(:) = qcg(idjd,:)
   write(6,*) 'qcg ',diag_temp
   diag_temp(:) = qlg(idjd,:)
   write(6,*) 'qlg ',diag_temp
   diag_temp(:) = qfg(idjd,:)
   write(6,*) 'qfg ',diag_temp
   diag_temp(:) = qtg(idjd,:)
   write(6,*) 'qtg ',diag_temp
end if
#endif

return
end subroutine newcloud

 
! This routine is part of the prognostic cloud scheme. It calculates rainfall
! and the evaporation of rain, and also does the frozen precipitation. It is
! called by progcld.
!
! INPUT/OUTPUT
!
! Input:
!
! from arguments
!      tdt - leapfrog timestep (seconds)
!      rhoa - air density (kg/m**3)
!      dz - layer thicknes (m)
!      prf - pressure at full levels (in hPa. NB: not SI units)
!
! In/Out:
!
! from arguments
!      ttg - temperature (K)
!      qlg - cloud liquid water mixing ratio (kg/kg)
!      qfg - cloud ice mixing ratio (kg/kg)
!      qrg - falling rain (kg/kg)
!      qsng - falling snow (kg/kg)
!      qgrg - falling graupel (kg/kg)
!      precs - amount of stratiform precipitation in timestep (mm)
!      qtg - water vapour mixing ratio (kg/kg) - called qg in C-CAM
!      stratcloud - stratiform cloud fraction
!      cfrain - falling rain fraction
!      cfsnow - falling snow fraction
!      cfgraupel - falling graupel fraction
!
! Output:
!
! from arguments
!      preci - amount of stratiform snowfall in timestep (mm)
!      precg - amount of stratiform graupel in timestep (mm)
!      qevap - evaporation of rainfall (kg/kg)
!      qsubl - sublimation of snowfall (kg/kg)
!      qauto - autoconversion of cloud liquid water (kg/kg)
!      qcoll - collection by rain of cloud liquid water (kg/kg)
!      qaccr - accretion by snow of cloud liquid water (kg/kg)
!
!**************************************************************************

subroutine newsnowrain(tdt_in,rhoa,dz,prf,cdrop,ttg,qlg,qfg,qrg,qsng,qgrg,precs,qtg,stratcloud,cfrain,    &
                       cfsnow,cfgraupel,preci,precg,qevap,qsubl,qauto,qcoll,qaccr,qaccf,fluxr,            &
                       fluxi,fluxs,fluxg,fluxm,fluxf,pfstayice,pfstayliq,pqfsedice,pslopes,prscav,        &
                       condx,ktsav,idjd,mydiag,ncloud,nevapls,ldr,rcm,imax,kl)
!$acc routine vector

use const_phys                    ! Physical constants
use estab                         ! Liquid saturation function
use parm_m, only : diag, nmr, nmaxpr
                                  ! Model configuration

implicit none

integer, intent(in) :: idjd, ncloud, nevapls, ldr
integer, intent(in) :: imax, kl
real, intent(in) :: tdt_in
real, dimension(imax,kl), intent(in) :: rhoa
real, dimension(imax,kl), intent(in) :: dz
real, dimension(imax,kl), intent(in) :: prf
real, dimension(imax,kl), intent(in) :: cdrop
real, dimension(imax,kl), intent(inout) :: ttg
real, dimension(imax,kl), intent(inout) :: qlg
real, dimension(imax,kl), intent(inout) :: qfg
real, dimension(imax,kl), intent(inout) :: qrg
real, dimension(imax,kl), intent(inout) :: qsng
real, dimension(imax,kl), intent(inout) :: qgrg
real, dimension(imax,kl), intent(inout) :: qtg
real, dimension(imax,kl), intent(inout) :: stratcloud
real, dimension(imax,kl), intent(inout) :: cfrain
real, dimension(imax,kl), intent(inout) :: cfsnow
real, dimension(imax,kl), intent(inout) :: cfgraupel
real, dimension(imax,kl), intent(out) :: qevap
real, dimension(imax,kl), intent(out) :: qsubl
real, dimension(imax,kl), intent(out) :: qauto
real, dimension(imax,kl), intent(out) :: qcoll
real, dimension(imax,kl), intent(out) :: qaccr
real, dimension(imax,kl), intent(out) :: qaccf
real, dimension(imax,kl), intent(out) :: pqfsedice
real, dimension(imax,kl), intent(out) :: pfstayice
real, dimension(imax,kl), intent(out) :: pfstayliq
real, dimension(imax,kl), intent(out) :: pslopes
real, dimension(imax,kl), intent(out) :: prscav
real, dimension(imax,kl), intent(out) :: fluxr
real, dimension(imax,kl), intent(out) :: fluxi
real, dimension(imax,kl), intent(out) :: fluxs
real, dimension(imax,kl), intent(out) :: fluxg
real, dimension(imax,kl), intent(out) :: fluxm
real, dimension(imax,kl), intent(out) :: fluxf
real, dimension(imax), intent(in) :: condx
real, dimension(imax), intent(inout) :: precs
real, dimension(imax), intent(inout) :: preci
real, dimension(imax), intent(inout) :: precg
real, intent(in) :: rcm
integer, dimension(imax), intent(in) :: ktsav
logical, intent(in) :: mydiag

real, dimension(imax,kl) :: fluxautorain, fluxautosnow, fluxautograupel
real, dimension(imax,kl) :: cfautorain, cfautosnow, cfautograupel
real, dimension(imax,kl) :: rhov, rhol, rhoi, rhos, rhog, rhor
real, dimension(imax,kl) :: clfr,cifr,qsatg
real, dimension(imax) :: fthruliq,foutliq,fthruice,foutice
real, dimension(imax) :: fthrusnow,foutsnow,fthrugraupel,foutgraupel
real, dimension(imax) :: vi2, vr2, vs2, vg2
real, dimension(imax) :: fluxice,fluxsnow,fluxgraupel,fluxrain
real rhoiin,rhoiout,rhorin,rhorout
real rhosin,rhosout,rhogin,rhogout
real cffluxout
real, dimension(imax) :: crfra,cifra,csfra,cgfra
real, dimension(imax) :: mxclfrrain,rdclfrrain,mxclfrice,rdclfrice
real, dimension(imax) :: mxclfrsnow,rdclfrsnow,mxclfrgraupel,rdclfrgraupel
real cffluxin
real rg, rl, rn, rf, rs
real, dimension(imax) :: rhodz,evap,clrevap,fr,sublflux
real, dimension(imax) :: fcol
real alph, alphaf, qpf
real, dimension(imax) :: pk, csbsav
real n0s, aprpr, bprpr, curly
real, dimension(imax) :: cfmelt, fluxmelt, fluxfreeze
real slopes_g, slopes_s, xwgt, qsl
real, dimension(imax) :: denfac
real, dimension(imax) :: xfrac_graupel, xfrac_snow, xfrac_ice
real, dimension(imax) :: rhototf
real, dimension(imax) :: gam1, deles
real, dimension(kl) :: diag_temp

integer k, n, njumps, iq
real scm3, tdt
real qcrit, qcic, ql, dqls, Crate, ql1, ql2
real Frb, cdts, selfcoll
real qfs, dqfs
real fsclr_g, fsclr_s, fsclr_i
real qvp, iflux, lflux
real drf, drl
real dqf, dqs, dql
real cdt, dttg, csb, bf
real qrn, qsn, qif, qf
real coll
real es
real cftmp, cltmp
real slopes_r, slopes_i
real esi, apr, bpr, cev
real dqsdt, bl, satevap

real, parameter :: n0r = 8.e6        ! intercept for rain
real, parameter :: n0g = 4.e6        ! intercept for graupel
real, parameter :: rho_r = 1.0e3     ! rain density
real, parameter :: rho_s = 0.1e3     ! snow density
real, parameter :: rho_g = 0.4e3     ! grauple density
real, parameter :: qr0_crt = 2.e-4   ! rain -> snow or graupel density threshold
real, parameter :: qi0_crt = 8.e-5   ! ice -> snow density threshold
real, parameter :: qs0_crt = 6.e-3   ! snow -> graupel density threshold
real, parameter :: c_piacr = 0.1     ! accretion rate of rain -> ice
real, parameter :: c_psaut = 1.e-3   ! autoconversion rate of ice -> snow
!real, parameter :: c_pgacs = 1.e-3  ! snow -> graupel "accretion" eff
real, parameter :: sfcrho = 1.2      ! reference density rho_0
real, parameter :: vdifu = 2.11e-5
real, parameter :: tcond = 2.36e-2
real, parameter :: visk = 1.259e-5
real, parameter :: gam263 = 1.456943 ! gamma function for 2.63
real, parameter :: gam275 = 1.608355 ! gamma function for 2.75
real, parameter :: gam325 = 2.54925  ! gamma function for 3.25
real, parameter :: gam350 = 3.323363 ! gamma function for 3.5
real, parameter :: gam380 = 4.694155 ! gamma function for 3.8
real, parameter :: alin = 842.
real, parameter :: clin = 4.8
real, parameter :: gcon = 44.628 ! = 40.74*sqrt(sfcrho)
!real, parameter :: tau_s = 90.   ! (sec) snow melt
!real, parameter :: tau_g = 180.  ! (sec) graupel melt

scm3 = (visk/vdifu)**(1./3.)

do concurrent (k = 1:kl)
  fluxr(:,k)           = 0.
  fluxi(:,k)           = 0.
  fluxs(:,k)           = 0.
  fluxg(:,k)           = 0. 
  fluxm(:,k)           = 0.  
  fluxf(:,k)           = 0.
  fluxautorain(:,k)    = 0.
  fluxautosnow(:,k)    = 0.
  fluxautograupel(:,k) = 0.
  qevap(:,k)           = 0.
  qauto(:,k)           = 0.
  qcoll(:,k)           = 0.
  qsubl(:,k)           = 0.
  qaccr(:,k)           = 0.
  qaccf(:,k)           = 0.
  pqfsedice(:,k)       = 0.
  prscav(:,k)          = 0.  
  pfstayice(:,k)       = 0.  
  pfstayliq(:,k)       = 0. 
  pslopes(:,k)         = 0.
  pk(:)                = 100.*prf(:,k)
  qsatg(:,k)           = qsati(pk(:),ttg(:,k),imax)
  cifr(:,k)            = qfg(:,k)*stratcloud(:,k)/max( qlg(:,k)+qfg(:,k), 1.e-30 )
  clfr(:,k)            = qlg(:,k)*stratcloud(:,k)/max( qlg(:,k)+qfg(:,k), 1.e-30 )
  cfautorain(:,k)      = 0.
  cfautosnow(:,k)      = 0.
  cfautograupel(:,k)   = 0.
end do

! Use full timestep for autoconversion
!njumps = 1
tdt = tdt_in

do concurrent (k = 1:kl-1)
  do concurrent (iq = 1:imax)
    if ( clfr(iq,k)>0. ) then
      qcrit = (4.*pi/3.)*rhow*Rcm**3*cdrop(iq,k)/rhoa(iq,k)
      qcic  = qlg(iq,k)/clfr(iq,k) !In cloud value
      if ( qcic>=qcrit ) then
        Crate    = Aurate*rhoa(iq,k)*(rhoa(iq,k)/(cdrop(iq,k)*rhow))**(1./3.)
        ql1      = 1./pow75(qcic**(-4./3.)+(4./3.)*Crate*tdt)
        ql1      = max( ql1, qcrit ) !Intermediate qlg after auto
        Frb      = dz(iq,k)*rhoa(iq,k)*(qcic-ql1)/tdt
        Frb      = min( Frb, 1.e10 ) ! prevent overflow
        cdts     = tdt*0.5*Ecol*0.24*pow75(Frb) ! old
        selfcoll = min( ql1, ql1*cdts )
        ql2      = ql1 - selfcoll
        ql       = clfr(iq,k)*ql2
        dqls     = max( qlg(iq,k)-ql, 0. )
        cfautorain(iq,k) = clfr(iq,k)
        qauto(iq,k)      = qauto(iq,k) + dqls
        qlg(iq,k)        = qlg(iq,k)   - dqls
        fluxautorain(iq,k) = dqls*rhoa(iq,k)*dz(iq,k)
      end if ! pcic>=qcrit
    end if   ! clfr>0.
  end do ! iq loop
end do   ! k loop

! calculate rate of precipitation of frozen cloud water to snow
if ( ncloud>=3 ) then

  do concurrent (k = 1:kl)
    do concurrent (iq = 1:imax)
      
      ! autoconversion of ice to snow (from Lin et al 1983)
      ! Threshold from WSM6 scheme, Hong et al 2004, Eq(13) : qi0_crt ~8.e-5
      if ( qfg(iq,k)*rhoa(iq,k)>qi0_crt ) then
        qfs  = max( qfg(iq,k)-qi0_crt/rhoa(iq,k), 0. )
        cdts = tdt*c_psaut*exp(0.025*(ttg(iq,k)-tfrz))
        dqfs = max( min( qfg(iq,k), qfs*cdts ), 0. )
        cfautosnow(iq,k)   = cifr(iq,k)
        qfg(iq,k)          = qfg(iq,k) - dqfs
        fluxautosnow(iq,k) = dqfs*rhoa(iq,k)*dz(iq,k)
      end if
    
      ! autoconversion of snow to graupel (from Lin et al 1983)
      if ( qsng(iq,k)*rhoa(iq,k)>qs0_crt ) then
        qfs  = max( qsng(iq,k)-qs0_crt/rhoa(iq,k), 0. )
        cdts = tdt*1.e-3*exp(0.09*(ttg(iq,k)-tfrz))
        dqfs = max( min( qsng(iq,k), qfs*cdts ), 0.) 
        cfautograupel(iq,k)   = cfsnow(iq,k)
        qsng(iq,k)            = qsng(iq,k) - dqfs
        fluxautograupel(iq,k) = dqfs*rhoa(iq,k)*dz(iq,k)
      end if

    end do ! iq loop 
  end do   ! k loop
  
end if ! ( ncloud>=3 )

! update density and area fractions
do concurrent (k = 1:kl)
  cifr(:,k) = stratcloud(:,k)*qfg(:,k)/max(qlg(:,k)+qfg(:,k),1.e-30 )
  clfr(:,k) = max( stratcloud(:,k)-cifr(:,k), 0. )
  rhov(:,k) = qtg(:,k)*rhoa(:,k)
  rhoi(:,k) = qfg(:,k)*rhoa(:,k)
  rhol(:,k) = qlg(:,k)*rhoa(:,k)
  rhor(:,k) = qrg(:,k)*rhoa(:,k)
  rhos(:,k) = qsng(:,k)*rhoa(:,k)
  rhog(:,k) = qgrg(:,k)*rhoa(:,k)
end do


#ifndef GPU
if ( diag .and. mydiag ) then
  diag_temp(:) = stratcloud(idjd,:)
  write(6,*) 'stratcloud',diag_temp
  diag_temp(:) = cifr(idjd,:)
  write(6,*) 'cifr      ',diag_temp
  diag_temp(:) = clfr(idjd,:)
  write(6,*) 'clfr      ',diag_temp
  diag_temp(:) = cfrain(idjd,:)
  write(6,*) 'cfrain    ',diag_temp
  diag_temp(:) = cfsnow(idjd,:)
  write(6,*) 'cfsnow    ',diag_temp
  diag_temp(:) = cfgraupel(idjd,:) 
  write(6,*) 'cfgraupel ',diag_temp
  diag_temp(:) = qlg(idjd,:) 
  write(6,*) 'qlg ',diag_temp
  diag_temp(:) = qfg(idjd,:)
  write(6,*) 'qfg ',diag_temp
  diag_temp(:) = qrg(idjd,:)
  write(6,*) 'qrg ',diag_temp
  diag_temp(:) = qsng(idjd,:)
  write(6,*) 'qsng',diag_temp
  diag_temp(:) = qgrg(idjd,:)
  write(6,*) 'qgrg',diag_temp
endif  ! (diag.and.mydiag)
#endif


! Use sub time-step if required
if ( ncloud>=3 ) then
  njumps = int(tdt_in/(maxlintime+0.01)) + 1
  tdt    = tdt_in/real(njumps)
else
  njumps = 1
  tdt = tdt_in
end if

do n = 1,njumps

  fluxgraupel(:)   = 0.
  mxclfrgraupel(:) = 0. ! max overlap graupel fraction
  rdclfrgraupel(:) = 0. ! rnd overlap graupel fraction
  cgfra(:)         = 0. ! total graupel fraction = mx+rd-mx*rd
  vg2(:)           = 0.1

  fluxsnow(:)   = 0.
  mxclfrsnow(:) = 0. ! max overlap snow fraction
  rdclfrsnow(:) = 0. ! rnd overlap snow fraction
  csfra(:)      = 0. ! total snow fraction = mx+rd-mx*rd
  vs2(:)        = 0.1

  fluxice(:)   = 0.
  mxclfrice(:) = 0. ! max overlap ice fraction
  rdclfrice(:) = 0. ! rnd overlap ice fraction
  cifra(:)     = 0. ! total ice fraction = mx+rd-mx*rd 
  vi2(:)       = 0.1 ! Assume no cloud at top level

  fluxrain(:)   = 0.
  mxclfrrain(:) = 0. ! max overlap rain fraction
  rdclfrrain(:) = 0. ! rnd overlap rain fraction
  crfra(:)      = 1.e-6 ! total rain fraction = mx+rd-mx*rd
  vr2(:)        = 0.


  ! Now work down through the levels...
  do k = kl-1,1,-1
  
    ! misc fields
    do concurrent (iq = 1:imax)
      pk(iq)     = 100.*prf(iq,k)
      rhodz(iq)  = rhoa(iq,k)*dz(iq,k)
      denfac(iq) = sqrt(sfcrho/rhoa(iq,k))
      fluxmelt(iq)   = 0.
      fluxfreeze(iq) = 0.
      cfmelt(iq)     = 0.
    end do
    
    if ( ncloud>=3 ) then
  
      ! Graupel ---------------------------------------------------------------------------
      do concurrent (iq = 1:imax)      

        sublflux(iq) = 0.
        fluxgraupel(iq) = fluxgraupel(iq) + fluxautograupel(iq,k)*tdt/tdt_in
      
        ! Detect max/random overlap clouds that are separated by a clear layer
        if ( (stratcloud(iq,k)>=1.e-10.and.stratcloud(iq,k+1)<1.e-10) .or. nmr==0 ) then
          rdclfrgraupel(iq) = rdclfrgraupel(iq) + mxclfrgraupel(iq) - rdclfrgraupel(iq)*mxclfrgraupel(iq)
          mxclfrgraupel(iq) = 0.
        end if
        cgfra(iq) = max( rdclfrgraupel(iq) + mxclfrgraupel(iq) - rdclfrgraupel(iq)*mxclfrgraupel(iq), 1.e-15 )
       
        ! graupel fall speed (from Lin et al 1983 - see GFDL AM3)
        rg = max( fluxgraupel(iq)/dz(iq,k), 0. )
        if ( cgfra(iq)>=1.e-10 ) then
          vg2(iq) = max( 0.1, 5.34623815*(rg/cgfra(iq))**0.125 )
        end if

        ! Set up the parameters for the flux-divergence calculation
        alph         = tdt*vg2(iq)/dz(iq,k)
        alph         = max( min( alph, 50. ), 0. )
        foutgraupel(iq)  = 1. - exp(-alph)        !analytical
        fthrugraupel(iq) = 1. - foutgraupel(iq)/alph  !analytical

        alphaf = hls*qsatg(iq,k)/(rvap*ttg(iq,k)**2)
        gam1(iq) = hlscp*alphaf !(L/cp)*dqsdt (HBG notation)
        
        if ( fluxgraupel(iq)>0. ) then        
      
          ! Melt falling graupel (based on Lin et al 83)
          rg = max(fluxgraupel(iq), 0.)/dz(iq,k)
          if ( ttg(iq,k)>tfrz .and. rg>1.e-15 ) then
            slopes_g          = ( max( fluxgraupel(iq), 0. )/dz(iq,k)/(pi*n0g*rho_g))**0.25
            qvp               = rhov(iq,k)/rhoa(iq,k)
            cdt               = tdt*2.*pi*n0g/hlf*(tcond*(ttg(iq,k)-tfrz)/rhoa(iq,k)-vdifu*hl*(qsatg(iq,k)-qvp))              &
                               *(0.78*slopes_g**2+0.31*scm3*gam275*sqrt(gcon/visk)*slopes_g**2.75*sqrt(denfac(iq)))
            drf               = max( min( rg, cdt ), 0. )
            iflux             = min( drf*dz(iq,k), fluxgraupel(iq) ) ! flux of graupel
            drf               = iflux/dz(iq,k)                   ! mass of graupel
            dqf               = drf/rhoa(iq,k)                   ! mixing ratio of graupel
            fluxmelt(iq)      = fluxmelt(iq)    + iflux
            fluxgraupel(iq)   = fluxgraupel(iq) - iflux
            dttg              = -hlfcp*dqf
            ttg(iq,k)         = ttg(iq,k) + dttg
            qsatg(iq,k)       = qsatg(iq,k) + gam1(iq)*dttg/hlscp
            rdclfrgraupel(iq) = rdclfrgraupel(iq)*(1.-drf/rg)
            mxclfrgraupel(iq) = mxclfrgraupel(iq)*(1.-drf/rg)
            cftmp             = mxclfrgraupel(iq) + rdclfrgraupel(iq) - mxclfrgraupel(iq)*rdclfrgraupel(iq)
            cfmelt(iq)        = max( cfmelt(iq), max(cgfra(iq)-cftmp,0.) )
            cgfra(iq)         = cftmp
          end if
        
          ! Sublimation of graupel is neglected in the UM and ACCESS 1.3.
          ! (Currently treated the same as LDR97 ice sublimation)
          qvp = rhov(iq,k)/rhoa(iq,k)
          if ( qvp<qsatg(iq,k) ) then ! sublime graupel
            slopes_g        = ( max(fluxgraupel(iq),0.)/dz(iq,k)/(pi*n0g*rho_g))**0.25
            fsclr_g         = max( (1.-cifr(iq,k)-clfr(iq,k))*fluxgraupel(iq), 0. )  
            cdt             = 2.*pi*vdifu*tcond*rvap*n0g*ttg(iq,k)**2                                                    &
                             *(0.78*slopes_g**2+0.31*scm3*gam275*sqrt(gcon/visk)*slopes_g**2.75*sqrt(denfac(iq))) &
                             /(tcond*rvap*ttg(iq,k)**2+hls**2*vdifu*qsatg(iq,k)*rhoa(iq,k))
            dqs             = tdt*cdt*(qsatg(iq,k)-qvp)
            dqs             = min( dqs, (qsatg(iq,k)-qvp)/(1.+gam1(iq)) ) !Don't supersat.
            sublflux(iq)    = min( dqs*rhodz(iq), fsclr_g ) ! flux of graupel
            drf             = sublflux(iq)/dz(iq,k)         ! mass of graupel
            dqs             = drf/rhoa(iq,k)                ! mixing ratio of graupel
            fluxgraupel(iq) = fluxgraupel(iq) - sublflux(iq)
            fsclr_g         = fsclr_g     - sublflux(iq)
            rhov(iq,k)      = rhov(iq,k)  + drf        
            qsubl(iq,k)     = qsubl(iq,k) + dqs
            dttg            = -hlscp*dqs
            ttg(iq,k)       = ttg(iq,k) + dttg
            qsatg(iq,k)     = qsatg(iq,k) + gam1(iq)*dttg/hlscp
          end if
        
          ! Accretion of cloud liquid by falling graupel (from Lin et al 1983 - pgacw)
          ! This calculation uses the incoming fluxgraupel without subtracting sublimation
          ! (since subl occurs only outside cloud), so add sublflux back to fluxgraupel.
          rl = rhol(iq,k)
          if ( fluxgraupel(iq)+sublflux(iq)>0. .and. rl>1.e-15 .and. ttg(iq,k)<tfrz ) then
            slopes_g        = ( max(fluxgraupel(iq)+sublflux(iq),0.)/dz(iq,k)/(pi*n0g*rho_g))**0.25
            cdt             = tdt*pi*n0g*gam350*gcon/4.0*slopes_g**3.5/sqrt(rhoa(iq,k))
            drl             = max( min( cgfra(iq)*rl, rl*cdt/(1.+0.5*cdt) ), 0. ) ! mass of liquid
            lflux           = drl*dz(iq,k)           ! flux of liquid
            dql             = drl/rhoa(iq,k)         ! mixing ratio of liquid
            fluxgraupel(iq) = fluxgraupel(iq) + lflux        
            rhol(iq,k)      = rhol(iq,k) - drl
            qaccr(iq,k)     = qaccr(iq,k) + dql
            dttg            = hlfcp*dql
            ttg(iq,k)       = ttg(iq,k) + dttg
            qsatg(iq,k)     = qsatg(iq,k) + gam1(iq)*dttg/hlscp
            cftmp           = clfr(iq,k)*drl/rl
            clfr(iq,k)      = clfr(iq,k) - cftmp
            mxclfrgraupel(iq) = max( mxclfrgraupel(iq), cftmp )
          end if
        
          ! Accretion of rain by falling graupel (from Lin et al 1983 - pgacr)
          ! (Neglected in UM and ACCESS 1.3)
          rn = rhor(iq,k)
          if ( fluxgraupel(iq)+sublflux(iq)>0. .and. rn>1.e-15 .and. ttg(iq,k)<tfrz ) then
            slopes_g        = ( max( fluxgraupel(iq)+sublflux(iq), 0. )/dz(iq,k)/(pi*n0g*rho_g))**0.25
            slopes_r        = (( max( rn*dz(iq,k), 0. )/max( crfra(iq),1.e-15 )/tdt)**0.22)/714.        
            qrn             = rn/rhoa(iq,k)            
            cdt             = tdt*pi*pi*n0g*n0r*abs(vg2(iq)-vr2(iq))*qrn*(rho_r/rhoa(iq,k))   &
                             *(5.*slopes_r**6*slopes_g+2.*slopes_r**5*slopes_g**2      &
                             +0.5*slopes_r**4*slopes_g**3)          
            drl             = max( min( cgfra(iq)*rn, rn*cdt/(1.+0.5*cdt) ), 0. ) ! mass of rain
            lflux           = drl*dz(iq,k)   ! flux of rain
            dql             = drl/rhoa(iq,k) ! mixing ratio of rain
            fluxgraupel(iq) = fluxgraupel(iq) + lflux
            rhor(iq,k)      = rhor(iq,k) - drl
            dttg            = hlfcp*dql
            ttg(iq,k)       = ttg(iq,k) + dttg
            qsatg(iq,k)     = qsatg(iq,k) + gam1(iq)*dttg/hlscp 
            cftmp           = cfrain(iq,k)*drl/rn
            cfrain(iq,k)    = cfrain(iq,k) - cftmp
            mxclfrgraupel(iq) = max( mxclfrgraupel(iq), cftmp )
          end if
        
          ! Accretion of cloud ice by falling graupel (from Lin et al 1983 - pgaci)
          ! (Neglected in UM and ACCESS 1.3)
          rf = rhoi(iq,k)
          if ( fluxgraupel(iq)+sublflux(iq)>0. .and. rf>1.e-15 .and. ttg(iq,k)<tfrz ) then
            slopes_g        = ( max(fluxgraupel(iq)+sublflux(iq),0.)/dz(iq,k)/(pi*n0g*rho_g))**0.25
            cdt             = tdt*0.1*pi*n0g*gam350*gcon/4.*slopes_g**3.5/sqrt(rhoa(iq,k))
            drf             = max( min( cgfra(iq)*rf, rf*cdt/(1.+0.5*cdt) ), 0. ) ! mass of ice
            iflux           = drf*dz(iq,k)    ! flux of ice
            dqf             = drf/rhoa(iq,k)  ! mixing ratio of ice
            fluxgraupel(iq) = fluxgraupel(iq) + iflux
            rhoi(iq,k)      = rhoi(iq,k) - drf
            qaccf(iq,k)     = qaccf(iq,k) + dqf      
            cftmp           = cifr(iq,k)*drf/rf
            cifr(iq,k)      = cifr(iq,k) - cftmp
            mxclfrgraupel(iq) = max( mxclfrgraupel(iq), cftmp )
          end if
        
          ! Accretion of snow by falling graupel (from Lin et al 1983 - pgacs )
          rs = rhos(iq,k)
          if ( fluxgraupel(iq)+sublflux(iq)>0. .and. rs>1.e-15 .and. ttg(iq,k)<tfrz ) then
            n0s             = 2.e6*exp(-0.12*max(ttg(iq,k)-tfrz,-200.))        
            slopes_s        = ( max(rs,0.)/(pi*rho_s*n0s))**0.25
            slopes_g        = ( max(fluxgraupel(iq)+sublflux(iq),0.)/dz(iq,k)/(pi*n0g*rho_g))**0.25
            qsn             = rs/rhoa(iq,k)  
            cdt             = tdt*pi*pi*n0g*n0s*abs(vg2(iq)-vs2(iq))*qsn*(rho_s/rhoa(iq,k))   &
                             *(5.*slopes_s**6*slopes_g+2.*slopes_s**5*slopes_g**2   &
                             +0.5*slopes_s**4*slopes_g**3)        
            drf             = max( min( cgfra(iq)*rs, rs*cdt/(1.+0.5*cdt) ), 0. ) ! mass of snow
            iflux           = drf*dz(iq,k)    ! flux of snow
            dqf             = drf/rhoa(iq,k)  ! mixing ratio of snow
            fluxgraupel(iq) = fluxgraupel(iq) + iflux
            rhos(iq,k)      = rhos(iq,k) - drf
            qaccf(iq,k)     = qaccf(iq,k) + dqf
            cftmp           = cfsnow(iq,k)*drf/rs
            cfsnow(iq,k)    = cfsnow(iq,k) - cftmp
            mxclfrgraupel(iq) = max( mxclfrgraupel(iq), cftmp )
          end if
       
        end if  ! fluxgraupel>0.

      end do ! iq loop  

     
      ! Snow ------------------------------------------------------------------------------
      do concurrent (iq = 1:imax)
      
        sublflux(iq) = 0.
        fluxsnow(iq) = fluxsnow(iq) + fluxautosnow(iq,k)*tdt/tdt_in
      
        ! Detect max/random overlap clouds that are separated by a clear layer
        if ( (stratcloud(iq,k)>=1.e-10.and.stratcloud(iq,k+1)<1.e-10) .or. nmr==0 ) then
          rdclfrsnow(iq) = rdclfrsnow(iq) + mxclfrsnow(iq) - rdclfrsnow(iq)*mxclfrsnow(iq)
          mxclfrsnow(iq) = 0.
        end if
        csfra(iq) = max( rdclfrsnow(iq) + mxclfrsnow(iq) - rdclfrsnow(iq)*mxclfrsnow(iq), 1.e-15 )
  
        ! Snow fall speed (from Lin et al 1983 - see GFDL AM3)
        rs = max( fluxsnow(iq)/dz(iq,k), 0. )
        if ( csfra(iq)>=1.e-10 ) then
          vs2(iq) = max( 0.1, 1.82*(rs/csfra(iq))**0.0625 )
        end if

        ! Set up the parameters for the flux-divergence calculation
        alph          = tdt*vs2(iq)/dz(iq,k)
        alph          = max( min( alph, 50. ), 0. )
        foutsnow(iq)  = 1. - exp(-alph)          !analytical
        fthrusnow(iq) = 1. - foutsnow(iq)/alph  !analytical

        alphaf = hls*qsatg(iq,k)/(rvap*ttg(iq,k)**2)
        gam1(iq) = hlscp*alphaf !(L/cp)*dqsdt (HBG notation)       

        if ( fluxsnow(iq)>0. ) then
          
          ! Melt falling snow if > 0 deg C due to rain accretion
          ! (based on Lin et al 83, but using 0.65 and 0.44 coeffs following the UM approach)
          rs = max(fluxsnow(iq), 0.)/dz(iq,k)
          if ( ttg(iq,k)>tfrz .and. rs>1.e-15 ) then
            n0s            = 2.e6*exp(-0.12*max(ttg(iq,k)-tfrz,-200.))        
            slopes_s       = ( max(fluxsnow(iq),0.)/dz(iq,k)/(pi*rho_s*n0s))**0.25
            qvp            = rhov(iq,k)/rhoa(iq,k)  
            cdt            = tdt*2.*pi*n0s/hlf*(tcond*(ttg(iq,k)-tfrz)/rhoa(iq,k)-vdifu*hl*(qsatg(iq,k)-qvp))          &
                                     *(0.65*slopes_s**2+0.44*scm3*gam263*sqrt(clin/visk)*slopes_s**2.63*sqrt(denfac(iq)))
            drf            = max( min( rs, cdt ), 0. ) 
            iflux          = min( drf*dz(iq,k), fluxsnow(iq) )    ! flux of snow
            drf            = iflux/dz(iq,k)                      ! mass of snow
            dqf            = drf/rhoa(iq,k)                      ! mixing ratio of snow
            fluxmelt(iq)   = fluxmelt(iq) + iflux
            fluxsnow(iq)   = fluxsnow(iq) - iflux
            dttg           = -hlfcp*dqf
            ttg(iq,k)      = ttg(iq,k) + dttg
            qsatg(iq,k)    = qsatg(iq,k) + gam1(iq)*dttg/hlscp
            rdclfrsnow(iq) = rdclfrsnow(iq)*(1.-drf/rs)
            mxclfrsnow(iq) = mxclfrsnow(iq)*(1.-drf/rs)
            cftmp          = mxclfrsnow(iq) + rdclfrsnow(iq) - mxclfrsnow(iq)*rdclfrsnow(iq)
            cfmelt(iq)     = max( cfmelt(iq), max( csfra(iq)-cftmp, 0. ) )
            csfra(iq)      = cftmp      
          end if
        
          ! Compute the sublimation of snow falling from level k+1 into level k
          ! (Currently treated the same as LDR97 ice sublimation - see UM and ACCESS 1.3)
          qvp = rhov(iq,k)/rhoa(iq,k)
          if ( qvp<qsatg(iq,k) ) then ! sublime snow
            n0s          = 2.e6*exp(-0.12*max(ttg(iq,k)-tfrz,-200.))        
            slopes_s     = ( max(fluxsnow(iq),0.)/dz(iq,k)/(pi*rho_s*n0s))**0.25
            fsclr_s      = max( (1.-cifr(iq,k)-clfr(iq,k))*fluxsnow(iq), 0. )  
            cdt          = 2.*pi*vdifu*tcond*rvap*n0s*ttg(iq,k)**2                                                 &
                               *(0.65*slopes_s**2+0.44*scm3*gam263*sqrt(clin/visk)*slopes_s**2.63*sqrt(denfac(iq))) &
                               /(tcond*rvap*ttg(iq,k)**2+hls**2*vdifu*qsatg(iq,k)*rhoa(iq,k))
            dqs          = tdt*cdt*(qsatg(iq,k)-qvp)
            dqs          = min( dqs, (qsatg(iq,k)-qvp)/(1.+gam1(iq)) ) !Don't supersat.
            sublflux(iq) = min( dqs*rhodz(iq), fsclr_s ) ! flux of snow
            drf          = sublflux(iq)/dz(iq,k)         ! mass of snow
            dqs          = drf/rhoa(iq,k)                ! mixing ratio of snow
            fluxsnow(iq) = fluxsnow(iq) - sublflux(iq)
            fsclr_s      = fsclr_s  - sublflux(iq)
            rhov(iq,k)   = rhov(iq,k)   + drf
            qsubl(iq,k)  = qsubl(iq,k)  + dqs
            dttg         = -hlscp*dqs
            ttg(iq,k)    = ttg(iq,k) + dttg
            qsatg(iq,k)  = qsatg(iq,k) + gam1(iq)*dttg/hlscp
          end if
        
          ! Accretion of cloud liquid by falling snow (from Lin et al 1983 - psacw)
          rl = rhol(iq,k)
          if ( fluxsnow(iq)+sublflux(iq)>0. .and. rl>1.e-15 .and. ttg(iq,k)<tfrz ) then
            n0s      = 2.e6*exp(-0.12*max(ttg(iq,k)-tfrz,-200.))        
            slopes_s = ( max(fluxsnow(iq)+sublflux(iq),0.)/dz(iq,k)/(pi*rho_s*n0s))**0.25
            cdt          = tdt*denfac(iq)*pi*clin*gam325*n0s/4.*slopes_s**3.25
            drl          = max( min( csfra(iq)*rl, rl*cdt/(1.+0.5*cdt) ), 0. ) ! mass of liquid
            lflux        = drl*dz(iq,k)                                        ! flux of liquid
            dql          = drl/rhoa(iq,k)                                      ! mixing ratio of liquid
            fluxsnow(iq) = fluxsnow(iq) + lflux
            rhol(iq,k)   = rhol(iq,k)   - drl
            qaccr(iq,k)  = qaccr(iq,k)  + dql
            dttg         = hlfcp*dql
            ttg(iq,k)    = ttg(iq,k) + dttg
            qsatg(iq,k)  = qsatg(iq,k) + gam1(iq)*dttg/hlscp
            cftmp        = clfr(iq,k)*drl/rl
            clfr(iq,k)   = clfr(iq,k) - cftmp
            mxclfrsnow(iq) = max( mxclfrsnow(iq), cftmp )
          end if
        
          ! Accretion of rain by falling snow to form snow (from Lin et al 1983 - psacr)
          rn = rhor(iq,k)
          if ( fluxsnow(iq)+sublflux(iq)>0. .and. rn>1.e-15 .and. ttg(iq,k)<tfrz ) then
            n0s          = 2.e6*exp(-0.12*max(ttg(iq,k)-tfrz,-200.))        
            slopes_s     = ( max(fluxsnow(iq)+sublflux(iq),0.)/dz(iq,k)/(pi*rho_s*n0s))**0.25
            slopes_r     = (( max(rn*dz(iq,k),0.)/max(crfra(iq),1.e-15)/tdt)**0.22)/714.
            qrn          = rn/rhoa(iq,k)  
            cdt          = tdt*pi*pi*n0r*n0s*abs(vs2(iq)-vr2(iq))*qrn*(rho_r/rhoa(iq,k))         &
                                *(5.*slopes_r**6*slopes_s+2.*slopes_r**5*slopes_s**2  &
                                 +0.5*slopes_r**4*slopes_s**3)
            drl          = max( min( crfra(iq)*rn, rn*cdt/(1.+0.5*cdt) ), 0. ) ! mass of rain
            lflux        = drl*dz(iq,k)                                                 ! flux of rain
            dql          = drl/rhoa(iq,k)                                               ! mixing ratio of rain
            fluxsnow(iq) = fluxsnow(iq) + lflux
            rhor(iq,k)   = rhor(iq,k)   - drl
            dttg         = hlfcp*dql
            ttg(iq,k)    = ttg(iq,k) + dttg
            qsatg(iq,k)  = qsatg(iq,k) + gam1(iq)*dttg/hlscp  
            cftmp        = cfrain(iq,k)*drl/rn
            cfrain(iq,k) = cfrain(iq,k) - cftmp
            mxclfrsnow(iq) = max( mxclfrsnow(iq), cftmp )
          end if
        
          ! Accretion of rain by falling snow to form graupel (neglected in Lin83 but included in UM)   
    
          ! Accretion of cloud ice by falling snow (from HDC 2004 - psaci)
          rf = rhoi(iq,k)
          if ( fluxsnow(iq)+sublflux(iq)>0. .and. rf>1.e-15 .and. ttg(iq,k)<tfrz ) then
            n0s          = 2.e6*exp(-0.12*max(ttg(iq,k)-tfrz,-200.))        
            slopes_s     = ( max(fluxsnow(iq)+sublflux(iq),0.)/dz(iq,k)/(pi*rho_s*n0s))**0.25
            esi          = exp(0.05*max(ttg(iq,k)-tfrz,-100.))       ! efficiency
            cdt          = tdt*denfac(iq)*27.737*n0s*esi*slopes_s**3.41
            drf          = max( min( csfra(iq)*rf, rf*cdt/(1.+0.5*cdt) ), 0. ) ! mass of ice
            iflux        = drf*dz(iq,k)                                        ! flux of ice
            dqf          = drf/rhoa(iq,k)                                      ! mixing ratio of ice
            fluxsnow(iq) = fluxsnow(iq) + iflux
            rhoi(iq,k)   = rhoi(iq,k)   - drf
            qaccf(iq,k)  = qaccf(iq,k)  + dqf
            cftmp        = cifr(iq,k)*drf/rf
            cifr(iq,k)   = cifr(iq,k) - cftmp
            mxclfrsnow(iq) = max( mxclfrsnow(iq), cftmp )
          end if
        
        end if  ! fluxsnow(iq)>0.

      end do ! iq loop  
      
    end if ! ncloud>=3

        
    ! Ice ---------------------------------------------------------------------------------
    do concurrent (iq = 1:imax)
        
      sublflux(iq) = 0.  
  
      ! Set up the rate constant for ice sublimation
      ! MJT notes - curly and Csbsav depend on vi2(:,k+1), so vi2(:,k) can be updated below
      slopes_i = 1.6e3*10**(-0.023*(ttg(iq,k)-tfrz))
      es = qsatg(iq,k)*pk(iq)/epsil
      Aprpr = (hls/(rKa*ttg(iq,k)))*(hls/(rvap*ttg(iq,k))-1.)
      Bprpr = rvap*ttg(iq,k)/((Dva/pk(iq))*es)
      if ( nevapls==-1 .or. (nevapls==-2.and.condx(iq)>0..and.k<=ktsav(iq)) ) then
        curly = 0.
      else
        curly = 0.65*slopes_i**2+0.493*slopes_i*sqrt(slopes_i*vi2(iq)*rhoa(iq,k)/um) !Factor in curly brackets
      end if
      ! Define the rate constant for sublimation of snow, omitting factor rhoi
      Csbsav(iq) = 4.*curly/(rhoa(iq,k)*qsatg(iq,k)*(Aprpr+Bprpr)*pi*vi2(iq)*rho_s)
    
      ! Detect max/random overlap clouds that are separated by a clear layer
      if ( (stratcloud(iq,k)>=1.e-10.and.stratcloud(iq,k+1)<1.e-10) .or. nmr==0 ) then
        rdclfrice(iq) = rdclfrice(iq) + mxclfrice(iq) - rdclfrice(iq)*mxclfrice(iq)
        mxclfrice(iq) = 0.
      end if
      cifra(iq) = max( rdclfrice(iq) + mxclfrice(iq) - rdclfrice(iq)*mxclfrice(iq), 1.e-15 )
      
    end do ! iq loop  
  
    ! Set up snow fall speed field
    select case(abs(ldr))
      case(1)
        where ( cifr(:,k)>=1.e-10 )
          vi2(:) = max( 0.1, 3.23*(max(rhoi(:,k),0.)/cifr(:,k))**0.17 )  ! Ice fall speed from LDR 1997
        end where
      case(2)
        where ( cifr(:,k)>=1.e-10 )
          vi2(:) = 0.9*3.23*(rhoi(:,k)/cifr(:,k))**0.17
        end where
      case(3)
        where ( cifr(:,k)>=1.e-10 )
          vi2(:) = max( 0.1, 2.05+0.35*log10(rhoi(:,k)/rhoa(:,k)/cifr(:,k)) )
        end where
      case(4)
        where ( cifr(:,k)>=1.e-10 )
          vi2(:) = 1.4*3.23*(rhoi(:,k)/cifr(:,k))**0.17
        end where
      case(5)
        where ( cifr(:,k)>=1.e-10 )  
          vi2(:) = max( 0.1, 3.29*(max( rhoi(:,k), 0. )/cifr(:,k))**0.16 ) ! from Lin et al 1983 
        end where  
      case(11)
        ! following are alternative slightly-different versions of above
        ! used for I runs from 29/4/05 till 30/8/05
        ! for given qfg, large cifr implies small ice crystals, 
        ! with a small fall speed. 
        ! Note that for very small qfg, cifr is small.
        ! But rhoi is like qfg, so ratio should also be small and OK.
        vi2(:) = max( vi2(:), 3.23*(rhoi(:,k)/max(cifr(:,k),1.e-30))**0.17 )
      case(22)
        vi2(:) = max( vi2(:), 0.9*3.23*(rhoi(:,k)/max(cifr(:,k),1.e-30))**0.17 )
      case(33)
        ! following max gives vi2=.1 for qfg=cifr=0
        vi2(:) = max( vi2(:), 2.05+0.35*log10(max(rhoi(:,k)/rhoa(:,k),2.68e-36)/max(cifr(:,k),1.e-30)) )
      case(55)
        vi2(:) = max( vi2(:), 3.29*(max(rhoi(:,k),0.)/cifr(:,k))**0.16 ) ! from Lin et al 1983   
    end select
    vi2 = max( vi2, 0.001 )  

    do concurrent (iq = 1:imax)
    
      ! Set up the parameters for the flux-divergence calculation
      alph         = tdt*vi2(iq)/dz(iq,k)
      alph         = max( min( alph, 50. ), 0. )
      foutice(iq)  = 1. - exp(-alph)    !analytical
      fthruice(iq) = 1. - foutice(iq)/alph  !analytical  

      alphaf = hls*qsatg(iq,k)/(rvap*ttg(iq,k)**2)
      gam1(iq) = hlscp*alphaf !(L/cp)*dqsdt (HBG notation)      
      
      if ( fluxice(iq)>0. ) then
        
        ! Melt falling ice if > 0 deg C
        if ( ttg(iq,k)>tfrz ) then
          qif           = fluxice(iq)/rhodz(iq)      !Mixing ratio of ice
          fluxmelt(iq)  = fluxmelt(iq) + fluxice(iq)
          dttg          = -hlfcp*qif
          ttg(iq,k)     = ttg(iq,k) + dttg
          qsatg(iq,k)   = qsatg(iq,k) + gam1(iq)*dttg/hlscp
          cfmelt(iq)    = max( cfmelt(iq), cifra(iq) )
          fluxice(iq)   = 0.
          cifra(iq)     = 0.
          rdclfrice(iq) = 0.
          mxclfrice(iq) = 0.
        end if
      
        ! Compute the sublimation of ice falling from level k+1 into level k
        qvp = rhov(iq,k)/rhoa(iq,k)
        if ( qvp<qsatg(iq,k) ) then ! sublime ice
          fsclr_i      = (1.-cifr(iq,k)-clfr(iq,k))*fluxice(iq)  
          Csb          = Csbsav(iq)*fluxice(iq)/tdt
          bf           = 1. + 0.5*Csb*tdt*(1.+gam1(iq))
          dqs          = max( 0., tdt*(Csb/bf)*(qsatg(iq,k)-qvp) )
          dqs          = min( dqs, (qsatg(iq,k)-qvp)/(1.+gam1(iq)) ) !Don't supersat.
          sublflux(iq) = min( dqs*rhodz(iq), fsclr_i ) ! flux of ice
          drf          = sublflux(iq)/dz(iq,k)         ! mass of ice
          dqs          = drf/rhoa(iq,k)                ! mixing ratio of ice     
          fluxice(iq)  = fluxice(iq) - sublflux(iq)
          fsclr_i      = fsclr_i - sublflux(iq)
          rhov(iq,k)   = rhov(iq,k)  + drf
          qsubl(iq,k)  = qsubl(iq,k) + dqs
          dttg         = -hlscp*dqs
          ttg(iq,k)    = ttg(iq,k) + dttg
          qsatg(iq,k)  = qsatg(iq,k) + gam1(iq)*dttg/hlscp
        end if
      
        ! Accretion of cloud liquid by falling ice (neglected in Lin et al 1983, but
        ! included in UM and ACCESS 1.3 as piacw)
        ! This calculation uses the incoming fluxice without subtracting sublimation
        ! (since subl occurs only outside cloud), so add sublflux back to fluxice.
        rl = rhol(iq,k)
        if ( fluxice(iq)+sublflux(iq)>0. .and. rl>1.e-15 ) then
          slopes_i    = 1.6e3*10**(-0.023*(ttg(iq,k)-tfrz))
          cdt         = Eac*slopes_i*(fluxice(iq)+sublflux(iq))/(2.*rhosno)
          drl         = max( min( cifra(iq)*rl, rl*cdt/(1.+0.5*cdt) ), 0. ) ! mass of liquid
          lflux       = drl*dz(iq,k)   ! flux of liquid
          dql         = drl/rhoa(iq,k) ! mixing ratio of liquid
          fluxice(iq) = fluxice(iq) + lflux
          rhol(iq,k)  = rhol(iq,k)  - drl
          qaccr(iq,k) = qaccr(iq,k) + dql
          dttg        = hlfcp*dql
          ttg(iq,k)   = ttg(iq,k) + dttg
          qsatg(iq,k) = qsatg(iq,k) + gam1(iq)*dttg/hlscp
          cftmp       = clfr(iq,k)*drl/rl
          clfr(iq,k)  = clfr(iq,k) - cftmp
          mxclfrice(iq) = max( mxclfrice(iq), cftmp )
        end if
      
      end if ! fluxice(iq)>0.
      
    end do ! iq loop   

    if ( ncloud>=3 ) then
      do concurrent (iq = 1:imax)
        if ( fluxice(iq)>0. ) then
            
          ! Accretion of rain by falling ice to produce ice (from Lin et al 1983 - piacr)
          ! (see UM and ACCESS 1.3 piacr-c for an alternate formulation)
          rn  = rhor(iq,k)
          if ( fluxice(iq)+sublflux(iq)>0. .and. rn>1.e-15 .and. ttg(iq,k)<tfrz ) then
            qf           = max(fluxice(iq)+sublflux(iq),0.)/rhodz(iq)
            cdt          = tdt*denfac(iq)*c_piacr*qf/sqrt(rhoa(iq,k))
            drl          = max( min( cifra(iq)*rn, rn*cdt/(1.+0.5*cdt) ), 0. ) ! mass of rain
            lflux        = drl*dz(iq,k)   ! flux of rain
            dql          = drl/rhoa(iq,k) ! mixing ratio of rain
            fluxice(iq)  = fluxice(iq) + lflux
            rhor(iq,k)   = rhor(iq,k)  - drl
            dttg         = hlfcp*dql
            ttg(iq,k)    = ttg(iq,k) + dttg
            qsatg(iq,k)  = qsatg(iq,k) + gam1(iq)*dttg/hlscp
            cftmp        = cfrain(iq,k)*drl/rn
            cfrain(iq,k) = cfrain(iq,k) - cftmp
            mxclfrice(iq) = max( mxclfrice(iq), cftmp )
          end if
      
          ! Accretion of rain by falling ice to produce graupel (Neglected in Lin et al 1983)
          ! (see UM and ACCESS 1.3 piacr-g for an alternate formulation)    
          
        end if ! fluxice>0.
      end do ! iq loop
    end if ! ncloud>=3  
    
    
    ! store slope for aerosols
    do concurrent (iq = 1:imax)
      slopes_i      = 1.6e3*10**(-0.023*(ttg(iq,k)-tfrz))
      pslopes(iq,k) = pslopes(iq,k) + slopes_i*tdt/tdt_in  
    end do
    
    
    ! Rain --------------------------------------------------------------------------------
    do concurrent (iq = 1:imax)
      evap(iq) = 0.

      ! Add flux of melted snow to fluxrain
      fluxrain(iq) = fluxrain(iq) + fluxmelt(iq) + fluxautorain(iq,k)*tdt/tdt_in
      mxclfrrain(iq) = max( mxclfrrain(iq), cfmelt(iq) )
    
      ! Detect maximum/random overlap clouds that are separated by a clear layer
      if ( (stratcloud(iq,k)>=1.e-10.and.stratcloud(iq,k+1)<1.e-10) .or. nmr==0 ) then
        rdclfrrain(iq) = rdclfrrain(iq) + mxclfrrain(iq) - rdclfrrain(iq)*mxclfrrain(iq)
        mxclfrrain(iq) = 0.
      end if
      crfra(iq) = max( rdclfrrain(iq) + mxclfrrain(iq) - rdclfrrain(iq)*mxclfrrain(iq), 1.e-15 )
      
    end do ! iq loop  
    
    ! Calculate rain fall speed (MJT suggestion)
    if ( ncloud>=2 ) then
      do concurrent (iq = 1:imax)
        Fr(iq)       = max( fluxrain(iq)/tdt/max(crfra(iq),1.e-15),0.)
        vr2(iq)      = max( 0.1, 11.3*Fr(iq)**(1./9.)/sqrt(rhoa(iq,k)) )  !Actual fall speed
        !vr2(iq)     = max( 0.1, 5./sqrt(rhoa(iq,k)) )                    !Nominal fall speed
        alph         = tdt*vr2(iq)/dz(iq,k)
        alph         = max( min( alph, 50. ), 0. )
        foutliq(iq)  = 1. - exp(-alph)
        fthruliq(iq) = 1. - foutliq(iq)/alph
      end do
    else
      vr2(:) = 9.e9   
      foutliq(:)  = 1.
      fthruliq(:) = 1.
    end if
    
    do concurrent (iq = 1:imax)

      alphaf = hls*qsatg(iq,k)/(rvap*ttg(iq,k)**2)
      gam1(iq)   = hlscp*alphaf !(L/cp)*dqsdt (HBG notation)    
      
    end do ! iq loop  
    
    if ( ncloud>=3 ) then
      do concurrent (iq = 1:imax)  
        if ( fluxrain(iq)>0. ) then
        
          ! Freezing rain to produce graupel (pgfr)
          ! (Neglected in UM and ACCESS 1.3)
          rn = max(fluxrain(iq),0.)/dz(iq,k)
          if ( rn>1.e-15 .and. ttg(iq,k)<tfrz ) then
            slopes_r        = (( max(fluxrain(iq),0.)/max(crfra(iq),1.e-15)/tdt)**0.22)/714.
            ! MJT notes - limit temperature to -100 C to avoid overflow with single precision
            cdt             = tdt*20.e2*pi**2*n0r*(rho_r/rhoa(iq,k))*slopes_r**7 &
                                   *(exp(-0.66*max(ttg(iq,k)-tfrz,-100.))-1.)
            drl             = max( min( rn, rn*cdt/(1.+0.5*cdt) ), 0. )
            lflux           = min( drl*dz(iq,k), fluxrain(iq) ) ! flux
            lflux           = min( lflux, rhodz(iq)*(tfrz-ttg(iq,k))/hlfcp ) ! do not overshoot tfrz
            drl             = lflux/dz(iq,k) ! mass
            dql             = drl/rhoa(iq,k) ! mixing ratio
            fluxrain(iq)    = fluxrain(iq)    - lflux
            fluxgraupel(iq) = fluxgraupel(iq) + lflux
            fluxfreeze(iq)  = fluxfreeze(iq)  + lflux
            dttg            = hlfcp*dql
            ttg(iq,k)       = ttg(iq,k) + dttg
            qsatg(iq,k)     = qsatg(iq,k) + gam1(iq)*dttg/hlscp
            rdclfrrain(iq)  = rdclfrrain(iq)*(1.-drl/rn)
            mxclfrrain(iq)  = mxclfrrain(iq)*(1.-drl/rn)
            cltmp           = mxclfrrain(iq) + rdclfrrain(iq) - mxclfrrain(iq)*rdclfrrain(iq)
            mxclfrgraupel(iq) = max( mxclfrgraupel(iq), max(crfra(iq)-cltmp, 0.) )
            crfra(iq)       = cltmp
          end if
          
        end if ! fluxrain>0.
      end do ! iq loop
    end if ! ncloud>=3  

    deles(:) = esdiffx(ttg(:,k),imax)
    do concurrent (iq = 1:imax)  
      if ( fluxrain(iq)>0. ) then
        
        ! Evaporation of rain
        qpf         = fluxrain(iq)/rhodz(iq) !Mix ratio of rain which falls into layer
        clrevap(iq) = (1.-clfr(iq,k)-cifr(iq,k))*qpf
        qsl         = qsatg(iq,k) + epsil*deles(iq)/pk(iq)
        qvp         = rhov(iq,k)/rhoa(iq,k)
        if ( crfra(iq)>0. ) then
          es          = qsl*pk(iq)/epsil
          Apr         = (hl/(rKa*ttg(iq,k)))*(hl/(rvap*ttg(iq,k))-1.)
          Bpr         = rvap*ttg(iq,k)/((Dva/pk(iq))*es)
          Fr(iq)      = fluxrain(iq)/tdt/max(crfra(iq), 1.e-15)
          Cev         = crfra(iq)*3.8e2*sqrt(Fr(iq)/rhoa(iq,k))/(qsl*(Apr+Bpr))
          dqsdt       = hl*qsl/(rvap*ttg(iq,k)**2)
          bl          = 1. + 0.5*Cev*tdt*(1.+hlcp*dqsdt)
          evap(iq)    = tdt*(Cev/bl)*(qsl-qvp)
          satevap     = (qsl-qvp)/(1.+hlcp*dqsdt)  !Evap to saturate
          ! vr2=11.3*Fr(iq)**(1./9.)/sqrt(rhoa(mg,k)) !Actual fall speed
          ! vr2=5./sqrt(rhoa(mg,k))               !Nominal fall speed
          evap(iq) = max( 0., min( evap(iq), satevap, clrevap(iq) ) )
        end if ! crfra>0.
        
      end if ! fluxrain>0.
    end do ! iq loop
    
    select case(nevapls)
      case(-1)  
        evap(:) = 0.
      case(-2)
        where ( k<=ktsav(:) .and. condx(:)>0. )
          evap(:) = 0.
        end where
      case(-3)
        evap(:) = 0.5*evap
      case(-4)
        where ( k<=ktsav(:) .and. condx(:)>0. )
          evap(:) = 0.5*evap ! usual
        end where
    end select
      
    fcol(:) = 0.
    Fr(:) = 0.
    do concurrent (iq = 1:imax)  
      if ( fluxrain(iq)>0. ) then
   
        drl        = evap(iq)*rhoa(iq,k) ! mass
        rhov(iq,k) = rhov(iq,k) + drl
        ttg(iq,k)  = ttg(iq,k) - hlcp*evap(iq)
        !frclr  = rhodz(iq)*(clrevap(iq)-evap(iq)) ! flux over tdt
      
        ! Now do the collection of liquid cloud by rain term (cf. pracc in Lin83).
        rl = rhol(iq,k)
        if ( rl>1.e-15 ) then
          Fr(iq)       = max(fluxrain(iq),0.)/tdt/max(crfra(iq),1.e-15)
          fcol(iq)     = crfra(iq)
          cdt          = tdt*Ecol*0.24*fcol(iq)*pow75(Fr(iq))
          coll         = max( min( rhol(iq,k), rhol(iq,k)*cdt/(1.+0.5*cdt) ), 0. ) ! mass
          lflux        = coll*dz(iq,k)                                            ! flux
          dql          = coll/rhoa(iq,k)                                          ! mixing ratio
          fluxrain(iq) = fluxrain(iq) + lflux
          rhol(iq,k)   = rhol(iq,k)   - coll
          qcoll(iq,k)  = qcoll(iq,k)  + dql
          cltmp        = clfr(iq,k)*coll/rl
          clfr(iq,k)   = clfr(iq,k) - cltmp
          mxclfrrain(iq) = max( mxclfrrain(iq), cltmp )
        end if
      
        ! subtract evaporated rain
        lflux        = evap(iq)*rhodz(iq)
        fluxrain(iq) = max( fluxrain(iq) - lflux, 0. ) !To avoid roundoff -ve's
        
      end if ! fluxrain>0.  
    end do ! iq loop
      
    if ( ncloud>=3 ) then
      do concurrent (iq = 1:imax)  
        if ( fluxrain(iq)>0. ) then
        
          ! Accretion of cloud snow by rain (from Lin et al 1983 - pracs)
          rs = max( rhos(iq,k), 0. )
          if ( rs>1.e-15 .and. ttg(iq,k)>tfrz+1. ) then
            n0s          = 2.e6*exp(-0.12*max(ttg(iq,k)-tfrz,-200.))        
            slopes_s     = ( max(rs,0.)/(pi*rho_s*n0s))**0.25
            slopes_r     = (( max(fluxrain(iq),0.)/max(crfra(iq),1.e-15)/tdt)**0.22)/714.  
            qsn          = max( rs/rhoa(iq,k), 0. )  
            cdt          = tdt*pi*pi*n0r*n0s*abs(vr2(iq)-vs2(iq))*qsn*(rho_s/rhoa(iq,k))        &
                                *(5.*slopes_s**6*slopes_r+2.*slopes_s**5*slopes_r**2  &
                                 +0.5*slopes_s**4*slopes_r**3)
            drf          = max( min( crfra(iq)*rs, rs*cdt/(1.+0.5*cdt) ), 0. ) ! mass
            lflux        = drf*dz(iq,k)                                     ! flux
            dqf          = drf/rhoa(iq,k)                                   ! mixing ratio
            fluxrain(iq) = fluxrain(iq) + lflux
            rhos(iq,k)   = rhos(iq,k)   - drf
            dttg         = hlfcp*dqf
            ttg(iq,k)    = ttg(iq,k) - dttg
            qsatg(iq,k)  = qsatg(iq,k) - gam1(iq)*dttg/hlscp      
            cftmp        = cfsnow(iq,k)*drf/rs
            cfsnow(iq,k) = cfsnow(iq,k) - cftmp
            mxclfrrain(iq) = max( mxclfrrain(iq), cftmp )
          end if
          
        end if ! fluxrain>0.
      end do ! iq loop
    end if ! ncloud>=3

    do concurrent (iq = 1:imax)  
    
      ! store for aerosols
      qevap(iq,k) = qevap(iq,k) + evap(iq)
      prscav(iq,k) = prscav(iq,k) + tdt*0.24*fcol(iq)*pow75(Fr(iq))   !Strat only
      
    end do ! iq loop
    
    
    ! Liquid ------------------------------------------------------------------------------
    ! (Currently cloud droplet settling is negected, although included in UM and ACCESS 1.3)


    ! Misc ------------------------------------------------------------------------------

    if ( ncloud>=3 ) then  
      do concurrent (iq = 1:imax)    
        if ( fluxrain(iq)>0. ) then    
            
          ! Accretion of cloud ice by rain to produce snow or grauple (from Lin et al 1983 - praci)
          ! (Neglected in UM and ACCESS 1.3)
          rf = rhoi(iq,k)
          rn = fluxrain(iq)/dz(iq,k)
          if ( rf>1.e-15 .and. ttg(iq,k)<tfrz ) then
            if ( rn>qr0_crt ) then
              xwgt = 1.
            else
              xwgt = 0.  
            end  if
            slopes_r        = (( max(fluxrain(iq),0.)/max(crfra(iq),1.e-15)/tdt)**0.22)/714.  
            cdt             = tdt*pi*n0r*alin*gam380/4.*slopes_r**3.8*denfac(iq)
            drf             = max( min( crfra(iq)*rf, rf*cdt/(1.+0.5*cdt) ), 0. ) ! mass
            iflux           = drf*dz(iq,k)                                                                          ! flux
            rhoi(iq,k)      = rhoi(iq,k)      - drf
            fluxgraupel(iq) = fluxgraupel(iq) + iflux*xwgt
            fluxsnow(iq)    = fluxsnow(iq)    + iflux*(1.-xwgt)
            qaccf(iq,k)     = qaccf(iq,k)  + drf
            cftmp           = cifr(iq,k)*drf/rf
            cifr(iq,k)      = cifr(iq,k) - cftmp
            mxclfrgraupel(iq) = max( mxclfrgraupel(iq), cftmp*xwgt )
            mxclfrsnow(iq)    = max( mxclfrsnow(iq), cftmp*(1.-xwgt) )
          end if
          
        end if ! fluxrain>0.  
      end do ! iq loop
    end if ! ncloud>=3
  
   
    ! Update fluxes and area fractions for graupel, snow, ice and rain

    rhototf(:)       = rhog(:,k) + rhos(:,k) + rhoi(:,k)
    xfrac_graupel(:) = rhog(:,k)/max(rhototf(:),1.e-20)
    xfrac_snow(:)    = rhos(:,k)/max(rhototf(:),1.e-20)
    xfrac_ice(:)     = max( 0., 1.-xfrac_graupel(:)-xfrac_snow(:) )
    
    ! Melting and freezing
    fluxm(:,k) = fluxm(:,k) + fluxmelt(:)
    fluxf(:,k) = fluxf(:,k) + fluxfreeze(:)

    
    if ( ncloud>=3 ) then
      do concurrent (iq = 1:imax)  
        
        ! Grauple
        ! calculate maximum and random overlap for falling graupel
        pfstayice(iq,k) = pfstayice(iq,k) + fluxgraupel(iq)*(1.-fthrugraupel(iq))/tdt_in ! Save flux for the wet deposition scheme.  
        pqfsedice(iq,k) = pqfsedice(iq,k) + xfrac_graupel(iq)*foutgraupel(iq)*tdt/tdt_in ! Save sedimentation rate for aerosol scheme
        if ( fluxgraupel(iq)<=0. ) then
          rdclfrgraupel(iq) = 0.
          mxclfrgraupel(iq) = 0.
        end if
        mxclfrgraupel(iq) = max( mxclfrgraupel(iq), cfgraupel(iq,k) ) ! for rhogout
        cgfra(iq) = max( 1.e-15, mxclfrgraupel(iq)+rdclfrgraupel(iq)-mxclfrgraupel(iq)*rdclfrgraupel(iq) ) ! rnd overlap
        ! Compute fluxes into the box
        cffluxin = cgfra(iq) - cfgraupel(iq,k)
        rhogin   = fluxgraupel(iq)/dz(iq,k)
        ! Compute the fluxes of snow leaving the box
        cffluxout = cfgraupel(iq,k)*foutgraupel(iq)
        rhogout   = rhog(iq,k)*foutgraupel(iq)
        ! Update the rhos and cfsnow fields
        cfgraupel(iq,k) = cfgraupel(iq,k) - cffluxout + cffluxin*(1.-fthrugraupel(iq))
        rhog(iq,k)      = rhog(iq,k) - rhogout + rhogin*(1.-fthrugraupel(iq))
        fluxgraupel(iq) = max( rhogout*dz(iq,k) + fluxgraupel(iq)*fthrugraupel(iq), 0. )
        if ( fluxgraupel(iq)<1.e-20 ) then
          rhog(iq,k) = rhog(iq,k) + fluxgraupel(iq)/dz(iq,k)
          fluxgraupel(iq) = 0.
        end if
        ! Now fluxgraupel is flux leaving layer k
        fluxg(iq,k) = fluxg(iq,k) + fluxgraupel(iq)
      
        ! Snow
        ! calculate maximum and random overlap for falling snow
        pfstayice(iq,k) = pfstayice(iq,k) + fluxsnow(iq)*(1.-fthrusnow(iq))/tdt_in ! Save flux for the wet deposition scheme.
        pqfsedice(iq,k) = pqfsedice(iq,k) + xfrac_snow(iq)*foutsnow(iq)*tdt/tdt_in ! Save sedimentation rate for aerosol scheme
        if ( fluxsnow(iq)<=0. ) then
          rdclfrsnow(iq) = 0.
          mxclfrsnow(iq) = 0.
        end if
        mxclfrsnow(iq) = max( mxclfrsnow(iq), cfsnow(iq,k) ) ! for rhosout
        csfra(iq) = max( 1.e-15, mxclfrsnow(iq)+rdclfrsnow(iq)-mxclfrsnow(iq)*rdclfrsnow(iq) ) 
        ! Compute fluxes into the box
        cffluxin = csfra(iq) - cfsnow(iq,k)
        rhosin   = fluxsnow(iq)/dz(iq,k)
        ! Compute the fluxes of snow leaving the box
        cffluxout = cfsnow(iq,k)*foutsnow(iq)
        rhosout   = rhos(iq,k)*foutsnow(iq)
        ! Update the rhos and cfsnow fields
        cfsnow(iq,k) = cfsnow(iq,k) - cffluxout + cffluxin*(1.-fthrusnow(iq))
        rhos(iq,k)   = rhos(iq,k) - rhosout + rhosin*(1.-fthrusnow(iq))
        fluxsnow(iq) = max( rhosout*dz(iq,k) + fluxsnow(iq)*fthrusnow(iq), 0. )
        if ( fluxsnow(iq)<1.e-20 ) then
          rhos(iq,k) = rhos(iq,k) + fluxsnow(iq)/dz(iq,k)
          fluxsnow(iq) = 0.
        end if
        ! Now fluxsnow is flux leaving layer k
        fluxs(iq,k) = fluxs(iq,k) + fluxsnow(iq)
        
      end do ! iq loop  
    end if ! ncloud>=3

    
    do concurrent (iq = 1:imax)
    
      ! Ice
      ! calculate maximum and random overlap for falling ice
      pfstayice(iq,k) = pfstayice(iq,k) + fluxice(iq)*(1.-fthruice(iq))/tdt_in ! Save flux for the wet deposition scheme.
      pqfsedice(iq,k) = pqfsedice(iq,k) + xfrac_ice(iq)*foutice(iq)*tdt/tdt_in ! Save sedimentation rate for aerosol scheme
      if ( fluxice(iq)<=0. ) then
        rdclfrice(iq) = 0.
        mxclfrice(iq) = 0.
      end if
      mxclfrice(iq) = max( mxclfrice(iq), cifr(iq,k) ) ! for rhoiout
      cifra(iq) = max( 1.e-15, mxclfrice(iq)+rdclfrice(iq)-mxclfrice(iq)*rdclfrice(iq) ) !rnd overlap the mx and rd ice fractions
      ! Compute fluxes into the box
      cffluxin = cifra(iq) - cifr(iq,k)
      rhoiin   = fluxice(iq)/dz(iq,k)
      ! Compute the fluxes of ice leaving the box
      cffluxout = cifr(iq,k)*foutice(iq)
      rhoiout = rhoi(iq,k)*foutice(iq)
      ! Update the rhoi and cifr fields
      cifr(iq,k)  = min( 1.-clfr(iq,k), cifr(iq,k)-cffluxout+cffluxin*(1.-fthruice(iq)) )
      rhoi(iq,k)  = rhoi(iq,k) - rhoiout + rhoiin*(1.-fthruice(iq))
      fluxice(iq) = max( rhoiout*dz(iq,k) + fluxice(iq)*fthruice(iq), 0. )
      if ( fluxice(iq)<1.e-20 ) then
        rhoi(iq,k) = rhoi(iq,k) + fluxice(iq)/dz(iq,k)
        fluxice(iq) = 0.
      end if
      ! Now fluxice is flux leaving layer k
      fluxi(iq,k) = fluxi(iq,k) + fluxice(iq)
  
      ! Rain
      ! Calculate the raining cloud cover down to this level, for stratiform (crfra).
      pfstayliq(iq,k) = pfstayliq(iq,k) + fluxrain(iq)*(1.-fthruliq(iq))/tdt_in ! store liquid flux for aerosols
      if ( fluxrain(iq)<=0. ) then
        rdclfrrain(iq) = 0.
        mxclfrrain(iq) = 0.
      end if
      mxclfrrain(iq) = max( mxclfrrain(iq), cfrain(iq,k) ) ! for rhorout    
      crfra(iq) = max( 1.e-15, rdclfrrain(iq)+mxclfrrain(iq)-rdclfrrain(iq)*mxclfrrain(iq) )
      ! Compute fluxes into the box
      cffluxin = crfra(iq) - cfrain(iq,k)
      rhorin   = fluxrain(iq)/dz(iq,k)
      ! Compute the fluxes of rain leaving the box
      ! Use the flux-divergent form as in Rotstayn (QJRMS, 1997)
      cffluxout = cfrain(iq,k)*foutliq(iq)
      rhorout   = rhor(iq,k)*foutliq(iq)
      ! Update the rhor and cfrain fields
      cfrain(iq,k) = cfrain(iq,k) - cffluxout + cffluxin*(1.-fthruliq(iq))
      rhor(iq,k)   = rhor(iq,k) - rhorout + rhorin*(1.-fthruliq(iq))
      fluxrain(iq) = max( rhorout*dz(iq,k) + fluxrain(iq)*fthruliq(iq), 0. )
      if ( fluxrain(iq)<1.e-20 ) then
        rhor(iq,k) = rhor(iq,k) + fluxrain(iq)/dz(iq,k)
        fluxrain(iq) = 0.
      end if
      ! Now fluxrain is flux leaving layer k
      fluxr(iq,k) = fluxr(iq,k) + fluxrain(iq)
      
    end do ! iq loop  
    
  end do ! k loop
  
end do   ! n


! store precip, snow and graupel
precs(:) = precs + fluxr(:,1) + fluxi(:,1) + fluxs(:,1) + fluxg(:,1)
preci(:) = preci + fluxi(:,1) + fluxs(:,1)
precg(:) = precg + fluxg(:,1)

do concurrent (k = 1:kl)
  ! Re-create qtg, qrg, qlg, qfg, qsng and qgrg fields
  qtg(:,k)  = rhov(:,k)/rhoa(:,k)
  qrg(:,k)  = rhor(:,k)/rhoa(:,k)
  qfg(:,k)  = rhoi(:,k)/rhoa(:,k)
  qlg(:,k)  = rhol(:,k)/rhoa(:,k)
  qsng(:,k) = rhos(:,k)/rhoa(:,k)
  qgrg(:,k) = rhog(:,k)/rhoa(:,k)

  ! Remove small amounts of cloud and precip
  where ( qlg(:,k)<1.e-10 )
    qtg(:,k)  = qtg(:,k) + qlg(:,k)
    ttg(:,k)  = ttg(:,k) - hlcp*qlg(:,k)
    qlg(:,k)  = 0.
    clfr(:,k) = 0.
  end where
  where ( qfg(:,k)<1.e-10 )
    qtg(:,k)  = qtg(:,k) + qfg(:,k)
    ttg(:,k)  = ttg(:,k) - hlscp*qfg(:,k)
    qfg(:,k)  = 0.
    cifr(:,k) = 0.
  end where
  where ( qrg(:,k)<1.e-10 )
    qtg(:,k)    = qtg(:,k) + qrg(:,k)
    ttg(:,k)    = ttg(:,k) - hlcp*qrg(:,k)
    qrg(:,k)    = 0.
    cfrain(:,k) = 0.
  end where
  where ( qsng(:,k)<1.e-10 )
    qtg(:,k)    = qtg(:,k) + qsng(:,k)
    ttg(:,k)    = ttg(:,k) - hlscp*qsng(:,k)
    qsng(:,k)   = 0.
    cfsnow(:,k) = 0.
  end where
  where ( qgrg(:,k)<1.e-10 )
    qtg(:,k)       = qtg(:,k) + qgrg(:,k)
    ttg(:,k)       = ttg(:,k) - hlscp*qgrg(:,k)
    qgrg(:,k)      = 0.
    cfgraupel(:,k) = 0.
  end where
  qtg(:,k) = max( qtg(:,k), 0. )
    
  stratcloud(:,k) = clfr(:,k) + cifr(:,k)
  
end do  

#ifndef GPU
!      Adjust cloud fraction (and cloud cover) after precipitation
if ( nmaxpr==1 .and. mydiag ) then
  write(6,*) 'diags from newrain for idjd ',idjd
  diag_temp(:) = stratcloud(idjd,:)
  write (6,"('stratcloud',9f8.3/6x,9f8.3)") diag_temp
  diag_temp(:) = cfrain(idjd,:)
  write (6,"('cfrain    ',9f8.3/6x,9f8.3)") diag_temp
  diag_temp(:) = cfsnow(idjd,:)
  write (6,"('cfsnow    ',9f8.3/6x,9f8.3)") diag_temp
  diag_temp(:) = cfgraupel(idjd,:)
  write (6,"('cfgraupel ',9f8.3/6x,9f8.3)") diag_temp
end if

! Diagnostics for debugging
if ( diag .and. mydiag ) then
  diag_temp(:) = stratcloud(idjd,:)
  write(6,*) 'stratcloud',diag_temp
  diag_temp(:) = cifr(idjd,:)
  write(6,*) 'cifr',diag_temp
  diag_temp(:) = clfr(idjd,:)
  write(6,*) 'clfr',diag_temp
  diag_temp(:) = ttg(idjd,:)
  write(6,*) 'ttg',diag_temp
  diag_temp(:) = qsatg(idjd,:)
  write(6,*) 'qsatg',diag_temp         
  diag_temp(:) = qlg(idjd,:)
  write(6,*) 'qlg',diag_temp
  diag_temp(:) = qfg(idjd,:)
  write(6,*) 'qfg',diag_temp
  diag_temp(:) = qrg(idjd,:)
  write(6,*) 'qrg',diag_temp
  diag_temp(:) = qsng(idjd,:)
  write(6,*) 'qsng',diag_temp
  diag_temp(:) = qgrg(idjd,:)
  write(6,*) 'qgrg',diag_temp
  diag_temp(:) = qsubl(idjd,:)
  write(6,*) 'qsubl',diag_temp
  diag_temp(:) = rhoa(idjd,:)
  write(6,*) 'rhoa',diag_temp
  diag_temp(:) = rhos(idjd,:)
  write(6,*) 'rhos',diag_temp
  diag_temp(:) = fluxs(idjd,:)
  write(6,*) 'fluxs ',diag_temp
  !diag_temp(1:kl-1) = foutice(idjd,1:kl-1)
  !write(6,*) 'foutice',diag_temp(1:kl-1)
  !diag_temp(1:kl-1) = fthruice(idjd,1:kl-1)
  !write(6,*) 'fthruice',diag_temp(1:kl-1)
  diag_temp(:) = pqfsedice(idjd,:)
  write(6,*) 'pqfsedice',diag_temp
  diag_temp(:) = fluxm(idjd,:)
  write(6,*) 'fluxm',diag_temp
  write(6,*) 'cifra,fluxsnow',cifra(idjd),fluxsnow(idjd)
end if  ! (diag.and.mydiag)
#endif

return
end subroutine newsnowrain
    
subroutine progcloud(dt,qc,qtot,press,rho,fice,qs,t,rhcrit, &
                     dpsldt,nettend,stratcloud,imax,kl)
!$acc routine vector

use const_phys                    ! Physical constants
use parm_m, only : qgmin          ! Model configuration

implicit none

integer, intent(in) :: imax, kl
integer k
real, dimension(imax,kl), intent(inout) :: qc ! condensate = qf + ql
real, dimension(imax,kl), intent(in) :: qtot, rho, fice, qs, t, rhcrit, press
real, dimension(imax,kl), intent(in) :: dpsldt
real, dimension(imax,kl), intent(inout) :: nettend
real, dimension(imax,kl), intent(inout) :: stratcloud
real, dimension(imax) :: aa, bb, cc, at, a_dt, b_dt, cf1, cfeq, cfbar
real, dimension(imax) :: qv, omega, hlrvap, dqsdT, gamma, xf, dqs
real, intent(in) :: dt
real erosion_scale
real, parameter :: u00ramp = 0.01

! background erosion scale in 1/secs
erosion_scale = 1.E-6

!if ( ncloud>=5 ) then
!  ! convert convective mass flux from half levels to full levels
!  do k = 1,kl-1
!    cmflx(:,k) = rathb(k)*fluxtot(:,k)+ratha(k)*fluxtot(:,k+1)
!  end do
!  cmflx(:,kl) = rathb(kl)*fluxtot(:,kl)
!else ! ncloud==4
!  ! use convective area fraction in leoncld.f, instead of convective mass flux
!  cmflx = 0.
!end if

! calculate dqs = ((omega + grav*Mc)/(cp*rho)+nettend)*dqsdT*dt
!                 -------------------------------------------------------
!                 1 + (stratcloud + 0.5*da)*gamma
! MJT notes - GFDL AM adds (stratcloud+0.5*at*da)*gamma term

! Change in saturated volume fraction
! da = -0.5*(1.-cf)^2*dqs/(qs-qv)
! MJT notes - Tiedtke 93 does not use 0.5

! gamma = L/cp*dqsdT

! Follow GFDL AM approach since da=da(dqs), hence need to solve the above
! quadratic equation for dqs if da/=0

! dqs*dqs*AA + dqs*BB + CC = 0
! AA = 0.25*gamma*(1-cf)^2/(qs-qv)
! BB = -(1+gamma*cf)
! CC = ((omega + grav*mflx)/(cp*rho)+netten)*dqsdT*dt

do k = 1,kl
  stratcloud(:,k) = max( min( stratcloud(:,k), 1. ), 0. )  
    
  qv = qtot(:,k) - qc(:,k)  
  ! calculate vertical velocity, dqs/dT and gamma 
  omega = press(:,k)*dpsldt(:,k)
  hlrvap = (hl+fice(:,k)*hlf)/rvap
  dqsdT = qs(:,k)*hlrvap/(t(:,k)**2)
  gamma = (hlcp+fice(:,k)*hlfcp)*dqsdT
  
  xf = max(min( (qv/qs(:,k) - rhcrit(:,k) - u00ramp ) / ( 2.*u00ramp ), 1. ), 0. ) ! MJT suggestion
  
  !cc = ((omega + grav*cmflx(:,k))/(cp*rho(:,k))+nettend(:,k))*dt*dqsdT
  cc = (omega/(cp*rho(:,k))+nettend(:,k))*dt*dqsdT ! neglect cmflx
  at = 1.-stratcloud(:,k)
  aa = 0.5*at*at/max( qs(:,k)-qv, 1.e-20 )
  bb = 1.+gamma*stratcloud(:,k)
  where ( cc<=0. .and. xf>0. )
    !dqs = ( bb - sqrt( bb*bb - 2.*gamma*xf*aa*cc ) ) / ( gamma*xf*aa ) ! GFDL style
    !dqs = min( dqs, cc/(1. + 0.5*bb) )                                 ! GFDL style
    dqs = 2.*cc/( bb + sqrt( bb*bb - 2.*gamma*xf*aa*cc ) ) ! alternative form of quadratic equation
                                                           ! note that aa and bb have been multipled by 2 and -1, respectively.
    ! Large scale cloud formation via condensation (A)
    a_dt = -xf*aa*dqs
  elsewhere
    ! da = 0, so dqs can be solved from a linear equation
    dqs = cc/bb
    ! Large scale cloud formation via condensation (A)
    a_dt = 0.
  end where

  ! Large scale cloud destruction via erosion (B)
  b_dt = stratcloud(:,k)*erosion_scale*dt*max(qs(:,k)-qv, 1.e-20)/max(qc(:,k), 1.e-20)

  ! Integrate
  !   dcf/dt = (1-cf)*A - cf*B
  ! to give (use cf' = A-cf*(A+B))
  !   cf(t=1) = cfeq + (cf(t=0) - cfeq)*exp(-(A+B)*dt)
  !   cfeq = A/(A+B)
  ! Average cloud fraction over the interval t=tau to t=tau+1
  !   cfbar = cfeq - (cf(t=1) - cf(t=0))/((A+B)*dt)
  ! cfeq is the equilibrum cloud fraction that is approached with
  ! a time scale of 1/(A+B)
  where ( a_dt>1.e-20 .or. b_dt>1.e-20 )
    cfeq  = a_dt/(a_dt+b_dt)
    cf1   = cfeq + (stratcloud(:,k) - cfeq)*exp(-a_dt-b_dt)
    cfbar = cfeq + (stratcloud(:,k) - cf1 )/(a_dt+b_dt)
  elsewhere
    cfeq  = stratcloud(:,k)
    cf1   = stratcloud(:,k)
    cfbar = stratcloud(:,k)
  end where

  ! Change in condensate
  ! dqc = -dqs*(stratcloud+0.5*da) = -dqs*cfbar
  ! MJT notes - missing erosion term -cfbar*erosion_scale*dt*(qs-qv)
  qc(:,k) = max(min( qc(:,k) - max(cfbar,1.e-20)*dqs, qtot(:,k)-qgmin ), 0. )

  ! Change in cloud fraction
  where ( qc(:,k)>1.e-20 )
    stratcloud(:,k) = max(min( cf1, 1.), 1.e-20 )
  elsewhere
    ! MJT notes - cloud fraction is maintained (da=0.) while condesate evaporates (dqc<0.) until
    ! the condesate dissipates
    stratcloud(:,k) = 0.
    qc(:,k) = 0.
  end where

  ! Reset tendency and mass flux for next time-step
  nettend(:,k) = 0.
  
end do

return
end subroutine progcloud    
    
pure function pow75_s(x) result(ans)
!$acc routine vector
implicit none
real, intent(in) :: x
real ans, y
y=sqrt(x)
ans=y*sqrt(y)
end function pow75_s

pure function pow75_v(x) result(ans)
!$acc routine vector
implicit none
real, dimension(:), intent(in) :: x
real, dimension(size(x)) :: ans, y
y=sqrt(x)
ans=y*sqrt(y)
end function pow75_v    
    
end module leoncld_mod
