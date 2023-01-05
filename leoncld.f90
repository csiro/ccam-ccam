! Conformal Cubic Atmospheric Model
    
! Copyright 2015-2022 Commonwealth Scientific Industrial Research Organisation (CSIRO)
    
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
! and intercept to be consistent with Rotstayn 97.

! ldr    = 0    Diagnosed cloud scheme (depreciated)
! ldr   /= 0    Prognostic cloud condensate (different ice fall speed options)
    
! ncloud = 0    Standard LDR cloud microphysics with water vapour, liquid cloud and ice cloud
! ncloud = 2    Same as ncloud=0, but with prognostic rain and modified cfrac
! ncloud = 3    Same as ncloud=2, but with prognostic graupel and snow, as well as modified cfrac
! ncloud = 4    Use prognostic cloud fraction based on Tiedtke from GFDL-CM3
! ncloud = 10   Same as ncloud=0 with Tiedtke from GFDL-CM3
! ncloud = 12   Same as ncloud=2 with Tiedtke from GFDL-CM3
! ncloud = 13   Same as ncloud=3 with Tiedtke from GFDL-CM3 (i.e., same as ncloud=4)
! ncloud = 20   Same as ncloud=3 with MG cloud
! ncloud = 21   Same as ncloud=3 with 2nd moment condensate
! ncloud = 22   Same as ncloud=3 with MG cloud and 2nd momement condensate
! ncloud = 100  Use Lin et al 2nd moment microphysics
! ncloud = 110  Same as ncloud=100 with Tiedtke from GFDL-CM3
! ncloud = 120  Same as ncloud=100 with MG cloud fraction
   
!                            Water vapour (qg)
!
!   Cloud water (qlg,cfrac)                      Cloud ice (qfg,cfrac)
!
!   Rain (qrg,rfrac)                             Snow (qsg,sfrac)         Graupel (qgrg,gfrac)

! qg, qlg, qfg, qrg, qsg and qgrg are mixing ratios (g/g) and cfrac, rfrac, sfrac, gfrac are area cover
! fractions
    
module leoncld_mod

implicit none

private
public leoncld_work
public rhow, rhoice
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

! The following are used in the Manton-Cotton rain parameterization
real, parameter :: Ec=0.55                    !Mean collection efficiency for cloud drops
real, parameter :: Aurate=0.104*9.80616*Ec/um !Part of rate constant

! Parameters related to snow
real, parameter :: rhosno=100. !Assumed density of snow in kg/m^3
real, parameter :: Eac=0.7     !Mean collection efficiency of ql by snow

! Parameters related to rain
real, parameter :: Ecol=0.7  !Mean collection efficiency of ql by rain

! Parameters related to cloud radiative properties
real, parameter :: aice=1.016 !Constant in Platt optical depth for ice (SI units)
real, parameter :: bice=0.68  !Constant in Platt optical depth for ice (SI units)


contains
    
! This subroutine is the interface for the LDR cloud microphysics
subroutine leoncld_work(condg,conds,condx,gfrac,ktsav,                                  &
#ifndef GPUPHYSICS    
                        ppfevap,ppfmelt,ppfprec,ppfsnow,ppfsubl,                        &
                        pplambs,ppmaccr,ppmrate,ppqfsedice,pprfreeze,pprscav,           &
#endif
                        precip,ps,qfg,qg,qgrg,qlg,qrg,qsng,rfrac,sfrac,t,               &
                        stratcloud,cdrop,fluxr,fluxm,fluxf,fluxi,fluxs,fluxg,qevap,     &
                        qsubl,qauto,qcoll,qaccr,vi,                                     &
                        idjd,mydiag,ncloud,nevapls,ldr,rcm,imax,kl)
!$acc routine vector

use const_phys                    ! Physical constants
use estab                         ! Liquid saturation function
use mgcloud_m , only : mg_2cond
                                  ! MG cloud microphysics            
use parm_m, only : iaero, nmaxpr, dt
                                  ! Model configuration
use sigs_m                        ! Atmosphere sigma levels

implicit none

integer, intent(in) :: idjd, ncloud, nevapls, ldr
integer, intent(in) :: imax, kl
integer, dimension(imax), intent(in) :: ktsav
real, dimension(imax,kl), intent(inout) :: gfrac, rfrac, sfrac
real, dimension(imax,kl), intent(inout) :: qg, qlg, qfg, qrg, qsng, qgrg
real, dimension(imax,kl), intent(inout) :: t
real, dimension(imax,kl), intent(inout) :: stratcloud, cdrop
#ifndef GPUPHYSICS
real, dimension(imax,kl), intent(out) :: ppfevap
real, dimension(imax,kl), intent(out) :: ppfmelt
real, dimension(imax,kl), intent(out) :: ppfprec
real, dimension(imax,kl), intent(out) :: ppfsnow
real, dimension(imax,kl), intent(out) :: ppfsubl
real, dimension(imax,kl), intent(out) :: pplambs
real, dimension(imax,kl), intent(out) :: ppmaccr
real, dimension(imax,kl), intent(out) :: ppmrate
real, dimension(imax,kl), intent(out) :: ppqfsedice
real, dimension(imax,kl), intent(out) :: pprfreeze
real, dimension(imax,kl), intent(out) :: pprscav
#endif
real, dimension(imax,kl), intent(out) :: fluxr
real, dimension(imax,kl), intent(out) :: fluxm
real, dimension(imax,kl), intent(out) :: fluxf
real, dimension(imax,kl), intent(out) :: fluxi
real, dimension(imax,kl), intent(out) :: fluxs
real, dimension(imax,kl), intent(out) :: fluxg
real, dimension(imax,kl), intent(out) :: qevap
real, dimension(imax,kl), intent(out) :: qsubl
real, dimension(imax,kl), intent(out) :: qauto
real, dimension(imax,kl), intent(out) :: qcoll
real, dimension(imax,kl), intent(out) :: qaccr
real, dimension(imax,kl), intent(out) :: vi
real, dimension(imax), intent(inout) :: condg
real, dimension(imax), intent(inout) :: conds
real, dimension(imax), intent(inout) :: condx
real, dimension(imax), intent(inout) :: precip
real, dimension(imax), intent(in) :: ps
real, intent(in) :: rcm
logical, intent(in) :: mydiag

integer, dimension(imax) :: kbase,ktop  !Bottom and top of convective cloud 
real, dimension(imax,kl) :: prf      !Pressure on full levels (hPa)
real, dimension(imax,kl) :: rhoa     !Air density (kg/m3)
real, dimension(imax,kl) :: dz       !Layer thickness (m)
real, dimension(imax,kl) :: ccov     !Cloud cover (may differ from cloud frac if vertically subgrid)
real, dimension(imax,kl) :: qcl      !Vapour mixing ratio inside convective cloud
real, dimension(imax,kl) :: qenv     !Vapour mixing ratio outside convective cloud
real, dimension(imax,kl) :: tenv     !Temperature outside convective cloud
real, dimension(imax) :: precs       !Amount of stratiform precipitation in timestep (mm)
real, dimension(imax) :: preci       !Amount of stratiform snowfall in timestep (mm)
real, dimension(imax) :: precg       !Amount of stratiform graupel in timestep (mm)
real, dimension(imax) :: wcon        !Convective cloud water content (in-cloud, prescribed)

integer k, iq
real, dimension(imax,kl) :: qaccf
#ifndef GPUPHYSICS
real, dimension(imax,kl) :: pqfsedice, pslopes, prscav
#endif
real, dimension(imax,kl) :: prf_temp
real, dimension(imax) :: fl, diag_temp
real invdt

! meterological fields
do k = 1,kl
  do iq = 1,imax
    prf_temp(iq,k) = ps(iq)*sig(k)
    prf(iq,k)      = 0.01*prf_temp(iq,k)    !ps is SI units
    rhoa(iq,k)     = prf_temp(iq,k)/(rdry*t(iq,k))        ! air density
    dz(iq,k)       = -rdry*dsig(k)*t(iq,k)/(grav*sig(k)) ! level thickness in metres 
    dz(iq,k)       = min( max(dz(iq,k), 1.), 2.e4 )
  end do
end do 

! default values
precs(:) = 0. ! rain
preci(:) = 0. ! snow
precg(:) = 0. ! graupel


!     Calculate precipitation and related processes
if ( ncloud==21 .or. ncloud==22 ) then
  call mg_2cond
else
  call newsnowrain(dt,rhoa,dz,prf,cdrop,t,qlg,qfg,qrg,qsng,qgrg,                    &
                   precs,qg,stratcloud,rfrac,sfrac,gfrac,preci,precg,qevap,qsubl,   &
                   qauto,qcoll,qaccr,qaccf,fluxr,fluxi,fluxs,fluxg,fluxm,           &
                   fluxf,                                                           &
#ifndef GPUPHYSICS
                   pqfsedice,pslopes,prscav,                                        &
#endif
                   vi,                                                              &
                   condx,ktsav,idjd,mydiag,ncloud,nevapls,ldr,rcm,imax,kl)
end if


#ifdef debug
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


#ifndef GPUPHYSICS
!--------------------------------------------------------------
! Store data needed by prognostic aerosol scheme
! MJT notes - invert levels for aerosol code
if ( abs(iaero)>=2 ) then
  invdt = 1./dt
  ppfprec(:,1) = 0.   !At TOA
  ppfmelt(:,1) = 0.   !At TOA
  ppfsnow(:,1) = 0.   !At TOA
  pprfreeze(:,1) = 0. !At TOA
  do k = 1,kl-1
    do iq = 1,imax
      ppfprec(iq,kl-k) = (fluxr(iq,k+1)+fluxm(iq,k)-fluxf(iq,k))*invdt     !flux *entering* layer k
      ppfmelt(iq,kl-k) = fluxm(iq,k)*invdt                                 !flux melting in layer k
      ppfsnow(iq,kl-k) = (fluxi(iq,k+1)+fluxs(iq,k+1)+fluxg(iq,k+1) &
                        -fluxm(iq,k)+fluxf(iq,k))*invdt                    !flux *entering* layer k
      pprfreeze(iq,kl-k) = fluxf(iq,k)*invdt                               !flux freezing in layer k
    end do
  end do
  do k = 1,kl
    do iq = 1,imax
      ppfevap(iq,kl-k+1)    = qevap(iq,k)*rhoa(iq,k)*dz(iq,k)*invdt
      ppfsubl(iq,kl-k+1)    = qsubl(iq,k)*rhoa(iq,k)*dz(iq,k)*invdt !flux sublimating or staying in k
      pplambs(iq,kl-k+1)    = pslopes(iq,k)
      ppmrate(iq,kl-k+1)    = (qauto(iq,k)+qcoll(iq,k))*invdt
      ppmaccr(iq,kl-k+1)    = qaccr(iq,k)*invdt
      ppqfsedice(iq,kl-k+1) = pqfsedice(iq,k)
      pprscav(iq,kl-k+1)    = prscav(iq,k)
    end do
  end do
end if
!--------------------------------------------------------------
#endif


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
                       fluxi,fluxs,fluxg,fluxm,fluxf,                                                     &
#ifndef GPUPHYSICS
                       pqfsedice,pslopes,prscav,                                                          &
#endif
                       vi,                                                                                &
                       condx,ktsav,idjd,mydiag,ncloud,nevapls,ldr,rcm,imax,kl)
!$acc routine vector

use const_phys                    ! Physical constants
use estab                         ! Liquid saturation function
use parm_m, only : diag, nmr, nmaxpr
                                  ! Model configuration

implicit none

integer, intent(in) :: idjd, ncloud, nevapls, ldr
integer, intent(in) :: imax, kl
real, intent(in) :: tdt_in, rcm
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
#ifndef GPUPHYSICS
real, dimension(imax,kl), intent(out) :: pqfsedice
real, dimension(imax,kl), intent(out) :: pslopes
real, dimension(imax,kl), intent(out) :: prscav
#endif
real, dimension(imax,kl), intent(out) :: fluxr
real, dimension(imax,kl), intent(out) :: fluxi
real, dimension(imax,kl), intent(out) :: fluxs
real, dimension(imax,kl), intent(out) :: fluxg
real, dimension(imax,kl), intent(out) :: fluxm
real, dimension(imax,kl), intent(out) :: fluxf
real, dimension(imax,kl), intent(out) :: vi
real, dimension(imax), intent(in) :: condx
real, dimension(imax), intent(inout) :: precs
real, dimension(imax), intent(inout) :: preci
real, dimension(imax), intent(inout) :: precg
integer, dimension(imax), intent(in) :: ktsav
logical, intent(in) :: mydiag

integer, parameter :: graupel = 1
integer, parameter :: snow = 2
integer, parameter :: ice = 3
integer, parameter :: rain = 4
integer, parameter :: ncldtr = 4

real, dimension(imax,kl) :: fluxautorain, fluxautosnow, fluxautograupel
real, dimension(imax,ncldtr,kl) :: flux3
real, dimension(imax,kl) :: cfautorain, cfautosnow, cfautograupel
real, dimension(imax,kl) :: rhov, rhol
real, dimension(imax,ncldtr,kl) :: rho
real, dimension(imax,kl) :: clfr,cifr,qsatg
real, dimension(imax,ncldtr,kl) :: cf3
real, dimension(imax,ncldtr) :: fthru, fout
real, dimension(imax,ncldtr) :: v2
real, dimension(imax,ncldtr) :: flux2
real, dimension(imax,ncldtr) :: rhoin, rhoout
real, dimension(imax,ncldtr) :: cffluxout
real, dimension(imax,ncldtr) :: cf2
real, dimension(imax,ncldtr) :: mxclfr, rdclfr
real, dimension(imax,ncldtr) :: cffluxin
real rg, rl, rn, rf, rs
real, dimension(imax) :: rhodz,evap,clrevap,fr,sublflux
real, dimension(imax) :: fcol
real alph, alphaf
real, dimension(imax) :: qpf
real, dimension(imax,kl) :: pk3
real, dimension(imax) :: pk
real, dimension(imax) :: csbsav
real n0s
real aprpr, bprpr, curly
real, dimension(imax) :: cfmelt, fluxmelt, fluxfreeze
real slopes_g, slopes_s, xwgt
real, dimension(imax) :: qsl
real, dimension(imax) :: denfac
real, dimension(imax) :: gam1, deles
real, dimension(kl) :: diag_temp

integer k, n, njumps, iq, nn
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

fluxr            = 0.
fluxi            = 0.
fluxs            = 0.
fluxg            = 0. 
fluxm            = 0.  
fluxf            = 0.
fluxautorain     = 0.
fluxautosnow     = 0.
fluxautograupel  = 0.
qevap            = 0.
qauto            = 0.
qcoll            = 0.
qsubl            = 0.
qaccr            = 0.
qaccf            = 0.
#ifndef GPUPHYSICS
pqfsedice        = 0.
prscav           = 0.  
pslopes          = 0.
#endif
pk3              = 100.*prf
qsatg            = qsati(pk3,ttg)
cifr             = qfg*stratcloud/max( qlg+qfg, 1.e-30 )
clfr             = qlg*stratcloud/max( qlg+qfg, 1.e-30 )
cfautorain       = 0.
cfautosnow       = 0.
cfautograupel    = 0.

! Use full timestep for autoconversion
!njumps = 1
tdt = tdt_in

do k = 1,kl-1
  do iq = 1,imax
    if ( clfr(iq,k)>0. ) then
      qcrit = (4.*pi/3.)*rhow*Rcm**3*cdrop(iq,k)/rhoa(iq,k)
      qcic  = qlg(iq,k)/clfr(iq,k) !In cloud value
      if ( qcic>=qcrit ) then
        Crate    = Aurate*rhoa(iq,k)*(rhoa(iq,k)/(cdrop(iq,k)*rhow))**(1./3.)
        ql1      = 1./(qcic**(-4./3.)+(4./3.)*Crate*tdt)**0.75
        ql1      = max( ql1, qcrit ) !Intermediate qlg after auto
        Frb      = dz(iq,k)*rhoa(iq,k)*(qcic-ql1)/tdt
        Frb      = min( Frb, 1.e10 ) ! prevent overflow
        cdts     = tdt*0.5*Ecol*0.24*Frb**0.75 ! old
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
if ( ncloud==3 .or. ncloud==4 .or. ncloud==13 ) then

  do k = 1,kl
    do iq = 1,imax
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
  
end if ! ( ncloud==3 .or. ncloud==4 .or. ncloud==13 )


! update density and area fractions
do k = 1,kl
  do iq = 1,imax
    cifr(iq,k) = stratcloud(iq,k)*qfg(iq,k)/max(qlg(iq,k)+qfg(iq,k),1.e-30 )
    clfr(iq,k) = max( stratcloud(iq,k)-cifr(iq,k), 0. )
    rhov(iq,k) = qtg(iq,k)*rhoa(iq,k)
    rhol(iq,k) = qlg(iq,k)*rhoa(iq,k)
    rho(iq,ice,k)     = qfg(iq,k)*rhoa(iq,k)
    rho(iq,rain,k)    = qrg(iq,k)*rhoa(iq,k)
    rho(iq,snow,k)    = qsng(iq,k)*rhoa(iq,k)
    rho(iq,graupel,k) = qgrg(iq,k)*rhoa(iq,k)
    cf3(iq,graupel,k) = cfgraupel(iq,k)
    cf3(iq,snow,k)    = cfsnow(iq,k)
    cf3(iq,rain,k)    = cfrain(iq,k)
  end do
end do

flux3 = 0.

#ifdef debug
if ( diag .and. mydiag ) then
  diag_temp(:) = stratcloud(idjd,:)
  write(6,*) 'stratcloud',diag_temp
  diag_temp(:) = cifr(idjd,:)
  write(6,*) 'cifr      ',diag_temp
  diag_temp(:) = clfr(idjd,:)
  write(6,*) 'clfr      ',diag_temp
  diag_temp(:) = cf3(idjd,:,rain)
  write(6,*) 'cfrain    ',diag_temp
  diag_temp(:) = cf3(idjd,:,snow)
  write(6,*) 'cfsnow    ',diag_temp
  diag_temp(:) = cf3(idjd,:,graupel) 
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
if ( ncloud==3 .or. ncloud==4 .or. ncloud==13 ) then
  njumps = int(tdt_in/(maxlintime+0.01)) + 1
  tdt    = tdt_in/real(njumps)
else
  njumps = 1
  tdt = tdt_in
end if

do n = 1,njumps

  do nn = 1,ncldtr
    do iq = 1,imax
      flux2(iq,nn)  = 0.
      mxclfr(iq,nn) = 0. ! max overlap fraction
      rdclfr(iq,nn) = 0. ! rnd overlap fraction
      cf2(iq,nn)    = 0. ! total fraction = mx+rd-mx*rd
    end do
  end do


  v2(:,graupel) = 0.1
  v2(:,snow)    = 0.1
  v2(:,ice)     = 0.1 ! Assume no cloud at top level
  vi(:,kl)      = v2(:,ice)
  v2(:,rain)    = 0.
  cf2(:,rain)   = 1.e-6


  ! Now work down through the levels...
  do k = kl-1,1,-1
  
    ! misc fields
    pk(:)     = 100.*prf(:,k)
    rhodz(:)  = rhoa(:,k)*dz(:,k)
    denfac(:) = sqrt(sfcrho/rhoa(:,k))
    fluxmelt(:)   = 0.
    fluxfreeze(:) = 0.
    cfmelt(:)     = 0.


    ! Detect max/random overlap clouds that are separated by a clear layer
    do nn = 1,ncldtr
      do iq = 1,imax
        if ( (stratcloud(iq,k)>=1.e-10.and.stratcloud(iq,k+1)<1.e-10) .or. nmr==0 ) then
          rdclfr(iq,nn) = rdclfr(iq,nn) + mxclfr(iq,nn) - rdclfr(iq,nn)*mxclfr(iq,nn)
          mxclfr(iq,nn) = 0.
        end if
      end do
    end do
    mxclfr(:,rain) = max( mxclfr(:,rain), cfmelt(:) )    
    cf2(:,:) = max( rdclfr(:,:) + mxclfr(:,:) - rdclfr(:,:)*mxclfr(:,:), 1.e-10 )
    

    if ( ncloud==3 .or. ncloud==4 .or. ncloud==13 ) then
  
      ! Graupel ---------------------------------------------------------------------------
      do iq = 1,imax

        sublflux(iq) = 0.
        flux2(iq,graupel) = flux2(iq,graupel) + fluxautograupel(iq,k)*tdt/tdt_in

        ! graupel fall speed (from Lin et al 1983 - see GFDL AM3)

        rg = max( flux2(iq,graupel)/dz(iq,k), 0. )
        if ( cf2(iq,graupel)>=1.e-10 ) then
          v2(iq,graupel) = max( 0.1, 5.34623815*(rg/cf2(iq,graupel))**0.125 )
        end if

        ! Set up the parameters for the flux-divergence calculation
        alph = tdt*v2(iq,graupel)/dz(iq,k)
        alph = max( min( alph, 50. ), 0. )
        fout(iq,graupel)  = 1. - exp(-alph)             !analytical
        fthru(iq,graupel) = 1. - fout(iq,graupel)/alph  !analytical

      end do
        
      if ( any(flux2(:,graupel)>0.) ) then        

        do iq = 1,imax

          alphaf   = hls*qsatg(iq,k)/(rvap*ttg(iq,k)**2)
          gam1(iq) = hlscp*alphaf !(L/cp)*dqsdt (HBG notation)
      
          ! Melt falling graupel (based on Lin et al 83)
          rg = max(flux2(iq,graupel), 0.)/dz(iq,k)
          if ( ttg(iq,k)>tfrz .and. rg>1.e-10 ) then
            slopes_g           = ( max( flux2(iq,graupel), 0. )/dz(iq,k)/(pi*n0g*rho_g))**0.25
            qvp                = rhov(iq,k)/rhoa(iq,k)
            cdt                = tdt*2.*pi*n0g/hlf*(tcond*(ttg(iq,k)-tfrz)/rhoa(iq,k)-vdifu*hl*(qsatg(iq,k)-qvp))              &
                                *(0.78*slopes_g**2+0.31*scm3*gam275*sqrt(gcon/visk)*slopes_g**2.75*sqrt(denfac(iq)))
            drf                = max( min( rg, cdt ), 0. )
            iflux              = min( drf*dz(iq,k), flux2(iq,graupel) ) ! flux of graupel
            drf                = iflux/dz(iq,k)                    ! mass of graupel
            dqf                = drf/rhoa(iq,k)                    ! mixing ratio of graupel
            fluxmelt(iq)       = fluxmelt(iq)    + iflux
            flux2(iq,graupel)  = flux2(iq,graupel)    - iflux
            dttg               = -hlfcp*dqf
            ttg(iq,k)          = ttg(iq,k) + dttg
            qsatg(iq,k)        = qsatg(iq,k) + gam1(iq)*dttg/hlscp
            rdclfr(iq,graupel) = rdclfr(iq,graupel)*(1.-drf/rg)
            mxclfr(iq,graupel) = mxclfr(iq,graupel)*(1.-drf/rg)
            cftmp              = mxclfr(iq,graupel) + rdclfr(iq,graupel) - mxclfr(iq,graupel)*rdclfr(iq,graupel)
            cfmelt(iq)         = max( cfmelt(iq), max(cf2(iq,graupel)-cftmp,0.) )
            cf2(iq,graupel)    = cftmp
          end if
        
          ! Sublimation of graupel is neglected in the UM and ACCESS 1.3.
          ! (Currently treated the same as LDR97 ice sublimation)
          qvp = rhov(iq,k)/rhoa(iq,k)
          if ( flux2(iq,graupel)>0. .and. qvp<qsatg(iq,k) ) then ! sublime graupel
            slopes_g        = ( max(flux2(iq,graupel),0.)/dz(iq,k)/(pi*n0g*rho_g))**0.25
            fsclr_g         = max( (1.-cifr(iq,k)-clfr(iq,k))*flux2(iq,graupel), 0. )  
            cdt             = 2.*pi*vdifu*tcond*rvap*n0g*ttg(iq,k)**2                                                    &
                             *(0.78*slopes_g**2+0.31*scm3*gam275*sqrt(gcon/visk)*slopes_g**2.75*sqrt(denfac(iq))) &
                             /(tcond*rvap*ttg(iq,k)**2+hls**2*vdifu*qsatg(iq,k)*rhoa(iq,k))
            dqs             = tdt*cdt*(qsatg(iq,k)-qvp)
            dqs             = min( dqs, (qsatg(iq,k)-qvp)/(1.+gam1(iq)) ) !Don't supersat.
            sublflux(iq)    = min( dqs*rhodz(iq), fsclr_g ) ! flux of graupel
            drf             = sublflux(iq)/dz(iq,k)         ! mass of graupel
            dqs             = drf/rhoa(iq,k)                ! mixing ratio of graupel
            flux2(iq,graupel)    = flux2(iq,graupel) - sublflux(iq)
            fsclr_g         = fsclr_g      - sublflux(iq)
            rhov(iq,k)      = rhov(iq,k)  + drf        
            qsubl(iq,k)     = qsubl(iq,k) + dqs
            dttg            = -hlscp*dqs
            ttg(iq,k)       = ttg(iq,k) + dttg
            qsatg(iq,k)     = qsatg(iq,k) + gam1(iq)*dttg/hlscp
          end if
        
          ! Accretion of cloud liquid by falling graupel (from Lin et al 1983 - pgacw)
          ! This calculation uses the incoming flux2(:,graupel) without subtracting sublimation
          ! (since subl occurs only outside cloud), so add sublflux back to flux2(:,graupel).
          rl = rhol(iq,k)
          if ( flux2(iq,graupel)+sublflux(iq)>0. .and. rl>1.e-10 .and. ttg(iq,k)<tfrz ) then
            slopes_g        = ( max(flux2(iq,graupel)+sublflux(iq),0.)/dz(iq,k)/(pi*n0g*rho_g))**0.25
            cdt             = tdt*pi*n0g*gam350*gcon/4.0*slopes_g**3.5/sqrt(rhoa(iq,k))
            drl             = max( min( cf2(iq,graupel)*rl, rl*cdt/(1.+0.5*cdt) ), 0. ) ! mass of liquid
            lflux           = drl*dz(iq,k)           ! flux of liquid
            dql             = drl/rhoa(iq,k)         ! mixing ratio of liquid
            flux2(iq,graupel)    = flux2(iq,graupel) + lflux        
            rhol(iq,k)      = rhol(iq,k) - drl
            qaccr(iq,k)     = qaccr(iq,k) + dql
            dttg            = hlfcp*dql
            ttg(iq,k)       = ttg(iq,k) + dttg
            qsatg(iq,k)     = qsatg(iq,k) + gam1(iq)*dttg/hlscp
            cftmp           = clfr(iq,k)*drl/rl
            clfr(iq,k)      = clfr(iq,k) - cftmp
            mxclfr(iq,graupel)   = max( mxclfr(iq,graupel), cftmp )
          end if
        
          ! Accretion of rain by falling graupel (from Lin et al 1983 - pgacr)
          ! (Neglected in UM and ACCESS 1.3)
          rn = rho(iq,rain,k)
          if ( flux2(iq,graupel)+sublflux(iq)>0. .and. rn>1.e-10 .and. ttg(iq,k)<tfrz ) then
            slopes_g           = ( max( flux2(iq,graupel)+sublflux(iq), 0. )/dz(iq,k)/(pi*n0g*rho_g))**0.25
            slopes_r           = (( max( rn*dz(iq,k), 0. )/max( cf2(iq,rain),1.e-10 )/tdt)**0.22)/714.        
            qrn                = rn/rhoa(iq,k)            
            cdt                = tdt*pi*pi*n0g*n0r*abs(v2(iq,graupel)-v2(iq,rain))*qrn*(rho_r/rhoa(iq,k))   &
                                *(5.*slopes_r**6*slopes_g+2.*slopes_r**5*slopes_g**2                        &
                                +0.5*slopes_r**4*slopes_g**3)          
            drl                = max( min( cf2(iq,graupel)*rn, rn*cdt/(1.+0.5*cdt) ), 0. ) ! mass of rain
            lflux              = drl*dz(iq,k)   ! flux of rain
            dql                = drl/rhoa(iq,k) ! mixing ratio of rain
            flux2(iq,graupel)  = flux2(iq,graupel) + lflux
            rho(iq,rain,k)     = rho(iq,rain,k) - drl
            dttg               = hlfcp*dql
            ttg(iq,k)          = ttg(iq,k) + dttg
            qsatg(iq,k)        = qsatg(iq,k) + gam1(iq)*dttg/hlscp 
            cftmp              = cf3(iq,rain,k)*drl/rn
            cf3(iq,rain,k)     = cf3(iq,rain,k) - cftmp
            mxclfr(iq,graupel) = max( mxclfr(iq,graupel), cftmp )
          end if
        
          ! Accretion of cloud ice by falling graupel (from Lin et al 1983 - pgaci)
          ! (Neglected in UM and ACCESS 1.3)
          rf = rho(iq,ice,k)
          if ( flux2(iq,graupel)+sublflux(iq)>0. .and. rf>1.e-10 .and. ttg(iq,k)<tfrz ) then
            slopes_g           = ( max(flux2(iq,graupel)+sublflux(iq),0.)/dz(iq,k)/(pi*n0g*rho_g))**0.25
            cdt                = tdt*0.1*pi*n0g*gam350*gcon/4.*slopes_g**3.5/sqrt(rhoa(iq,k))
            drf                = max( min( cf2(iq,graupel)*rf, rf*cdt/(1.+0.5*cdt) ), 0. ) ! mass of ice
            iflux              = drf*dz(iq,k)    ! flux of ice
            dqf                = drf/rhoa(iq,k)  ! mixing ratio of ice
            flux2(iq,graupel)  = flux2(iq,graupel) + iflux
            rho(iq,ice,k)      = rho(iq,ice,k) - drf
            qaccf(iq,k)        = qaccf(iq,k) + dqf      
            cftmp              = cifr(iq,k)*drf/rf
            cifr(iq,k)         = cifr(iq,k) - cftmp
            mxclfr(iq,graupel) = max( mxclfr(iq,graupel), cftmp )
          end if
        
          ! Accretion of snow by falling graupel (from Lin et al 1983 - pgacs )
          rs = rho(iq,snow,k)
          if ( flux2(iq,graupel)+sublflux(iq)>0. .and. rs>1.e-10 .and. ttg(iq,k)<tfrz ) then
            n0s                = 2.e6*exp(-0.12*max(ttg(iq,k)-tfrz,-200.))        
            slopes_s           = ( max(rs,0.)/(pi*rho_s*n0s))**0.25
            slopes_g           = ( max(flux2(iq,graupel)+sublflux(iq),0.)/dz(iq,k)/(pi*n0g*rho_g))**0.25
            qsn                = rs/rhoa(iq,k)  
            cdt                = tdt*pi*pi*n0g*n0s*abs(v2(iq,graupel)-v2(iq,snow))*qsn*(rho_s/rhoa(iq,k))   &
                                *(5.*slopes_s**6*slopes_g+2.*slopes_s**5*slopes_g**2                        &
                                +0.5*slopes_s**4*slopes_g**3)        
            drf                = max( min( cf2(iq,graupel)*rs, rs*cdt/(1.+0.5*cdt) ), 0. ) ! mass of snow
            iflux              = drf*dz(iq,k)    ! flux of snow
            dqf                = drf/rhoa(iq,k)  ! mixing ratio of snow
            flux2(iq,graupel)  = flux2(iq,graupel) + iflux
            rho(iq,snow,k)     = rho(iq,snow,k) - drf
            qaccf(iq,k)        = qaccf(iq,k) + dqf
            cftmp              = cf3(iq,snow,k)*drf/rs
            cf3(iq,snow,k)     = cf3(iq,snow,k) - cftmp
            mxclfr(iq,graupel) = max( mxclfr(iq,graupel), cftmp )
          end if

        end do
       
      end if  ! flux2(:,graupel)>0.

     
      ! Snow ------------------------------------------------------------------------------
      do iq = 1,imax

        sublflux(iq) = 0.
        flux2(iq,snow) = flux2(iq,snow) + fluxautosnow(iq,k)*tdt/tdt_in
      
        ! Snow fall speed (from Lin et al 1983 - see GFDL AM3)
        rs = max( flux2(iq,snow)/dz(iq,k), 0. )
        if ( cf2(iq,snow)>=1.e-10 ) then
          v2(iq,snow) = max( 0.1, 1.82*(rs/cf2(iq,snow))**0.0625 )
        end if

        ! Set up the parameters for the flux-divergence calculation
        alph          = tdt*v2(iq,snow)/dz(iq,k)
        alph         = max( min( alph, 50. ), 0. )
        fout(iq,snow)  = 1. - exp(-alph)          !analytical
        fthru(iq,snow) = 1. - fout(iq,snow)/alph  !analytical
 
     end do

      if ( any( flux2(:,snow)>0. ) ) then

        do iq = 1,imax

          alphaf = hls*qsatg(iq,k)/(rvap*ttg(iq,k)**2)
          gam1(iq) = hlscp*alphaf !(L/cp)*dqsdt (HBG notation)
          
          ! Melt falling snow if > 0 deg C due to rain accretion
          ! (based on Lin et al 83, but using 0.65 and 0.44 coeffs following the UM approach)
          rs = max(flux2(iq,snow), 0.)/dz(iq,k)
          if ( ttg(iq,k)>tfrz .and. rs>1.e-10 ) then
            n0s             = 2.e6*exp(-0.12*max(ttg(iq,k)-tfrz,-200.))        
            slopes_s        = ( max(flux2(iq,snow),0.)/dz(iq,k)/(pi*rho_s*n0s))**0.25
            qvp             = rhov(iq,k)/rhoa(iq,k)  
            cdt             = tdt*2.*pi*n0s/hlf*(tcond*(ttg(iq,k)-tfrz)/rhoa(iq,k)-vdifu*hl*(qsatg(iq,k)-qvp))          &
                              *(0.65*slopes_s**2+0.44*scm3*gam263*sqrt(clin/visk)*slopes_s**2.63*sqrt(denfac(iq)))
            drf             = max( min( rs, cdt ), 0. ) 
            iflux           = min( drf*dz(iq,k), flux2(iq,snow) )   ! flux of snow
            drf             = iflux/dz(iq,k)                      ! mass of snow
            dqf             = drf/rhoa(iq,k)                      ! mixing ratio of snow
            fluxmelt(iq)    = fluxmelt(iq) + iflux
            flux2(iq,snow)  = flux2(iq,snow) - iflux
            dttg            = -hlfcp*dqf
            ttg(iq,k)       = ttg(iq,k) + dttg
            qsatg(iq,k)     = qsatg(iq,k) + gam1(iq)*dttg/hlscp
            rdclfr(iq,snow) = rdclfr(iq,snow)*(1.-drf/rs)
            mxclfr(iq,snow) = mxclfr(iq,snow)*(1.-drf/rs)
            cftmp           = mxclfr(iq,snow) + rdclfr(iq,snow) - mxclfr(iq,snow)*rdclfr(iq,snow)
            cfmelt(iq)      = max( cfmelt(iq), max( cf2(iq,snow)-cftmp, 0. ) )
            cf2(iq,snow)    = cftmp      
          end if
        
          ! Compute the sublimation of snow falling from level k+1 into level k
          ! (Currently treated the same as LDR97 ice sublimation - see UM and ACCESS 1.3)
          qvp = rhov(iq,k)/rhoa(iq,k)
          if ( flux2(iq,snow)>0. .and. qvp<qsatg(iq,k) ) then ! sublime snow
            n0s          = 2.e6*exp(-0.12*max(ttg(iq,k)-tfrz,-200.))        
            slopes_s     = ( max(flux2(iq,snow),0.)/dz(iq,k)/(pi*rho_s*n0s))**0.25
            fsclr_s      = max( (1.-cifr(iq,k)-clfr(iq,k))*flux2(iq,snow), 0. )  
            cdt          = 2.*pi*vdifu*tcond*rvap*n0s*ttg(iq,k)**2                                                 &
                               *(0.65*slopes_s**2+0.44*scm3*gam263*sqrt(clin/visk)*slopes_s**2.63*sqrt(denfac(iq))) &
                               /(tcond*rvap*ttg(iq,k)**2+hls**2*vdifu*qsatg(iq,k)*rhoa(iq,k))
            dqs          = tdt*cdt*(qsatg(iq,k)-qvp)
            dqs          = min( dqs, (qsatg(iq,k)-qvp)/(1.+gam1(iq)) ) !Don't supersat.
            sublflux(iq) = min( dqs*rhodz(iq), fsclr_s ) ! flux of snow
            drf          = sublflux(iq)/dz(iq,k)         ! mass of snow
            dqs          = drf/rhoa(iq,k)                ! mixing ratio of snow
            flux2(iq,snow) = flux2(iq,snow) - sublflux(iq)
            fsclr_s      = fsclr_s  - sublflux(iq)
            rhov(iq,k)   = rhov(iq,k)   + drf
            qsubl(iq,k)  = qsubl(iq,k)  + dqs
            dttg         = -hlscp*dqs
            ttg(iq,k)    = ttg(iq,k) + dttg
            qsatg(iq,k)  = qsatg(iq,k) + gam1(iq)*dttg/hlscp
          end if
        
          ! Accretion of cloud liquid by falling snow (from Lin et al 1983 - psacw)
          rl = rhol(iq,k)
          if ( flux2(iq,snow)+sublflux(iq)>0. .and. rl>1.e-10 .and. ttg(iq,k)<tfrz ) then
            n0s           = 2.e6*exp(-0.12*max(ttg(iq,k)-tfrz,-200.))        
            slopes_s      = ( max(flux2(iq,snow)+sublflux(iq),0.)/dz(iq,k)/(pi*rho_s*n0s))**0.25
            cdt           = tdt*denfac(iq)*pi*clin*gam325*n0s/4.*slopes_s**3.25
            drl           = max( min( cf2(iq,snow)*rl, rl*cdt/(1.+0.5*cdt) ), 0. ) ! mass of liquid
            lflux         = drl*dz(iq,k)                                                 ! flux of liquid
            dql           = drl/rhoa(iq,k)                                               ! mixing ratio of liquid
            flux2(iq,snow)  = flux2(iq,snow) + lflux
            rhol(iq,k)    = rhol(iq,k)   - drl
            qaccr(iq,k)   = qaccr(iq,k)  + dql
            dttg          = hlfcp*dql
            ttg(iq,k)     = ttg(iq,k) + dttg
            qsatg(iq,k)   = qsatg(iq,k) + gam1(iq)*dttg/hlscp
            cftmp         = clfr(iq,k)*drl/rl
            clfr(iq,k)    = clfr(iq,k) - cftmp
            mxclfr(iq,snow) = max( mxclfr(iq,snow), cftmp )
          end if
        
          ! Accretion of rain by falling snow to form snow (from Lin et al 1983 - psacr)
          rn = rho(iq,rain,k)
          if ( flux2(iq,snow)+sublflux(iq)>0. .and. rn>1.e-10 .and. ttg(iq,k)<tfrz ) then
            n0s             = 2.e6*exp(-0.12*max(ttg(iq,k)-tfrz,-200.))        
            slopes_s        = ( max(flux2(iq,snow)+sublflux(iq),0.)/dz(iq,k)/(pi*rho_s*n0s))**0.25
            slopes_r        = (( max(rn*dz(iq,k),0.)/max(cf2(iq,rain),1.e-10)/tdt)**0.22)/714.
            qrn             = rn/rhoa(iq,k)  
            cdt             = tdt*pi*pi*n0r*n0s*abs(v2(iq,snow)-v2(iq,rain))*qrn*(rho_r/rhoa(iq,k)) &
                              *(5.*slopes_r**6*slopes_s+2.*slopes_r**5*slopes_s**2                  &
                              +0.5*slopes_r**4*slopes_s**3)
            drl             = max( min( cf2(iq,rain)*rn, rn*cdt/(1.+0.5*cdt) ), 0. ) ! mass of rain
            lflux           = drl*dz(iq,k)                                                 ! flux of rain
            dql             = drl/rhoa(iq,k)                                               ! mixing ratio of rain
            flux2(iq,snow)  = flux2(iq,snow) + lflux
            rho(iq,rain,k)  = rho(iq,rain,k)   - drl
            dttg            = hlfcp*dql
            ttg(iq,k)       = ttg(iq,k) + dttg
            qsatg(iq,k)     = qsatg(iq,k) + gam1(iq)*dttg/hlscp  
            cftmp           = cf3(iq,rain,k)*drl/rn
            cf3(iq,rain,k)  = cf3(iq,rain,k) - cftmp
            mxclfr(iq,snow) = max( mxclfr(iq,snow), cftmp )
          end if
        
          ! Accretion of rain by falling snow to form graupel (neglected in Lin83 but included in UM)   
    
          ! Accretion of cloud ice by falling snow (from HDC 2004 - psaci)
          rf = rho(iq,ice,k)
          if ( flux2(iq,snow)+sublflux(iq)>0. .and. rf>1.e-10 .and. ttg(iq,k)<tfrz ) then
            n0s             = 2.e6*exp(-0.12*max(ttg(iq,k)-tfrz,-200.))        
            slopes_s        = ( max(flux2(iq,snow)+sublflux(iq),0.)/dz(iq,k)/(pi*rho_s*n0s))**0.25
            esi             = exp(0.05*max(ttg(iq,k)-tfrz,-100.))       ! efficiency
            cdt             = tdt*denfac(iq)*27.737*n0s*esi*slopes_s**3.41
            drf             = max( min( cf2(iq,snow)*rf, rf*cdt/(1.+0.5*cdt) ), 0. ) ! mass of ice
            iflux           = drf*dz(iq,k)                                                 ! flux of ice
            dqf             = drf/rhoa(iq,k)                                               ! mixing ratio of ice
            flux2(iq,snow)  = flux2(iq,snow) + iflux
            rho(iq,ice,k)   = rho(iq,ice,k)   - drf
            qaccf(iq,k)     = qaccf(iq,k)  + dqf
            cftmp           = cifr(iq,k)*drf/rf
            cifr(iq,k)      = cifr(iq,k) - cftmp
            mxclfr(iq,snow) = max( mxclfr(iq,snow), cftmp )
          end if

        end do
        
      end if  ! flux2(iq,snow)>0.
     
    end if ! ncloud==3 .or. ncloud==4 .or. ncloud==13

        
    ! Ice ---------------------------------------------------------------------------------
    sublflux(:) = 0.
   
    ! Set up the rate constant for ice sublimation
    ! MJT notes - curly and Csbsav depend on v2(:,ice,k+1), so v2(:,ice,k) can be updated below
    do iq = 1,imax
      if ( nevapls==-1 .or. (nevapls==-2.and.condx(iq)>0..and.k<=ktsav(iq)) ) then
        ! Define the rate constant for sublimation of snow, omitting factor rho(:,ice,:)
        Csbsav(iq) = 0.
      else
        slopes_i = 1.6e3*10**(-0.023*(ttg(iq,k)-tfrz))
        es = qsatg(iq,k)*pk(iq)/epsil
        Aprpr = (hls/(rKa*ttg(iq,k)))*(hls/(rvap*ttg(iq,k))-1.)
        Bprpr = rvap*ttg(iq,k)/((Dva/pk(iq))*es)
        curly = 0.65*slopes_i**2+0.493*slopes_i*sqrt(slopes_i*v2(iq,ice)*rhoa(iq,k)/um) !Factor in curly brackets
       ! Define the rate constant for sublimation of snow, omitting factor rho(:,ice,:)
        Csbsav(iq) = 4.*curly/(rhoa(iq,k)*qsatg(iq,k)*(Aprpr+Bprpr)*pi*v2(iq,ice)*rho_s)
      end if
    end do
     
    ! Set up snow fall speed field
    select case(abs(ldr))
      case(1)
        where ( cifr(:,k)>=1.e-10 )
          v2(:,ice) = max( 0.1, 3.23*(max(rho(:,ice,k),0.)/cifr(:,k))**0.17 )  ! Ice fall speed from LDR 1997
        end where
      case(2)
        where ( cifr(:,k)>=1.e-10 )
          v2(:,ice) = 0.9*3.23*(rho(:,ice,k)/cifr(:,k))**0.17
        end where
      case(3)
        where ( cifr(:,k)>=1.e-10 )
          v2(:,ice) = max( 0.1, 2.05+0.35*log10(rho(:,ice,k)/rhoa(:,k)/cifr(:,k)) )
        end where
      case(4)
        where ( cifr(:,k)>=1.e-10 )
          v2(:,ice) = 1.4*3.23*(rho(:,ice,k)/cifr(:,k))**0.17
        end where
      case(5)
        where ( cifr(:,k)>=1.e-10 )  
          v2(:,ice) = max( 0.1, 3.29*(max( rho(:,ice,k), 0. )/cifr(:,k))**0.16 ) ! from Lin et al 1983 
        end where  
      case(11)
        ! following are alternative slightly-different versions of above
        ! used for I runs from 29/4/05 till 30/8/05
        ! for given qfg, large cifr implies small ice crystals, 
        ! with a small fall speed. 
        ! Note that for very small qfg, cifr is small.
        ! But rho(:,ice,:) is like qfg, so ratio should also be small and OK.
        v2(:,ice) = max( v2(:,ice), 3.23*(rho(:,ice,k)/max(cifr(:,k),1.e-30))**0.17 )
      case(22)
        v2(:,ice) = max( v2(:,ice), 0.9*3.23*(rho(:,ice,k)/max(cifr(:,k),1.e-30))**0.17 )
      case(33)
        ! following max gives v2(:,ice)=.1 for qfg=cifr=0
        v2(:,ice) = max( v2(:,ice), 2.05+0.35*log10(max(rho(:,ice,k)/rhoa(:,k),2.68e-36)/max(cifr(:,k),1.e-30)) )
      case(55)
        v2(:,ice) = max( v2(:,ice), 3.29*(max(rho(:,ice,k),0.)/cifr(:,k))**0.16 ) ! from Lin et al 1983   
    end select
    v2(:,ice) = max( v2(:,ice), 0.001 )  

    ! Set up the parameters for the flux-divergence calculation
    do iq = 1,imax
      alph          = tdt*v2(iq,ice)/dz(iq,k)
      alph          = max( min( alph, 50. ), 0. )
      fout(iq,ice)  = 1. - exp(-alph)         !analytical
      fthru(iq,ice) = 1. - fout(iq,ice)/alph  !analytical  
    end do 
      
    if ( any( flux2(:,ice)>0. ) ) then

      do iq = 1,imax

        alphaf = hls*qsatg(iq,k)/(rvap*ttg(iq,k)**2)
        gam1(iq) = hlscp*alphaf !(L/cp)*dqsdt (HBG notation)   
        
        ! Melt falling ice if > 0 deg C
        if ( ttg(iq,k)>tfrz .and. flux2(iq,ice)>0. ) then
          qif           = flux2(iq,ice)/rhodz(iq)      !Mixing ratio of ice
          fluxmelt(iq)  = fluxmelt(iq) + flux2(iq,ice)
          dttg          = -hlfcp*qif
          ttg(iq,k)     = ttg(iq,k) + dttg
          qsatg(iq,k)   = qsatg(iq,k) + gam1(iq)*dttg/hlscp
          cfmelt(iq)    = max( cfmelt(iq), cf2(iq,ice) )
          flux2(iq,ice)  = 0.
          cf2(iq,ice)    = 0.
          rdclfr(iq,ice) = 0.
          mxclfr(iq,ice) = 0.
        end if
      
        ! Compute the sublimation of ice falling from level k+1 into level k
        qvp = rhov(iq,k)/rhoa(iq,k)
        if ( flux2(iq,ice)>0. .and. qvp<qsatg(iq,k) ) then ! sublime ice
          fsclr_i      = (1.-cifr(iq,k)-clfr(iq,k))*flux2(iq,ice)  
          Csb          = Csbsav(iq)*flux2(iq,ice)/tdt
          bf           = 1. + 0.5*Csb*tdt*(1.+gam1(iq))
          dqs          = max( 0., tdt*(Csb/bf)*(qsatg(iq,k)-qvp) )
          dqs          = min( dqs, (qsatg(iq,k)-qvp)/(1.+gam1(iq)) ) !Don't supersat.
          sublflux(iq) = min( dqs*rhodz(iq), fsclr_i ) ! flux of ice
          drf          = sublflux(iq)/dz(iq,k)         ! mass of ice
          dqs          = drf/rhoa(iq,k)                ! mixing ratio of ice     
          flux2(iq,ice) = flux2(iq,ice) - sublflux(iq)
          fsclr_i      = fsclr_i - sublflux(iq)
          rhov(iq,k)   = rhov(iq,k)  + drf
          qsubl(iq,k)  = qsubl(iq,k) + dqs
          dttg         = -hlscp*dqs
          ttg(iq,k)    = ttg(iq,k) + dttg
          qsatg(iq,k)  = qsatg(iq,k) + gam1(iq)*dttg/hlscp
        end if
      
        ! Accretion of cloud liquid by falling ice (neglected in Lin et al 1983, but
        ! included in UM and ACCESS 1.3 as piacw)
        ! This calculation uses the incoming flux2(:,ice) without subtracting sublimation
        ! (since subl occurs only outside cloud), so add sublflux back to flux2(:,ice).
        rl = rhol(iq,k)
        if ( flux2(iq,ice)+sublflux(iq)>0. .and. rl>1.e-10 ) then
          slopes_i      = 1.6e3*10**(-0.023*(ttg(iq,k)-tfrz))
          cdt           = Eac*slopes_i*(flux2(iq,ice)+sublflux(iq))/(2.*rhosno)
          drl           = max( min( cf2(iq,ice)*rl, rl*cdt/(1.+0.5*cdt) ), 0. ) ! mass of liquid
          lflux         = drl*dz(iq,k)   ! flux of liquid
          dql           = drl/rhoa(iq,k) ! mixing ratio of liquid
          flux2(iq,ice)  = flux2(iq,ice) + lflux
          rhol(iq,k)    = rhol(iq,k)  - drl
          qaccr(iq,k)   = qaccr(iq,k) + dql
          dttg          = hlfcp*dql
          ttg(iq,k)     = ttg(iq,k) + dttg
          qsatg(iq,k)   = qsatg(iq,k) + gam1(iq)*dttg/hlscp
          cftmp         = clfr(iq,k)*drl/rl
          clfr(iq,k)    = clfr(iq,k) - cftmp
          mxclfr(iq,2)  = max( mxclfr(iq,2), cftmp )
        end if

      end do
      
      if ( ncloud==3 .or. ncloud==4 .or. ncloud==13 ) then
        ! Accretion of rain by falling ice to produce ice (from Lin et al 1983 - piacr)
        ! (see UM and ACCESS 1.3 piacr-c for an alternate formulation)
        do iq = 1,imax
          rn  = rho(iq,rain,k)
          if ( flux2(iq,ice)+sublflux(iq)>0. .and. rn>1.e-10 .and. ttg(iq,k)<tfrz ) then
            qf             = max(flux2(iq,ice)+sublflux(iq),0.)/rhodz(iq)  
            cdt            = tdt*denfac(iq)*c_piacr*qf/sqrt(rhoa(iq,k))
            drl            = max( min( cf2(iq,ice)*rn, rn*cdt/(1.+0.5*cdt) ), 0. ) ! mass of rain
            lflux          = drl*dz(iq,k)   ! flux of rain
            dql            = drl/rhoa(iq,k) ! mixing ratio of rain
            flux2(iq,ice)  = flux2(iq,ice) + lflux
            rho(iq,rain,k) = rho(iq,rain,k)  - drl
            dttg           = hlfcp*dql
            ttg(iq,k)      = ttg(iq,k) + dttg
            qsatg(iq,k)    = qsatg(iq,k) + gam1(iq)*dttg/hlscp
            cftmp          = cf3(iq,rain,k)*drl/rn
            cf3(iq,rain,k) = cf3(iq,rain,k) - cftmp
            mxclfr(iq,2)   = max( mxclfr(iq,2), cftmp )
          end if
        end do
      end if 
      
      ! Accretion of rain by falling ice to produce graupel (Neglected in Lin et al 1983)
      ! (see UM and ACCESS 1.3 piacr-g for an alternate formulation)
      
    end if ! flux2(iq,ice)>0.
    
#ifndef GPUPHYSICS
    ! store slope for aerosols
    do iq = 1,imax
      slopes_i      = 1.6e3*10**(-0.023*(ttg(iq,k)-tfrz))
      pslopes(iq,k) = pslopes(iq,k) + slopes_i*tdt/tdt_in  
    end do
#endif
    
    
    ! Rain --------------------------------------------------------------------------------
    evap(:) = 0.

    ! Add flux of melted snow to flux2(:,rain)
    flux2(:,rain) = flux2(:,rain) + fluxmelt(:) + fluxautorain(:,k)*tdt/tdt_in
      
    ! Calculate rain fall speed (MJT suggestion)
    if ( ncloud==2 .or. ncloud==3 .or. ncloud==4 .or. ncloud==12 .or. ncloud==13 ) then
      do iq = 1,imax
        Fr(iq)         = max( flux2(iq,rain)/tdt/max(cf2(iq,rain),1.e-10),0.)
        v2(iq,rain)    = max( 0.1, 11.3*Fr(iq)**(1./9.)/sqrt(rhoa(iq,k)) )  !Actual fall speed
        !v2(iq,rain)   = max( 0.1, 5./sqrt(rhoa(iq,k)) )                    !Nominal fall speed
        alph           = tdt*v2(iq,rain)/dz(iq,k)
        alph           = max( min( alph, 50. ), 0. )
        fout(iq,rain)  = 1. - exp(-alph)
        fthru(iq,rain) = 1. - fout(iq,rain)/alph
      end do
    else
      v2(:,rain) = 9.e9   
      fout(:,rain)  = 1.
      fthru(:,rain) = 1.
    end if
        
    if ( any( flux2(:,rain)>0. ) ) then    
    
      do iq = 1,imax
        alphaf = hls*qsatg(iq,k)/(rvap*ttg(iq,k)**2)
        gam1(iq) = hlscp*alphaf !(L/cp)*dqsdt (HBG notation)
      end do

      if ( ncloud==3 .or. ncloud==4 .or. ncloud==13 ) then
      
        do iq = 1,imax
          rn = max(flux2(iq,rain),0.)/dz(iq,k)
          if ( rn>1.e-10 .and. ttg(iq,k)<tfrz ) then
            slopes_r           = (( max(flux2(iq,rain),0.)/max(cf2(iq,rain),1.e-10)/tdt)**0.22)/714.
            ! MJT notes - limit temperature to -100 C to avoid overflow with single precision
            cdt                = tdt*20.e2*pi**2*n0r*(rho_r/rhoa(iq,k))*slopes_r**7 &
                                  *(exp(-0.66*max(ttg(iq,k)-tfrz,-100.))-1.)
            drl                = max( min( rn, rn*cdt/(1.+0.5*cdt) ), 0. )
            lflux              = min( drl*dz(iq,k), flux2(iq,rain) ) ! flux
            lflux              = min( lflux, rhodz(iq)*(tfrz-ttg(iq,k))/hlfcp ) ! do not overshoot tfrz
            drl                = lflux/dz(iq,k) ! mass
            dql                = drl/rhoa(iq,k) ! mixing ratio
            flux2(iq,rain)     = flux2(iq,rain)    - lflux
            flux2(iq,graupel)  = flux2(iq,graupel)    + lflux
            fluxfreeze(iq)     = fluxfreeze(iq)  + lflux
            dttg               = hlfcp*dql
            ttg(iq,k)          = ttg(iq,k) + dttg
            qsatg(iq,k)        = qsatg(iq,k) + gam1(iq)*dttg/hlscp
            rdclfr(iq,rain)    = rdclfr(iq,rain)*(1.-drl/rn)
            mxclfr(iq,rain)    = mxclfr(iq,rain)*(1.-drl/rn)
            cltmp              = mxclfr(iq,rain) + rdclfr(iq,rain) - mxclfr(iq,rain)*rdclfr(iq,rain)
            mxclfr(iq,graupel) = max( mxclfr(iq,graupel), max(cf2(iq,rain)-cltmp, 0.) )
            cf2(iq,rain)       = cltmp
          end if
        end do

      end if ! ncloud==3 .or. ncloud==4 .or. ncloud==13
          
      ! Evaporation of rain
      qpf(:)     = flux2(:,rain)/rhodz !Mix ratio of rain which falls into layer
      clrevap(:) = (1.-clfr(:,k)-cifr(:,k))*qpf
      where ( ttg(:,k)<tfrz .and. ttg(:,k)>=tice )
        qsl(:)   = qsatg(:,k) + epsil*esdiffx(ttg(:,k))/pk
      elsewhere
        qsl(:)   = qsatg(:,k)
      end where
      do iq = 1,imax
        qvp     = rhov(iq,k)/rhoa(iq,k)
        if ( flux2(iq,rain)>0. .and. cf2(iq,rain)>0. ) then
          es       = qsl(iq)*pk(iq)/epsil
          Apr      = (hl/(rKa*ttg(iq,k)))*(hl/(rvap*ttg(iq,k))-1.)
          Bpr      = rvap*ttg(iq,k)/((Dva/pk(iq))*es)
          Fr(iq)   = flux2(iq,rain)/tdt/max(cf2(iq,rain), 1.e-10)
          Cev      = cf2(iq,rain)*3.8e2*sqrt(Fr(iq)/rhoa(iq,k))/(qsl(iq)*(Apr+Bpr))
          dqsdt    = hl*qsl(iq)/(rvap*ttg(iq,k)**2)
          bl       = 1. + 0.5*Cev*tdt*(1.+hlcp*dqsdt)
          evap(iq) = tdt*(Cev/bl)*(qsl(iq)-qvp)
          satevap  = (qsl(iq)-qvp)/(1.+hlcp*dqsdt)  !Evap to saturate
          ! v2(:,rain)=11.3*Fr(iq)**(1./9.)/sqrt(rhoa(mg,k)) !Actual fall speed
          ! v2(:,rain)=5./sqrt(rhoa(mg,k))                   !Nominal fall speed
          evap(iq) = max( 0., min( evap(iq), satevap, clrevap(iq) ) )
        end if
      end do
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

      do iq = 1,imax

        drl        = evap(iq)*rhoa(iq,k) ! mass
        rhov(iq,k) = rhov(iq,k) + drl
        ttg(iq,k)  = ttg(iq,k) - hlcp*evap(iq)
        !frclr(iq)  = rhodz(iq)*(clrevap(iq)-evap(iq)) ! flux over tdt

        ! Now do the collection of liquid cloud by rain term (cf. pracc in Lin83).
        fcol(iq) = 0.
        Fr(iq) = 0.
        rl = rhol(iq,k)
        if ( flux2(iq,rain)>0. .and. rl>1.e-10 ) then
          Fr(iq)        = max(flux2(iq,rain),0.)/tdt/max(cf2(iq,rain),1.e-10)
          fcol(iq)      = cf2(iq,rain)
          cdt           = tdt*Ecol*0.24*fcol(iq)*Fr(iq)**0.75
          coll          = max( min( rhol(iq,k), rhol(iq,k)*cdt/(1.+0.5*cdt) ), 0. ) ! mass
          lflux         = coll*dz(iq,k)                                            ! flux
          dql           = coll/rhoa(iq,k)                                          ! mixing ratio
          flux2(iq,rain)  = flux2(iq,rain) + lflux
          rhol(iq,k)    = rhol(iq,k)   - coll
          qcoll(iq,k)   = qcoll(iq,k)  + dql
          cltmp         = clfr(iq,k)*coll/rl
          clfr(iq,k)    = clfr(iq,k) - cltmp
          mxclfr(iq,rain) = max( mxclfr(iq,rain), cltmp )
        end if
      
        ! subtract evaporated rain
        lflux        = evap(iq)*rhodz(iq)
        flux2(iq,rain) = max( flux2(iq,rain) - lflux, 0. ) !To avoid roundoff -ve's

      end do
      
      if ( ncloud==3 .or. ncloud==4 .or. ncloud==13 ) then
          
        ! Accretion of cloud snow by rain (from Lin et al 1983 - pracs)
        do iq = 1,imax
          rs = max( rho(iq,snow,k), 0. )
          if ( flux2(iq,rain)>0. .and. rs>1.e-10 .and. ttg(iq,k)>tfrz+1. ) then
            n0s             = 2.e6*exp(-0.12*max(ttg(iq,k)-tfrz,-200.))        
            slopes_s        = ( max(rs,0.)/(pi*rho_s*n0s))**0.25
            slopes_r        = (( max(flux2(iq,rain),0.)/max(cf2(iq,rain),1.e-10)/tdt)**0.22)/714.  
            qsn             = max( rs/rhoa(iq,k), 0. )  
            cdt             = tdt*pi*pi*n0r*n0s*abs(v2(iq,rain)-v2(iq,snow))*qsn*(rho_s/rhoa(iq,k))  &
                              *(5.*slopes_s**6*slopes_r+2.*slopes_s**5*slopes_r**2                   &
                              +0.5*slopes_s**4*slopes_r**3)
            drf             = max( min( cf2(iq,rain)*rs, rs*cdt/(1.+0.5*cdt) ), 0. ) ! mass
            lflux           = drf*dz(iq,k)                                     ! flux
            dqf             = drf/rhoa(iq,k)                                   ! mixing ratio
            flux2(iq,rain)  = flux2(iq,rain) + lflux
            rho(iq,snow,k)  = rho(iq,snow,k)   - drf
            dttg            = hlfcp*dqf
            ttg(iq,k)       = ttg(iq,k) - dttg
            qsatg(iq,k)     = qsatg(iq,k) - gam1(iq)*dttg/hlscp      
            cftmp           = cf3(iq,snow,k)*drf/rs
            cf3(iq,snow,k)  = cf3(iq,snow,k) - cftmp
            mxclfr(iq,rain) = max( mxclfr(iq,rain), cftmp )
          end if
        end do

      end if ! ncloud==3 .or. ncloud==4 .or. ncloud==13        
        
      ! store for aerosols
      qevap(:,k) = qevap(:,k) + evap
#ifndef GPUPHYSICS
      prscav(:,k) = prscav(:,k) + tdt*0.24*fcol*Fr**0.75   !Strat only
#endif
      
    end if ! flux2(:,rain)>0.

  
    ! Liquid ------------------------------------------------------------------------------
    ! (Currently cloud droplet settling is negected, although included in UM and ACCESS 1.3)


    ! Misc ------------------------------------------------------------------------------

    if ( any( flux2(:,rain)>0. ) ) then
    
      if ( ncloud==3 .or. ncloud==4 .or. ncloud==13 ) then  
        ! Accretion of cloud ice by rain to produce snow or grauple (from Lin et al 1983 - praci)
        ! (Neglected in UM and ACCESS 1.3)
        do iq = 1,imax
          rf = rho(iq,ice,k)
          rn = flux2(iq,rain)/dz(iq,k)
          if ( rn>qr0_crt ) then
            xwgt = 1.
          else
            xwgt = 0.  
          end  if
          if ( flux2(iq,rain)>0. .and. rf>1.e-10 .and. ttg(iq,k)<tfrz ) then
            slopes_r           = (( max(flux2(iq,rain),0.)/max(cf2(iq,rain),1.e-10)/tdt)**0.22)/714.  
            cdt                = tdt*pi*n0r*alin*gam380/4.*slopes_r**3.8*denfac(iq)
            drf                = max( min( cf2(iq,rain)*rf, rf*cdt/(1.+0.5*cdt) ), 0. ) ! mass
            iflux              = drf*dz(iq,k)                                                                          ! flux
            rho(iq,ice,k)      = rho(iq,ice,k)      - drf
            flux2(iq,graupel)  = flux2(iq,graupel) + iflux*xwgt
            flux2(iq,snow)     = flux2(iq,snow)    + iflux*(1.-xwgt)
            qaccf(iq,k)        = qaccf(iq,k)  + drf
            cftmp              = cifr(iq,k)*drf/rf
            cifr(iq,k)         = cifr(iq,k) - cftmp
            mxclfr(iq,graupel) = max( mxclfr(iq,graupel), cftmp*xwgt )
            mxclfr(iq,snow)    = max( mxclfr(iq,snow), cftmp*(1.-xwgt) )
          end if
        end do
        
      end if  
      
    end if  
  
   
    ! Update melting and freezing fluxes
    fluxm(:,k) = fluxm(:,k) + fluxmelt(:)
    fluxf(:,k) = fluxf(:,k) + fluxfreeze(:)
    
    ! Save velocity
    vi(:,k) = v2(:,ice)
    
#ifndef GPUPHYSICS    
    pqfsedice(:,k) = pqfsedice(:,k) + fout(:,ice)*tdt/tdt_in ! Save sedimentation rate for aerosol scheme
#endif
    cf3(:,ice,k) = cifr(:,k)


    ! Update graupel, snow, ice and rain    
    ! calculate maximum and random overlap
    do nn = 1,ncldtr
      do iq = 1,imax
        if ( flux2(iq,nn)<=0. ) then
          rdclfr(iq,nn) = 0.
          mxclfr(iq,nn) = 0.
        end if
        mxclfr(iq,nn) = max( mxclfr(iq,nn), cf3(iq,nn,k) ) ! for rhoout
        cf2(iq,nn) = max( 1.e-10, mxclfr(iq,nn)+rdclfr(iq,nn)-mxclfr(iq,nn)*rdclfr(iq,nn) ) ! rnd overlap
        cffluxin(iq,nn) = cf2(iq,nn) - cf3(iq,nn,k)
        ! Compute fluxes into the box
        rhoin(iq,nn)   = flux2(iq,nn)/dz(iq,k)
        ! Compute the fluxes of snow leaving the box
        ! Use the flux-divergent form as in Rotstayn (QJRMS, 1997)
        cffluxout(iq,nn) = cf3(iq,nn,k)*fout(iq,nn)
        rhoout(iq,nn) = rho(iq,nn,k)*fout(iq,nn)
        ! Update the rho and cf3 fields
        cf3(iq,nn,k) = cf3(iq,nn,k) - cffluxout(iq,nn) + cffluxin(iq,nn)*(1.-fthru(iq,nn))
        rho(iq,nn,k) = rho(iq,nn,k) - rhoout(iq,nn) + rhoin(iq,nn)*(1.-fthru(iq,nn))
        flux2(iq,nn) = max( rhoout(iq,nn)*dz(iq,k) + flux2(iq,nn)*fthru(iq,nn), 0. )
        if ( flux2(iq,nn)<1.e-6 ) then
          rho(iq,nn,k) = rho(iq,nn,k) + flux2(iq,nn)/dz(iq,k)
          flux2(iq,nn) = 0.
        end if
        ! Now flux2 is flux leaving layer k
        flux3(iq,nn,k) = flux3(iq,nn,k) + flux2(iq,nn)
      end do  
    end do
      

    cifr(:,k) = min( cf3(:,ice,k), 1.-clfr(:,k) )
    
  end do ! k loop
  
end do   ! n


do k = 1,kl
  do iq = 1,imax
    cfgraupel(iq,k) = cf3(iq,graupel,k) 
    cfsnow(iq,k)    = cf3(iq,snow,k) 
    cfrain(iq,k)    = cf3(iq,rain,k)
    fluxg(iq,k)     = flux3(iq,graupel,k)
    fluxs(iq,k)     = flux3(iq,snow,k)
    fluxi(iq,k)     = flux3(iq,ice,k)
    fluxr(iq,k)     = flux3(iq,rain,k)
  end do
end do 

! store precip, snow and graupel
precs(:) = precs + fluxr(:,1) + fluxi(:,1) + fluxs(:,1) + fluxg(:,1)
preci(:) = preci + fluxi(:,1) + fluxs(:,1)
precg(:) = precg + fluxg(:,1)

do k = 1,kl
  do iq = 1,imax

    ! Re-create qtg, qrg, qlg, qfg, qsng and qgrg fields
    qtg(iq,k)  = rhov(iq,k)/rhoa(iq,k)
    qlg(iq,k)  = rhol(iq,k)/rhoa(iq,k)
    qrg(iq,k)  = rho(iq,rain,k)/rhoa(iq,k)
    qfg(iq,k)  = rho(iq,ice,k)/rhoa(iq,k)
    qsng(iq,k) = rho(iq,snow,k)/rhoa(iq,k)
    qgrg(iq,k) = rho(iq,graupel,k)/rhoa(iq,k)

    ! Remove small amounts of cloud and precip
    if ( qlg(iq,k)<1.e-10 ) then
      qtg(iq,k)  = qtg(iq,k) + qlg(iq,k)
      ttg(iq,k)  = ttg(iq,k) - hlcp*qlg(iq,k)
      qlg(iq,k)  = 0.
      clfr(iq,k) = 0.
    end if
    if ( qfg(iq,k)<1.e-10 ) then
      qtg(iq,k)  = qtg(iq,k) + qfg(iq,k)
      ttg(iq,k)  = ttg(iq,k) - hlscp*qfg(iq,k)
      qfg(iq,k)  = 0.
      cifr(iq,k) = 0.
    end if
    if ( qrg(iq,k)<1.e-10 ) then
      qtg(iq,k)    = qtg(iq,k) + qrg(iq,k)
      ttg(iq,k)    = ttg(iq,k) - hlcp*qrg(iq,k)
      qrg(iq,k)    = 0.
      cfrain(iq,k) = 0.
    end if
    if ( qsng(iq,k)<1.e-10 ) then
      qtg(iq,k)    = qtg(iq,k) + qsng(iq,k)
      ttg(iq,k)    = ttg(iq,k) - hlscp*qsng(iq,k)
      qsng(iq,k)   = 0.
      cfsnow(iq,k) = 0.
    end if
    if ( qgrg(iq,k)<1.e-10 ) then
      qtg(iq,k)       = qtg(iq,k) + qgrg(iq,k)
      ttg(iq,k)       = ttg(iq,k) - hlscp*qgrg(iq,k)
      qgrg(iq,k)      = 0.
      cfgraupel(iq,k) = 0.
    end if
    qtg(iq,k) = max( qtg(iq,k), 0. )
    
    stratcloud(iq,k) = clfr(iq,k) + cifr(iq,k)

  end do
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
#endif

#ifdef debug
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
  diag_temp(:) = rho(idjd,snow,:)
  write(6,*) 'rhos',diag_temp
  diag_temp(:) = fluxs(idjd,:)
  write(6,*) 'fluxs ',diag_temp
  !diag_temp(1:kl-1) = foutice(idjd,1:kl-1)
  !write(6,*) 'foutice',diag_temp(1:kl-1)
  !diag_temp(1:kl-1) = fthruice(idjd,1:kl-1)
  !write(6,*) 'fthruice',diag_temp(1:kl-1)
#ifndef GPUPHYSICS
  diag_temp(:) = pqfsedice(idjd,:)
  write(6,*) 'pqfsedice',diag_temp
#endif
  diag_temp(:) = fluxm(idjd,:)
  write(6,*) 'fluxm',diag_temp
  write(6,*) 'cifra,fluxsnow',cf2(idjd,ice),fluxsnow(idjd)
end if  ! (diag.and.mydiag)
#endif

return
end subroutine newsnowrain
    
end module leoncld_mod
