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
                        ppfevap,ppfmelt,ppfprec,ppfsnow,ppfsubl,                        &
                        pplambs,ppmaccr,ppmrate,ppqfsedice,pprfreeze,pprscav,precip,    &
                        ps,qfg,qg,qgrg,qlg,qrg,qsng,rfrac,sfrac,t,                      &
                        stratcloud,cdrop,fluxr,fluxm,fluxf,fluxi,fluxs,fluxg,qevap,     &
                        qsubl,qauto,qcoll,qaccr,vi,vs,vg,                               &
                        idjd,mydiag,ncloud,nevapls,ldr,rcm,imax,kl,                     &
                        ppleo_pcaut,ppleo_psaut,ppleo_pgaut,ppleo_pgmlt,                &
                        ppleo_pgsub,ppleo_pgacw,ppleo_pgacr,ppleo_pgaci,                &
                        ppleo_pgacs,ppleo_psmlt,ppleo_pssub,ppleo_psacw,                &
                        ppleo_psacr,ppleo_psaci,ppleo_pimlt,ppleo_pisub,                &
                        ppleo_piacw,ppleo_piacr,ppleo_psure,ppleo_prevp,                &
                        ppleo_pracc,ppleo_pracs,ppleo_praci)

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
real, dimension(imax,kl), intent(out) :: vs
real, dimension(imax,kl), intent(out) :: vg
real, dimension(imax), intent(inout) :: condg
real, dimension(imax), intent(inout) :: conds
real, dimension(imax), intent(inout) :: condx
real, dimension(imax), intent(inout) :: precip
real, dimension(imax,kl), intent(inout) :: ppleo_pcaut,ppleo_psaut,ppleo_pgaut,ppleo_pgmlt,&
                                           ppleo_pgsub,ppleo_pgacw,ppleo_pgacr,ppleo_pgaci,&
                                           ppleo_pgacs,ppleo_psmlt,ppleo_pssub,ppleo_psacw,&
                                           ppleo_psacr,ppleo_psaci,ppleo_pimlt,ppleo_pisub,&
                                           ppleo_piacw,ppleo_piacr,ppleo_psure,ppleo_prevp,&
                                           ppleo_pracc,ppleo_pracs,ppleo_praci
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
real, dimension(imax,kl) :: pqfsedice, pslopes, prscav
real, dimension(imax) :: prf_temp, fl
real, dimension(imax) :: diag_temp
real invdt

! meterological fields
do k = 1,kl
  prf_temp(:) = ps*sig(k)
  prf(:,k)    = 0.01*prf_temp    !ps is SI units
  rhoa(:,k)   = prf_temp/(rdry*t(:,k))             ! air density
  dz(:,k)     = -rdry*dsig(k)*t(:,k)/(grav*sig(k)) ! level thickness in metres 
  dz(:,k)     = min( max(dz(:,k), 1.), 2.e4 )
end do
 
! default values
precs(:) = 0. ! rain
preci(:) = 0. ! snow
precg(:) = 0. ! graupel


!     Calculate precipitation and related processes
if ( ncloud==21 .or. ncloud==22 ) then
  call mg_2cond
else
  call newsnowrain(dt,rhoa,dz,prf,cdrop,t,qlg,qfg,qrg,qsng,qgrg,                                       &
                   precs,qg,stratcloud,rfrac,sfrac,gfrac,preci,precg,qevap,qsubl,                      &
                   qauto,qcoll,qaccr,qaccf,fluxr,fluxi,fluxs,fluxg,fluxm,                              &
                   fluxf,pqfsedice,pslopes,prscav,vi,vs,vg,                                            &
                   condx,ktsav,idjd,mydiag,ncloud,nevapls,ldr,rcm,imax,kl,                             &
                   ppleo_pcaut,ppleo_psaut,ppleo_pgaut,ppleo_pgmlt,ppleo_pgsub,ppleo_pgacw,ppleo_pgacr,&
                   ppleo_pgaci,ppleo_pgacs,ppleo_psmlt,ppleo_pssub,ppleo_psacw,ppleo_psacr,ppleo_psaci,&
                   ppleo_pimlt,ppleo_pisub,ppleo_piacw,ppleo_piacr,ppleo_psure,ppleo_prevp,ppleo_pracc,&
                   ppleo_pracs,ppleo_praci)
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
    ppfprec(:,kl+1-k) = (fluxr(:,k+1)+fluxm(:,k)-fluxf(:,k))*invdt     !flux *entering* layer k
    ppfmelt(:,kl+1-k) = fluxm(:,k)*invdt                               !flux melting in layer k
    ppfsnow(:,kl+1-k) = (fluxi(:,k+1)+fluxs(:,k+1)+fluxg(:,k+1) &
                        -fluxm(:,k)+fluxf(:,k))*invdt                  !flux *entering* layer k
    pprfreeze(:,kl+1-k) = fluxf(:,k)*invdt                             !flux freezing in layer k
  end do
  do k = 1,kl
    ppfevap(:,kl+1-k)    = qevap(:,k)*rhoa(:,k)*dz(:,k)*invdt
    ppfsubl(:,kl+1-k)    = qsubl(:,k)*rhoa(:,k)*dz(:,k)*invdt !flux sublimating or staying in k
    pplambs(:,kl+1-k)    = pslopes(:,k)
    ppmrate(:,kl+1-k)    = (qauto(:,k)+qcoll(:,k))*invdt
    ppmaccr(:,kl+1-k)    = qaccr(:,k)*invdt
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
                       fluxi,fluxs,fluxg,fluxm,fluxf,pqfsedice,pslopes,prscav,vi,vs,vg,                   &
                       condx,ktsav,idjd,mydiag,ncloud,nevapls,ldr,rcm,imax,kl,                            &
                       pleo_pcaut,pleo_psaut,pleo_pgaut,pleo_pgmlt,pleo_pgsub,pleo_pgacw,pleo_pgacr,      &
                       pleo_pgaci,pleo_pgacs,pleo_psmlt,pleo_pssub,pleo_psacw,pleo_psacr,pleo_psaci,      &
                       pleo_pimlt,pleo_pisub,pleo_piacw,pleo_piacr,pleo_psure,pleo_prevp,pleo_pracc,      &
                       pleo_pracs,pleo_praci)

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
real, dimension(imax,kl), intent(inout) :: pleo_pcaut,pleo_psaut,pleo_pgaut,pleo_pgmlt,&
                                           pleo_pgsub,pleo_pgacw,pleo_pgacr,pleo_pgaci,&
                                           pleo_pgacs,pleo_psmlt,pleo_pssub,pleo_psacw,&
                                           pleo_psacr,pleo_psaci,pleo_pimlt,pleo_pisub,&
                                           pleo_piacw,pleo_piacr,pleo_psure,pleo_prevp,&
                                           pleo_pracc,pleo_pracs,pleo_praci
real, dimension(imax,kl), intent(out) :: qevap
real, dimension(imax,kl), intent(out) :: qsubl
real, dimension(imax,kl), intent(out) :: qauto
real, dimension(imax,kl), intent(out) :: qcoll
real, dimension(imax,kl), intent(out) :: qaccr
real, dimension(imax,kl), intent(out) :: qaccf
real, dimension(imax,kl), intent(out) :: pqfsedice
real, dimension(imax,kl), intent(out) :: pslopes
real, dimension(imax,kl), intent(out) :: prscav
real, dimension(imax,kl), intent(out) :: fluxr
real, dimension(imax,kl), intent(out) :: fluxi
real, dimension(imax,kl), intent(out) :: fluxs
real, dimension(imax,kl), intent(out) :: fluxg
real, dimension(imax,kl), intent(out) :: fluxm
real, dimension(imax,kl), intent(out) :: fluxf
real, dimension(imax,kl), intent(out) :: vi
real, dimension(imax,kl), intent(out) :: vs
real, dimension(imax,kl), intent(out) :: vg
real, dimension(imax), intent(in) :: condx
real, dimension(imax), intent(inout) :: precs
real, dimension(imax), intent(inout) :: preci
real, dimension(imax), intent(inout) :: precg
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
real, dimension(imax) :: rhoiin,rhoiout,rhorin,rhorout
real, dimension(imax) :: rhosin,rhosout,rhogin,rhogout
real, dimension(imax) :: cffluxout
real, dimension(imax) :: crfra,cifra,csfra,cgfra
real, dimension(imax) :: mxclfrrain,rdclfrrain,mxclfrice,rdclfrice
real, dimension(imax) :: mxclfrsnow,rdclfrsnow,mxclfrgraupel,rdclfrgraupel
real, dimension(imax) :: cffluxin
real rg, rl, rn, rf, rs
real, dimension(imax) :: rhodz,evap,clrevap,fr,sublflux
real, dimension(imax) :: fcol
real alph
real, dimension(imax) :: alphaf, qpf
real, dimension(imax) :: pk, csbsav
real n0s
real, dimension(imax) :: aprpr, bprpr, curly
real, dimension(imax) :: cfmelt, fluxmelt, fluxfreeze
real slopes_g, slopes_s, xwgt
real, dimension(imax) :: qsl
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

fluxr           = 0.
fluxi           = 0.
fluxs           = 0.
fluxg           = 0. 
fluxm           = 0.  
fluxf           = 0.
fluxautorain    = 0.
fluxautosnow    = 0.
fluxautograupel = 0.
qevap           = 0.
qauto           = 0.
qcoll           = 0.
qsubl           = 0.
qaccr           = 0.
qaccf           = 0.
pqfsedice       = 0.
prscav          = 0.  
pslopes         = 0.
do k = 1,kl
  pk(:)         = 100.*prf(:,k)
  qsatg(:,k)    = qsati(pk(:),ttg(:,k),imax)
end do
cifr            = qfg*stratcloud/max( qlg+qfg, 1.e-30 )
clfr            = qlg*stratcloud/max( qlg+qfg, 1.e-30 )
cfautorain      = 0.
cfautosnow      = 0.
cfautograupel   = 0.

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
        pleo_pcaut(iq,k) = pleo_pcaut(iq,k) + dqls / tdt            ! sny 01: autoconversion to cloud water 
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
        pleo_psaut(iq,k)   = pleo_psaut(iq,k) + dqfs / tdt         ! sny 02: autoconversion of ice to snow
      end if
    
      ! autoconversion of snow to graupel (from Lin et al 1983)
      if ( qsng(iq,k)*rhoa(iq,k)>qs0_crt ) then
        qfs  = max( qsng(iq,k)-qs0_crt/rhoa(iq,k), 0. )
        cdts = tdt*1.e-3*exp(0.09*(ttg(iq,k)-tfrz))
        dqfs = max( min( qsng(iq,k), qfs*cdts ), 0.) 
        cfautograupel(iq,k)   = cfsnow(iq,k)
        qsng(iq,k)            = qsng(iq,k) - dqfs
        fluxautograupel(iq,k) = dqfs*rhoa(iq,k)*dz(iq,k)
        pleo_pgaut(iq,k)   = pleo_pgaut(iq,k) + dqfs / tdt         ! sny 03: autoconversion of snow to graupel
      end if

    end do ! iq loop 
  end do   ! k loop
  
end if ! ( ncloud==3 .or. ncloud==4 .or. ncloud==13 )

! update density and area fractions
do k = 1,kl
  cifr(:,k) = stratcloud(:,k)*qfg(:,k)/max(qlg(:,k)+qfg(:,k),1.e-30 )
  clfr(:,k) = max( stratcloud(:,k)-cifr(:,k), 0. )
  rhov(:,k) = qtg(:,k)*rhoa(:,k)
  rhoi(:,k) = qfg(:,k)*rhoa(:,k)
  rhol(:,k) = qlg(:,k)*rhoa(:,k)
  rhor(:,k) = qrg(:,k)*rhoa(:,k)
  rhos(:,k) = qsng(:,k)*rhoa(:,k)
  rhog(:,k) = qgrg(:,k)*rhoa(:,k)
end do


#ifdef debug
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
if ( ncloud==3 .or. ncloud==4 .or. ncloud==13 ) then
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
  vg(:,kl)         = vg2(:)

  fluxsnow(:)   = 0.
  mxclfrsnow(:) = 0. ! max overlap snow fraction
  rdclfrsnow(:) = 0. ! rnd overlap snow fraction
  csfra(:)      = 0. ! total snow fraction = mx+rd-mx*rd
  vs2(:)        = 0.1
  vs(:,kl)      = vs2(:)

  fluxice(:)   = 0.
  mxclfrice(:) = 0. ! max overlap ice fraction
  rdclfrice(:) = 0. ! rnd overlap ice fraction
  cifra(:)     = 0. ! total ice fraction = mx+rd-mx*rd 
  vi2(:)       = 0.1 ! Assume no cloud at top level
  vi(:,kl)     = vi2(:)

  fluxrain(:)   = 0.
  mxclfrrain(:) = 0. ! max overlap rain fraction
  rdclfrrain(:) = 0. ! rnd overlap rain fraction
  crfra(:)      = 1.e-6 ! total rain fraction = mx+rd-mx*rd
  vr2(:)        = 0.


  ! Now work down through the levels...
  do k = kl-1,1,-1
  
    ! misc fields
    do iq = 1,imax
      pk(iq)     = 100.*prf(iq,k)
      rhodz(iq)  = rhoa(iq,k)*dz(iq,k)
      denfac(iq) = sqrt(sfcrho/rhoa(iq,k))
      fluxmelt(iq)   = 0.
      fluxfreeze(iq) = 0.
      cfmelt(iq)     = 0.
    end do
    
    if ( ncloud==3 .or. ncloud==4 .or. ncloud==13 ) then
  
      ! Graupel ---------------------------------------------------------------------------
      sublflux(:) = 0.
      fluxgraupel(:) = fluxgraupel + fluxautograupel(:,k)*tdt/tdt_in
      
      ! Detect max/random overlap clouds that are separated by a clear layer
      where ( (stratcloud(:,k)>=1.e-10.and.stratcloud(:,k+1)<1.e-10) .or. nmr==0 )
        rdclfrgraupel(:) = rdclfrgraupel + mxclfrgraupel - rdclfrgraupel*mxclfrgraupel
        mxclfrgraupel(:) = 0.
      end where
      cgfra(:) = max( rdclfrgraupel + mxclfrgraupel - rdclfrgraupel*mxclfrgraupel, 1.e-10 )
       
      ! graupel fall speed (from Lin et al 1983 - see GFDL AM3)
      do iq = 1,imax
        rg = max( fluxgraupel(iq)/dz(iq,k), 0. )
        if ( cgfra(iq)>=1.e-10 ) then
          vg2(iq) = max( 0.1, 5.34623815*(rg/cgfra(iq))**0.125 )
        end if
      end do

      ! Set up the parameters for the flux-divergence calculation
      do iq = 1,imax
        alph         = tdt*vg2(iq)/dz(iq,k)
        alph         = max( min( alph, 50. ), 0. )
        foutgraupel(iq)  = 1. - exp(-alph)        !analytical
        fthrugraupel(iq) = 1. - foutgraupel(iq)/alph  !analytical
      end do

      alphaf(:) = hls*qsatg(:,k)/(rvap*ttg(:,k)**2)
      gam1(:)   = hlscp*alphaf !(L/cp)*dqsdt (HBG notation)
        
      if ( any(fluxgraupel>0.) ) then        
      
        ! Melt falling graupel (based on Lin et al 83)
        do iq = 1,imax
          rg = max(fluxgraupel(iq), 0.)/dz(iq,k)
          if ( ttg(iq,k)>tfrz .and. rg>1.e-10 ) then
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
            pleo_pgmlt(iq,k)  = pleo_pgmlt(iq,k) + dqf             ! sny 04: Melt falling graupel                                 
          end if
        end do
        
        ! Sublimation of graupel is neglected in the UM and ACCESS 1.3.
        ! (Currently treated the same as LDR97 ice sublimation)
        do iq = 1,imax
          qvp = rhov(iq,k)/rhoa(iq,k)
          if ( fluxgraupel(iq)>0. .and. qvp<qsatg(iq,k) ) then ! sublime graupel
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
            pleo_pgsub(iq,k)  = pleo_pgsub(iq,k) + dqs              ! sny 05: Sublimation of graupel 
          end if
        end do
        
        ! Accretion of cloud liquid by falling graupel (from Lin et al 1983 - pgacw)
        ! This calculation uses the incoming fluxgraupel without subtracting sublimation
        ! (since subl occurs only outside cloud), so add sublflux back to fluxgraupel.
        do iq = 1,imax
          rl = rhol(iq,k)
          if ( fluxgraupel(iq)+sublflux(iq)>0. .and. rl>1.e-10 .and. ttg(iq,k)<tfrz ) then
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
            pleo_pgacw(iq,k)  = pleo_pgacw(iq,k) + dql                ! sny 06: Accretion of cloud liquid by falling graupel
          end if
        end do
        
        ! Accretion of rain by falling graupel (from Lin et al 1983 - pgacr)
        ! (Neglected in UM and ACCESS 1.3)
        do iq = 1,imax
          rn = rhor(iq,k)
          if ( fluxgraupel(iq)+sublflux(iq)>0. .and. rn>1.e-10 .and. ttg(iq,k)<tfrz ) then
            slopes_g        = ( max( fluxgraupel(iq)+sublflux(iq), 0. )/dz(iq,k)/(pi*n0g*rho_g))**0.25
            slopes_r        = (( max( rn*dz(iq,k), 0. )/max( crfra(iq),1.e-10 )/tdt)**0.22)/714.        
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
            pleo_pgacr(iq,k)  = pleo_pgacr(iq,k) + dql              ! sny 07: Accretion of rain by falling graupel
          end if
        end do 
        
        ! Accretion of cloud ice by falling graupel (from Lin et al 1983 - pgaci)
        ! (Neglected in UM and ACCESS 1.3)
        do iq = 1,imax
          rf = rhoi(iq,k)
          if ( fluxgraupel(iq)+sublflux(iq)>0. .and. rf>1.e-10 .and. ttg(iq,k)<tfrz ) then
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
            pleo_pgaci(iq,k)  = pleo_pgaci(iq,k) + dqf             ! sny 08:  Accretion of cloud ice by falling graupel
          end if
        end do
        
        ! Accretion of snow by falling graupel (from Lin et al 1983 - pgacs )
        do iq = 1,imax
          rs = rhos(iq,k)
          if ( fluxgraupel(iq)+sublflux(iq)>0. .and. rs>1.e-10 .and. ttg(iq,k)<tfrz ) then
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
            pleo_pgacs(iq,k)  = pleo_pgacs(iq,k) + dqf             ! sny 09: Accretion of snow by falling graupel
          end if
        end do
       
      end if  ! fluxgraupel>0.

     
      ! Snow ------------------------------------------------------------------------------
      sublflux(:) = 0.
      fluxsnow(:) = fluxsnow + fluxautosnow(:,k)*tdt/tdt_in
      
      ! Detect max/random overlap clouds that are separated by a clear layer
      where ( (stratcloud(:,k)>=1.e-10.and.stratcloud(:,k+1)<1.e-10) .or. nmr==0 )
        rdclfrsnow(:) = rdclfrsnow + mxclfrsnow - rdclfrsnow*mxclfrsnow
        mxclfrsnow(:) = 0.
      end where
      csfra(:) = max( rdclfrsnow + mxclfrsnow - rdclfrsnow*mxclfrsnow, 1.e-10 )

      ! Snow fall speed (from Lin et al 1983 - see GFDL AM3)
      do iq = 1,imax
        rs = max( fluxsnow(iq)/dz(iq,k), 0. )
        if ( csfra(iq)>=1.e-10 ) then
          vs2(iq) = max( 0.1, 1.82*(rs/csfra(iq))**0.0625 )
        end if
      end do

      ! Set up the parameters for the flux-divergence calculation
      do iq = 1,imax
        alph          = tdt*vs2(iq)/dz(iq,k)
        alph         = max( min( alph, 50. ), 0. )
        foutsnow(iq)  = 1. - exp(-alph)          !analytical
        fthrusnow(iq) = 1. - foutsnow(iq)/alph  !analytical
      end do

      alphaf(:) = hls*qsatg(:,k)/(rvap*ttg(:,k)**2)
      gam1(:)   = hlscp*alphaf(:) !(L/cp)*dqsdt (HBG notation)     

      if ( any( fluxsnow>0. ) ) then
          
        ! Melt falling snow if > 0 deg C due to rain accretion
        ! (based on Lin et al 83, but using 0.65 and 0.44 coeffs following the UM approach)
        do iq = 1,imax
          rs = max(fluxsnow(iq), 0.)/dz(iq,k)
          if ( ttg(iq,k)>tfrz .and. rs>1.e-10 ) then
            n0s            = 2.e6*exp(-0.12*max(ttg(iq,k)-tfrz,-200.))        
            slopes_s       = ( max(fluxsnow(iq),0.)/dz(iq,k)/(pi*rho_s*n0s))**0.25
            qvp            = rhov(iq,k)/rhoa(iq,k)  
            cdt            = tdt*2.*pi*n0s/hlf*(tcond*(ttg(iq,k)-tfrz)/rhoa(iq,k)-vdifu*hl*(qsatg(iq,k)-qvp))          &
                                     *(0.65*slopes_s**2+0.44*scm3*gam263*sqrt(clin/visk)*slopes_s**2.63*sqrt(denfac(iq)))
            drf            = max( min( rs, cdt ), 0. ) 
            iflux          = min( drf*dz(iq,k), fluxsnow(iq) )   ! flux of snow
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
            pleo_psmlt(iq,k)  = pleo_psmlt(iq,k) / dqf         ! sny 10: Melt falling snow if > 0 deg C due to rain accretion
          end if
        end do 
        
        ! Compute the sublimation of snow falling from level k+1 into level k
        ! (Currently treated the same as LDR97 ice sublimation - see UM and ACCESS 1.3)
        do iq = 1,imax
          qvp = rhov(iq,k)/rhoa(iq,k)
          if ( fluxsnow(iq)>0. .and. qvp<qsatg(iq,k) ) then ! sublime snow
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
            pleo_pssub(iq,k)  = pleo_pssub(iq,k) + dqs           ! sny 11: sublimation of snow falling from level k+1 into level k 
          end if
        end do   
        
        ! Accretion of cloud liquid by falling snow (from Lin et al 1983 - psacw)
        do iq = 1,imax
          rl = rhol(iq,k)
          if ( fluxsnow(iq)+sublflux(iq)>0. .and. rl>1.e-10 .and. ttg(iq,k)<tfrz ) then
            n0s      = 2.e6*exp(-0.12*max(ttg(iq,k)-tfrz,-200.))        
            slopes_s = ( max(fluxsnow(iq)+sublflux(iq),0.)/dz(iq,k)/(pi*rho_s*n0s))**0.25
            cdt          = tdt*denfac(iq)*pi*clin*gam325*n0s/4.*slopes_s**3.25
            drl          = max( min( csfra(iq)*rl, rl*cdt/(1.+0.5*cdt) ), 0. ) ! mass of liquid
            lflux        = drl*dz(iq,k)                                                 ! flux of liquid
            dql          = drl/rhoa(iq,k)                                               ! mixing ratio of liquid
            fluxsnow(iq) = fluxsnow(iq) + lflux
            rhol(iq,k)   = rhol(iq,k)   - drl
            qaccr(iq,k)  = qaccr(iq,k)  + dql
            dttg         = hlfcp*dql
            ttg(iq,k)    = ttg(iq,k) + dttg
            qsatg(iq,k)  = qsatg(iq,k) + gam1(iq)*dttg/hlscp
            cftmp        = clfr(iq,k)*drl/rl
            clfr(iq,k)   = clfr(iq,k) - cftmp
            mxclfrsnow(iq) = max( mxclfrsnow(iq), cftmp )
            pleo_psacw(iq,k)  = pleo_psacw(iq,k) + dql           ! sny 12: Accretion of cloud liquid by falling snow
          end if
        end do 
        
        ! Accretion of rain by falling snow to form snow (from Lin et al 1983 - psacr)
        do iq = 1,imax
          rn = rhor(iq,k)
          if ( fluxsnow(iq)+sublflux(iq)>0. .and. rn>1.e-10 .and. ttg(iq,k)<tfrz ) then
            n0s          = 2.e6*exp(-0.12*max(ttg(iq,k)-tfrz,-200.))        
            slopes_s     = ( max(fluxsnow(iq)+sublflux(iq),0.)/dz(iq,k)/(pi*rho_s*n0s))**0.25
            slopes_r     = (( max(rn*dz(iq,k),0.)/max(crfra(iq),1.e-10)/tdt)**0.22)/714.
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
            pleo_psacr(iq,k)  = pleo_psacr(iq,k) + dql            ! sny 13:  Accretion of rain by falling snow to form snow
          end if
        end do
        
        ! Accretion of rain by falling snow to form graupel (neglected in Lin83 but included in UM)   
    
        ! Accretion of cloud ice by falling snow (from HDC 2004 - psaci)
        do iq = 1,imax
          rf = rhoi(iq,k)
          if ( fluxsnow(iq)+sublflux(iq)>0. .and. rf>1.e-10 .and. ttg(iq,k)<tfrz ) then
            n0s          = 2.e6*exp(-0.12*max(ttg(iq,k)-tfrz,-200.))        
            slopes_s     = ( max(fluxsnow(iq)+sublflux(iq),0.)/dz(iq,k)/(pi*rho_s*n0s))**0.25
            esi          = exp(0.05*max(ttg(iq,k)-tfrz,-100.))       ! efficiency
            cdt          = tdt*denfac(iq)*27.737*n0s*esi*slopes_s**3.41
            drf          = max( min( csfra(iq)*rf, rf*cdt/(1.+0.5*cdt) ), 0. ) ! mass of ice
            iflux        = drf*dz(iq,k)                                                 ! flux of ice
            dqf          = drf/rhoa(iq,k)                                               ! mixing ratio of ice
            fluxsnow(iq) = fluxsnow(iq) + iflux
            rhoi(iq,k)   = rhoi(iq,k)   - drf
            qaccf(iq,k)  = qaccf(iq,k)  + dqf
            cftmp        = cifr(iq,k)*drf/rf
            cifr(iq,k)   = cifr(iq,k) - cftmp
            mxclfrsnow(iq) = max( mxclfrsnow(iq), cftmp )
            pleo_psaci(iq,k)  = pleo_psaci(iq,k) + dqf            ! sny 14: Accretion of cloud ice by falling snow
          end if
        end do
        
      end if  ! fluxsnow(iq)>0.
     
    end if ! ncloud==3 .or. ncloud==4 .or. ncloud==13

        
    ! Ice ---------------------------------------------------------------------------------
    sublflux(:) = 0.
   
    ! Set up the rate constant for ice sublimation
    ! MJT notes - curly and Csbsav depend on vi2(:,k+1), so vi2(:,k) can be updated below
    do iq = 1,imax
      slopes_i = 1.6e3*10**(-0.023*(ttg(iq,k)-tfrz))
      es = qsatg(iq,k)*pk(iq)/epsil
      Aprpr(iq) = (hls/(rKa*ttg(iq,k)))*(hls/(rvap*ttg(iq,k))-1.)
      Bprpr(iq) = rvap*ttg(iq,k)/((Dva/pk(iq))*es)
      if ( nevapls==-1 .or. (nevapls==-2.and.condx(iq)>0..and.k<=ktsav(iq)) ) then
        curly(iq) = 0.
      else
        curly(iq) = 0.65*slopes_i**2+0.493*slopes_i*sqrt(slopes_i*vi2(iq)*rhoa(iq,k)/um) !Factor in curly brackets
      end if
    end do
    ! Define the rate constant for sublimation of snow, omitting factor rhoi
    Csbsav(:) = 4.*curly(:)/(rhoa(:,k)*qsatg(:,k)*(Aprpr(:)+Bprpr(:))*pi*vi2*rho_s)
    
    ! Detect max/random overlap clouds that are separated by a clear layer
    where ( (stratcloud(:,k)>=1.e-10.and.stratcloud(:,k+1)<1.e-10) .or. nmr==0 )
      rdclfrice(:) = rdclfrice + mxclfrice - rdclfrice*mxclfrice
      mxclfrice(:) = 0.
    end where
    cifra(:) = max( rdclfrice + mxclfrice - rdclfrice*mxclfrice, 1.e-10 )
 
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

    ! Set up the parameters for the flux-divergence calculation
    do iq = 1,imax
      alph         = tdt*vi2(iq)/dz(iq,k)
      alph         = max( min( alph, 50. ), 0. )
      foutice(iq)  = 1. - exp(-alph)    !analytical
      fthruice(iq) = 1. - foutice(iq)/alph  !analytical  
    end do 

    alphaf(:) = hls*qsatg(:,k)/(rvap*ttg(:,k)**2)
    gam1(:)   = hlscp*alphaf !(L/cp)*dqsdt (HBG notation)   
      
    if ( any( fluxice>0. ) ) then
        
      ! Melt falling ice if > 0 deg C
      do iq = 1,imax
        if ( ttg(iq,k)>tfrz .and. fluxice(iq)>0. ) then
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
          pleo_pimlt(iq,k)  = pleo_pimlt(iq,k) + qif          ! sny 15: Melt falling ice if > 0 deg C to form cloud water
        end if
      end do
      
      ! Compute the sublimation of ice falling from level k+1 into level k
      do iq = 1,imax
        qvp = rhov(iq,k)/rhoa(iq,k)
        if ( fluxice(iq)>0. .and. qvp<qsatg(iq,k) ) then ! sublime ice
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
          pleo_pisub(iq,k)  = pleo_pisub(iq,k) + dqs                         ! sny 16 sublimation of ice
        end if
      end do 
      
      ! Accretion of cloud liquid by falling ice (neglected in Lin et al 1983, but
      ! included in UM and ACCESS 1.3 as piacw)
      ! This calculation uses the incoming fluxice without subtracting sublimation
      ! (since subl occurs only outside cloud), so add sublflux back to fluxice.
      do iq = 1,imax
        rl = rhol(iq,k)
        if ( fluxice(iq)+sublflux(iq)>0. .and. rl>1.e-10 ) then
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
          pleo_piacw(iq,k)  = pleo_piacw(iq,k) + dql                        ! sny 17: Accretion of cloud liquid by falling ice
        end if
      end do
      
      if ( ncloud==3 .or. ncloud==4 .or. ncloud==13 ) then
        ! Accretion of rain by falling ice to produce ice (from Lin et al 1983 - piacr)
        ! (see UM and ACCESS 1.3 piacr-c for an alternate formulation)
        do iq = 1,imax
          rn  = rhor(iq,k)
          if ( fluxice(iq)+sublflux(iq)>0. .and. rn>1.e-10 .and. ttg(iq,k)<tfrz ) then
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
            pleo_piacr(iq,k)  = pleo_piacr(iq,k) + dql                     ! sny 18: Accretion of rain by falling ice to produce ice
          end if
        end do
      end if 
      
      ! Accretion of rain by falling ice to produce graupel (Neglected in Lin et al 1983)
      ! (see UM and ACCESS 1.3 piacr-g for an alternate formulation)
      
    end if ! fluxice(iq)>0.
    
    ! store slope for aerosols
    do iq = 1,imax
      slopes_i      = 1.6e3*10**(-0.023*(ttg(iq,k)-tfrz))
      pslopes(iq,k) = pslopes(iq,k) + slopes_i*tdt/tdt_in  
    end do
    
    
    ! Rain --------------------------------------------------------------------------------
    evap(:) = 0.

    ! Add flux of melted snow to fluxrain
    fluxrain(:) = fluxrain(:) + fluxmelt(:) + fluxautorain(:,k)*tdt/tdt_in
    mxclfrrain(:) = max( mxclfrrain(:), cfmelt(:) )
    
    ! Detect maximum/random overlap clouds that are separated by a clear layer
    where ( (stratcloud(:,k)>=1.e-10.and.stratcloud(:,k+1)<1.e-10) .or. nmr==0 )
      rdclfrrain(:) = rdclfrrain + mxclfrrain - rdclfrrain*mxclfrrain
      mxclfrrain(:) = 0.
    end where
    crfra(:) = max( rdclfrrain + mxclfrrain - rdclfrrain*mxclfrrain, 1.e-10 )
      
    ! Calculate rain fall speed (MJT suggestion)
    if ( ncloud==2 .or. ncloud==3 .or. ncloud==4 .or. ncloud==12 .or. ncloud==13 ) then
      do iq = 1,imax
        Fr(iq)       = max( fluxrain(iq)/tdt/max(crfra(iq),1.e-10),0.)
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
    
    alphaf(:) = hls*qsatg(:,k)/(rvap*ttg(:,k)**2)
    gam1(:)   = hlscp*alphaf !(L/cp)*dqsdt (HBG notation)
    
    if ( any( fluxrain>0. ) ) then    
    
      if ( ncloud==3 .or. ncloud==4 .or. ncloud==13 ) then
      
        do iq = 1,imax
          rn = max(fluxrain(iq),0.)/dz(iq,k)
          if ( rn>1.e-10 .and. ttg(iq,k)<tfrz ) then
            slopes_r        = (( max(fluxrain(iq),0.)/max(crfra(iq),1.e-10)/tdt)**0.22)/714.
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
            pleo_psure(iq,k)  = pleo_psure(iq,k) + dql                ! sny 19: NOT SURE what process here
          end if
        end do

      end if ! ncloud==3 .or. ncloud==4 .or. ncloud==13
          
      ! Evaporation of rain
      qpf(:)     = fluxrain/rhodz !Mix ratio of rain which falls into layer
      clrevap(:) = (1.-clfr(:,k)-cifr(:,k))*qpf
      where ( ttg(:,k)<tfrz .and. ttg(:,k)>=tice )
        qsl(:)   = qsatg(:,k) + epsil*esdiffx(ttg(:,k),imax)/pk
      elsewhere
        qsl(:)   = qsatg(:,k)
      end where
      do iq = 1,imax
        qvp     = rhov(iq,k)/rhoa(iq,k)
        if ( fluxrain(iq)>0. .and. crfra(iq)>0. ) then
          es       = qsl(iq)*pk(iq)/epsil
          Apr      = (hl/(rKa*ttg(iq,k)))*(hl/(rvap*ttg(iq,k))-1.)
          Bpr      = rvap*ttg(iq,k)/((Dva/pk(iq))*es)
          Fr(iq)   = fluxrain(iq)/tdt/max(crfra(iq), 1.e-10)
          Cev      = crfra(iq)*3.8e2*sqrt(Fr(iq)/rhoa(iq,k))/(qsl(iq)*(Apr+Bpr))
          dqsdt    = hl*qsl(iq)/(rvap*ttg(iq,k)**2)
          bl       = 1. + 0.5*Cev*tdt*(1.+hlcp*dqsdt)
          evap(iq) = tdt*(Cev/bl)*(qsl(iq)-qvp)
          satevap  = (qsl(iq)-qvp)/(1.+hlcp*dqsdt)  !Evap to saturate
          ! vr2=11.3*Fr(iq)**(1./9.)/sqrt(rhoa(mg,k)) !Actual fall speed
          ! vr2=5./sqrt(rhoa(mg,k))               !Nominal fall speed
          evap(iq) = max( 0., min( evap(iq), satevap, clrevap(iq) ) )
          pleo_prevp(iq,k)  = pleo_prevp(iq,k) + evap(iq)                  ! sny 20: Evaporation of rain 
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
      end do  
      
      ! Now do the collection of liquid cloud by rain term (cf. pracc in Lin83).
      fcol(:) = 0.
      Fr(:) = 0.
      do iq = 1,imax
        rl = rhol(iq,k)
        if ( fluxrain(iq)>0. .and. rl>1.e-10 ) then
          Fr(iq)       = max(fluxrain(iq),0.)/tdt/max(crfra(iq),1.e-10)
          fcol(iq)     = crfra(iq)
          cdt          = tdt*Ecol*0.24*fcol(iq)*Fr(iq)**0.75
          coll         = max( min( rhol(iq,k), rhol(iq,k)*cdt/(1.+0.5*cdt) ), 0. ) ! mass
          lflux        = coll*dz(iq,k)                                            ! flux
          dql          = coll/rhoa(iq,k)                                          ! mixing ratio
          fluxrain(iq) = fluxrain(iq) + lflux
          rhol(iq,k)   = rhol(iq,k)   - coll
          qcoll(iq,k)  = qcoll(iq,k)  + dql
          cltmp        = clfr(iq,k)*coll/rl
          clfr(iq,k)   = clfr(iq,k) - cltmp
          mxclfrrain(iq) = max( mxclfrrain(iq), cltmp )
          pleo_pracc(iq,k)  = pleo_pracc(iq,k) + dql              ! sny 21: collection of liquid cloud by rain
        end if
      end do
      
      ! subtract evaporated rain
      do iq = 1,imax
        lflux        = evap(iq)*rhodz(iq)
        fluxrain(iq) = max( fluxrain(iq) - lflux, 0. ) !To avoid roundoff -ve's
      end do
      
      if ( ncloud==3 .or. ncloud==4 .or. ncloud==13 ) then
          
        ! Accretion of cloud snow by rain (from Lin et al 1983 - pracs)
        do iq = 1,imax
          rs = max( rhos(iq,k), 0. )
          if ( fluxrain(iq)>0. .and. rs>1.e-10 .and. ttg(iq,k)>tfrz+1. ) then
            n0s          = 2.e6*exp(-0.12*max(ttg(iq,k)-tfrz,-200.))        
            slopes_s     = ( max(rs,0.)/(pi*rho_s*n0s))**0.25
            slopes_r     = (( max(fluxrain(iq),0.)/max(crfra(iq),1.e-10)/tdt)**0.22)/714.  
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
            pleo_pracs(iq,k)  = pleo_pracs(iq,k) + dqf             ! sny 22: Accretion of cloud snow by rain
          end if
        end do

      end if ! ncloud==3 .or. ncloud==4 .or. ncloud==13        
        
      ! store for aerosols
      qevap(:,k) = qevap(:,k) + evap
      prscav(:,k) = prscav(:,k) + tdt*0.24*fcol*Fr**0.75   !Strat only
      
    end if ! fluxrain>0.

  
    ! Liquid ------------------------------------------------------------------------------
    ! (Currently cloud droplet settling is negected, although included in UM and ACCESS 1.3)


    ! Misc ------------------------------------------------------------------------------

    if ( any( fluxrain>0. ) ) then
    
      if ( ncloud==3 .or. ncloud==4 .or. ncloud==13 ) then  
        ! Accretion of cloud ice by rain to produce snow or grauple (from Lin et al 1983 - praci)
        ! (Neglected in UM and ACCESS 1.3)
        do iq = 1,imax
          rf = rhoi(iq,k)
          rn = fluxrain(iq)/dz(iq,k)
          if ( rn>qr0_crt ) then
            xwgt = 1.
          else
            xwgt = 0.  
          end  if
          if ( fluxrain(iq)>0. .and. rf>1.e-10 .and. ttg(iq,k)<tfrz ) then
            slopes_r        = (( max(fluxrain(iq),0.)/max(crfra(iq),1.e-10)/tdt)**0.22)/714.  
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
            pleo_praci(iq,k)  = pleo_praci(iq,k) + drf/rhoa(iq,k)       ! sny 23: Accretion of cloud ice by rain to produce snow or grauple
          end if
        end do
        
      end if  
      
    end if  
  
   
    ! Update fluxes and area fractions for graupel, snow, ice and rain

    rhototf(:)       = rhog(:,k) + rhos(:,k) + rhoi(:,k)
    xfrac_graupel(:) = rhog(:,k)/max(rhototf(:),1.e-10)
    xfrac_snow(:)    = rhos(:,k)/max(rhototf(:),1.e-10)
    xfrac_ice(:)     = max( 0., 1.-xfrac_graupel(:)-xfrac_snow(:) )
    
    ! Melting and freezing
    fluxm(:,k) = fluxm(:,k) + fluxmelt(:)
    fluxf(:,k) = fluxf(:,k) + fluxfreeze(:)
    
    ! Save velocity
    vg(:,k) = vg2(:)
    vs(:,k) = vs2(:)
    vi(:,k) = vi2(:)
    
    if ( ncloud==3 .or. ncloud==4 .or. ncloud==13 ) then
        
      ! Grauple
      ! calculate maximum and random overlap for falling graupel
      pqfsedice(:,k) = pqfsedice(:,k) + xfrac_graupel(:)*foutgraupel(:)*tdt/tdt_in ! Save sedimentation rate for aerosol scheme
      where ( fluxgraupel(:)<=0. )
        rdclfrgraupel(:) = 0.
        mxclfrgraupel(:) = 0.
      end where
      mxclfrgraupel(:) = max( mxclfrgraupel, cfgraupel(:,k) ) ! for rhogout
      cgfra(:) = max( 1.e-10, mxclfrgraupel+rdclfrgraupel-mxclfrgraupel*rdclfrgraupel ) ! rnd overlap
      ! Compute fluxes into the box
      cffluxin(:) = cgfra - cfgraupel(:,k)
      rhogin(:)   = fluxgraupel(:)/dz(:,k)
      ! Compute the fluxes of snow leaving the box
      cffluxout(:) = cfgraupel(:,k)*foutgraupel
      rhogout(:)   = rhog(:,k)*foutgraupel
      ! Update the rhos and cfsnow fields
      cfgraupel(:,k) = cfgraupel(:,k) - cffluxout(:) + cffluxin(:)*(1.-fthrugraupel)
      rhog(:,k)      = rhog(:,k) - rhogout + rhogin*(1.-fthrugraupel)
      fluxgraupel(:) = max( rhogout*dz(:,k) + fluxgraupel*fthrugraupel, 0. )
      where ( fluxgraupel<1.e-6 )
        rhog(:,k) = rhog(:,k) + fluxgraupel/dz(:,k)
        fluxgraupel = 0.
      end where  
      ! Now fluxgraupel is flux leaving layer k
      fluxg(:,k) = fluxg(:,k) + fluxgraupel
      
      ! Snow
      ! calculate maximum and random overlap for falling snow
      pqfsedice(:,k) = pqfsedice(:,k) + xfrac_snow(:)*foutsnow(:)*tdt/tdt_in ! Save sedimentation rate for aerosol scheme
      where ( fluxsnow(:)<=0. )
        rdclfrsnow(:) = 0.
        mxclfrsnow(:) = 0.
      end where
      mxclfrsnow(:) = max( mxclfrsnow, cfsnow(:,k) ) ! for rhosout
      csfra(:) = max( 1.e-10, mxclfrsnow+rdclfrsnow-mxclfrsnow*rdclfrsnow ) 
      ! Compute fluxes into the box
      cffluxin(:) = csfra - cfsnow(:,k)
      rhosin(:)   = fluxsnow(:)/dz(:,k)
      ! Compute the fluxes of snow leaving the box
      cffluxout(:) = cfsnow(:,k)*foutsnow
      rhosout(:)   = rhos(:,k)*foutsnow
      ! Update the rhos and cfsnow fields
      cfsnow(:,k) = cfsnow(:,k) - cffluxout + cffluxin*(1.-fthrusnow)
      rhos(:,k)   = rhos(:,k) - rhosout + rhosin*(1.-fthrusnow)
      fluxsnow(:) = max( rhosout*dz(:,k) + fluxsnow*fthrusnow, 0. )
      where ( fluxsnow<1.e-6 )
        rhos(:,k) = rhos(:,k) + fluxsnow/dz(:,k)
        fluxsnow = 0.
      end where  
      ! Now fluxsnow is flux leaving layer k
      fluxs(:,k) = fluxs(:,k) + fluxsnow
        
    end if ! ncloud==3 .or. ncloud==4 .or. ncloud==13

    
    ! Ice
    ! calculate maximum and random overlap for falling ice
    pqfsedice(:,k) = pqfsedice(:,k) + xfrac_ice(:)*foutice(:)*tdt/tdt_in ! Save sedimentation rate for aerosol scheme
    where ( fluxice(:)<=0. )
      rdclfrice(:) = 0.
      mxclfrice(:) = 0.
    end where
    mxclfrice(:) = max( mxclfrice, cifr(:,k) ) ! for rhoiout
    cifra(:) = max( 1.e-10, mxclfrice+rdclfrice-mxclfrice*rdclfrice ) !rnd overlap the mx and rd ice fractions
    ! Compute fluxes into the box
    cffluxin(:) = cifra - cifr(:,k)
    rhoiin(:)   = fluxice/dz(:,k)
    ! Compute the fluxes of ice leaving the box
    cffluxout(:) = cifr(:,k)*foutice
    rhoiout(:)   = rhoi(:,k)*foutice
    ! Update the rhoi and cifr fields
    cifr(:,k)  = min( 1.-clfr(:,k), cifr(:,k)-cffluxout(:)+cffluxin*(1.-fthruice) )
    rhoi(:,k)  = rhoi(:,k) - rhoiout + rhoiin*(1.-fthruice)
    fluxice(:) = max( rhoiout*dz(:,k) + fluxice*fthruice, 0. )
    where ( fluxice<1.e-6 )
      rhoi(:,k) = rhoi(:,k) + fluxice/dz(:,k)
      fluxice = 0.
    end where  
    ! Now fluxice is flux leaving layer k
    fluxi(:,k) = fluxi(:,k) + fluxice  
  
    ! Rain
    ! Calculate the raining cloud cover down to this level, for stratiform (crfra).
    where ( fluxrain(:)<=0. )
      rdclfrrain(:) = 0.
      mxclfrrain(:) = 0.
    end where
    mxclfrrain(:) = max( mxclfrrain, cfrain(:,k) ) ! for rhorout    
    crfra(:) = max( 1.e-10, rdclfrrain+mxclfrrain-rdclfrrain*mxclfrrain )
    ! Compute fluxes into the box
    cffluxin(:) = crfra - cfrain(:,k)
    rhorin(:)   = fluxrain(:)/dz(:,k)
    ! Compute the fluxes of rain leaving the box
    ! Use the flux-divergent form as in Rotstayn (QJRMS, 1997)
    cffluxout(:) = cfrain(:,k)*foutliq
    rhorout(:)   = rhor(:,k)*foutliq
    ! Update the rhor and cfrain fields
    cfrain(:,k) = cfrain(:,k) - cffluxout + cffluxin*(1.-fthruliq)
    rhor(:,k)   = rhor(:,k) - rhorout + rhorin*(1.-fthruliq)
    fluxrain(:) = max( rhorout*dz(:,k) + fluxrain*fthruliq, 0. )

    !SONNY fluxrain_save(:,k) = fluxrain(:)
    where ( fluxrain<1.e-6 )
      rhor(:,k) = rhor(:,k) + fluxrain/dz(:,k)
      fluxrain = 0.
    end where
    ! Now fluxrain is flux leaving layer k
    fluxr(:,k) = fluxr(:,k) + fluxrain
    
  end do ! k loop
  
end do   ! n

! process rate: devide to total (large)time step
pleo_pcaut = pleo_pcaut 
pleo_psaut = pleo_psaut 
pleo_pgaut = pleo_pgaut 
pleo_pgmlt = pleo_pgmlt / tdt_in
pleo_pgsub = pleo_pgsub / tdt_in
pleo_pgacw = pleo_pgacw / tdt_in
pleo_pgacr = pleo_pgacr / tdt_in
pleo_pgaci = pleo_pgaci / tdt_in
pleo_pgacs = pleo_pgacs / tdt_in
pleo_psmlt = pleo_psmlt / tdt_in
pleo_pssub = pleo_pssub / tdt_in
pleo_psacw = pleo_psacw / tdt_in
pleo_psacr = pleo_psacr / tdt_in
pleo_psaci = pleo_psaci / tdt_in
pleo_pimlt = pleo_pimlt / tdt_in
pleo_pisub = pleo_pisub / tdt_in
pleo_piacw = pleo_piacw / tdt_in
pleo_piacr = pleo_piacr / tdt_in
pleo_psure = pleo_psure / tdt_in
pleo_prevp = pleo_prevp / tdt_in
pleo_pracc = pleo_pracc / tdt_in
pleo_pracs = pleo_pracs / tdt_in
pleo_praci = pleo_praci / tdt_in

! store precip, snow and graupel
precs(:) = precs + fluxr(:,1) + fluxi(:,1) + fluxs(:,1) + fluxg(:,1)
preci(:) = preci + fluxi(:,1) + fluxs(:,1)
precg(:) = precg + fluxg(:,1)

do k = 1,kl
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
    
end module leoncld_mod
