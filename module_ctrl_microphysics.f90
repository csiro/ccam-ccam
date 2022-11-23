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
    
module module_ctrl_microphysics

implicit none

private
public ctrl_microphysics
public cloud_aerosol_mode

integer, save :: cloud_aerosol_mode = 0     ! 0=original, 1=standard feedback to aerosols

contains
    
!====================================================================================================
! SUBROUTINE ctrl_microphysics
! subroutine to call cloud microphysics
! The current available option are LEO and LIN microphysics
!====================================================================================================
subroutine ctrl_microphysics

use aerointerface                 ! Aerosol interface
use arrays_m                      ! Atmosphere dyamics prognostic arrays
use cc_mpi                        ! CC MPI routines
use cc_omp                        ! CC OpenMP routines
use cfrac_m                       ! Cloud fraction
use cloudmod                      ! Prognostic cloud fraction
use const_phys                    ! Physical constants
use estab                         ! Liquid saturation function
use filnames_m                    ! Filenames
use kuocom_m                      ! Convection parameters
use kuocomb_m                     ! JLM convection
use latlong_m                     ! Lat/lon coordinates
use leoncld_mod                   ! Prognostic cloud condensate
use liqwpar_m                     ! Cloud water mixing ratios
use map_m                         ! Grid map arrays
use module_aux_rad                ! Additional cloud and radiation routines
use morepbl_m                     ! Additional boundary layer diagnostics
use newmpar_m                     ! Grid parameters
use nharrs_m                      ! Non-hydrostatic atmosphere arrays
use parm_m, only : dt,idjd,     &
      iaero,irest,ktau,nwt        ! Model configuration
use pbl_m                         ! Boundary layer arrays
use prec_m                        ! Precipitation
use raddiag_m                     ! Radiation diagnostic
use screen_m                      ! Screen level diagnostics
use sflux_m                       ! Surface flux routines
use sigs_m                        ! Atmosphere sigma levels
use soil_m                        ! Soil and surface data
use soilsnow_m                    ! Soil, snow and surface data
use work3f_m                      ! Grid work arrays
use vvel_m                        ! Additional vertical velocity

implicit none

integer :: tile, is, ie, k, ij, n, iq
integer :: idjd_t
real, dimension(imax,kl) :: lcfrac, lgfrac
real, dimension(imax,kl) :: lqg, lqgrg, lqlg, lqfg, lqlrad, lqfrad, lqrg, lqsng, lrfrac, lsfrac, lt
real, dimension(imax,kl) :: ldpsldt, lnettend, lstratcloud, lclcon, lcdrop, lrhoa
real, dimension(imax,kl) :: lqccon
real, dimension(imax,kl) :: lrkmsave, lrkhsave
real, dimension(imax,kl) :: lfluxr, lfluxm, lfluxf, lfluxi, lfluxs, lfluxg
real, dimension(imax,kl) :: lqevap, lqsubl, lqauto, lqcoll, lqaccr, lqaccf
real, dimension(imax,kl) :: lvi, lvs, lvg
real, dimension(imax,kl) :: lppfevap, lppfmelt, lppfprec, lppfsnow, lppfsubl
real, dimension(imax,kl) :: lpplambs, lppmaccr, lppmrate, lppqfsedice, lpprfreeze, lpprscav
real, dimension(ifull,kl) :: clcon, cdrop
real, dimension(ifull,kl) :: fluxr, fluxm, fluxf, fluxi, fluxs, fluxg
real, dimension(ifull,kl) :: fevap, fsubl, fauto, fcoll, faccr, faccf
real, dimension(ifull,kl) :: vi, vs, vg
real, dimension(ifull,kl) :: dz, rhoa
real fcol, fr, qtot, xic, xsn, xgr, vave, alph
logical :: mydiag_t

!----------------------------------------------------------------------------
! Prepare inputs for /loud microphysics

!$acc enter data create(dz,rhoa,cdrop,clcon)
!$acc enter data create(vi,vs,vg)
!$acc enter data create(fluxr,fluxm,fluxf,fluxi,fluxs,fluxg,fevap,fsubl,fauto,fcoll,faccr)
!$acc enter data create(ppfevap,ppfmelt,ppfprec,ppfsnow)
!$acc enter data create(ppfsubl,pplambs,ppmaccr,ppmrate,ppqfsedice,pprfreeze,pprscav)

!$omp do schedule(static) private(lrhoa,lcdrop,lclcon)
!$acc parallel present(sig,dsig,dz,rhoa,cdrop,clcon) &
!$acc present(kbsav,ktsav,ps,t,condc)
!$acc loop gang private(lrhoa,lcdrop,lclcon)
do tile = 1,ntiles
  is = (tile-1)*imax + 1
  ie = tile*imax
  
  do k = 1,kl
    lrhoa(:,k) = ps(is:ie)*sig(k)/(rdry*t(is:ie,k))
    rhoa(is:ie,k) = lrhoa(:,k)
    dz(is:ie,k) = -rdry*dsig(k)*t(is:ie,k)/(grav*sig(k)) 
  end do
  
  ! Calculate droplet concentration from aerosols (for non-convective faction of grid-box)  
  call aerodrop(is,lcdrop,lrhoa,outconv=.TRUE.)
  cdrop(is:ie,:) = lcdrop(:,:)

  ! Calculate convective cloud fraction
  call convectivecloudfrac(lclcon,kbsav(is:ie),ktsav(is:ie),condc(is:ie))
  clcon(is:ie,:) = lclcon(:,:)
  
end do
!$acc end parallel
!$omp end do nowait

!----------------------------------------------------------------------------
! Update cloud fraction

!$acc enter data create(qlrad,qfrad,stratcloud,nettend)
!$acc enter data create(rkmsave,rkhsave)
!$acc enter data create(qccon)
!$acc update device(qlrad,qfrad,stratcloud,nettend)
!$acc update device(rkmsave,rkhsave)

!$omp do schedule(static) private(is,ie),                                      &
!$omp private(lcfrac),                                                         &
!$omp private(lqccon,lqfg,lqfrad,lqg,lqlg,lqlrad,lt),                          &
!$omp private(ldpsldt,lnettend,lstratcloud,lclcon,lcdrop,lrkmsave,lrkhsave),   &
!$omp private(idjd_t,mydiag_t)
!$acc parallel present(qlrad,qfrad,stratcloud,nettend) &
!$acc present(rkmsave,rkhsave) present(em,dpsldt) &
!$acc present(qccon) present(land,kbsav,ktsav,ps,pblh,qg,qlg,qfg,t,cdrop,clcon,cfrac)
!$acc loop gang private(lcfrac,lqg,lqlg,lqfg,lqlrad,lqfrad,lt,ldpsldt,lclcon) &
!$acc private(lcdrop,lstratcloud,lnettend,lrkmsave,lrkhsave,lqccon)
do tile = 1,ntiles
  is = (tile-1)*imax + 1
  ie = tile*imax

  idjd_t = mod(idjd-1,imax) + 1
  mydiag_t = ((idjd-1)/imax==tile-1).AND.mydiag

  lcfrac   = cfrac(is:ie,:)
  lqg      = qg(is:ie,:)
  lqlg     = qlg(is:ie,:)
  lqfg     = qfg(is:ie,:)
  lqlrad   = qlrad(is:ie,:)
  lqfrad   = qfrad(is:ie,:)
  lt       = t(is:ie,:)
  ldpsldt  = dpsldt(is:ie,:)
  lclcon   = clcon(is:ie,:)
  lcdrop   = cdrop(is:ie,:)
  lstratcloud = stratcloud(is:ie,:)
  if ( ncloud==4 .or. (ncloud>=10.and.ncloud<=13) .or. ncloud==110 ) then
    lnettend = nettend(is:ie,:)
    lrkmsave = rkmsave(is:ie,:)
    lrkhsave = rkhsave(is:ie,:)
  end if

  call update_cloud_fraction(lcfrac,kbsav(is:ie),ktsav(is:ie),land(is:ie),             &
              ps(is:ie),lqccon,lqfg,lqfrad,lqg,lqlg,lqlrad,lt,                         &
              ldpsldt,lnettend,lstratcloud,lclcon,lcdrop,em(is:ie),pblh(is:ie),idjd_t, &
              mydiag_t,ncloud,nclddia,ldr,rcrit_l,rcrit_s,rcm,cld_decay,               &
              vdeposition_mode,tiedtke_form,lrkmsave,lrkhsave,imax,kl)

  cfrac(is:ie,:) = lcfrac
  qccon(is:ie,:) = lqccon
  qg(is:ie,:)    = lqg
  qlg(is:ie,:)   = lqlg
  qfg(is:ie,:)   = lqfg
  qlrad(is:ie,:) = lqlrad
  qfrad(is:ie,:) = lqfrad
  t(is:ie,:)     = lt
  stratcloud(is:ie,:) = lstratcloud
  if ( ncloud==4 .OR. (ncloud>=10.AND.ncloud<=13) .or. ncloud==110 ) then
    nettend(is:ie,:) = lnettend
  end if
end do
!$acc end parallel
!$omp end do nowait

!$acc exit data copyout(qlrad,qfrad,nettend)
!$acc exit data delete(rkmsave,rkhsave)
!$acc exit data copyout(qccon)

!----------------------------------------------------------------------------
! Update cloud condensate
select case ( interp_ncloud(ldr,ncloud) )
  case("LEON")

!$acc enter data create(gfrac,sfrac,qgrg)
!$acc enter data create(qrg,qsng,rfrac)
!$acc update device(gfrac,sfrac,qgrg)
!$acc update device(qrg,qsng,rfrac)
  
    !$omp do schedule(static) private(is,ie),                                     &
    !$omp private(lgfrac,lrfrac,lsfrac),                                          &
    !$omp private(lppfevap,lppfmelt,lppfprec,lppfsnowlppfsubl),                   &
    !$omp private(lpplambs,lppmaccr,lppmrate,lppqfsedice,lpprfreeze,lpprscav),    &
    !$omp private(lqfg,lqg,lqgrg,lqlg,lqrg,lqsng,lt),                             &
    !$omp private(lstratcloud,lclcon,lcdrop),                                     &
    !$omp private(idjd_t,mydiag_t)
    !$acc parallel present(gfrac,sfrac,qgrg)             &
    !$acc present(qrg,qsng,stratcloud,rfrac) &
    !$acc present(ppfevap,ppfmelt,ppfprec,ppfsnow)                                   &
    !$acc present(ppfsubl,pplambs,ppmaccr,ppmrate,ppqfsedice,pprfreeze,pprscav)       &
    !$acc present(fluxr,fluxm,fluxf,fluxi,fluxs,fluxg,fevap,fsubl,fauto,fcoll,faccr) &
    !$acc present(vi,vs,vg) present(dz,rhoa) &
    !$acc present(ktsav,ps,qg,qlg,qfg,t,precip,condx,conds,condg,cdrop)
    !$acc loop gang private(lgfrac,lrfrac,lsfrac,lqg,lqgrg,lqlg,lqfg,lqrg,lqsng)   &
    !$acc private(lt,lcdrop,lstratcloud,lfluxr,lfluxm,lfluxf,lfluxi,lfluxs,lfluxg) &
    !$acc private(lqevap,lqsubl,lqauto,lqcoll,lqaccr,lvi,lvs,lvg)                  &
    !$acc private(lppfevap,lppfmelt,lppfprec,lppfsnow,lppfsubl,lpplambs,lppmaccr)  &
    !$acc private(lppmrate,lppqfsedice,lpprfreeze,lpprscav)
    do tile = 1,ntiles
      is = (tile-1)*imax + 1
      ie = tile*imax

      idjd_t = mod(idjd-1,imax) + 1
      mydiag_t = ((idjd-1)/imax==tile-1).AND.mydiag

      lgfrac   = gfrac(is:ie,:)
      lrfrac   = rfrac(is:ie,:)
      lsfrac   = sfrac(is:ie,:)
      lqg      = qg(is:ie,:)
      lqgrg    = qgrg(is:ie,:)
      lqlg     = qlg(is:ie,:)
      lqfg     = qfg(is:ie,:)
      lqrg     = qrg(is:ie,:)
      lqsng    = qsng(is:ie,:)
      lt       = t(is:ie,:)
      lcdrop   = cdrop(is:ie,:)
      lstratcloud = stratcloud(is:ie,:)

      call leoncld_work(condg(is:ie),conds(is:ie),condx(is:ie),lgfrac,ktsav(is:ie),           &
              lppfevap,lppfmelt,lppfprec,lppfsnow,lppfsubl,                                   &
              lpplambs,lppmaccr,lppmrate,lppqfsedice,lpprfreeze,lpprscav,precip(is:ie),       &
              ps(is:ie),lqfg,lqg,lqgrg,lqlg,lqrg,lqsng,lrfrac,lsfrac,lt,                      &
              lstratcloud,lcdrop,lfluxr,lfluxm,lfluxf,lfluxi,lfluxs,lfluxg,lqevap,lqsubl,     &
              lqauto,lqcoll,lqaccr,lvi,lvs,lvg,                                               &
              idjd_t,mydiag_t,ncloud,nevapls,ldr,rcm,imax,kl)

      gfrac(is:ie,:) = lgfrac
      rfrac(is:ie,:) = lrfrac
      sfrac(is:ie,:) = lsfrac
      qg(is:ie,:)    = lqg
      qlg(is:ie,:)   = lqlg
      qfg(is:ie,:)   = lqfg
      qrg(is:ie,:)   = lqrg
      qsng(is:ie,:)  = lqsng
      qgrg(is:ie,:)  = lqgrg
      t(is:ie,:)     = lt
      stratcloud(is:ie,:) = lstratcloud
      fluxr(is:ie,:) = lfluxr/dt
      fluxm(is:ie,:) = lfluxm/dt
      fluxf(is:ie,:) = lfluxf/dt
      fluxi(is:ie,:) = lfluxi/dt
      fluxs(is:ie,:) = lfluxs/dt
      fluxg(is:ie,:) = lfluxg/dt
      fevap(is:ie,:) = lqevap(:,:)*rhoa(is:ie,:)*dz(is:ie,:)/dt
      fsubl(is:ie,:) = lqsubl(:,:)*rhoa(is:ie,:)*dz(is:ie,:)/dt
      fauto(is:ie,:) = lqauto(:,:)*rhoa(is:ie,:)*dz(is:ie,:)/dt
      fcoll(is:ie,:) = lqcoll(:,:)*rhoa(is:ie,:)*dz(is:ie,:)/dt
      faccr(is:ie,:) = lqaccr(:,:)*rhoa(is:ie,:)*dz(is:ie,:)/dt
      vi(is:ie,:) = lvi
      vs(is:ie,:) = lvs
      vg(is:ie,:) = lvg   
      ! backwards compatible data for aerosols
      if ( abs(iaero)>=2 ) then
#ifdef GPU
        call set_for_pgi_bug(ppfevap,ppfmelt,ppfprec,ppfsnow,ppfsubl,pplambs, &
                             ppmaccr,ppmrate,ppqfsedice,pprfreeze,pprscav,    &
                             lppfevap,lppfmelt,lppfprec,lppfsnow,lppfsubl,lpplambs, &
                             lppmaccr,lppmrate,lppqfsedice,lpprfreeze,lpprscav, &
                             is,ie,ifull,imax,kl)
#else
        ppfevap(is:ie,:)    = lppfevap
        ppfmelt(is:ie,:)    = lppfmelt
        ppfprec(is:ie,:)    = lppfprec
        ppfsnow(is:ie,:)    = lppfsnow
        ppfsubl(is:ie,:)    = lppfsubl
        pplambs(is:ie,:)    = lpplambs
        ppmaccr(is:ie,:)    = lppmaccr
        ppmrate(is:ie,:)    = lppmrate
        ppqfsedice(is:ie,:) = lppqfsedice
        pprfreeze(is:ie,:)  = lpprfreeze
        pprscav(is:ie,:)    = lpprscav
#endif
      end if

    end do
    !$acc end parallel
    !$omp end do nowait

!$acc exit data copyout(gfrac,sfrac,qgrg)
!$acc exit data copyout(qrg,qsng,rfrac)

  case("LIN")
    if ( myid==0 ) then
      write(6,*) "LIN microphysics ",ncloud
      write(6,*) "ERROR: Not supported"
    end if
    call ccmpi_abort(-1)
      
  case default
    write(6,*) "ERROR: unknown mp_physics option "
    call ccmpi_abort(-1)
      
end select
  
!$acc update self(t,qfg)
!$acc exit data copyout(dz,rhoa,cdrop,clcon)
!$acc exit data copyout(stratcloud)

! Aerosol feedbacks
if ( abs(iaero)>=2 .and. (interp_ncloud(ldr,ncloud)/="LEON".or.cloud_aerosol_mode>0)  ) then
  !$omp do schedule(static) private(is,ie,iq,k,fcol,fr,qtot,xic,xsn,xgr,vave,alph)
  do tile = 1,ntiles
    is = (tile-1)*imax + 1
    ie = tile*imax
    
    ! fluxi - ice flux leaving layer k to k-1 (kg/m2/s)
    ! fluxs - snow flux leaving layer k to k-1 (kg/m2/s)
    ! fluxg - graupel flux leving layer k to k-1 (kg/m2/s)
    ! fluxr - rain flux leaving layer k to k-1 (kg/m2/s)
    ! fluxm - ice melting flux in layer k (kg/m2/s)
    ! fluxf - liquid freezing flux in layer k (kg/m2/s)
    
    ! fevap - evaporation of rainfall flux (kg/m2/s)
    ! fsubl - sublimation of snow, ice and graupel flux (kg/m2/s)
    ! fauto - autoconversion flux for rainfall (kg/m2/s)
    ! fcoll - collection of cloud liquid water by rain (kg/m2/s)
    ! faccr - accretion of cloud liquid water by snow, ice and graupel (kg/m2/s)
    
    ! vi - icefall velocity (m/s)
    ! vs - snowfall velocity (m/s)
    ! vg - graupelfall velocity (m/s)

    ppfprec(is:ie,1)   = 0. !At TOA
    ppfmelt(is:ie,1)   = 0. !At TOA
    ppfsnow(is:ie,1)   = 0. !At TOA
    pprfreeze(is:ie,1) = 0. !At TOA
    do k = 1,kl-1
      ! rainfall flux (entering from above) (kg/m2/s)  
      ppfprec(is:ie,kl+1-k) = fluxr(is:ie,k+1)+fluxm(is:ie,k)-fluxf(is:ie,k)
      ! snowfall flux (entering from above) (kg/m2/s)
      ppfsnow(is:ie,kl+1-k) = fluxi(is:ie,k+1)+fluxs(is:ie,k+1)+fluxg(is:ie,k+1) &
                              -fluxm(is:ie,k)+fluxf(is:ie,k)
      ! snowfall flux melting in layer k (kg/m2/s)
      ppfmelt(is:ie,kl+1-k) = fluxm(is:ie,k)
      ! rainfall flux freezing in layer k (kg/m2/s)
      pprfreeze(is:ie,kl+1-k) = fluxf(is:ie,k) 
    end do
    do k = 1,kl
      ! rainall flux evaporating in layer k  
      ppfevap(is:ie,kl+1-k) = fevap(is:ie,k)
      ! snowfall flux evaporating in layer k
      ppfsubl(is:ie,kl+1-k) = fsubl(is:ie,k)
      ! precipitation formation rate (kg/kg/s)
      ppmrate(is:ie,kl+1-k) = (fauto(is:ie,k)+fcoll(is:ie,k))/(rhoa(is:ie,k)*dz(is:ie,k))
      ! liquid accertion rate (kg/kg/s)
      ppmaccr(is:ie,kl+1-k) = faccr(is:ie,k)/(rhoa(is:ie,k)*dz(is:ie,k))
      ! slope (lambda) for snow crystal size distribution (m**-1)
      pplambs(is:ie,kl+1-k) = 1.6e3*10**(-0.023*(t(iq,k)-tfrz))          
      ! Fraction rain scavenging rate in time-step
      do iq = is,ie
        fcol = rfrac(iq,k)
        Fr = fluxr(iq,k)/max(rfrac(iq,k),1.e-10)
        pprscav(iq,kl+1-k) = dt*0.24*fcol*Fr**0.75
      end do  
      ! Fractional ice sedimentation in time-step
      do iq = is,ie
        qtot = max( qfg(iq,k) + qsng(iq,k) + qgrg(iq,k), 1.e-10 )
        xsn = qsng(iq,k)/qtot
        xgr = qgrg(iq,k)/qtot
        xic = 1. - xsn - xgr
        vave = xic*vi(iq,k) + xsn*vs(iq,k) + xgr*vg(iq,k)
        alph = dt*vave/dz(iq,k)
        ppqfsedice(iq,kl+1-k) = 1. - exp(-alph)
      end do  
    end do
  end do ! tile
  !$omp end do nowait  
end if   ! if abs(iaero)>=2 .and. (interp_ncloud(ldr,ncloud)/="LEON".or.cloud_aerosol_mode>0)
  
!$acc exit data delete(vi,vs,vg)
!$acc exit data delete(fluxr,fluxm,fluxf,fluxi,fluxs,fluxg,fevap,fsubl,fauto,fcoll,faccr)
!$acc exit data copyout(ppfevap,ppfmelt,ppfprec,ppfsnow)
!$acc exit data copyout(ppfsubl,pplambs,ppmaccr,ppmrate,ppqfsedice,pprfreeze,pprscav)
  
  !! Estimate cloud droplet size
  !call cloud3(lrdrop,lrice,lconl,lconi,lcfrac,lqlrad,lqfrad,lpress,lt,lcdrop,imax,kl)

  ! cloud optical depth and emissivity ----------------------------
  ! Bands based on Slingo      
  !            BAND               SLINGO                 EBERT AND CURRY
  !
  !             1               0.25-0.69                0.25 - 0.7
  !             2               0.69-1.19                0.7 - 1.3
  !             3               1.19-2.38                1.3 - 2.5
  !             4               2.38-4.00                2.5 - 3.5    
  !do k = 1,kl
  !  rdrop = real(lrdrop(:,k))
  !  rice = real(lrice(:,k))
  !  lwp(:) = -dsig(k)*lqlrad(:,k)*ps(js:je)/grav
  !  iwp(:) = -dsig(k)*lqfrad(:,k)*ps(js:je)/grav
  !  tau_liq(:,1) = lwp(:)*1000.*(0.02817 + (1.305/rdrop))
  !  tau_liq(:,2) = lwp(:)*1000.*(0.02682 + (1.346/rdrop))
  !  tau_liq(:,3) = lwp(:)*1000.*(0.02264 + (1.454/rdrop))
  !  tau_liq(:,4) = lwp(:)*1000.*(0.01281 + (1.641/rdrop))
  !  tau_ice(:) = iwp(:)*1000.*(0.003448 + (2.431/rice))
  !  cloud_tau(js:je,k) = sum(tau_liq,dim=2)/4. + tau_ice ! 4 bands
  !  kliq(:) = 140.
  !  kice(:) = 4.83591 + 1758.511/rice
  !  cloud_emiss(js:je,k) = 1. - exp(-min(kliq(:)*lwp(:) + kice(:)*iwp(:),20.))
  !end do 
  
return
end subroutine ctrl_microphysics

!====================================================================================================
! SUBROUTINE interp_ncloud
!   
! subroutine to select the cloud microphysics scheme for CCAM
!====================================================================================================
function interp_ncloud(ldr, ncloud) result(mp_physics)

implicit none

integer, intent(in) :: ldr
integer, intent(in) :: ncloud
character(len=10) :: mp_physics

mp_physics = "ERROR"

if ( ldr /= 0 ) then
  select case(ncloud)
    case(0,2,3,4,10,12,13,20,21,22)
      mp_physics = "LEON"
    case(100,110,120)
      mp_physics = "LIN"
    end select
end if

return
end function interp_ncloud  
  
end module module_ctrl_microphysics
