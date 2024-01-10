! Conformal Cubic Atmospheric Model
    
! Copyright 2015-2023 Commonwealth Scientific Industrial Research Organisation (CSIRO)
    
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
public cloud_aerosol_mode, lin_aerosolmode, maxlintime

integer, save :: cloud_aerosol_mode = 0     ! 0=original, 1=standard feedback to aerosols
integer, save :: lin_aerosolmode    = 0     ! 0=off, 1=aerosol indirect effects for Lin microphysics
real, save :: maxlintime            = 120.  ! time-step for Lin microphysics

! ldr  = 0      Diagnosed cloud scheme (depreciated)
! ldr /= 0      Prognostic cloud condensate (different ice fall speed options)
  
! ncloud = 0    Standard LDR cloud microphysics with water vapour, liquid cloud and ice cloud
! ncloud = 2    Same as ncloud=0, but with prognostic rain and modified cfrac
! ncloud = 3    Same as ncloud=2, but with prognostic graupel and snow, as well as modified cfrac
! ncloud = 4    Use prognostic cloud fraction based on Tiedtke from GFDL-CM3
! ncloud = 10   Same as ncloud=0 with Tiedtke from GFDL-CM3
! ncloud = 12   Same as ncloud=2 with Tiedtke from GFDL-CM3
! ncloud = 13   Same as ncloud=3 with Tiedtke from GFDL-CM3 (i.e., same as ncloud=4)
! ncloud = 100  Use Lin et al 2nd moment microphysics
! ncloud = 110  Same as ncloud=100 with Tiedtke from GFDL-CM3

!                            Water vapour (qg)
!
!   Cloud water (qlg,cfrac)                      Cloud ice (qfg,cfrac)
!
!   Rain (qrg,rfrac)                             Snow (qsg,sfrac)         Graupel (qgrg,gfrac)

! qg, qlg, qfg, qrg, qsg and qgrg are mixing ratios (g/g) and cfrac, rfrac, sfrac, gfrac are area cover
! fractions

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
use kuocomb_m                     ! JLM convection
use latlong_m                     ! Lat/lon coordinates
use leoncld_mod                   ! Prognostic cloud condensate
use liqwpar_m                     ! Cloud water mixing ratios
use map_m                         ! Grid map arrays
use module_aux_cosp               ! COSP cloud sat simulator
use module_aux_rad                ! Additional cloud and radiation routines
use module_mp_sbu_ylin            ! Lin cloud microphysics
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

include 'kuocom.h'                ! Convection parameters
  
integer :: tile, js, je, k, n, iq
integer :: njumps, m, idjd_t
real, dimension(imax,kl) :: lcfrac, lgfrac
real, dimension(imax,kl) :: lqg, lqgrg, lqlg, lqfg, lqlrad, lqfrad, lqrg, lqsng, lrfrac, lsfrac, lt
real, dimension(imax,kl) :: ldpsldt, lnettend, lstratcloud, lclcon, lcdrop, lrhoa
real, dimension(imax,kl) :: lqccon
real, dimension(imax,kl) :: lrkmsave, lrkhsave
real, dimension(imax,kl) :: lfluxr, lfluxm, lfluxf, lfluxi, lfluxs, lfluxg
real, dimension(imax,kl) :: lqevap, lqsubl, lqauto, lqcoll, lqaccr
real, dimension(imax,kl) :: lvi
real, dimension(imax,kl) :: l_rliq, l_rice, l_cliq, l_cice, lp
real, dimension(imax,kl) :: l_rliq_in, l_rice_in, l_rsno_in
real, dimension(imax,kl) :: lppfevap, lppfmelt, lppfprec, lppfsnow, lppfsubl
real, dimension(imax,kl) :: lpplambs, lppmaccr, lppmrate, lppqfsedice, lpprfreeze, lpprscav
real, dimension(ifull,kl) :: clcon, cdrop
real, dimension(ifull,kl) :: fluxr, fluxm, fluxf, fluxi, fluxs, fluxg
real, dimension(ifull,kl) :: fevap, fsubl, fauto, fcoll, faccr, faccf
real, dimension(ifull,kl) :: vi
real, dimension(ifull,kl) :: dz, rhoa

real(kind=8), dimension(imax,0:kl) :: zlevv
real(kind=8), dimension(imax,kl) :: riz, tothz, thz, dzw
real(kind=8), dimension(imax,kl) :: znc, znr, zni, zns
real(kind=8), dimension(imax,kl) :: zpres, zrhoa, zcdrop, zvi
real(kind=8), dimension(imax,kl) :: zEFFC1D, zEFFI1D, zEFFS1D, zEFFR1D
real(kind=8), dimension(imax,kl) :: zqg, zqlg, zqfg, zqrg, zqsng
real(kind=8), dimension(imax,kl) :: zfluxr, zfluxi, zfluxs, zfluxm, zfluxf
real(kind=8), dimension(imax,kl) :: zqevap, zqsubl, zqauto, zqcoll, zqaccr
#ifdef debug
real(kind=8), dimension(imax,kl) :: zpsnow, zpsaut, zpsfw, zpsfi
real(kind=8), dimension(imax,kl) :: zpraci, zpiacr, zpsaci, zpsacw
real(kind=8), dimension(imax,kl) :: zpsdep, zpssub, zpracs, zpsacr
real(kind=8), dimension(imax,kl) :: zpsmlt, zpsmltevp, zprain
real(kind=8), dimension(imax,kl) :: zpraut, zpracw, zprevp, zpgfr
real(kind=8), dimension(imax,kl) :: zpvapor, zpclw, zpladj, zpcli
real(kind=8), dimension(imax,kl) :: zpimlt, zpihom, zpidw, zpiadj
real(kind=8), dimension(imax,kl) :: zqschg
#endif
!real, dimension(imax,kl) :: precrz, preciz, precsz
real(kind=8), dimension(imax) :: pptrain, pptsnow, pptice
real(kind=8) tdt
real prf_temp, prf, fcol, fr, alph
logical mydiag_t


!----------------------------------------------------------------------------
! Prepare inputs for cloud microphysics


!$omp do schedule(static) private(js,je,k,lrhoa,lcdrop)
do tile = 1,ntiles
  js = (tile-1)*imax + 1
  je = tile*imax
  
  do k = 1,kl
    lrhoa(:,k) = ps(js:je)*sig(k)/(rdry*t(js:je,k))
    rhoa(js:je,k) = lrhoa(:,k)
    dz(js:je,k) = -rdry*dsig(k)*t(js:je,k)/(grav*sig(k)) 
  end do
  
  ! Calculate droplet concentration from aerosols (for non-convective faction of grid-box)
  ! xtg is unchanged since updating GPU
  call aerodrop(js,lcdrop,lrhoa,outconv=.true.)
  cdrop(js:je,:) = lcdrop(:,:)
end do
!$omp end do nowait


!$omp do schedule(static) private(js,je,lclcon)
do tile = 1,ntiles
  js = (tile-1)*imax + 1
  je = tile*imax
  ! Calculate convective cloud fraction
  call convectivecloudfrac(lclcon,kbsav(js:je),ktsav(js:je),condc(js:je),acon,bcon)
  clcon(js:je,:) = lclcon(:,:)
end do
!$omp end do nowait


!----------------------------------------------------------------------------
! Update cloud fraction

!$omp do schedule(static) private(js,je,idjd_t,mydiag_t,lcfrac),               &
!$omp private(lqccon,lqfg,lqfrad,lqg,lqlg,lqlrad,lt),                          &
!$omp private(ldpsldt,lnettend,lstratcloud,lclcon,lcdrop,lrkmsave,lrkhsave)
do tile = 1,ntiles
  js = (tile-1)*imax + 1
  je = tile*imax

  idjd_t = mod(idjd-1,imax) + 1
  mydiag_t = ((idjd-1)/imax==tile-1).AND.mydiag

  lqg      = qg(js:je,:)
  lqlg     = qlg(js:je,:)
  lqfg     = qfg(js:je,:)
  lt       = t(js:je,:)
  ldpsldt  = dpsldt(js:je,:)
  lclcon   = clcon(js:je,:)
  lcdrop   = cdrop(js:je,:)
  lstratcloud = stratcloud(js:je,:)
  if ( ncloud==4 .or. (ncloud>=10.and.ncloud<=13) .or. ncloud==110 ) then
    lnettend = nettend(js:je,:)
    lrkmsave = rkmsave(js:je,:)
    lrkhsave = rkhsave(js:je,:)
  end if

  call update_cloud_fraction(lcfrac,land(js:je),                                       &
              ps(js:je),lqccon,lqfg,lqfrad,lqg,lqlg,lqlrad,lt,                         &
              ldpsldt,lnettend,lstratcloud,lclcon,lcdrop,em(js:je),pblh(js:je),idjd_t, &
              mydiag_t,ncloud,nclddia,ldr,rcrit_l,rcrit_s,rcm,cld_decay,               &
              vdeposition_mode,tiedtke_form,lrkmsave,lrkhsave)

  cfrac(js:je,:) = lcfrac
  qccon(js:je,:) = lqccon
  qg(js:je,:)    = lqg
  qlg(js:je,:)   = lqlg
  qfg(js:je,:)   = lqfg
  qlrad(js:je,:) = lqlrad
  qfrad(js:je,:) = lqfrad
  t(js:je,:)     = lt
  stratcloud(js:je,:) = lstratcloud
  if ( ncloud==4 .OR. (ncloud>=10.AND.ncloud<=13) .or. ncloud==110 ) then
    nettend(js:je,:) = lnettend
  end if
end do
!$omp end do nowait


!----------------------------------------------------------------------------
! Update cloud condensate
select case ( interp_ncloud(ldr,ncloud) )
  case("LEON")
  
    !$omp do schedule(static) private(js,je,idjd_t,mydiag_t),                     &
    !$omp private(lgfrac,lrfrac,lsfrac),                                          &
    !$omp private(lppfevap,lppfmelt,lppfprec,lppfsnow,lppfsubl),                  &
    !$omp private(lpplambs,lppmaccr,lppmrate,lppqfsedice,lpprfreeze,lpprscav),    &
    !$omp private(lqfg,lqg,lqgrg,lqlg,lqrg,lqsng,lt),                             &
    !$omp private(lstratcloud,lcdrop),                                            &
    !$omp private(lfluxr,lfluxm,lfluxf,lfluxi,lfluxs,lfluxg),                     &
    !$omp private(lqevap,lqsubl,lqauto,lqcoll,lqaccr,lvi)
    do tile = 1,ntiles
      js = (tile-1)*imax + 1
      je = tile*imax

      idjd_t = mod(idjd-1,imax) + 1
      mydiag_t = ((idjd-1)/imax==tile-1).AND.mydiag

      lgfrac   = gfrac(js:je,:)
      lrfrac   = rfrac(js:je,:)
      lsfrac   = sfrac(js:je,:)
      lqg      = qg(js:je,:)
      lqgrg    = qgrg(js:je,:)
      lqlg     = qlg(js:je,:)
      lqfg     = qfg(js:je,:)
      lqrg     = qrg(js:je,:)
      lqsng    = qsng(js:je,:)
      lt       = t(js:je,:)
      lcdrop   = cdrop(js:je,:)
      lstratcloud = stratcloud(js:je,:)

      call leoncld_work(condg(js:je),conds(js:je),condx(js:je),lgfrac,ktsav(js:je),           &
              lppfevap,lppfmelt,lppfprec,lppfsnow,lppfsubl,                                   &
              lpplambs,lppmaccr,lppmrate,lppqfsedice,lpprfreeze,lpprscav,                     &
              precip(js:je),ps(js:je),lqfg,lqg,lqgrg,lqlg,lqrg,lqsng,lrfrac,lsfrac,lt,        &
              lstratcloud,lcdrop,lfluxr,lfluxm,lfluxf,lfluxi,lfluxs,lfluxg,lqevap,lqsubl,     &
              lqauto,lqcoll,lqaccr,lvi,                                                       &
              idjd_t,mydiag_t,ncloud,nevapls,ldr,rcm,imax,kl)

      gfrac(js:je,:) = lgfrac
      rfrac(js:je,:) = lrfrac
      sfrac(js:je,:) = lsfrac
      qg(js:je,:)    = lqg
      qlg(js:je,:)   = lqlg
      qfg(js:je,:)   = lqfg
      qrg(js:je,:)   = lqrg
      qsng(js:je,:)  = lqsng
      qgrg(js:je,:)  = lqgrg
      t(js:je,:)     = lt
      stratcloud(js:je,:) = lstratcloud
      fluxr(js:je,:) = lfluxr/dt
      fluxm(js:je,:) = lfluxm/dt
      fluxf(js:je,:) = lfluxf/dt
      fluxi(js:je,:) = lfluxi/dt
      fluxs(js:je,:) = lfluxs/dt
      fluxg(js:je,:) = lfluxg/dt
      fevap(js:je,:) = lqevap(:,:)*rhoa(js:je,:)*dz(js:je,:)/dt
      fsubl(js:je,:) = lqsubl(:,:)*rhoa(js:je,:)*dz(js:je,:)/dt
      fauto(js:je,:) = lqauto(:,:)*rhoa(js:je,:)*dz(js:je,:)/dt
      fcoll(js:je,:) = lqcoll(:,:)*rhoa(js:je,:)*dz(js:je,:)/dt
      faccr(js:je,:) = lqaccr(:,:)*rhoa(js:je,:)*dz(js:je,:)/dt
      vi(js:je,:) = lvi
      ! backwards compatible data for aerosols
      if ( abs(iaero)>=2 ) then
        ppfevap(js:je,:)    = lppfevap
        ppfmelt(js:je,:)    = lppfmelt
        ppfprec(js:je,:)    = lppfprec
        ppfsnow(js:je,:)    = lppfsnow
        ppfsubl(js:je,:)    = lppfsubl
        pplambs(js:je,:)    = lpplambs
        ppmaccr(js:je,:)    = lppmaccr
        ppmrate(js:je,:)    = lppmrate
        ppqfsedice(js:je,:) = lppqfsedice
        pprfreeze(js:je,:)  = lpprfreeze
        pprscav(js:je,:)    = lpprscav
      end if
    end do
    !$omp end do nowait

    !$omp do schedule(static) private(js,je,iq,k,lcfrac,lqlg,lqfg,lt,lcdrop,lp) &
    !$omp private(l_rliq,l_rice,l_cliq,l_cice)
    do tile = 1,ntiles
      js = (tile-1)*imax + 1
      je = tile*imax

      lcfrac = cfrac(js:je,:)
      lqlg = qlrad(js:je,:)
      lqfg = qfrad(js:je,:)
      lt = t(js:je,:)
      lcdrop = cdrop(js:je,:)
      do k = 1,kl
        lp(:,k) = ps(js:je)*sig(k)
      end do
      call cloud3(l_rliq,l_rice,l_cliq,l_cice,lcfrac,lqlg,lqfg,lp,lt,lcdrop,imax,kl)
      stras_rliq(js:je,:) = l_rliq   ! save effective radius for cosp
      stras_rice(js:je,:) = l_rice
      stras_rsno(js:je,:) = 0.
      stras_rrai(js:je,:) = 0.
      stras_cliq(js:je,:) = l_cliq
      stras_cice(js:je,:) = l_cice
    end do
    !$omp end do nowait

  case("LIN")
      
    !$omp do schedule(static) private(js,je,riz,zlevv,m,zqg,zqlg,zqrg,zqfg)   &
    !$omp private(zqsng,k,iq,prf_temp,prf,tothz,thz,zrhoa,zpres,dzw,znc)      &
    !$omp private(zcdrop,znr,zni,zns,pptrain,pptsnow,pptice,njumps,tdt,n)     &
    !$omp private(zeffc1d,zeffi1d,zeffs1d,zeffr1d,zfluxr,zfluxi,zfluxs)       &
    !$omp private(zfluxm,zfluxf,zqevap,zqsubl,zqauto,zqcoll,zqaccr)           &
#ifdef debug
    !$omp private(zpsnow,zpsaut)                                              &
    !$omp private(zpsfw,zpsfi,zpraci,zpiacr,zpsaci,zpsacw,zpsdep,zpssub)      &
    !$omp private(zpracs,zpsacr,zpsmlt,zpsmltevp,zprain,zpraut,zpracw)        &
    !$omp private(zprevp,zpgfr,zpvapor,zpclw,zpladj,zpcli,zpimlt,zpihom)      &
    !$omp private(zpidw,zpiadj,zqschg)                                        &
#endif
    !$omp private(zvi)
    do tile = 1, ntiles
      js = (tile-1)*imax + 1 ! js:je inside 1:ifull
      je = tile*imax         ! len(js:je) = imax

      riz(1:imax,:) = 0._8 ! partition between snow and graupel
        
      ! pack data from ifull into imax
      zlevv(1:imax,0)   = real( zs(js:je)/grav, 8 )
      zlevv(1:imax,1)   = zlevv(1:imax,0) + real( bet(1)*t(js:je,1)/grav, 8 )
      do m = 2,kl
        zlevv(1:imax,m) = zlevv(1:imax,m-1) + real( (bet(m)*t(js:je,m)+betm(m)*t(js:je,m-1))/grav, 8 )
      end do
        
      zqg(1:imax,:) = real( qg(js:je,:), 8 )
      zqlg(1:imax,:) = real( qlg(js:je,:), 8 )
      zqrg(1:imax,:) = real( qrg(js:je,:), 8 )
      zqfg(1:imax,:) = real( qfg(js:je,:), 8 )
      zqsng(1:imax,:) = real( qsng(js:je,:) + qgrg(js:je,:), 8 )
      ! ----------------
      do k = 1,kl
        do iq = 1,imax
          prf_temp      = ps(iq+js-1)*sig(k)
          prf           = 0.01*prf_temp                        ! ps is SI units
          tothz(iq,k)   = real( (prf/1000.)**(rdry/cp), 8 )
          thz(iq,k)     = real( t(iq+js-1,k)/tothz(iq,k), 8 )
          zrhoa(iq,k)   = real( rhoa(iq+js-1,k), 8 )
          !zorhoa(iq,k)  = 1._8/zrhoa(iq,k)
          zpres(iq,k)   = real( prf_temp, 8 )
          !sqrhoz(iq,k)  = 1._8
          dzw(iq,k)     = real( dz(iq+js-1,k), 8 )
        end do
      end do

      znc(1:imax,:) = 0._8
      zcdrop(1:imax,:) = real( cdrop(js:je,:), 8 ) ! aerosol
      znr(1:imax,:)  = real( nr(js:je,:), 8 )
      zni(1:imax,:)  = real( ni(js:je,:), 8 )
      zns(1:imax,:)  = real( ns(js:je,:), 8 )

      pptrain(1:imax) = 0._8
      pptsnow(1:imax) = 0._8
      pptice(1:imax)  = 0._8
     

      ! Use sub time-step if required
      njumps = int(dt/(maxlintime+0.01)) + 1
      tdt    = real( dt/real(njumps), 8 )
      do n = 1,njumps
        call clphy1d_ylin(tdt, imax,                       &
                       zqg, zqlg, zqrg, zqfg, zqsng,       &
                       thz, tothz, zrhoa,                  &
                       zpres, zlevv, dzw,                  &
                       !precrz, preciz, precsz,            & !zdc 20220116
                       zEFFC1D, zEFFI1D, zEFFS1D, zEFFR1D, & !zdc 20220208
                       pptrain, pptsnow, pptice,           &
                       1, kl, riz,                         &
                       znc, znr, zni, zns,                 &
                       zfluxr,zfluxi,zfluxs,zfluxm,        &
                       zfluxf,zqevap,zqsubl,zqauto,zqcoll, &
                       zqaccr,zvi,                         & !aerosol scheme
#ifdef diagnostic
                       zpsnow,zpsaut,zpsfw,zpsfi,zpraci,   & !process rate cloud microphysics
                       zpiacr,zpsaci,zpsacw,zpsdep,        &
                       zpssub,zpracs,zpsacr,zpsmlt,        &
                       zpsmltevp,zprain,zpraut,zpracw,     &
                       zprevp,zpgfr,zpvapor,zpclw,         &
                       zpladj,zpcli,zpimlt,zpihom,         &
                       zpidw,zpiadj,zqschg,                &
#endif
                       zcdrop, lin_aerosolmode)              !aerosol feedback
      end do


      t(js:je,:) = real( thz(1:imax,:)*tothz(1:imax,:) )

      !unpack data from imax to ifull.

      qg(js:je,:)         = real( zqg(1:imax,:) )                        ! qv mixing ratio
      qlg(js:je,:)        = real( zqlg(1:imax,:) )                       ! ql mixing ratio
      qfg(js:je,:)        = real( zqfg(1:imax,:) )                       ! qf mixing ratio (ice)
      qrg(js:je,:)        = real( zqrg(1:imax,:) )                       ! qr mixing ratio (rain)
      qsng(js:je,:)       = real( zqsng(1:imax,:)*(1._8-riz(1:imax,:)) ) ! qs mixing ratio (snow)
      qgrg(js:je,:)       = real( zqsng(1:imax,:)*riz(1:imax,:) )        ! qg mixing ration (graupel)

      where ( qrg(js:je,:)>0. )
         rfrac(js:je,:) = 1.
      elsewhere
         rfrac(js:je,:) = 0.
      end where
      where ( qsng(js:je,:)>0. )
         sfrac(js:je,:) = 1.
      elsewhere
         sfrac(js:je,:) = 0.
      end where
      where ( qgrg(js:je,:)>0. )
         gfrac(js:je,:) = 1.
      elsewhere
         gfrac(js:je,:) = 0.
      end where
      where ( qlg(js:je,:)+qfg(js:je,:)>1.e-12 )
        stratcloud(js:je,:) = max( stratcloud(js:je,:), 1.e-8 )
      end where 
      
      !nc(js:je,:)        = real( znc(1:imax,:) )
      nr(js:je,:)         = real( znr(1:imax,:) )
      ni(js:je,:)         = real( zni(1:imax,:) )
      ns(js:je,:)         = real( zns(1:imax,:) )
      stras_rliq(js:je,:) = real( zEFFC1D(1:imax,:) )   ! save effective radius for cosp
      stras_rice(js:je,:) = real( zEFFI1D(1:imax,:) )
      stras_rsno(js:je,:) = real( zEFFS1D(1:imax,:) )
      stras_rrai(js:je,:) = real( zEFFR1D(1:imax,:) )

      fluxr(js:je,1)      = real( pptrain(1:imax) )
      fluxr(js:je,2:kl)   = real( zfluxr(1:imax,1:kl-1) )  ! flux for aerosol calculation
      fluxi(js:je,1)      = real( pptice(1:imax) )
      fluxi(js:je,2:kl)   = real( zfluxi(1:imax,1:kl-1) )
      fluxs(js:je,1)      = real( pptsnow(1:imax)*(1._8-riz(1:imax,1)) )
      fluxs(js:je,2:kl)   = real( zfluxs(1:imax,1:kl-1)*(1._8-riz(1:imax,1:kl-1)) )
      fluxg(js:je,1)      = real( pptsnow(1:imax)*riz(1:imax,1) )
      fluxg(js:je,2:kl)   = real( zfluxs(1:imax,1:kl-1)*riz(1:imax,1:kl-1) )
      fluxm(js:je,1:kl-1) = real( zfluxm(1:imax,1:kl-1) )
      fluxm(js:je,kl)     = 0.
      fluxf(js:je,1:kl-1) = real( zfluxf(1:imax,1:kl-1) )
      fluxf(js:je,kl)     = 0.
      fevap(js:je,:) = real( zqevap(1:imax,:) )
      fsubl(js:je,:) = real( zqsubl(1:imax,:) )
      fauto(js:je,:) = real( zqauto(1:imax,:) )
      fcoll(js:je,:) = real( zqcoll(1:imax,:) )
      faccr(js:je,:) = real( zqaccr(1:imax,:) )
      vi(js:je,:) = real( zvi(1:imax,:) )

      !if (process_rate_mode == 2) then
      !  psnow(js:je,:)   = real( zpsnow(1:imax,:) ) !process rate to understand cloud microphysics
      !  psaut(js:je,:)   = real( zpsaut(1:imax,:) )
      !  psfw(js:je,:)    = real( zpsfw(1:imax,:) )
      !  psfi(js:je,:)    = real( zpsfi(1:imax,:) )
      !  praci(js:je,:)   = real( zpraci(1:imax,:) )
      !  piacr(js:je,:)   = real( zpiacr(1:imax,:) )
      !  psaci(js:je,:)   = real( zpsaci(1:imax,:) )
      !  psacw(js:je,:)   = real( zpsacw(1:imax,:) )
      !  psdep(js:je,:)   = real( zpsdep(1:imax,:) )
      !  pssub(js:je,:)   = real( zpssub(1:imax,:) )
      !  pracs(js:je,:)   = real( zpracs(1:imax,:) )
      !  psacr(js:je,:)   = real( zpsacr(1:imax,:) )
      !  psmlt(js:je,:)   = real( zpsmlt(1:imax,:) )
      !  psmltevp(js:je,:)= real( zpsmltevp(1:imax,:) )
      !  prain(js:je,:)   = real( zprain(1:imax,:) )
      !  praut(js:je,:)   = real( zpraut(1:imax,:) )
      !  pracw(js:je,:)   = real( zpracw(1:imax,:) )
      !  prevp(js:je,:)   = real( zprevp(1:imax,:) )
      !  pgfr(js:je,:)    = real( zpgfr(1:imax,:) )
      !  pvapor(js:je,:)  = real( zpvapor(1:imax,:) )
      !  pclw(js:je,:)    = real( zpclw(1:imax,:) )
      !  pladj(js:je,:)   = real( zpladj(1:imax,:) )
      !  pcli(js:je,:)    = real( zpcli(1:imax,:) )
      !  pimlt(js:je,:)   = real( zpimlt(1:imax,:) )
      !  pihom(js:je,:)   = real( zpihom(1:imax,:) )
      !  pidw(js:je,:)    = real( zpidw(1:imax,:) )
      !  piadj(js:je,:)   = real( zpiadj(1:imax,:) )
      !  qschg(js:je,:)   = real( zqschg(1:imax,:) )
      !end if
          
      condx(js:je)  = condx(js:je) + real( pptrain(1:imax) + pptsnow(1:imax) + pptice(1:imax) )
      conds(js:je)  = conds(js:je) + real( pptsnow(1:imax)*(1._8-riz(1:imax,1)) + pptice(1:imax) )
      condg(js:je)  = condg(js:je) + real( pptsnow(1:imax)*riz(1:imax,1) ) ! for graupel
      precip(js:je) = precip(js:je) + real( pptrain(1:imax) + pptsnow(1:imax) + pptice(1:imax) )

    end do     !tile loop
    !$omp end do nowait

    !$omp do schedule(static) private(js,je,iq,k,lcfrac,lqlg,lqfg,lt,lcdrop,lp) &
    !$omp private(l_rliq,l_rice,l_cliq,l_cice,l_rliq_in,l_rice_in,l_rsno_in)
    do tile = 1,ntiles
      js = (tile-1)*imax + 1
      je = tile*imax

      lcfrac = cfrac(js:je,:)
      lqlg = qlrad(js:je,:)
      lqfg = qfrad(js:je,:)
      lt = t(js:je,:)
      lcdrop = cdrop(js:je,:)
      do k = 1,kl
        lp(:,k) = ps(js:je)*sig(k)
      end do
      l_rliq_in = stras_rliq(js:je,:)
      l_rice_in = stras_rice(js:je,:)
      l_rsno_in = stras_rsno(js:je,:)
      call cloud3(l_rliq,l_rice,l_cliq,l_cice,lcfrac,lqlg,lqfg,lp,lt,lcdrop,imax,kl, &
                  stras_rliq=l_rliq_in,stras_rice=l_rice_in,stras_rsno=l_rsno_in)
      stras_rliq(js:je,:) = l_rliq   ! save effective radius for cosp
      stras_rice(js:je,:) = l_rice
      stras_cliq(js:je,:) = l_cliq
      stras_cice(js:je,:) = l_cice
    end do
    !$omp end do nowait
      
  case default
    write(6,*) "ERROR: unknown mp_physics option "
    call ccmpi_abort(-1)
      
end select


! Aerosol feedbacks
if ( abs(iaero)>=2 ) then
  if ( interp_ncloud(ldr,ncloud)/="LEON".or.cloud_aerosol_mode>0  ) then
    !$omp do schedule(static) private(js,je,iq,k,fcol,fr,alph)
    do tile = 1,ntiles
      js = (tile-1)*imax + 1
      je = tile*imax
    
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

      ppfprec(js:je,1)   = 0. !At TOA
      ppfmelt(js:je,1)   = 0. !At TOA
      ppfsnow(js:je,1)   = 0. !At TOA
      pprfreeze(js:je,1) = 0. !At TOA
      do k = 1,kl-1
        do iq = js,je
          ! rainfall flux (entering from above) (kg/m2/s)  
          ppfprec(iq,kl+1-k) = fluxr(iq,k+1)+fluxm(iq,k)-fluxf(iq,k)
          ! snowfall flux (entering from above) (kg/m2/s)
          ppfsnow(iq,kl+1-k) = fluxi(iq,k+1)+fluxs(iq,k+1)+fluxg(iq,k+1) &
                                 -fluxm(iq,k)+fluxf(iq,k)
          ! snowfall flux melting in layer k (kg/m2/s)
          ppfmelt(iq,kl+1-k) = fluxm(iq,k)
          ! rainfall flux freezing in layer k (kg/m2/s)
          pprfreeze(iq,kl+1-k) = fluxf(iq,k) 
        end do
      end do
      do k = 1,kl
        do iq = js,je
          ! rainall flux evaporating in layer k  
          ppfevap(iq,kl+1-k) = fevap(iq,k)
          ! snowfall flux evaporating in layer k
          ppfsubl(iq,kl+1-k) = fsubl(iq,k)
          ! precipitation formation rate (kg/kg/s)
          ppmrate(iq,kl+1-k) = (fauto(iq,k)+fcoll(iq,k))/(rhoa(iq,k)*dz(iq,k))
          ! liquid accertion rate (kg/kg/s)
          ppmaccr(iq,kl+1-k) = faccr(iq,k)/(rhoa(iq,k)*dz(iq,k))
          ! slope (lambda) for snow crystal size djstribution (m**-1)
          pplambs(iq,kl+1-k) = 1.6e3*10**(-0.023*max(t(iq,k)-tfrz,-200.))          
          ! Fraction rain scavenging rate in time-step  
          fcol = rfrac(iq,k)
          Fr = fluxr(iq,k)/max(rfrac(iq,k),1.e-10)
          pprscav(iq,kl+1-k) = dt*0.24*fcol*Fr**0.75
          ! Fractional ice sedimentation in time-step
          alph = dt*vi(iq,k)/dz(iq,k)
          ppqfsedice(iq,kl+1-k) = 1. - exp(-alph)
        end do  
      end do
    end do ! tile
    !$omp end do nowait
  end if   ! interp_ncloud(ldr,ncloud)/="LEON".or.cloud_aerosol_mode>0
end if     ! abs(iaero)>=2

! update COSP cloud sat simulator if avaliable
! (must be downloaded and compiled independently of CCAM)
call cloud_simulator

return
end subroutine ctrl_microphysics

!====================================================================================================
! SUBROUTINE interp_ncloud
!   
! subroutine to select the cloud microphysics scheme for CCAM
!====================================================================================================
pure function interp_ncloud(ldr, ncloud) result(mp_physics)

implicit none

integer, intent(in) :: ldr
integer, intent(in) :: ncloud
character(len=10) :: mp_physics

mp_physics = "ERROR"

if ( ldr /= 0 ) then
  select case(ncloud)
    case(0,2,3,4,10,12,13)
      mp_physics = "LEON"
    case(100,110)
      mp_physics = "LIN"
  end select
end if

return
end function interp_ncloud  
  
end module module_ctrl_microphysics