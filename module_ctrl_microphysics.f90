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
public cloud_aerosol_mode, process_rate_mode
public lin_aerosolmode
public maxlintime

integer, save :: cloud_aerosol_mode = 0     ! 0=original, 1=standard feedback to aerosols
integer, save :: process_rate_mode  = 0     ! process rate for cloud microphysics (0=off)
integer, save :: lin_aerosolmode    = 0     ! aerosol in lin microphysics (0=off)
real, save :: maxlintime = 120.             ! Maximum time-step for Lin 2nd microphysics (sec)

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
use module_aux_cosp
use module_aux_rad                ! Additional cloud and radiation routines
use module_mp_sbu_ylin            ! Lin 2022 Scheme
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
real, dimension(imax,kl) :: ppleo_pcaut,ppleo_psaut,ppleo_pgaut,ppleo_pgmlt,ppleo_pgsub,ppleo_pgacw,&
                            ppleo_pgacr,ppleo_pgaci,ppleo_pgacs,ppleo_psmlt,ppleo_pssub,ppleo_psacw,&
                            ppleo_psacr,ppleo_psaci,ppleo_pimlt,ppleo_pisub,ppleo_piacw,ppleo_piacr,&
                            ppleo_psure,ppleo_prevp,ppleo_pracc,ppleo_pracs,ppleo_praci
real, dimension(ifull,kl) :: clcon, cdrop
real, dimension(ifull,kl) :: fluxm, fluxf
real, dimension(ifull,kl) :: fevap, fsubl, fauto, fcoll, faccr, faccf
real, dimension(ifull,kl) :: vi, vs, vg
real, dimension(ifull,kl) :: dz, rhoa
real, dimension(imax,kl) :: r_cfrac, r_qlrad, r_qfrad
real, dimension(imax,kl) :: r_cdrop
real, dimension(imax,kl) :: ptemp,ttemp
real(kind=8), dimension(imax,kl) :: Rdrop, Rice
real(kind=8), dimension(imax,kl) :: conl, coni
real fcol, fr, qtot, xic, xsn, xgr, vave, alph
logical :: mydiag_t

!====================================================================================================
!variable declaration for LIN cloud microphysics
!====================================================================================================
!intent(in)
integer                     :: ids,ide,jds,jde,kds,kde,ims,ime,jms,jme,kms,kme, &
                               its,ite,jts,jte,kts,kte
real                        :: dt_in
real                        :: ccn0 !1.0E8
real(kind=8), dimension(imax) :: ht
real, dimension(imax,kl)    :: w
real, dimension(imax,kl)         :: zlevv
real, dimension(imax,kl)    :: rho, pii, z, p_lin, dz8w
!intent in out
real, dimension(imax)       :: RAINNC,RAINNCV
real, dimension(imax)       :: SNOWNC,SNOWNCV
real, dimension(imax,kl)    :: th,qv_lin,qi,ql,qs,qr
real, dimension(imax,kl)    :: nnr,nni,nns
real, dimension(imax,kl)    :: Ri3D
real, dimension(imax,kl)    :: nn !nc,nr,ni,ns,nn
! intent(out)
real, dimension(imax,kl)    :: precr,preci,precs,eradc,eradi,erads,eradr
! local variables
integer                     :: min_q, max_q
real, dimension(imax)       :: rain, snow,ice
!real, dimension(kl)         :: qgz
real(kind=8), dimension(imax,kl) :: qvz,qlz,qrz,qiz,qsz,thz,tothz,rhoz,orhoz,sqrhoz
real(kind=8), dimension(imax,kl) :: zfluxr,zfluxi,zfluxs,zfluxg,zfluxm,    &
                               zfluxf,zfevap,zfsubl,zfauto,zfcoll,    &
                               zfaccr,zvi,zvs,zvg                         !for aerosol scheme
real(kind=8), dimension(imax,kl) :: zpsnow,zpsaut,zpsfw,zpsfi,zpraci,      &   !process rate to understand cloud microphysics
                               zpiacr,zpsaci,zpsacw,zpsdep,           &
                               zpssub,zpracs,zpsacr,zpsmlt,           &
                               zpsmltevp,zprain,zpraut,zpracw,        &
                               zprevp,zpgfr,zpvapor,zpclw,            &
                               zpladj,zpcli,zpimlt,zpihom,            &
                               zpidw,zpiadj,zpmidep
real, dimension(imax,kl)   ::  zzpsnow,zzpsaut,zzpsfw,zzpsfi,zzpraci, &   !process rate to understand cloud microphysics
                               zzpiacr,zzpsaci,zzpsacw,zzpsdep,       &
                               zzpssub,zzpracs,zzpsacr,zzpsmlt,       &
                               zzpsmltevp,zzprain,zzpraut,zzpracw,    &
                               zzprevp,zzpgfr,zzpvapor,zzpclw,        &
                               zzpladj,zzpcli,zzpimlt,zzpihom,        &
                               zzpidw,zzpiadj,zzpmidep

real(kind=8), dimension(imax,kl) :: zz,dzw,precrz,preciz,precsz
real(kind=8), dimension(imax,kl) :: prez
real(kind=8), dimension(imax,kl) :: EFFC1D,EFFI1D,EFFS1D,EFFR1D
real, dimension(imax,kl) :: r_effc1d, r_effi1d, r_effs1d
real(kind=8), dimension(imax,kl) :: riz
real                        :: rhoe_s
real(kind=8), dimension(imax)    :: pptice, pptrain, pptsnow
real, dimension(kl)         :: nnz
real(kind=8), dimension(imax,kl) :: ncz,nrz,niz,nsz,zcdrop
real(kind=8) :: tdt
integer :: njumps
integer i, j, m, cnt_sny, kr
real, dimension(imax)       :: prf
real, dimension(imax)       :: prf_temp
real, dimension(ifull_g,kl) :: data_g
real, dimension(ifull_g)    :: dsurf_g
real, dimension(ifull, kl)  :: data_kl1,data_kl2,data_kl3,data_kl4,    &
                               data_kl5,data_kl6,data_kl7,data_kl8,    &
                               data_kl9,data_kl10,data_kl11,data_kl12, &
                               data_kl13,data_kl14,data_kl15,data_kl16,&
                               data_kl17
!====================================================================================================

!----------------------------------------------------------------------------
! Prepare inputs for cloud microphysics

!$omp do schedule(static) private(is,ie,k,lrhoa,lcdrop,lclcon)
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
!$omp end do nowait


!----------------------------------------------------------------------------
! Update cloud fraction

!$omp do schedule(static) private(is,ie),                                      &
!$omp private(lcfrac),                                                         &
!$omp private(lqccon,lqfg,lqfrad,lqg,lqlg,lqlrad,lt),                          &
!$omp private(lqrg,lqsng,lqgrg),                                               &
!$omp private(ldpsldt,lnettend,lstratcloud,lclcon,lcdrop,lrkmsave,lrkhsave),   &
!$omp private(idjd_t,mydiag_t)
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
  lqrg     = qrg(is:ie,:)
  lqsng    = qsng(is:ie,:) 
  lqgrg    = qgrg(is:ie,:)
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
              lqrg,lqsng,lqgrg,                                                        &
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
  qrg(is:ie,:)   = lqrg
  qsng(is:ie,:)  = lqsng
  qgrg(is:ie,:)  = lqgrg 
  stratcloud(is:ie,:) = lstratcloud
  if ( ncloud==4 .OR. (ncloud>=10.AND.ncloud<=13) .or. ncloud==110 ) then
    nettend(is:ie,:) = lnettend
  end if
end do
!$omp end do nowait


!----------------------------------------------------------------------------
! Update cloud condensate
select case ( interp_ncloud(ldr,ncloud) )
  case("LEON")
  
    !$omp do schedule(static) private(is,ie),                                     &
    !$omp private(lgfrac,lrfrac,lsfrac),                                          &
    !$omp private(lppfevap,lppfmelt,lppfprec,lppfsnow,lppfsubl),                  &
    !$omp private(lpplambs,lppmaccr,lppmrate,lppqfsedice,lpprfreeze,lpprscav),    &
    !$omp private(lqfg,lqg,lqgrg,lqlg,lqrg,lqsng,lt),                             &
    !$omp private(lstratcloud,lclcon,lcdrop),                                     &
    !$omp private(idjd_t,mydiag_t)
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

      ! re-initialized the process rate vars for each time step
      ppleo_pcaut = 0.
      ppleo_psaut = 0.
      ppleo_pgaut = 0.
      ppleo_pgmlt = 0.
      ppleo_pgsub = 0.
      ppleo_pgacw = 0.
      ppleo_pgacr = 0.
      ppleo_pgaci = 0.
      ppleo_pgacs = 0.
      ppleo_psmlt = 0.
      ppleo_pssub = 0.
      ppleo_psacw = 0.
      ppleo_psacr = 0.
      ppleo_psaci = 0.
      ppleo_pimlt = 0.
      ppleo_pisub = 0.
      ppleo_piacw = 0.
      ppleo_piacr = 0.
      ppleo_psure = 0.
      ppleo_prevp = 0.
      ppleo_pracc = 0.
      ppleo_pracs = 0.
      ppleo_praci = 0.

      call leoncld_work(condg(is:ie),conds(is:ie),condx(is:ie),lgfrac,ktsav(is:ie),           &
              lppfevap,lppfmelt,lppfprec,lppfsnow,lppfsubl,                                   &
              lpplambs,lppmaccr,lppmrate,lppqfsedice,lpprfreeze,lpprscav,precip(is:ie),       &
              ps(is:ie),lqfg,lqg,lqgrg,lqlg,lqrg,lqsng,lrfrac,lsfrac,lt,                      &
              lstratcloud,lcdrop,lfluxr,lfluxm,lfluxf,lfluxi,lfluxs,lfluxg,lqevap,lqsubl,     &
              lqauto,lqcoll,lqaccr,lvi,lvs,lvg,                                               &
              idjd_t,mydiag_t,ncloud,nevapls,ldr,rcm,imax,kl,                                 &
              ppleo_pcaut,ppleo_psaut,ppleo_pgaut,ppleo_pgmlt,ppleo_pgsub,ppleo_pgacw,        &
              ppleo_pgacr,ppleo_pgaci,ppleo_pgacs,ppleo_psmlt,ppleo_pssub,ppleo_psacw,        &
              ppleo_psacr,ppleo_psaci,ppleo_pimlt,ppleo_pisub,ppleo_piacw,ppleo_piacr,        &
              ppleo_psure,ppleo_prevp,ppleo_pracc,ppleo_pracs,ppleo_praci)

      ! output LEO process rate
      if (process_rate_mode == 1) then
        leo_pcaut(is:ie,:) = ppleo_pcaut
        leo_psaut(is:ie,:) = ppleo_psaut
        leo_pgaut(is:ie,:) = ppleo_pgaut
        leo_pgmlt(is:ie,:) = ppleo_pgmlt
        leo_pgsub(is:ie,:) = ppleo_pgsub
        leo_pgacw(is:ie,:) = ppleo_pgacw
        leo_pgacr(is:ie,:) = ppleo_pgacr
        leo_pgaci(is:ie,:) = ppleo_pgaci
        leo_pgacs(is:ie,:) = ppleo_pgacs
        leo_psmlt(is:ie,:) = ppleo_psmlt
        leo_pssub(is:ie,:) = ppleo_pssub
        leo_psacw(is:ie,:) = ppleo_psacw
        leo_psacr(is:ie,:) = ppleo_psacr
        leo_psaci(is:ie,:) = ppleo_psaci
        leo_pimlt(is:ie,:) = ppleo_pimlt
        leo_pisub(is:ie,:) = ppleo_pisub
        leo_piacw(is:ie,:) = ppleo_piacw
        leo_piacr(is:ie,:) = ppleo_piacr
        leo_psure(is:ie,:) = ppleo_psure
        leo_prevp(is:ie,:) = ppleo_prevp
        leo_pracc(is:ie,:) = ppleo_pracc
        leo_pracs(is:ie,:) = ppleo_pracs
        leo_praci(is:ie,:) = ppleo_praci
      end if

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

      ! output for cosp
      r_cfrac = cfrac(is:ie,1:kl)
      r_qlrad = qlrad(is:ie,1:kl)
      r_qfrad = qfrad(is:ie,1:kl)
      r_cdrop = cdrop_aerosol(is:ie,1:kl)
      do k=1,kl
        ptemp(:,k)= sig(k)*ps(is:ie)
      end do
      ttemp       = t(is:ie,1:kl)
      call cloud3(Rdrop,Rice,conl,coni,r_cfrac,r_qlrad,r_qfrad,ptemp,ttemp,r_cdrop,imax,kl)

      its = is
      ite = ie
      kts = 1
      kte = kl
      do k = kts, kte
        kr = kte + kts - k
        do iq = its, ite
          i = iq - its + 1 ! i must be smaller than 97
          stras_rliq(iq,kr) = real(Rdrop(i,k))/1.E6 ! unit meters
          stras_rice(iq,kr) = real(Rice(i,k))/1.E6
          stras_rsno(iq,kr) = 0.
          stras_rrai(iq,kr) = 0.
        end do
      end do

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
      end if
    end do
    !$omp end do nowait

  case("LIN")
    cnt_sny = 0
    !if ( myid==0 ) then
    !  write(6,*) "LIN microphysics ",ncloud
    !end if
 
      ccn0 = 250
      riz  = 0
      rhoe_s=1.29

      njumps = int(dt/(maxlintime+0.01)) + 1
      tdt    = real(dt,8)/real(njumps,8)

      !$omp do schedule(static) private(is,ie),                                 &
      !$omp private(k,kts,kte,zlevv,m,ht),                                      &
      !$omp private(qvz,qlz,qrz,qiz,qsz),                                       &
      !$omp private(prf_temp,prf,tothz,thz,rhoz,orhoz,prez,sqrhoz,zz,dzw),      &
      !$omp private(ncz,zcdrop,nrz,niz,nsz,pptrain,pptsnow,pptice)
      do tile = 1, ntiles
        is = (tile-1)*imax + 1 ! is:ie inside 1:ifull
        ie = tile*imax       ! len(is:ie) = imax

        kts = 1
        kte = kl
        ! pack data from ifull into imax
        zlevv(1:imax,1)   = bet(1)*t(is:ie,1)/grav + zs(is:ie)/grav  
        do m=2,kl
          zlevv(1:imax,m) = zlevv(1:imax,m-1) + (bet(m)*t(is:ie,m)+betm(m)*t(is:ie,m-1))/grav
        end do
        ht(1:imax)        = real( zs(is:ie)/grav, 8 )
        
        qvz(1:imax,:) = real(qg(is:ie,:),8)
        qlz(1:imax,:) = real(qlg(is:ie,:),8)  
        qrz(1:imax,:) = real(qrg(is:ie,:),8)
        qiz(1:imax,:) = real(qfg(is:ie,:),8)  
        qsz(1:imax,:) = real(qsng(is:ie,:),8)
        ! ----------------
        do k = 1,kl
          prf_temp(1:imax)  = ps(is:ie)*sig(k)
          prf(1:imax)       = 0.01*prf_temp(1:imax)                        ! ps is SI units
          tothz(1:imax,k)   = real( (prf(1:imax)/1000.)**(rdry/cp), 8 )
          thz(1:imax,k)     = real(t(is:ie,k),8) / tothz(1:imax,k)
          rhoz(1:imax,k)    = real(rhoa(is:ie,k),8)
          orhoz(1:imax,k)   = 1._8/rhoz(1:imax,k)
          prez(1:imax,k)    = real( sig(k)*ps(is:ie), 8 )
          sqrhoz(1:imax,k)  = 1.0_8
          zz(1:imax,k)      = real( zlevv(1:imax,k), 8 )
          dzw(1:imax,k)     = real( dz(is:ie,k), 8 )       
        end do

        select case( lin_aerosolmode )
          case(0)
            ncz(1:imax,:)    = 0.                   
            zcdrop(1:imax,:) = 0.                   
          case(1)
            ncz(1:imax,:)    = 0.                   
            zcdrop(1:imax,:) = real(cdrop(is:ie,:)) ! aerosol
        end select

        nrz(1:imax,:)  = real( nr(is:ie,:), 8 )
        niz(1:imax,:)  = real( ni(is:ie,:), 8 )
        nsz(1:imax,:)  = real( ns(is:ie,:), 8 )


        pptrain(1:imax) = 0._8
        pptsnow(1:imax) = 0._8
        pptice(1:imax)  = 0._8

        ! Use sub time-step if required
        do n = 1,njumps
          CALL clphy1d_ylin(tdt, imax,                      &
                         qvz, qlz, qrz, qiz, qsz,           &
                         thz, tothz, rhoz, orhoz, sqrhoz,   &
                         prez, zz, dzw, ht,                 &
                         precrz, preciz, precsz,            & !zdc 20220116
                         EFFC1D, EFFI1D, EFFS1D, EFFR1D,    & !zdc 20220208
                         pptrain, pptsnow, pptice,          &
                         kts, kte, riz,                     &
                         ncz, nrz, niz, nsz,                &
                         zfluxr,zfluxi,zfluxs,zfluxg,zfluxm,&
                         zfluxf,zfevap,zfsubl,zfauto,zfcoll,&
                         zfaccr,zvi,zvs,zvg,                & !aerosol scheme
                         zpsnow,zpsaut,zpsfw,zpsfi,zpraci,  & !process rate cloud microphysics
                         zpiacr,zpsaci,zpsacw,zpsdep,       &
                         zpssub,zpracs,zpsacr,zpsmlt,       &
                         zpsmltevp,zprain,zpraut,zpracw,    &
                         zprevp,zpgfr,zpvapor,zpclw,        &
                         zpladj,zpcli,zpimlt,zpihom,        &
                         zpidw,zpiadj,zpmidep,               &
                         zcdrop, lin_aerosolmode)              !aerosol feedback
        end do

        t(is:ie,:) = real(thz(1:imax,:) * tothz(1:imax,:))
        gfrac(is:ie,:) = 0.        ! graupel area fraction

        where ( qrz(1:imax,:)>0. )
           rfrac(is:ie,:) = 1.
        elsewhere
           rfrac(is:ie,:) = 0.
        end where
        where ( qsz(1:imax,:)>0. )
           sfrac(is:ie,:) = 1.
        elsewhere
           sfrac(is:ie,:) = 0.
        end where

        !unpack data from imax to ifull.

        qg(is:ie,:)         = real(qvz(1:imax,:))                             ! qv mixing ratio
        qlg(is:ie,:)        = real(qlz(1:imax,:))                             ! ql mixing ratio
        qfg(is:ie,:)        = real(qiz(1:imax,:))                             ! qf mixing ratio (ice)
        qrg(is:ie,:)        = real(qrz(1:imax,:))                             ! qr mixing ratio (rain)
        qsng(is:ie,:)       = real(qsz(1:imax,:))*(1.-real(riz(1:imax,:)))    ! qs mixing ratio (snow)
        qgrg(is:ie,:)       = real(qsz(1:imax,:))*real(riz(1:imax,:))         ! qg mixing ration (graupel)
        stratcloud(is:ie,:) = cfrac(is:ie,:)

        !nc(is:ie,:)        = ncz(:)
        nr(is:ie,:)         = nrz(1:imax,:)
        ni(is:ie,:)         = niz(1:imax,:)
        ns(is:ie,:)         = nsz(1:imax,:)    
        !stras_rliq(is:ie,:) = real(EFFC1D(1:imax,:))             ! save efflective radius for cosp
        !stras_rice(is:ie,:) = real(EFFI1D(1:imax,:))
        stras_rsno(is:ie,:) = real(EFFS1D(1:imax,:))
        stras_rrai(is:ie,:) = real(EFFR1D(1:imax,:))

        fluxr(is:ie,:) = zfluxr(1:imax,:)                        ! flux for aerosol calculation
        fluxm(is:ie,:) = zfluxm(1:imax,:)
        fluxf(is:ie,:) = zfluxf(1:imax,:)
        fluxi(is:ie,:) = zfluxi(1:imax,:)
        fluxs(is:ie,:) = zfluxs(1:imax,:)
        fluxg(is:ie,:) = zfluxg(1:imax,:)
        fevap(is:ie,:) = zfevap(1:imax,:)
        fsubl(is:ie,:) = zfsubl(1:imax,:)
        fauto(is:ie,:) = zfauto(1:imax,:)
        fcoll(is:ie,:) = zfcoll(1:imax,:)
        faccr(is:ie,:) = zfaccr(1:imax,:)
        vi(is:ie,:) = zvi(1:imax,:)
        vs(is:ie,:) = zvs(1:imax,:)
        vg(is:ie,:) = zvg(1:imax,:)

        if (process_rate_mode == 2) then
          psnow(is:ie,:)   = real(zpsnow(1:imax,:))  !process rate to understand cloud microphysics
          psaut(is:ie,:)   = real(zpsaut(1:imax,:))
          psfw(is:ie,:)    = real(zpsfw(1:imax,:))
          psfi(is:ie,:)    = real(zpsfi(1:imax,:))
          praci(is:ie,:)   = real(zpraci(1:imax,:))
          piacr(is:ie,:)   = real(zpiacr(1:imax,:))
          psaci(is:ie,:)   = real(zpsaci(1:imax,:))
          psacw(is:ie,:)   = real(zpsacw(1:imax,:))
          psdep(is:ie,:)   = real(zpsdep(1:imax,:))
          pssub(is:ie,:)   = real(zpssub(1:imax,:))
          pracs(is:ie,:)   = real(zpracs(1:imax,:))
          psacr(is:ie,:)   = real(zpsacr(1:imax,:))
          psmlt(is:ie,:)   = real(zpsmlt(1:imax,:))
          psmltevp(is:ie,:)= real(zpsmltevp(1:imax,:))
          prain(is:ie,:)   = real(zprain(1:imax,:))
          praut(is:ie,:)   = real(zpraut(1:imax,:))
          pracw(is:ie,:)   = real(zpracw(1:imax,:))
          prevp(is:ie,:)   = real(zprevp(1:imax,:))
          pgfr(is:ie,:)    = real(zpgfr(1:imax,:))
          pvapor(is:ie,:)  = real(zpvapor(1:imax,:))
          pclw(is:ie,:)    = real(zpclw(1:imax,:))
          pladj(is:ie,:)   = real(zpladj(1:imax,:))
          pcli(is:ie,:)    = real(zpcli(1:imax,:))
          pimlt(is:ie,:)   = real(zpimlt(1:imax,:))
          pihom(is:ie,:)   = real(zpihom(1:imax,:))
          pidw(is:ie,:)    = real(zpidw(1:imax,:))
          piadj(is:ie,:)   = real(zpiadj(1:imax,:))
          pmidep(is:ie,:)  = real(zpmidep(1:imax,:))
        end if
          
        condx(is:ie)  = condx(is:ie) + real( pptrain(1:imax) + pptsnow(1:imax) + pptice(1:imax) )
        conds(is:ie)  = conds(is:ie) + real( pptsnow(1:imax) + pptice(1:imax) )
        condg(is:ie)  = 0.0          !condg(is:ie) + 0. ! for graupel
        precip(is:ie) = precip(is:ie)+ real( pptrain(1:imax) + pptsnow(1:imax) + pptice(1:imax) )


        r_cfrac = cfrac(is:ie,1:kl)
        r_qlrad = qlrad(is:ie,1:kl)
        r_qfrad = qfrad(is:ie,1:kl)
        r_cdrop = cdrop_aerosol(is:ie,1:kl)
        r_effc1d = real(effc1d)
        r_effi1d = real(effi1d)
        r_effs1d = real(effs1d)
        do k=1,kl
          ptemp(:,k)= sig(k)*ps(is:ie)
        end do
        ttemp       = t(is:ie,1:kl)
        call cloud3(Rdrop,Rice,conl,coni,r_cfrac,r_qlrad,r_qfrad,ptemp,ttemp,r_cdrop,imax,kl, &
                    stras_rliq=r_effc1d,stras_rice=r_effi1d,stras_rsno=r_effs1d)
        its = is
        ite = ie
        kts = 1
        kte = kl
        do k = kts, kte
          kr = kte + kts - k
          do iq = its, ite
            i = iq - its + 1 ! i must be smaller than 97
            stras_rliq(iq,kr) = real(Rdrop(i,k))/1.E6 ! unit meters
            stras_rice(iq,kr) = real(Rice(i,k))/1.E6
            stras_rsno(iq,k) = real(EFFS1D(i,k))
            stras_rrai(iq,k) = real(EFFR1D(i,k))
          end do
        end do
        
      end do     !tile loop
      !$omp end do nowait
      
      !if ( myid==0 ) then
      !  write(6,*) "DONE LIN microphysics ",ncloud
      !end if


  case default
    write(6,*) "ERROR: unknown mp_physics option "
    call ccmpi_abort(-1)
      
end select
  

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
      pplambs(is:ie,kl+1-k) = 1.6e3*10**(-0.023*(t(is:ie,k)-tfrz))          
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
  

#ifdef COSPP
  call cloud_simulator()
#endif 
  
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
