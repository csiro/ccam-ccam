! Conformal Cubic Atmospheric Model
    
! Copyright 2015-2026 Commonwealth Scientific Industrial Research Organisation (CSIRO)
    
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
    
! CCAM boundary layer turbulent mixing routines

! Currently local Ri and prognostic k-e schemes are supported.
! Also, options for non-local counter gradient terms are included.
! Local Ri supports Gelyen and Tiedtke schemes for shallow
! convection, whereas the k-e follows the EDMF approach where
! shallow convection is represented in the mass flux terms.

! nvmix=3  Local Ri mixing
! nvmix=6  Prognostic k-e tubulence closure
! nvmix=7  Jing Huang's local Ri scheme (with axmlsq=9.)
! nvmix=9  Coupled atmosphere-ocean mixing with k-e closure

! nlocal=0 No counter gradient term
! nlocal=6 Holtslag and Boville non-local term
! nlocal=7 Mass flux based counter gradient (requires nvmix=6 or nvmix=9)
    
module module_ctrl_turbmix

implicit none

private
public turbmix_init, turbmix

contains

subroutine turbmix_init

use arrays_m                        ! Atmosphere dyamics prognostic arrays
use cc_mpi                          ! CC MPI routines
use cfrac_m                         ! Cloud fraction
use const_phys                      ! Physical constants
use kuocom_m                        ! JLM convection
use liqwpar_m                       ! Cloud water mixing ratios
use newmpar_m                       ! Grid parameters
use nharrs_m                        ! Non-hydrostatic atmosphere arrays
use parm_m                          ! Model configuration
use sigs_m                          ! Atmosphere sigma levels
use tkeeps                          ! TKE-EPS boundary layer
use vertmix_m                       ! Local RI boundary layer

implicit none

integer k
real, dimension(kl) :: sigkap

select case ( interp_nvmix(nvmix) )
  case("localRI")
    ! set ksctop for shallow convection
    ksctop = 1    ! ksctop will be first level below sigkcst
    do while( sig(ksctop+1)>sigksct .and. sigksct>0. )  !  e.g. sigksct=.75
      ksctop = ksctop + 1
    end do
    kscbase = 1  ! kscbase will be first level above sigkcsb
    do while( sig(kscbase)>sigkscb .and. sigkscb>0. ) ! e.g. sigkscb=.99
      kscbase = kscbase + 1
    end do
    if ( ksc/=0 .and. nvmix/=6 ) then
      if ( myid==0 ) then
        write(6,*)'For shallow convection:'
        write(6,*)'ksc,kscbase,ksctop,kscsea ',ksc,kscbase,ksctop,kscsea
        write(6,"(' sigkscb,sigksct,tied_con,tied_over,tied_rh:',5f8.3)") sigkscb,sigksct,tied_con,tied_over,tied_rh
      end if
    end if  
    
  case("TKE-EPS")  
    ! reset time averaging at start of simulation
    if ( nvmix==6 .or. nvmix==9 ) then
      if ( .not.lrestart ) then
        if ( myid==0 ) then
          write(6,*) "Reset EMA for tke-eps"
        end if  
        do k = 1,kl
          sigkap(k)       = sig(k)**(-rdry/cp)
          thetal_ema(:,k) = (t(1:ifull,k)-(hl*qlg(1:ifull,k)+hls*qfg(1:ifull,k))/cp)*sigkap(k)
          qv_ema(:,k)     = qg(1:ifull,k)
          ql_ema(:,k)     = qlg(1:ifull,k)
          qf_ema(:,k)     = qfg(1:ifull,k)
          cf_ema(:,k)     = stratcloud(1:ifull,k)
          u_ema(:,k)      = u(1:ifull,k)
          v_ema(:,k)      = v(1:ifull,k)
          w_ema(:,k)      = 0.
          tke_ema(:,k)    = tke(1:ifull,k)
        end do
      end if
    end if  

  case default
    write(6,*) "ERROR: Unknown turbulent mixing option"
    
end select
    
return
end subroutine turbmix_init

subroutine turbmix

use arrays_m                        ! Atmosphere dyamics prognostic arrays
use cc_mpi                          ! CC MPI routines
use cfrac_m                         ! Cloud fraction
use const_phys                      ! Physical constants
use diag_m                          ! Diagnostic routines
use extraout_m                      ! Additional diagnostics
use kuocom_m                        ! JLM convection
use liqwpar_m                       ! Cloud water mixing ratios
use map_m                           ! Grid map arrays
use morepbl_m                       ! Additional boundary layer diagnostics
use newmpar_m                       ! Grid parameters
use nharrs_m                        ! Non-hydrostatic atmosphere arrays
use nsibd_m                         ! Land-surface arrays
use parm_m, only : idjd, nmlo, nvmix, dt
                                    ! Model configuration
use pbl_m                           ! Boundary layer arrays
use savuvt_m                        ! Saved dynamic arrays
use screen_m                        ! Screen level diagnostics
use sigs_m                          ! Atmosphere sigma levels
use soil_m, only : land             ! Soil and surface data
use soilsnow_m, only : fracice      ! Soil, snow and surface data
use tkeeps                          ! TKE-EPS boundary layer
use trimmix_m                       ! Tridiagonal solver for turbulent mixing
use work2_m                         ! Diagnostic arrays
use vertmix_m                       ! Local RI boundary layer

implicit none

integer is, ie, tile, k, iq
integer idjd_t
real, dimension(imax,kl) :: lt, lqg, lqfg, lqlg, lni
real, dimension(imax,kl) :: lu, lv, lstratcloud
real, dimension(imax,kl) :: lsavu, lsavv, ltke, leps, lshear
real, dimension(imax,kl) :: lthetal_ema, lqv_ema, lql_ema, lqf_ema, lcf_ema
real, dimension(imax,kl) :: ltke_ema
real, dimension(imax,kl) :: lrkmsave, lrkhsave
real, dimension(imax,kl) :: lat, lct, lqtr
real, dimension(imax) :: lou, lov, liu, liv, lrho
real tmnht, dzz, gt, rlogs1, rlogs2, rlogh1, rlog12, rong
logical mydiag_t

select case ( interp_nvmix(nvmix) )
  case("TKE-EPS")  
    ! k-e + MF closure scheme
    do tile = 1,ntiles
      is = (tile-1)*imax + 1
      ie = tile*imax
      idjd_t = mod(idjd-1,imax)+1
      mydiag_t = ((idjd-1)/imax==tile-1).and.mydiag
      lt = t(is:ie,:)
      lqg = qg(is:ie,:)
      lqfg = qfg(is:ie,:)
      lqlg = qlg(is:ie,:)
      lni = ni(is:ie,:)
      lstratcloud = stratcloud(is:ie,:)
      ltke   = tke(is:ie,:)
      leps   = eps(is:ie,:)
      lshear = shear(is:ie,:)
      lu     = u(is:ie,:)
      lv     = v(is:ie,:)
      lthetal_ema = thetal_ema(is:ie,:)
      lqv_ema     = qv_ema(is:ie,:)
      lql_ema     = ql_ema(is:ie,:)
      lqf_ema     = qf_ema(is:ie,:)
      lcf_ema     = cf_ema(is:ie,:)
      ltke_ema    = tke_ema(is:ie,:)
      call tkeeps_work(lt,em(is:ie),tss(is:ie),eg(is:ie),fg(is:ie),ps(is:ie),lqg,lqfg,lqlg,   &
                       lni,lstratcloud,cduv(is:ie),lu,lv,pblh(is:ie),                         &
                       ustar(is:ie),ltke,leps,lshear,land(is:ie),lthetal_ema,lqv_ema,lql_ema, &
                       lqf_ema,lcf_ema,ltke_ema,lrkmsave,lrkhsave,f(is:ie),imax,kl,tile)      
      t(is:ie,:)          = lt
      qg(is:ie,:)         = lqg
      qfg(is:ie,:)        = lqfg
      qlg(is:ie,:)        = lqlg
      ni(is:ie,:)         = lni
      stratcloud(is:ie,:) = lstratcloud
      tke(is:ie,:)        = ltke
      eps(is:ie,:)        = leps
      u(is:ie,:)          = lu
      v(is:ie,:)          = lv
      thetal_ema(is:ie,:) = lthetal_ema
      qv_ema(is:ie,:)     = lqv_ema
      ql_ema(is:ie,:)     = lql_ema
      qf_ema(is:ie,:)     = lqf_ema
      cf_ema(is:ie,:)     = lcf_ema      
      tke_ema(is:ie,:)    = ltke_ema
      rkmsave(is:ie,:)    = lrkmsave
      rkhsave(is:ie,:)    = lrkhsave  
      lrho(1:imax) = ps(is:ie)/(rdry*tss(is:ie))
      taux(is:ie) = lrho(1:imax)*cduv(is:ie)*u(is:ie,1)
      tauy(is:ie) = lrho(1:imax)*cduv(is:ie)*v(is:ie,1)
    end do ! tile = 1,ntiles

    
  case("localRI")
    ! JLM's local Ri scheme
    do tile = 1,ntiles
      is = (tile-1)*imax + 1
      ie = tile*imax
      idjd_t = mod(idjd-1,imax)+1
      mydiag_t = ((idjd-1)/imax==tile-1).and.mydiag
      lt = t(is:ie,:)
      lqg = qg(is:ie,:)
      lqfg = qfg(is:ie,:)
      lqlg = qlg(is:ie,:)
      lni = ni(is:ie,:)
      lstratcloud = stratcloud(is:ie,:)
      lu = u(is:ie,:)
      lv = v(is:ie,:)
      lsavu = savu(is:ie,:)
      lsavv = savv(is:ie,:)
      call vertmix(lt,tss(is:ie),eg(is:ie),fg(is:ie),kbsav(is:ie),ktsav(is:ie),convpsav(is:ie),  &
                   ps(is:ie),lqg,lqfg,lqlg,lni,lstratcloud,                                      &
                   condc(is:ie),cduv(is:ie),lu,lv,pblh(is:ie),lsavu,lsavv,land(is:ie),           &
                   tscrn(is:ie),qgscrn(is:ie),ustar(is:ie),f(is:ie),condx(is:ie),zs(is:ie),      &
                   lrkmsave,lrkhsave,                                                            &
                   idjd_t,mydiag_t)
      t(is:ie,:)          = lt
      qg(is:ie,:)         = lqg
      qfg(is:ie,:)        = lqfg
      qlg(is:ie,:)        = lqlg
      ni(is:ie,:)         = lni
      stratcloud(is:ie,:) = lstratcloud
      rkmsave(is:ie,:)    = lrkmsave
      rkhsave(is:ie,:)    = lrkhsave  
      u(is:ie,:)          = lu
      v(is:ie,:)          = lv
    end do ! tile = 1,ntiles

end select

return
end subroutine turbmix
    
!====================================================================================================
! SUBROUTINE interp_nvmix
!   
! subroutine to select the cloud microphysics scheme for CCAM
!====================================================================================================
pure function interp_nvmix(nvmix) result(mp_physics)

implicit none

integer, intent(in) :: nvmix
character(len=10) :: mp_physics

mp_physics = "ERROR"

select case(nvmix)
  case(6,9)
    mp_physics = "TKE-EPS"
  case default
    mp_physics = "localRI"
end select

return
end function interp_nvmix  


subroutine tkeeps_work(t,em,tss,eg,fg,ps,qg,qfg,qlg,ni,stratcloud,                &
                       cduv,u,v,pblh,ustar,tke,eps,shear,land,thetal_ema,qv_ema,  &
                       ql_ema,qf_ema,cf_ema,tke_ema,rkmsave,rkhsave,f,imax,kl,    &
                       tile)

use const_phys                   ! Physical constants
use parm_m, only : ds, nlocal, dt, cqmix, nvmix
                                 ! Model configuration
use sigs_m                       ! Atmosphere sigma levels
use tkeeps, only : tkemix        ! TKE-EPS boundary layer

implicit none

integer, intent(in) :: imax, kl, tile
integer k, nlocal_mode
real, dimension(imax,kl), intent(inout) :: t, qg, qfg, qlg, ni
real, dimension(imax,kl), intent(inout) :: stratcloud, u, v
real, dimension(imax,kl), intent(inout) :: tke, eps
real, dimension(imax,kl), intent(in) :: shear
real, dimension(imax,kl), intent(inout) :: thetal_ema, qv_ema, ql_ema, qf_ema, cf_ema
real, dimension(imax,kl), intent(inout) :: tke_ema
real, dimension(imax,kl), intent(out) :: rkmsave, rkhsave
real, dimension(imax,kl) :: zh
real, dimension(imax,kl) :: rhs, zg
real, dimension(imax,kl) :: rkm
real, dimension(imax), intent(inout) :: pblh, ustar, eg, fg, tss
real, dimension(imax), intent(in) :: em, ps, f
real, dimension(imax), intent(in) :: cduv
real, dimension(imax) :: rhos, dx
real, dimension(kl) :: sigkap, delh
real rong
logical, dimension(imax), intent(in) :: land


! estimate grid spacing for scale aware MF
dx(:) = ds/em(:)

! Set-up potential temperature transforms
rong = rdry/grav
do k = 1,kl
  delh(k)   = -rong*dsig(k)/sig(k)  ! sign of delh defined so always +ve
  sigkap(k) = sig(k)**(-rdry/cp)
end do      ! k loop
      
! note ksc/=0 options are clobbered when nvmix=6 or nvmix=9
! However, nvmix=6 or nvmix=9 supports its own shallow
! convection options with nlocal=7

! calculate height on full levels AGL
zg(:,1) = bet(1)*t(:,1)/grav
zh(:,1) = t(:,1)*delh(1)
do k = 2,kl
  zg(:,k) = zg(:,k-1) + (bet(k)*t(:,k)+betm(k)*t(:,k-1))/grav
  zh(:,k) = zh(:,k-1) + t(:,k)*delh(k)
end do ! k  loop
       
! near surface air density (see sflux.f and cable_ccam2.f90)
rhos(:) = ps(:)/(rdry*tss(:))
    
! transform to ocean reference frame and temp to theta
do k = 1,kl
  rhs(:,k) = t(:,k)*sigkap(k) ! theta
end do

! Evaluate EDMF scheme
select case(nvmix)
  case(6)  
    select case(nlocal)
      case(0) ! atm only, no counter gradient
        nlocal_mode = 1
      case(7) ! atm only, mass-flux counter gradient
        nlocal_mode = 0
      case default
        write(6,*) "ERROR: Invalid nlocal = ",nlocal
        stop -1
    end select
  case(9)
    select case(nlocal)
      case(0) ! coupled atm-ocn, no counter gradient
        nlocal_mode = 3
      case(7) ! coupled atm-ocn, mass-flux counter gradient
        nlocal_mode = 2
      case default
        write(6,*) "ERROR: Invalid nlocal = ",nlocal
        stop -1
    end select
end select
    
call tkemix(rkm,rhs,qg,qlg,qfg,ni,stratcloud,u,v,pblh,fg,eg,tss,               &
            cduv,ps,zg,zh,sig,rhos,ustar,dt,nlocal_mode,tke,eps,shear,         &
            dx,thetal_ema,qv_ema,ql_ema,qf_ema,cf_ema,tke_ema,land,tile,imax,  &
            kl)  

do k = 1,kl
  ! replace counter gradient term
  ! save Km and Kh for output
  rkmsave(:,k) = rkm(:,k)*cqmix
  rkhsave(:,k) = rkm(:,k)*cqmix
  ! transform winds back to Earth reference frame and theta to temp
  t(:,k) = rhs(:,k)/sigkap(k)
enddo    !  k loop

return
end subroutine tkeeps_work

end module module_ctrl_turbmix