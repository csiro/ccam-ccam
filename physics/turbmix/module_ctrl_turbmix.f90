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
    ! do nothing

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
use tkeeps                          ! TKE-EPS boundary layer
use trimmix_m                       ! Tridiagonal solver for turbulent mixing
use work2_m                         ! Diagnostic arrays
use vertmix_m                       ! Local RI boundary layer

implicit none

integer is, ie, tile, k, iq
integer idjd_t
real, dimension(imax,kl) :: lt, lqg, lqfg, lqlg, lni
real, dimension(imax,kl) :: lu, lv, lstratcloud
real, dimension(imax,kl) :: lsavu, lsavv
real, dimension(imax,kl) :: lrkmsave, lrkhsave
real tmnht, dzz, gt, rlogs1, rlogs2, rlogh1, rlog12, rong
logical mydiag_t

select case ( interp_nvmix(nvmix) )
  case("TKE-EPS")  
    ! k-e + MF closure scheme
    call tkeeps_work      
    
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
                   lrkmsave,lrkhsave,idjd_t,mydiag_t)
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


subroutine tkeeps_work

use arrays_m                        ! Atmosphere dyamics prognostic arrays
use cfrac_m                         ! Cloud fraction
use const_phys                      ! Physical constants
use extraout_m                      ! Additional diagnostics
use liqwpar_m                       ! Cloud water mixing ratios
use map_m                           ! Grid map arrays
use mlo_ctrl                        ! Ocean physics control layer
use morepbl_m                       ! Additional boundary layer diagnostics
use newmpar_m                       ! Grid parameters
use parm_m, only : ds, nlocal, dt, cqmix, nvmix
                                    ! Model configuration
use pbl_m                           ! Boundary layer arrays
use sigs_m                          ! Atmosphere sigma levels
use soil_m, only : land             ! Soil and surface data
use soilsnow_m                      ! Soil, snow and surface data
use tkeeps, only : tkemix, tke, eps, shear
                                    ! TKE-EPS boundary layer

implicit none

integer k, nlocal_mode
integer tile, js, je
#ifdef GPU
integer, dimension(ifull) :: ibot, mixind_o
real, dimension(ifull,kl) :: zh, zg
real, dimension(ifull,kl) :: rkm
real, dimension(ifull,kl) :: lqg, lqfg, lqlg, lstratcloud, lni
real, dimension(ifull,kl) :: ltke, leps, lshear, lu, lv
real, dimension(ifull,kl) :: rhs
real, dimension(ifull) :: rhos, dx
real, dimension(ifull) :: lrho
real, dimension(kl) :: sigkap, delh
real, dimension(ifull,ol) :: w_t, w_s, w_u, w_v
real, dimension(ifull,ol) :: w_tke, w_eps
real, dimension(ifull,ol) :: dwdx_o, dwdy_o
real, dimension(ifull,ol) :: deptho_depth
real, dimension(ifull,ol+1) :: deptho_depth_hl
real, dimension(ifull,ol) :: deptho_dz
real, dimension(ifull,2:ol) :: deptho_dz_hl
real, dimension(ifull) :: w_e, d_zcr, deptho_max, rbot
real, dimension(ifull) :: ubot_o, vbot_o, utop_o, vtop_o
real, dimension(ifull) :: cd_water, cdh_water, cdbot_water
real, dimension(ifull) :: wt0_o, wt0rad_o, wt0melt_o, wt0eg_o, wt0fb_o
real, dimension(ifull) :: ws0_o, ws0subsurf_o
real, dimension(ifull) :: wu0_o, wv0_o, zo_o
real, dimension(ifull) :: bf_o, mixdepth_o
real, dimension(ifull,ol) :: rad_o
real, dimension(ifull) :: i_u, i_v
real, dimension(ifull) :: icefg_a, imass, cd_ice, cdbot_ice
#else
integer, dimension(imax) :: ibot, mixind_o
real, dimension(imax,kl) :: zh, zg
real, dimension(imax,kl) :: rkm
real, dimension(imax,kl) :: lqg, lqfg, lqlg, lstratcloud, lni
real, dimension(imax,kl) :: ltke, leps, lshear, lu, lv
real, dimension(imax,kl) :: rhs
real, dimension(imax) :: rhos, dx
real, dimension(imax) :: lrho
real, dimension(kl) :: sigkap, delh
real, dimension(imax,ol) :: w_t, w_s, w_u, w_v
real, dimension(imax,ol) :: w_tke, w_eps
real, dimension(imax,ol) :: dwdx_o, dwdy_o
real, dimension(imax,ol) :: deptho_depth
real, dimension(imax,ol+1) :: deptho_depth_hl
real, dimension(imax,ol) :: deptho_dz
real, dimension(imax,2:ol) :: deptho_dz_hl
real, dimension(imax) :: w_e, d_zcr, deptho_max, rbot
real, dimension(imax) :: ubot_o, vbot_o, utop_o, vtop_o
real, dimension(imax) :: cd_water, cdh_water, cdbot_water
real, dimension(imax) :: wt0_o, wt0rad_o, wt0melt_o, wt0eg_o, wt0fb_o
real, dimension(imax) :: ws0_o, ws0subsurf_o
real, dimension(imax) :: wu0_o, wv0_o, zo_o
real, dimension(imax) :: bf_o, mixdepth_o
real, dimension(imax,ol) :: rad_o
real, dimension(imax) :: i_u, i_v
real, dimension(imax) :: icefg_a, imass, cd_ice, cdbot_ice
#endif
real rong


! Set-up potential temperature transforms
rong = rdry/grav
do k = 1,kl
  delh(k)   = -rong*dsig(k)/sig(k)  ! sign of delh defined so always +ve
  sigkap(k) = sig(k)**(-rdry/cp)
end do      ! k loop

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

#ifdef GPU  
! GPU version

! estimate grid spacing for scale aware MF
dx(1:ifull) = ds/em(1:ifull)
     
! note ksc/=0 options are clobbered when nvmix=6 or nvmix=9
! However, nvmix=6 or nvmix=9 supports its own shallow
! convection options with nlocal=7

! calculate height on full levels AGL
zg(1:ifull,1) = bet(1)*t(1:ifull,1)/grav
zh(1:ifull,1) = t(1:ifull,1)*delh(1)
do k = 2,kl
  zg(1:ifull,k) = zg(1:ifull,k-1) + (bet(k)*t(1:ifull,k)+betm(k)*t(1:ifull,k-1))/grav
  zh(1:ifull,k) = zh(1:ifull,k-1) + t(1:ifull,k)*delh(k)
end do ! k  loop
       
! near surface air density (see sflux.f and cable_ccam2.f90)
rhos(1:ifull) = ps(1:ifull)/(rdry*tss(1:ifull))
    
lqg  = qg(1:ifull,:)
lqfg = qfg(1:ifull,:)
lqlg = qlg(1:ifull,:)
lni  = ni(1:ifull,:)
lstratcloud = stratcloud(1:ifull,:)
ltke   = tke(1:ifull,:)
leps   = eps(1:ifull,:)
lshear = shear(1:ifull,:)
lu     = u(1:ifull,:)
lv     = v(1:ifull,:)
    
! transform to ocean reference frame and temp to theta
do k = 1,kl
  rhs(1:ifull,k) = t(1:ifull,k)*sigkap(k) ! theta
end do
  
w_e = 0.
deptho_max = 0.
d_zcr = 1.
deptho_depth = 0.
deptho_depth_hl = 0.
deptho_dz = 0.
deptho_dz_hl = 0.
w_t = 0.
w_s = 0.
w_u = 0.
w_v = 0.
dwdx_o = 0.
dwdy_o = 0.
ibot = ol
ubot_o = 0.
vbot_o = 0.
utop_o = 0.
vtop_o = 0.
cd_water = 0.
cdh_water = 0.
cdbot_water = 0.
wt0_o = 0.
wt0rad_o = 0.
wt0melt_o = 0.
wt0eg_o = 0.
wt0fb_o = 0.
ws0_o = 0.
ws0subsurf_o = 0.
wu0_o = 0.
wv0_o = 0.
zo_o = 0. 
bf_o = 0.
mixdepth_o = 0.
mixind_o = ol
rad_o = 0.
w_tke = 1.e-8
w_eps = 1.e-11
i_u = 0.
i_v = 0.
icefg_a = 0.
imass = 100.
cd_ice = 0.
cdbot_ice = 0.
if ( nlocal_mode==2 .or. nlocal_mode==3 ) then
  ! export ocean data before mixing  
  call mloexport("eta",w_e,0,0,ifull,imax)
  call mloexpdep("depth_hl",deptho_max,ol+1,0,ifull,imax)
  d_zcr(:) = max( 1. + max(w_e(:),-delwater)/max(deptho_max(:),1.e-4), minwater/max(deptho_max(:),1.e-4) )
  do k = 1,ol
    call mloexpdep("depth_fl",deptho_depth(:,k),k,0,ifull,imax)
    deptho_depth(:,k) = deptho_depth(:,k)*d_zcr(:)
  end do
  do k = 1,ol+1
    call mloexpdep("depth_hl",deptho_depth_hl(:,k),k,0,ifull,imax)
    deptho_depth_hl(:,k) = deptho_depth_hl(:,k)*d_zcr(:)
  end do
  do k = 1,ol
    call mloexpdep("depth_dz_fl",deptho_dz(:,k),k,0,ifull,imax)
    deptho_dz(:,k) = deptho_dz(:,k)*d_zcr(:)
  end do
  do k = 2,ol
    call mloexpdep("depth_dz_hl",deptho_dz_hl(:,k),k,0,ifull,imax)
    deptho_dz_hl(:,k) = deptho_dz_hl(:,k)*d_zcr(:)
  end do
  do k = 1,ol
    call mloexport("temp",w_t(:,k),k,0,ifull,imax)
    call mloexport("sal",w_s(:,k),k,0,ifull,imax)
    call mloexport("u",w_u(:,k),k,0,ifull,imax)
    call mloexport("v",w_v(:,k),k,0,ifull,imax)
    call mloexport("dwdx",dwdx_o(:,k),k,0,ifull,imax)
    call mloexport("dwdy",dwdy_o(:,k),k,0,ifull,imax)  
  end do
  rbot(:) = real(ibot(:))
  call mloexport("ibot",rbot(:),0,0,ifull,imax)
  ibot(:) = nint(rbot(:))
  call mloexport("ubot",ubot_o(:),0,0,ifull,imax)
  call mloexport("vbot",vbot_o(:),0,0,ifull,imax)
  call mloexport("utop",utop_o(:),0,0,ifull,imax)
  call mloexport("vtop",vtop_o(:),0,0,ifull,imax)
  call mlodiag("cd",cd_water(:),0,0,ifull,imax)
  call mlodiag("cdh",cdh_water(:),0,0,ifull,imax)
  call mlodiag("cd_bot",cdbot_water(:),0,0,ifull,imax)
  call mlodiag("wt0",wt0_o(:),0,0,ifull,imax)
  call mlodiag("wt0_rad",wt0rad_o(:),0,0,ifull,imax)
  call mlodiag("wt0_melt",wt0melt_o(:),0,0,ifull,imax)
  call mlodiag("wt0_eg",wt0eg_o(:),0,0,ifull,imax)
  call mlodiag("wt0_fb",wt0fb_o(:),0,0,ifull,imax)
  call mlodiag("ws0",ws0_o(:),0,0,ifull,imax)
  call mlodiag("ws0_subsurf",ws0subsurf_o(:),0,0,ifull,imax)
  call mlodiag("wu0",wu0_o(:),0,0,ifull,imax)
  call mlodiag("wv0",wv0_o(;),0,0,ifull,imax)
  call mlodiag("zo",zo_o(:),0,0,ifull,imax)   
  call mlodiag("bf",bf_o(:),0,0,ifull,imax)
  call mlodiag("mixdepth",mixdepth_o(:),0,0,ifull,imax)
  rbot(:) = real(mixind_o(:))
  call mlodiag("mixind",rbot(:),0,0,ifull,imax)
  mixind_o(:) = nint(rbot(:))
  do k = 1,ol
    call mlodiag("rad",rad_o(:,k),k,0,ifull,imax)
  end do    
  do k = 1,ol
    call mloexpturb("tke",w_tke(:,k),k,0,ifull,imax)
    call mloexpturb("eps",w_eps(:,k),k,0,ifull,imax)
  end do
  call mloexpice("u",i_u,0,ifull,imax)
  call mloexpice("v",i_v,0,ifull,imax)
  call mlodiagice("fg",ifull,imax)
  call mlodiagice("mass",imass,0,ifull,imax)
  call mlodiagice("cd",cd_ice,0,ifull,imax)
  call mlodiagice("cd_bot",cdbot_ice,0,ifull,imax)    
end if
   
call tkemix(rkm,rhs,lqg,lqlg,lqfg,lni,lstratcloud,lu,lv,pblh,                  &
            fg,eg,cduv,ps,zg,zh,sig,                                           &
            rhos,ustar,dt,nlocal_mode,ltke,leps,lshear,dx,land,                &
            w_t,w_s,w_u,w_v,w_tke,w_eps,deptho_depth,deptho_depth_hl,          &
            deptho_dz,deptho_dz_hl,dwdx_o,dwdy_o,ibot,ubot_o,vbot_o,           &
            utop_o,vtop_o,cd_water,cdh_water,cdbot_water,wt0_o,wt0rad_o,       &
            wt0melt_o,wt0eg_o,wt0fb_o,ws0_o,ws0subsurf_o,wu0_o,wv0_o,zo_o,     &
            bf_o,mixdepth_o,mixind_o,rad_o,                                    &
            i_u,i_v,fracice,icefg_a,imass,cd_ice,cdbot_ice,                    &
            ifull,kl,ol)  
  
qg(1:ifull,:)         = lqg
qfg(1:ifull,:)        = lqfg
qlg(1:ifull,:)        = lqlg
ni(1:ifull,:)         = lni
stratcloud(1:ifull,:) = lstratcloud
tke(1:ifull,:)        = ltke
eps(1:ifull,:)        = leps
u(1:ifull,:)          = lu
v(1:ifull,:)          = lv
  
do k = 1,kl
  ! replace counter gradient term with cqmix
  ! save Km and Kh for output
  rkmsave(:,k) = rkm(:,k)*cqmix
  rkhsave(:,k) = rkm(:,k)*cqmix
  ! transform winds back to Earth reference frame and theta to temp
  t(1:ifull,k) = rhs(1:ifull,k)/sigkap(k)
enddo    !  k loop
  
lrho(1:ifull) = ps(1:ifull)/(rdry*tss(1:ifull))
taux(1:ifull) = lrho(1:ifull)*cduv(1:ifull)*u(1:ifull,1)
tauy(1:ifull) = lrho(1:ifull)*cduv(1:ifull)*v(1:ifull,1)
  
if ( nlocal_mode==2 .or. nlocal_mode==3 ) then
  ! import ocean data after mixing  
  do k = 1,ol
    call mloimport("temp",w_t(:,k),k,0,ifull,imax)
    call mloimport("sal",w_s(:,k),k,0,ifull,imax)
    call mloimport("u",w_u(:,k),k,0,ifull,imax)
    call mloimport("v",w_v(:,k),k,0,ifull,imax)
  end do
  call mloimport("ubot",ubot_o,0,0,ifull,imax)
  call mloimport("vbot",vbot_o,0,0,ifull,imax)
  call mloimport("utop",utop_o,0,0,ifull,imax)
  call mloimport("vtop",vtop_o,0,0,ifull,imax)
  call mlo_updatediag("wu0",wu0_o,0,ifull,imax)
  call mlo_updatediag("wv0",wv0_o,0,ifull,imax)
  call mlo_updatediag("wt0",wt0_o,0,ifull,imax)    
  do k = 1,ol
    call mloimpturb("tke",w_tke(:,k),k,0,ifull,imax)
    call mloimpturb("eps",w_eps(:,k),k,0,ifull,imax)
  end do  
  call mloimpice("u",i_u,0,ifull,imax)
  call mloimpice("v",i_v,0,ifull,imax)    
  call mloimpice("fracice",fracice(js:je),0,ifull,imax)        
  ! export ocean SST
  call mlosurf("sst",tss(js:je),0,ifull,imax)
end if

#else

! CPU version
do tile = 1,ntiles
  js = (tile-1)*imax + 1
  je = tile*imax

  ! estimate grid spacing for scale aware MF
  dx(:) = ds/em(js:je)
     
  ! note ksc/=0 options are clobbered when nvmix=6 or nvmix=9
  ! However, nvmix=6 or nvmix=9 supports its own shallow
  ! convection options with nlocal=7

  ! calculate height on full levels AGL
  zg(:,1) = bet(1)*t(js:je,1)/grav
  zh(:,1) = t(js:je,1)*delh(1)
  do k = 2,kl
    zg(:,k) = zg(:,k-1) + (bet(k)*t(js:je,k)+betm(k)*t(js:je,k-1))/grav
    zh(:,k) = zh(:,k-1) + t(js:je,k)*delh(k)
  end do ! k  loop
       
  ! near surface air density (see sflux.f and cable_ccam2.f90)
  rhos(:) = ps(js:je)/(rdry*tss(js:je))
    
  lqg  = qg(js:je,:)
  lqfg = qfg(js:je,:)
  lqlg = qlg(js:je,:)
  lni  = ni(js:je,:)
  lstratcloud = stratcloud(js:je,:)
  
  ltke   = tke(js:je,:)
  leps   = eps(js:je,:)
  lshear = shear(js:je,:)
  
  lu     = u(js:je,:)
  lv     = v(js:je,:)
    
  ! transform to ocean reference frame and temp to theta
  do k = 1,kl
    rhs(:,k) = t(js:je,k)*sigkap(k) ! theta
  end do
  
  w_e = 0.
  deptho_max = 0.
  d_zcr = 1.
  deptho_depth = 0.
  deptho_depth_hl = 0.
  deptho_dz = 0.
  deptho_dz_hl = 0.
  w_t = 0.
  w_s = 0.
  w_u = 0.
  w_v = 0.
  dwdx_o = 0.
  dwdy_o = 0.
  ibot = ol
  ubot_o = 0.
  vbot_o = 0.
  utop_o = 0.
  vtop_o = 0.
  cd_water = 0.
  cdh_water = 0.
  cdbot_water = 0.
  wt0_o = 0.
  wt0rad_o = 0.
  wt0melt_o = 0.
  wt0eg_o = 0.
  wt0fb_o = 0.
  ws0_o = 0.
  ws0subsurf_o = 0.
  wu0_o = 0.
  wv0_o = 0.
  zo_o = 0. 
  bf_o = 0.
  mixdepth_o = 0.
  mixind_o = ol
  rad_o = 0.
  w_tke = 1.e-8
  w_eps = 1.e-11
  i_u = 0.
  i_v = 0.
  icefg_a = 0.
  imass = 100.
  cd_ice = 0.
  cdbot_ice = 0.
  if ( nlocal_mode==2 .or. nlocal_mode==3 ) then
    ! export ocean data before mixing  
    call mloexport("eta",w_e,0,0,water_g(tile),depth_g(tile),imax)
    call mloexpdep("depth_hl",deptho_max,ol+1,0,depth_g(tile),imax)
    d_zcr = max( 1. + max(w_e,-delwater)/max(deptho_max,1.e-4), minwater/max(deptho_max,1.e-4) )
    do k = 1,ol
      call mloexpdep("depth_fl",deptho_depth(:,k),k,0,depth_g(tile),imax)
      deptho_depth(:,k) = deptho_depth(:,k)*d_zcr
    end do
    do k = 1,ol+1
      call mloexpdep("depth_hl",deptho_depth_hl(:,k),k,0,depth_g(tile),imax)
      deptho_depth_hl(:,k) = deptho_depth_hl(:,k)*d_zcr
    end do
    do k = 1,ol
      call mloexpdep("depth_dz_fl",deptho_dz(:,k),k,0,depth_g(tile),imax)
      deptho_dz(:,k) = deptho_dz(:,k)*d_zcr
    end do
    do k = 2,ol
      call mloexpdep("depth_dz_hl",deptho_dz_hl(:,k),k,0,depth_g(tile),imax)
      deptho_dz_hl(:,k) = deptho_dz_hl(:,k)*d_zcr
    end do
    do k = 1,ol
      call mloexport("temp",w_t(:,k),k,0,water_g(tile),depth_g(tile),imax)
      call mloexport("sal",w_s(:,k),k,0,water_g(tile),depth_g(tile),imax)
      call mloexport("u",w_u(:,k),k,0,water_g(tile),depth_g(tile),imax)
      call mloexport("v",w_v(:,k),k,0,water_g(tile),depth_g(tile),imax)
      call mloexport("dwdx",dwdx_o(:,k),k,0,water_g(tile),depth_g(tile),imax)
      call mloexport("dwdy",dwdy_o(:,k),k,0,water_g(tile),depth_g(tile),imax)  
    end do
    rbot = real(ibot)
    call mloexport("ibot",rbot,0,0,water_g(tile),depth_g(tile),imax)
    ibot = nint(rbot)
    call mloexport("ubot",ubot_o,0,0,water_g(tile),depth_g(tile),imax)
    call mloexport("vbot",vbot_o,0,0,water_g(tile),depth_g(tile),imax)
    call mloexport("utop",utop_o,0,0,water_g(tile),depth_g(tile),imax)
    call mloexport("vtop",vtop_o,0,0,water_g(tile),depth_g(tile),imax)
    call mlodiag("cd",cd_water,0,0,dgwater_g(tile),depth_g(tile),imax)
    call mlodiag("cdh",cdh_water,0,0,dgwater_g(tile),depth_g(tile),imax)
    call mlodiag("cd_bot",cdbot_water,0,0,dgwater_g(tile),depth_g(tile),imax)
    call mlodiag("wt0",wt0_o,0,0,dgwater_g(tile),depth_g(tile),imax)
    call mlodiag("wt0_rad",wt0rad_o,0,0,dgwater_g(tile),depth_g(tile),imax)
    call mlodiag("wt0_melt",wt0melt_o,0,0,dgwater_g(tile),depth_g(tile),imax)
    call mlodiag("wt0_eg",wt0eg_o,0,0,dgwater_g(tile),depth_g(tile),imax)
    call mlodiag("wt0_fb",wt0fb_o,0,0,dgwater_g(tile),depth_g(tile),imax)
    call mlodiag("ws0",ws0_o,0,0,dgwater_g(tile),depth_g(tile),imax)
    call mlodiag("ws0_subsurf",ws0subsurf_o,0,0,dgwater_g(tile),depth_g(tile),imax)
    call mlodiag("wu0",wu0_o,0,0,dgwater_g(tile),depth_g(tile),imax)
    call mlodiag("wv0",wv0_o,0,0,dgwater_g(tile),depth_g(tile),imax)
    call mlodiag("zo",zo_o,0,0,dgwater_g(tile),depth_g(tile),imax)   
    call mlodiag("bf",bf_o,0,0,dgwater_g(tile),depth_g(tile),imax)
    call mlodiag("mixdepth",mixdepth_o,0,0,dgwater_g(tile),depth_g(tile),imax)
    rbot = real(mixind_o)
    call mlodiag("mixind",rbot,0,0,dgwater_g(tile),depth_g(tile),imax)
    mixind_o = nint(rbot)
    do k = 1,ol
      call mlodiag("rad",rad_o(:,k),k,0,dgwater_g(tile),depth_g(tile),imax)
    end do    
    do k = 1,ol
      call mloexpturb("tke",w_tke(:,k),k,0,turb_g(tile),depth_g(tile),imax)
      call mloexpturb("eps",w_eps(:,k),k,0,turb_g(tile),depth_g(tile),imax)
    end do
    call mloexpice("u",i_u,0,ice_g(tile),depth_g(tile),imax)
    call mloexpice("v",i_v,0,ice_g(tile),depth_g(tile),imax)
    call mlodiagice("fg",icefg_a,0,dgice_g(tile),depth_g(tile),imax)
    call mlodiagice("mass",imass,0,dgice_g(tile),depth_g(tile),imax)
    call mlodiagice("cd",cd_ice,0,dgice_g(tile),depth_g(tile),imax)
    call mlodiagice("cd_bot",cdbot_ice,0,dgice_g(tile),depth_g(tile),imax)    
  end if
   
  call tkemix(rkm,rhs,lqg,lqlg,lqfg,lni,lstratcloud,lu,lv,pblh(js:je),           &
              fg(js:je),eg(js:je),cduv(js:je),ps(js:je),zg,zh,sig,               &
              rhos,ustar(js:je),dt,nlocal_mode,ltke,leps,lshear,dx,land(js:je),  &
              w_t,w_s,w_u,w_v,w_tke,w_eps,deptho_depth,deptho_depth_hl,          &
              deptho_dz,deptho_dz_hl,dwdx_o,dwdy_o,ibot,ubot_o,vbot_o,           &
              utop_o,vtop_o,cd_water,cdh_water,cdbot_water,wt0_o,wt0rad_o,       &
              wt0melt_o,wt0eg_o,wt0fb_o,ws0_o,ws0subsurf_o,wu0_o,wv0_o,zo_o,     &
              bf_o,mixdepth_o,mixind_o,rad_o,                                    &
              i_u,i_v,fracice(js:je),icefg_a,imass,cd_ice,cdbot_ice,             &
              imax,kl,ol)  
  
  qg(js:je,:)         = lqg
  qfg(js:je,:)        = lqfg
  qlg(js:je,:)        = lqlg
  ni(js:je,:)         = lni
  stratcloud(js:je,:) = lstratcloud

  tke(js:je,:)        = ltke
  eps(js:je,:)        = leps

  u(js:je,:)          = lu
  v(js:je,:)          = lv
  
  do k = 1,kl
    ! replace counter gradient term with cqmix
    ! save Km and Kh for output
    rkmsave(js:je,k) = rkm(:,k)*cqmix
    rkhsave(js:je,k) = rkm(:,k)*cqmix
    ! transform winds back to Earth reference frame and theta to temp
    t(js:je,k) = rhs(:,k)/sigkap(k)
  enddo    !  k loop
  
  lrho(1:imax) = ps(js:je)/(rdry*tss(js:je))
  taux(js:je)  = lrho(1:imax)*cduv(js:je)*u(js:je,1)
  tauy(js:je)  = lrho(1:imax)*cduv(js:je)*v(js:je,1)
  
  if ( nlocal_mode==2 .or. nlocal_mode==3 ) then
    ! import ocean data after mixing  
    do k = 1,ol
      call mloimport("temp",w_t(:,k),k,0,water_g(tile),depth_g(tile),imax)
      call mloimport("sal",w_s(:,k),k,0,water_g(tile),depth_g(tile),imax)
      call mloimport("u",w_u(:,k),k,0,water_g(tile),depth_g(tile),imax)
      call mloimport("v",w_v(:,k),k,0,water_g(tile),depth_g(tile),imax)
    end do
    call mloimport("ubot",ubot_o,0,0,water_g(tile),depth_g(tile),imax)
    call mloimport("vbot",vbot_o,0,0,water_g(tile),depth_g(tile),imax)
    call mloimport("utop",utop_o,0,0,water_g(tile),depth_g(tile),imax)
    call mloimport("vtop",vtop_o,0,0,water_g(tile),depth_g(tile),imax)
    call mlo_updatediag("wu0",wu0_o,0,dgwater_g(tile),depth_g(tile),imax)
    call mlo_updatediag("wv0",wv0_o,0,dgwater_g(tile),depth_g(tile),imax)
    call mlo_updatediag("wt0",wt0_o,0,dgwater_g(tile),depth_g(tile),imax)    
    do k = 1,ol
      call mloimpturb("tke",w_tke(:,k),k,0,turb_g(tile),depth_g(tile),imax)
      call mloimpturb("eps",w_eps(:,k),k,0,turb_g(tile),depth_g(tile),imax)
    end do  
    call mloimpice("u",i_u,0,ice_g(tile),depth_g(tile),imax)
    call mloimpice("v",i_v,0,ice_g(tile),depth_g(tile),imax)    
    call mloimpice("fracice",fracice(js:je),0,ice_g(tile),depth_g(tile),imax)        
    ! export ocean SST
    call mlosurf("sst",tss(js:je),0,water_g(tile),ice_g(tile),depth_g(tile),imax)
  end if
  
end do ! tile = 1,ntiles
#endif

return
end subroutine tkeeps_work

end module module_ctrl_turbmix