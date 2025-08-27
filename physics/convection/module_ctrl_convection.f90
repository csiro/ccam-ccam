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
    
module module_ctrl_convection

implicit none

private
public ctrl_convection, ctrl_convection_init

contains
    
subroutine ctrl_convection

use arrays_m                      ! Atmosphere dyamics prognostic arrays
use cc_mpi                        ! CC MPI routines
use convjlm_m                     ! Convection
use convjlm22_m                   ! Convection v2
use kuocom_m                      ! JLM convection
use nlin_m                        ! Atmosphere non-linear dynamics
use soil_m                        ! Soil and surface data

implicit none

select case ( interp_convection(nkuo) )
  case("betts_conv")
    call betts(t,qg,tn,land,ps) ! not called these days
  case("john_conv22")
    call convjlm22              ! split convjlm
  case("john_conv")
    call convjlm                ! split convjlm
  case("grell_conv")
    call grell_ccam             ! grell convection    
  case("disable")
    ! do nothing
  case default
    write(6,*) "ERROR: unknown convection option nkuo=",nkuo
    call ccmpi_abort(-1) 
end select

return
end subroutine ctrl_convection

subroutine grell_ccam

! specify module to use global variables
use aerointerface                 ! Aerosol interface
use aerosol_arrays                ! Aerosol arrays
use arrays_m                      ! Atmosphere dyamics prognostic arrays
use const_phys                    ! Physical constants
use cu_gf_deep                    ! Grell convection
use kuocom_m                      ! JLM convection
use liqwpar_m                     ! Cloud water mixing ratios
use map_m                         ! Grid map arrays
use morepbl_m                     ! Additional boundary layer diagnostics
use newmpar_m                     ! Grid parameters
use parm_m                        ! Model configuration
use prec_m                        ! Precipitation
use sigs_m                        ! Atmosphere sigma levels
use soil_m                        ! Soil and surface data
use vvel_m                        ! Additional vertical velocity

implicit none

! declare all the variables being use here, remove those are not used
! first those in Grell
integer :: tile, js, je, k, n, iq
integer :: njumps, m, idjd_t
real, dimension(imax,kl) :: lqg
real, dimension(imax,kl) :: rhoa

integer :: itf,ktf,its,ite, kts,kte, kdt   ! check itf and ktf
integer :: dicycle                         ! flag, set = 1 will turn on if
integer :: ichoice                         ! flag, usually set to 0
integer :: nchem
integer :: ipr                             ! flag to turn on debug ...CHECK

real, dimension(imax)    :: ccn, mconv
real                     :: ccnclean, dtime
integer                  :: imid                 ! flag to turn on mid level convection
integer, dimension(imax) :: kpbl, tropics        ! level of boundary layer height
real, dimension(imax,kl) :: dhdt,g_rho,g_t,g_q,po,g_us,g_vs,tn
real, dimension(imax)    :: dx,z1,psur,xland,dq
real, dimension(imax,kl) :: zo
real, dimension(imax,10) :: forcing
real, dimension(imax,kl) :: q,qo,zuo,zdo,zdm, g_qfg, g_qlg
real, dimension(imax)    :: hfx,qfx,xmbm_in,xmbs_in, xmb_out,g_pre
real, dimension(imax)    :: pre_mid, pre_deep
real, dimension(imax,kl) :: omeg
integer, dimension(imax) :: csum
real, dimension(imax,kl) :: cnvwt,cupclw
real, dimension(imax,kl) :: outu_mid,outv_mid,outt_mid,outq_mid,outqc_mid,outliqice_mid
real, dimension(imax,kl) :: outu_deep,outv_deep,outt_deep,outq_deep,outqc_deep,outliqice_deep
real, dimension(imax)    :: edto,edtm,hkbo,xhkb,xmb,pwavo,ccnloss,             &
                            pwevo,bu,bud,cap_max,                              &
                            cap_max_increment,closure_n,psum,psumh,sigd
integer, dimension(imax) :: kbcon,ktop
integer, dimension(imax) :: ktop_mid, ktop_deep
real, dimension(imax)    :: frh_out
integer, dimension(imax) :: ierr
character(len=50), dimension(imax) :: ierrc                ! CHECK THIS
real, dimension(imax,kl,naero) :: chem3d                   ! CHECK THIS
real, dimension(imax,naero)    :: wetdpc_deep              ! CHECK THIS
real, dimension(imax,naero)    :: wetdpc_mid
logical                  :: do_smoke_transport             ! CHECK THIS
real, dimension(imax)    :: rand_mom,rand_vmas
real, dimension(imax,4)  :: rand_clos
integer                  :: nranflag
integer                  :: do_capsuppress
integer, dimension(imax) :: k22,jmin
! declare vars locally
real, dimension(imax,kl) :: thz, tothz
real prf_temp, prf
real, dimension(imax,kl) :: zpres, zcdrop
integer                  :: kk, i
real, dimension(naero)   :: fscav
real, dimension(imax)    :: cap_suppress_j
integer, dimension(2) :: posmin, posmax
integer, dimension(3) :: posmin3, posmax3

real :: maxconvtime = 120.  ! time-step for convection
real :: tdt

do tile = 1,ntiles
  js = (tile-1)*imax+1      ! js:je inside 1:ifull
  je = tile*imax            ! len(js:je) = imax

  qamin = qgmin
  
  ! working in ifull (js:je), need to unpack/devide to imax, e.g., t1(1:imax,:)=t(js:je,:)

  do k = 1,kl
    do iq = 1,imax
      prf_temp      = ps(iq+js-1)*sig(k)
      prf           = 0.01*prf_temp           ! ps is SI units
      tothz(iq,k)   = (prf/1000.)**(rdry/cp)
      thz(iq,k)     = t(iq+js-1,k)/tothz(iq,k)
      rhoa(iq,k)    = prf_temp/(rdry*t(iq+js-1,k))
      zpres(iq,k)   = prf_temp
      !dzw(iq,k)     = dz(iq+js-1,k)
    end do
  end do

  ! add input for convection for each imax section of the tile
  ! RHS: js to je, from ifull ! LHS: a section of imax (1:imax) --> goes to cumulus subroutine
  itf       = imax            ! not sure why sometime it goes from its --> itf
  ktf       = kl-1            ! MJT suggestion to avoid bug on line 1874 of cu_gf_deep.F90
  its       = 1               ! this for its:ite in dims
  ite       = imax
  kts       = 1               ! this for kts:kte in dims
  kte       = kl
  dicycle   = 0               ! diurnal cycle flag                                  ! affect xmb calculations
  ichoice   = 0               ! choice of closure, use "0" for ensemble average     ! affect xmb_ave calculations
  ipr       = 0               ! this flag can be used for debugging prints
  ccn       = 150.            ! not well tested yet
  ccnclean  = 0.
  dtime     = dt              ! dt over which forcing is applied
  do k = 1,kl
    ! boundary layer forcing (one closure for shallow)      
    dhdt(1:imax,k) = dmsedt_adv(js:je,k) + dmsedt_rad(js:je,k) + dmsedt_pbl(js:je,k)    
  end do
  
  where ( land(js:je) )       ! land mask
    xland(1:imax) = 1.
  elsewhere
    xland(1:imax) = 0.
  end where
  po(1:imax,:) = zpres(1:imax,:)*0.01   ! pressure mb
  zo(1:imax,1) = bet(1)*t(js:je,1)/grav ! heights above surface
  do k = 2,kl
    zo(1:imax,k) = zo(1:imax,k-1) + (bet(k)*t(js:je,k)+betm(k)*t(js:je,k-1))/grav ! heights above surface
  end do
  kpbl(1:imax) = 1            ! default value
  do k = 1,kl
    where ( zo(1:imax,k)<=pblh(js:je) )
      kpbl(1:imax) = k
    end where
  end do
  forcing    = 0.             ! only diagnostic
  g_t        = t(js:je,:)     ! t before forcing   ! check whether abs t or potential (t/sigkap(k)) or take tothz above
  g_q        = qg(js:je,:)    ! q before forcing
  z1         = zs(js:je)/grav ! terrain                                             ! elevation
  psur       = ps(js:je)*0.01 ! surface pressure (mb)
  g_us       = u(js:je,:)     ! u on mass points
  g_vs       = v(js:je,:)     ! v on mass points
  g_rho      = rhoa(1:imax,:)  ! density
  hfx        = fg(js:je)      ! w/m2, positive upward
  qfx        = eg(js:je)      ! w/m2, positive upward
  dx         = ds/em(js:je)   ! dx is grid point dependent here     ! CHECK IF THIS KM OR NOT, ds/em(js:je) IS IN METER
  do k = 1,kl                 
    omeg(1:imax,k) =  ps(js:je)*dpsldt(js:je,k)         ! omega (pa/s)
  end do
  mconv(:) = 0.                ! integrated vertical advection of moisture
  do k = 1,kl-1
    dq(1:imax) = qg(js:je,k+1) - qg(js:je,k)
    mconv(1:imax) = mconv(1:imax) + omeg(1:imax,k)*dq(1:imax)/grav
  end do
  mconv(1:imax) = max( mconv(1:imax), 0. )
  csum       = 0.              ! used to implement memory, set to zero if not avail
  cnvwt      = 0.              ! gfs needs this
  zuo        = 0.              ! nomalized updraft mass flux
  zdo        = 0.              ! nomalized downdraft mass flux
  zdm        = 0.              ! nomalized downdraft mass flux from mid scheme
  edto       = 0.              !
  edtm       = 0.              !
  xmb_out    = 0.              ! the xmb's may be needed for dicycle
  xmbm_in    = 0.              !
  xmbs_in    = 0.              !
  kbcon      = 0.              ! lfc of parcel from k22
  ktop       = 0.              ! cloud top
  cupclw     = 0.              ! used for direct coupling to radiation, but with tuning factors
  frh_out    = 0.              ! fractional coverage
  if ( abs(iaero)>=2 ) then
    nchem      = naero
    fscav(1:naero) = scav_effl(1:naero)
    chem3d(1:imax,1:kl,1:naero) = xtg(js:je,1:kl,1:naero)
    wetdpc_mid(1:imax,1:naero)= 0.
    wetdpc_deep(1:imax,1:naero)= 0.
    do_smoke_transport = .false.
  else
    nchem       = 0
    fscav       = 0.
    chem3d      = 0.
    wetdpc_deep = 0.
    do_smoke_transport = .false.
  end if
  rand_mom            = 0     ! for stochastics mom, if temporal and spatial patterns exist
  rand_vmas           = 0     ! for stochastics vertmass, if temporal and spatial patterns exist
  rand_clos           = 0     ! for stochastics closures, if temporal and spatial patterns exist
  nranflag            = 0     ! flag to what you want perturbed
  do_capsuppress      = 0
  cap_suppress_j      = 0
  k22                 = 0
  jmin                = 0
  kdt                 = 0
  tropics             = 0
  g_pre(1:imax) = 0.
  g_qfg(1:imax,1:kl) = qfg(js:je,1:kl) ! MJT suggestion
  g_qlg(1:imax,1:kl) = qlg(js:je,1:kl)
  
  ! Use sub time-step if required
  njumps = int(dtime/(maxconvtime+0.01)) + 1
  tdt    = real(dtime/real(njumps))
  do n = 1,njumps
      
    ierr       = 0                ! ierr flags are error flags, used for debugging
    ierrc      = ' '              ! the following should be set to zero if not available
      
    tn = g_t
    qo = g_q

    pre_mid   = 0. ! mm/s?
    outt_mid  = 0.
    outq_mid  = 0.
    outqc_mid = 0.
    outu_mid  = 0.
    outv_mid  = 0.
    outliqice_mid = 0.

    pre_deep   = 0. ! mm/s?
    outt_deep  = 0.
    outq_deep  = 0.
    outqc_deep = 0.
    outu_deep  = 0.
    outv_deep  = 0.
    outliqice_deep = 0.

    !! imid=1
    !call cu_gf_deep_run(   &
    !         itf           &
    !        ,ktf           &
    !        ,its           &
    !       ,ite           &
    !        ,kts           &
    !        ,kte           &
    !        ,dicycle       &  ! diurnal cycle flag
    !        ,ichoice       &  ! choice of closure, use "0" for ensemble average
    !        ,ipr           &  ! this flag can be used for debugging prints
    !        ,ccn           &  ! not well tested yet
    !        ,ccnclean      &
    !        ,tdt           &  ! dt over which forcing is applied
    !        ,1             &  ! flag to turn on mid level convection
    !        ,kpbl          &  ! level of boundary layer height
    !        ,dhdt          &  ! boundary layer forcing (one closure for shallow)
    !        ,xland         &  ! land mask
    !        ,zo            &  ! heights above surface
    !        ,forcing       &  ! only diagnostic
    !        ,g_t           &  ! t before forcing
    !        ,g_q           &  ! q before forcing
    !        ,z1            &  ! terrain
    !        ,tn            &  ! t including forcing
    !        ,qo            &  ! q including forcing
    !        ,po            &  ! pressure (mb)
    !        ,psur          &  ! surface pressure (mb)
    !        ,g_us          &  ! u on mass points
    !        ,g_vs          &  ! v on mass points
    !        ,g_rho         &  ! density
    !        ,hfx           &  ! w/m2, positive upward
    !        ,qfx           &  ! w/m2, positive upward
    !        ,dx            &  ! dx is grid point dependent here
    !        ,mconv         &  ! integrated vertical advection of moisture
    !        ,omeg          &  ! omega (pa/s)
    !        ,csum          &  ! used to implement memory, set to zero if not avail
    !        ,cnvwt         &  ! gfs needs this
    !        ,zuo           &  ! nomalized updraft mass flux
    !        ,zdo           &  ! nomalized downdraft mass flux
    !        ,zdm           &  ! nomalized downdraft mass flux from mid scheme
    !        ,edto          &  !
    !        ,edtm          &  !
    !        ,xmb_out       &  ! the xmb's may be needed for dicycle
    !        ,xmbm_in       &  !
    !        ,xmbs_in       &  !
    !        ,pre_mid       &  !
    !        ,outu_mid      &  ! momentum tendencies at mass points
    !        ,outv_mid      &  !
    !        ,outt_mid      &  ! temperature tendencies
    !        ,outq_mid      &  ! q tendencies
    !        ,outqc_mid     &  ! ql/qice tendencies
    !        ,kbcon         &  ! lfc of parcel from k22
    !        ,ktop          &  ! cloud top
    !        ,cupclw        &  ! used for direct coupling to radiation, but with tuning factors
    !        ,frh_out       &  ! fractional coverage
    !        ,ierr          &  ! ierr flags are error flags, used for debugging
    !        ,ierrc         &  ! the following should be set to zero if not available
    !        ,nchem         &
    !        ,fscav         &
    !        ,chem3d        &
    !        ,wetdpc_mid    &
    !        ,do_smoke_transport   &
    !        ,rand_mom      &  ! for stochastics mom, if temporal and spatial patterns exist
    !        ,rand_vmas     &  ! for stochastics vertmass, if temporal and spatial patterns exist
    !        ,rand_clos     &  ! for stochastics closures, if temporal and spatial patterns exist
    !        ,nranflag      &  ! flag to what you want perturbed
    !                          !! 1 = momentum transport
    !                          !! 2 = normalized vertical mass flux profile
    !                          !! 3 = closures
    !                          !! more is possible, talk to developer or
    !                          !! implement yourself. pattern is expected to be
    !                          !! betwee -1 and +1
    !        ,do_capsuppress,cap_suppress_j    &    !
    !        ,k22                              &    !
    !        ,jmin,kdt,tropics                 &
    !        ,outliqice_mid)
    !
    !ktop_mid = ktop
    !call neg_check('mid',1,tdt,g_q,outq_mid,outt_mid,outu_mid,outv_mid,      &
    !                    outqc_mid,pre_mid,its,ite,kts,kte,itf,ktf,ktop_mid)
    
    ! imid=0
    call cu_gf_deep_run(   &
             itf           &
            ,ktf           &
            ,its           &
            ,ite           &
            ,kts           &
            ,kte           &
            ,dicycle       &  ! diurnal cycle flag
            ,ichoice       &  ! choice of closure, use "0" for ensemble average
            ,ipr           &  ! this flag can be used for debugging prints
            ,ccn           &  ! not well tested yet
            ,ccnclean      &
            ,tdt           &  ! dt over which forcing is applied
            ,0             &  ! flag to turn on mid level convection
            ,kpbl          &  ! level of boundary layer height
            ,dhdt          &  ! boundary layer forcing (one closure for shallow)
            ,xland         &  ! land mask
            ,zo            &  ! heights above surface
            ,forcing       &  ! only diagnostic
            ,g_t           &  ! t before forcing
            ,g_q           &  ! q before forcing
            ,z1            &  ! terrain
            ,tn            &  ! t including forcing
            ,qo            &  ! q including forcing
            ,po            &  ! pressure (mb)
            ,psur          &  ! surface pressure (mb)
            ,g_us          &  ! u on mass points
            ,g_vs          &  ! v on mass points
            ,g_rho         &  ! density
            ,hfx           &  ! w/m2, positive upward
            ,qfx           &  ! w/m2, positive upward
            ,dx            &  ! dx is grid point dependent here
            ,mconv         &  ! integrated vertical advection of moisture
            ,omeg          &  ! omega (pa/s)
            ,csum          &  ! used to implement memory, set to zero if not avail
            ,cnvwt         &  ! gfs needs this
            ,zuo           &  ! nomalized updraft mass flux
            ,zdo           &  ! nomalized downdraft mass flux
            ,zdm           &  ! nomalized downdraft mass flux from mid scheme
            ,edto          &  !
            ,edtm          &  !
            ,xmb_out       &  ! the xmb's may be needed for dicycle
            ,xmbm_in       &  !
            ,xmbs_in       &  !
            ,pre_deep      &  !
            ,outu_deep     &  ! momentum tendencies at mass points
            ,outv_deep     &  !
            ,outt_deep     &  ! temperature tendencies
            ,outq_deep     &  ! q tendencies
            ,outqc_deep    &  ! ql/qice tendencies
            ,kbcon         &  ! lfc of parcel from k22
            ,ktop          &  ! cloud top
            ,cupclw        &  ! used for direct coupling to radiation, but with tuning factors
            ,frh_out       &  ! fractional coverage
            ,ierr          &  ! ierr flags are error flags, used for debugging
            ,ierrc         &  ! the following should be set to zero if not available
            ,nchem         &
            ,fscav         &
            ,chem3d        &
            ,wetdpc_deep   &
            ,do_smoke_transport   &
            ,rand_mom      &  ! for stochastics mom, if temporal and spatial patterns exist
            ,rand_vmas     &  ! for stochastics vertmass, if temporal and spatial patterns exist
            ,rand_clos     &  ! for stochastics closures, if temporal and spatial patterns exist
            ,nranflag      &  ! flag to what you want perturbed
                              !! 1 = momentum transport
                              !! 2 = normalized vertical mass flux profile
                              !! 3 = closures
                              !! more is possible, talk to developer or
                              !! implement yourself. pattern is expected to be
                              !! betwee -1 and +1
            ,do_capsuppress,cap_suppress_j    &    !
            ,k22                              &    !
            ,jmin,kdt,tropics                 &
            ,outliqice_deep)
    
    ktop_deep=ktop
    call neg_check('deep',1,tdt,g_q,outq_deep,outt_deep,outu_deep,outv_deep,   &
                        outqc_deep,pre_deep,its,ite,kts,kte,itf,ktf,ktop_deep)
    
    g_t(1:imax,:)   = g_t(1:imax,:)   + tdt*(outt_mid(1:imax,:)+outt_deep(1:imax,:))
    g_q(1:imax,:)   = g_q(1:imax,:)   + tdt*(outq_mid(1:imax,:)+outq_deep(1:imax,:))
    g_qlg(1:imax,:) = g_qlg(1:imax,:) + tdt*(outqc_mid(1:imax,:)*outliqice_mid(1:imax,:)           &
                                            +outqc_deep(1:imax,:)*outliqice_deep(1:imax,:))
    g_qfg(1:imax,:) = g_qfg(1:imax,:) + tdt*(outqc_mid(1:imax,:)*(1. - outliqice_mid(1:imax,:))    &
                                            +outqc_deep(1:imax,:)*(1. - outliqice_deep(1:imax,:)))
    g_us(1:imax,:)  = g_us(1:imax,:)  + tdt*(outu_mid(1:imax,:)+outu_deep(1:imax,:))
    g_vs(1:imax,:)  = g_vs(1:imax,:)  + tdt*(outv_mid(1:imax,:)+outv_deep(1:imax,:))
    g_pre(1:imax)   = g_pre(1:imax)   + tdt*(pre_mid(1:imax)+pre_deep(1:imax))

  end do ! smaller time step ( n=1,njumps )

  t(js:je,:)       = g_t(1:imax,:)
  qg(js:je,:)      = g_q(1:imax,:)
  fluxtot(js:je,:) = 0.                 !fluxtot(js:je,:)  !pre(1:imax,:)
  qlg(js:je,:)     = g_qlg(1:imax,:)
  qfg(js:je,:)     = g_qfg(1:imax,:)
  u(js:je,:)       = g_us(1:imax,:)
  v(js:je,:)       = g_vs(1:imax,:)

  condc(js:je) = g_pre(1:imax)          ! check unit,tendencies
  conds(js:je) = 0.
  condg(js:je) = 0.
  precc(js:je) = precc(js:je) + real(condc(js:je),8)
  condx(js:je) = condx(js:je) + condc(js:je)
        
  if ( abs(iaero)>=2 ) then
    xtg(js:je,1:kl,1:naero) = chem3d(1:imax,1:kl,1:naero)
    so2wd(js:je) = so2wd(js:je) + wetdpc_mid(1:imax,itracso2)  &
                                + wetdpc_deep(1:imax,itracso2)
    so4wd(js:je) = so4wd(js:je) + wetdpc_mid(1:imax,itracso4)  &
                                + wetdpc_deep(1:imax,itracso4)
    bcwd(js:je) = bcwd(js:je) + wetdpc_mid(1:imax,itracbc) + wetdpc_mid(1:imax,itracbc+1)   &
                              + wetdpc_deep(1:imax,itracbc) + wetdpc_deep(1:imax,itracbc+1)
    ocwd(js:je) = ocwd(js:je) + wetdpc_mid(1:imax,itracoc) + wetdpc_mid(1:imax,itracoc+1)   &
                              + wetdpc_deep(1:imax,itracoc) + wetdpc_deep(1:imax,itracoc+1)        
    do i = 1,ndust
      dustwd(js:je,i) = dustwd(js:je,i) + wetdpc_mid(1:imax,itracdu+i-1)  &
                                        + wetdpc_deep(1:imax,itracdu+i-1)
    end do
    saltwd(js:je) = saltwd(js:je) + wetdpc_mid(1:imax,itracsa) + wetdpc_mid(1:imax,itracsa+1)   &
                                  + wetdpc_deep(1:imax,itracsa) + wetdpc_deep(1:imax,itracsa+1)        
  end if

end do ! end tile loop (tile=1,ntiles)

return
end subroutine grell_ccam


! Initialise convection if required (called from indata)

subroutine ctrl_convection_init

use cc_mpi                        ! CC MPI routines
use convjlm_m                     ! Convection
use convjlm22_m                   ! Convection v2
use kuocom_m                      ! JLM convection

implicit none

select case ( interp_convection(nkuo) )
  case("betts_conv")
    ! do nothing
  case("john_conv22")
    call convjlm22_init
  case("john_conv")
    call convjlm_init
  case("grell_conv")
    ! do nothing    
  case("disable")
    ! do nothing
  case default
    write(6,*) "ERROR: Cannot initialise convection in indata with unknown option for nkuo = ",nkuo
    call ccmpi_abort(-1)
end select

return
end subroutine ctrl_convection_init


!====================================================================================================
! SUBROUTINE interp_nconvection
!   
! subroutine to select the cloud convection scheme for CCAM
!====================================================================================================
pure function interp_convection(nconvection) result(cv_physics)

implicit none

integer, intent(in) :: nconvection
character(len=20) :: cv_physics

cv_physics = "ERROR"

select case(nconvection)
  case(0)
    cv_physics = "disable"
  case(5)
    cv_physics = "betts_conv"
  case(21,22,29)
    cv_physics = "john_conv22"
  case(23,24)
    cv_physics = "john_conv"
  case(31)
    cv_physics = "grell_conv"
end select

return
end function interp_convection  
  
end module module_ctrl_convection
