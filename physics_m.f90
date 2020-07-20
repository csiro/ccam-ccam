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
module physics_m

implicit none

private
public physics

contains

subroutine physics

use aerosolldr                             ! LDR prognostic aerosols
use arrays_m                               ! Atmosphere dyamics prognostic arrays
use cc_mpi                                 ! CC MPI routines
use cc_omp                                 ! CC OpenMP routines
use cfrac_m                                ! Cloud fraction
use convjlm_m                              ! Convection
use convjlm22_m , alfin22 => alfin,      & ! Convection v2
                  timeconv22 => timeconv
use extraout_m                             ! Additional diagnostics
use gdrag_m                                ! Gravity wave drag
use histave_m                              ! Time average arrays
use kuocomb_m                              ! JLM convection
use liqwpar_m                              ! Cloud water mixing ratios
use map_m                                  ! Grid map arrays
use morepbl_m                              ! Additional boundary layer diagnostics
use newmpar_m                              ! Grid parameters
use nharrs_m                               ! Non-hydrostatic atmosphere arrays
use parm_m                                 ! Model configuration
use pbl_m                                  ! Boundary layer arrays
use prec_m                                 ! Precipitation
use sigs_m                                 ! Atmosphere sigma levels
use soil_m                                 ! Soil and surface data
use tracers_m                              ! Tracer data
use vvel_m                                 ! Additional vertical velocity
use work2_m                                ! Diagnostic arrays

implicit none

include 'kuocom.h'   ! kbsav,ktsav,convfact,convpsav,ndavconv

integer tile, is, ie
integer idjd_t
integer, dimension(imax)          :: lkbsav, lktsav
real, dimension(imax,kl,naero)    :: lxtg
real, dimension(imax,kl,ntrac)    :: ltr
real, dimension(imax,kl)          :: ldpsldt, lt, lqg
real, dimension(imax,kl)          :: lfluxtot
real, dimension(imax,kl)          :: lqlg, lu, lv, lqfg
real, dimension(imax,kl)          :: lcfrac
real, dimension(imax,ndust)       :: ldustwd
real, dimension(imax)             :: lso2wd, lso4wd
real, dimension(imax)             :: lbcwd, locwd, lsaltwd
real, dimension(imax,kl)          :: lqrg, lqsng, lqgrg, lstratcloud
logical mydiag_t

!$omp do schedule(static) private(is,ie),                                           &
!$omp private(lt,lu,lv,idjd_t,mydiag_t),                                            &
!$omp private(ldpsldt,lqg,lfluxtot),                                                &
!$omp private(lxtg,lso2wd,lso4wd,lbcwd,locwd,ldustwd,lsaltwd),                      &
!$omp private(lqlg,lqfg,lcfrac,ltr),                                                &
!$omp private(lqrg,lqsng,lqgrg,lstratcloud)
!$acc parallel copy(u,v) copy(eg,fg,cdtq,cduv,xtosav)                               &
!$acc  copy(t,qg,qlg,qfg,xtg,dustwd,so2wd,so4wd,bcwd)                               &
!$acc  copy(ocwd,saltwd,tr,precc,precip,timeconv,kbsav,ktsav,cfrac)                 &
!$acc  copy(cape,condc,condx,conds,condg,fluxtot,timeconv22)                        &
!$acc  copy(qrg,qsng,qgrg,stratcloud,kt_saved,kb_saved,aug,convpsav)                &
!$acc  copyin(tss,he)                                                               &
!$acc  copyin(dpsldt,alfin,alfin22,ps,pblh,wetfac,land,entrainn)                    &
!$acc  copyin(em,sgsave)
!$acc loop gang                                                                     &
!$acc  private(ldpsldt,lt,lqg,lqlg,lqfg,lcfrac,lu,lv)                               &
!$acc  private(lxtg,ldustwd,lso2wd,lso4wd,lbcwd,locwd,lsaltwd)                      &
!$acc  private(ltr,lfluxtot,lqrg,lqsng,lqgrg,lstratcloud)
do tile = 1,ntiles
  is = (tile-1)*imax + 1
  ie = tile*imax
  
  idjd_t = mod(idjd-1,imax) + 1
  mydiag_t = ((idjd-1)/imax==tile-1).and.mydiag
  
  lt = t(is:ie,:)
  lu = u(is:ie,:)
  lv = v(is:ie,:)

  ldpsldt   = dpsldt(is:ie,:)
  lqg       = qg(is:ie,:)
  lqlg      = qlg(is:ie,:)
  lqfg      = qfg(is:ie,:)
  lcfrac    = cfrac(is:ie,:)
  if ( abs(iaero)>=2 ) then
    lxtg    = xtg(is:ie,:,:)
    ldustwd = dustwd(is:ie,:)
    lso2wd  = so2wd(is:ie)
    lso4wd  = so4wd(is:ie)
    lbcwd   = bcwd(is:ie)
    locwd   = ocwd(is:ie)
    lsaltwd = saltwd(is:ie)
  end if
  if ( ngas>0 ) then
    ltr = tr(is:ie,:,:)
  end if

  lqrg        = qrg(is:ie,:)
  lqsng       = qsng(is:ie,:)
  lqgrg       = qgrg(is:ie,:)
  lstratcloud = stratcloud(is:ie,:)

  ! MISC (PARALLEL) -------------------------------------------------------
  ! initialse surface rainfall to zero
  condc(is:ie) = 0. ! default convective rainfall (assumed to be rain)
  condx(is:ie) = 0. ! default total precip = rain + ice + snow + graupel (convection and large scale)
  conds(is:ie) = 0. ! default total ice + snow (convection and large scale)
  condg(is:ie) = 0. ! default total graupel (convection and large scale)
  ! Held & Suarez or no surf fluxes
  if ( ntsur<=1 .or. nhstest==2 ) then 
    eg(is:ie)   = 0.
    fg(is:ie)   = 0.
    cdtq(is:ie) = 0.
    cduv(is:ie) = 0.
  end if     ! (ntsur<=1.or.nhstest==2) 
  ! Save aerosol concentrations for outside convective fraction of grid box
  if ( abs(iaero)>=2 ) then
    xtosav(is:ie,1:kl,1:naero) = xtg(is:ie,1:kl,1:naero) ! Aerosol mixing ratio outside convective cloud
  end if
#ifndef GPU
  is = (tile-1)*imax + 1
  ie = tile*imax
  call nantest("start of physics",is,ie)
#endif


  ! GWDRAG ----------------------------------------------------------------
#ifndef GPU
  call START_LOG(gwdrag_begin)
#endif
  if ( ngwd<0 ) then
    call gwdrag(lt,lu,lv,tss(is:ie),he(is:ie),idjd_t,mydiag_t,imax,kl) ! <0 for split - only one now allowed
  end if
#ifndef GPU
  call nantest("after gravity wave drag",is,ie)
  call END_LOG(gwdrag_end)
#endif


  ! CONVECTION ------------------------------------------------------------
#ifndef GPU
  call START_LOG(convection_begin)
#endif
  convh_ave(is:ie,1:kl) = convh_ave(is:ie,1:kl) - lt(:,1:kl)*real(nperday)/real(nperavg)        
  ! Select convection scheme
  select case ( nkuo )
    !case(5)
    !  call betts(t,qg,tn,land,ps) ! not called these days
    case(21,22)
      call convjlm22(alfin22(is:ie),ldpsldt,lt,lqg,                  &
           ps(is:ie),lfluxtot,convpsav(is:ie),cape(is:ie),           &
           lxtg,lso2wd,lso4wd,lbcwd,locwd,ldustwd,lsaltwd,           &
           lqlg,condc(is:ie),precc(is:ie),condx(is:ie),conds(is:ie), &
           condg(is:ie),precip(is:ie),                               &
           pblh(is:ie),fg(is:ie),wetfac(is:ie),land(is:ie),lu,lv,    &
           timeconv22(is:ie),em(is:ie),                              &
           kbsav(is:ie),ktsav(is:ie),ltr,lqfg,lcfrac,sgsave(is:ie),  &
           kt_saved(is:ie),kb_saved(is:ie),aug(is:ie),               &
           idjd_t,mydiag_t,entrain,detrain,mbase,iterconv,           &
           nuvconv,alfsea,methdetr,methprec,fldown,alflnd,rhcv,      &
           convtime,nkuo,rhsat,nevapls,                              &
           tied_con,mdelay,convfact,ncvcloud,ldr,rhmois,imax,kl)     ! split convjlm
    case(23,24)
      call convjlm(alfin(is:ie),ldpsldt,lt,lqg,                     &
           ps(is:ie),lfluxtot,convpsav(is:ie),cape(is:ie),          &
           lxtg,lso2wd,lso4wd,lbcwd,locwd,ldustwd,lsaltwd,          &
           lqlg,condc(is:ie),precc(is:ie),condx(is:ie),             &
           conds(is:ie),condg(is:ie),precip(is:ie),                 &
           pblh(is:ie),fg(is:ie),wetfac(is:ie),land(is:ie),         &
           entrainn(is:ie),lu,lv,timeconv(is:ie),em(is:ie),         &
           kbsav(is:ie),ktsav(is:ie),ltr,lqfg,lcfrac,sgsave(is:ie), &
           idjd_t,mydiag_t,entrain,detrain,mbase,nbase,iterconv,    &
           nuvconv,alfsea,methdetr,methprec,fldown,alflnd,detrainx, &
           sigkscb,dsig2,sigksct,rhcv,sig_ct,convtime,tied_con,     &
           mdelay,nevapcc,convfact,ncvcloud,ldr,rhmois,imax,kl)     ! split convjlm
  end select

  call fixqg(lt,lqg,lqlg,lqfg,lqrg,lqsng,lqgrg,lstratcloud,imax,kl)

  t(is:ie,:)       = lt
  qg(is:ie,:)      = lqg
  qlg(is:ie,:)     = lqlg
  qfg(is:ie,:)     = lqfg
  fluxtot(is:ie,:) = lfluxtot
  u(is:ie,:)       = lu
  v(is:ie,:)       = lv
  cfrac(is:ie,:)   = lcfrac
  if ( abs(iaero)>=2 ) then
    xtg(is:ie,:,:)  = lxtg
    dustwd(is:ie,:) = ldustwd
    so2wd(is:ie) = lso2wd
    so4wd(is:ie) = lso4wd
    bcwd(is:ie) = lbcwd
    ocwd(is:ie) = locwd
    saltwd(is:ie) = lsaltwd
  end if
  if ( ngas>0 ) then
    tr(is:ie,:,:) = ltr
  end if
       
  qrg(is:ie,:)        = lqrg
  qsng(is:ie,:)       = lqsng
  qgrg(is:ie,:)       = lqgrg
  stratcloud(is:ie,:) = lstratcloud

#ifndef GPU
  call nantest("after convection",is,ie)
  call END_LOG(convection_end)
#endif
end do
!$acc end parallel
!$omp end do

return
end subroutine physics

subroutine fixqg(t,qg,qlg,qfg,qrg,qsng,qgrg,stratcloud,imax,kl)
!$acc routine vector

use const_phys, only : hlcp, hlscp                   ! Physical constants
use parm_m, only : qg_fix                            ! Model configuration

implicit none

integer k
integer, intent(in) :: imax, kl
real, dimension(imax,kl), intent(inout) :: t
real, dimension(imax,kl), intent(inout) :: qg
real, dimension(imax,kl), intent(inout) :: qlg
real, dimension(imax,kl), intent(inout) :: qfg
real, dimension(imax,kl), intent(inout) :: qrg
real, dimension(imax,kl), intent(inout) :: qsng
real, dimension(imax,kl), intent(inout) :: qgrg
real, dimension(imax,kl), intent(inout) :: stratcloud

real, dimension(1:imax) :: qtot, tliq

! requires qg_fix>=1
if ( qg_fix<=0 ) return

do k = 1,kl
  qtot(:) = max( qg(:,k) + qlg(:,k) + qfg(:,k), 0. )
  tliq(:) = t(:,k) - hlcp*qlg(:,k) - hlscp*qfg(:,k)
  
  qfg(:,k)   = max( qfg(:,k), 0. ) 
  qlg(:,k)   = max( qlg(:,k), 0. )
  qrg(:,k)   = max( qrg(:,k), 0. )
  qsng(:,k)  = max( qsng(:,k), 0. )
  qgrg(:,k)  = max( qgrg(:,k), 0. )
  
  qg(:,k) = max( qtot(:) - qlg(:,k) - qfg(:,k), 0. )
  t(:,k)  = tliq(:) + hlcp*qlg(:,k) + hlscp*qfg(:,k)
  where ( qlg(:,k)+qfg(:,k)>1.E-8 )
    stratcloud(:,k) = max( stratcloud(:,k), 1.E-8 )
  elsewhere
    stratcloud(:,k) = 0.  
  end where
end do

return
end subroutine fixqg
    
end module physics_m
