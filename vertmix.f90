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
    
! CCAM boundary layer turbulent mixing routines

! Currently local Ri and prognostic k-e schemes are supported.
! Also, options for non-local counter gradient terms are included.
! Local Ri supports Gelyen and Tiedtke schemes for shallow
! convection, whereas the k-e follows the EDMF approach where
! shallow convection is represented in the mass flux terms.

! nvmix=3  Local Ri mixing
! nvmix=6  Prognostic k-e tubulence closure
! nvmix=7  Jing Huang's local Ri scheme (with axmlsq=9.)

! nlocal=0 No counter gradient term
! nlocal=6 Holtslag and Boville non-local term
! nlocal=7 Mass flux based counter gradient (requires nvmix=6)
      
! kscmom   0 off, 1 turns on shal_conv momentum (usual) (requires nvmix<6)
    
module  vertmix_m

use mlo, only : waterdata,icedata   ! Ocean physics and prognostic arrays

implicit none

private

public vertmix,vertmix_init

integer, save :: kscbase=-1, ksctop=-1

contains

subroutine vertmix_init

use cc_mpi                              ! CC MPI routines
use parm_m                              ! Model configuration
use sigs_m                              ! Atmosphere sigma levels

implicit none

include 'kuocom.h'                  ! Convection parameters

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

return
end subroutine vertmix_init

subroutine vertmix

use aerosolldr                      ! LDR prognostic aerosols
use arrays_m                        ! Atmosphere dyamics prognostic arrays
use carbpools_m                     ! Carbon pools
use cc_mpi                          ! CC MPI routines
use cc_omp                          ! CC OpenMP routines
use cfrac_m                         ! Cloud fraction
use diag_m                          ! Diagnostic routines
use extraout_m                      ! Additional diagnostics
use kuocomb_m                       ! JLM convection
use liqwpar_m                       ! Cloud water mixing ratios
use map_m                           ! Grid map arrays
use mlo                             ! Ocean physics and prognostic arrays
use morepbl_m                       ! Additional boundary layer diagnostics
use newmpar_m                       ! Grid parameters
use nharrs_m                        ! Non-hydrostatic atmosphere arrays
use nsibd_m                         ! Land-surface arrays
use parm_m, only : idjd, nmlo, iaero, nvmix
                                    ! Model configuration
use pbl_m                           ! Boundary layer arrays
use savuvt_m                        ! Saved dynamic arrays
use screen_m                        ! Screen level diagnostics
use soil_m, only : land             ! Soil and surface data
use soilsnow_m, only : fracice      ! Soil, snow and surface data
use tkeeps                          ! TKE-EPS boundary layer
#ifndef scm
use trvmix, only : tracervmix       ! Tracer mixing routines
use tracermodule                    ! Tracer routines
use tracers_m                       ! Tracer data
#endif
use work2_m                         ! Diagnostic arrays

implicit none

include 'kuocom.h'                  ! Convection parameters

integer :: is, ie, tile, k
integer :: idjd_t
real, dimension(imax,kl,naero) :: lxtg
real, dimension(ifull,kl) :: at_save, ct_save
real, dimension(imax,kl) :: lt, lqg, lqfg,  lqlg
real, dimension(imax,kl) :: lcfrac, lu, lv, lstratcloud
real, dimension(imax,kl) :: lsavu, lsavv, ltke, leps, lshear
real, dimension(imax,kl) :: lat, lct
real, dimension(ifull) :: uadj, vadj
real, dimension(imax) :: lou, lov, liu, liv
logical :: mydiag_t
#ifdef scm
real, dimension(imax,kl) :: lwth_flux, lwq_flux, luw_flux, lvw_flux
real, dimension(imax,kl) :: ltkesave, lepssave, lrkmsave, lrkhsave
real, dimension(imax,kl) :: lbuoyproduction, lshearproduction, ltotaltransport
real, dimension(imax,kl-1) :: lmfsave
#else
real, dimension(imax,numtracer) :: lco2em
real, dimension(imax,kl,ntrac) :: ltr
real, dimension(imax,kl) :: loh, lstrloss, ljmcf
#endif
   
if ( nmlo/=0 ) then
!$omp do schedule(static) private(is,ie,k),             &
!$omp private(lou,lov,liu,liv)
  do tile = 1,ntiles
    is = (tile-1)*imax + 1
    ie = tile*imax
    lou = 0.
    lov = 0.
    liu = 0.
    liv = 0.
    call mloexport(2,lou,1,0,water_g(tile),wpack_g(:,tile),wfull_g(tile))
    call mloexport(3,lov,1,0,water_g(tile),wpack_g(:,tile),wfull_g(tile))
    call mloexpice(liu, 9,0,ice_g(tile),wpack_g(:,tile),wfull_g(tile))
    call mloexpice(liv,10,0,ice_g(tile),wpack_g(:,tile),wfull_g(tile))
    uadj(is:ie) = (1.-fracice(is:ie))*lou + fracice(is:ie)*liu
    vadj(is:ie) = (1.-fracice(is:ie))*lov + fracice(is:ie)*liv
  end do
!$omp end do nowait
else
!$omp do schedule(static) private(is,ie,k)
  do tile = 1,ntiles
    is = (tile-1)*imax + 1
    ie = tile*imax
    uadj(is:ie) = 0.
    vadj(is:ie) = 0.
  end do
!$omp end do nowait  
end if

select case(nvmix)
  case(6)  
    ! k-e + MF closure scheme
    
!$omp do schedule(static) private(is,ie,k),             &
!$omp private(lt,lqg,lqfg,lqlg),                        &
!$omp private(lstratcloud,lxtg,lu,lv,ltke,leps,lshear), &
!$omp private(lat,lct,lsavu,lsavv,idjd_t,mydiag_t)
#ifdef scm
!!$acc parallel copy(t,qg,qlg,qfg,stratcloud,xtg,tke,eps,u,v, &
!!$acc   pblh,ustar)                                          &
!!$acc copyin(shear,uadj,vadj,em,tss,eg,fg,ps,cduv)           &
!!$acc copyout(at_save,ct_save,wth_flux,wq_flux,uw_flux,      &
!!$acc   vw_flux,mfsave,tkesave,epssave,rkmsave,rkhsave,      &
!!$acc   buoyproduction,shearproduction,totaltransport)       &
!!$acc loop gang private(lt,lqg,lqfg,lqlg,lstratcloud,lxtg,   &
!!$acc   ltke,leps,lshear,lat,lct,lu,lv,lwth_flux,lwq_flux,   &
!!$acc   luw_flux,lvw_flux,lmfsave,ltkesave,lepssave,         &
!!$acc   lrkmsave,lrkhsave,lbuoyproduction,lshearproduction,  &
!!$acc   ltotaltransport)
#else
!!$acc parallel copy(t,qg,qlg,qfg,stratcloud,xtg,tke,eps,u,v, &
!!$acc   pblh,ustar)                                          &
!!$acc copyin(shear,uadj,vadj,em,tss,eg,fg,ps,cduv)           &
!!$acc copyout(at_save,ct_save)                               &
!!$acc loop gang private(lt,lqg,lqfg,lqlg,lstratcloud,lxtg,   &
!!$acc   ltke,leps,lshear,lat,lct,lu,lv)
#endif
    do tile = 1,ntiles
      is = (tile-1)*imax + 1
      ie = tile*imax

      idjd_t = mod(idjd-1,imax)+1
      mydiag_t = ((idjd-1)/imax==tile-1).and.mydiag
  
      lt = t(is:ie,:)
      lqg = qg(is:ie,:)
      lqfg = qfg(is:ie,:)
      lqlg = qlg(is:ie,:)
      lstratcloud = stratcloud(is:ie,:)
      if ( abs(iaero)>=2 ) then
        lxtg = xtg(is:ie,:,:)
      end if
      ltke   = tke(is:ie,:)
      leps   = eps(is:ie,:)
      lshear = shear(is:ie,:)
      ! Adjustment for moving ocean surface
      do k = 1,kl
        lu(:,k) = u(is:ie,k) - uadj(is:ie)
        lv(:,k) = v(is:ie,k) - vadj(is:ie)
      end do  
    
      call tkeeps_work(lt,em(is:ie),tss(is:ie),eg(is:ie),fg(is:ie),                                      &
                       ps(is:ie),lqg,lqfg,lqlg,lstratcloud,lxtg,cduv(is:ie),lu,lv,pblh(is:ie),           &
                       ustar(is:ie),ltke,leps,lshear,lat,lct,                                            &
#ifdef scm
                       lwth_flux,lwq_flux,luw_flux,lvw_flux,lmfsave,ltkesave,lepssave,lrkmsave,lrkhsave, &
                       lbuoyproduction,lshearproduction,ltotaltransport,                                 &
#endif
                       imax,kl,naero)      
                       
      t(is:ie,:)          = lt
      qg(is:ie,:)         = lqg
      qfg(is:ie,:)        = lqfg
      qlg(is:ie,:)        = lqlg
      stratcloud(is:ie,:) = lstratcloud
      at_save(is:ie,:)    = lat
      ct_save(is:ie,:)    = lct
      if ( abs(iaero)>=2 ) then
        xtg(is:ie,:,:) = lxtg
      end if
      tke(is:ie,:) = ltke
      eps(is:ie,:) = leps
      do k = 1,kl  
        u(is:ie,k) = lu(:,k) + uadj(is:ie)
        v(is:ie,k) = lv(:,k) + vadj(is:ie)
      end do  
#ifdef scm
#ifdef _OPENMP
      write(6,*) "ERROR: scm requires OMP is disabled"
      stop
#endif
      rkmsave(is:ie,:) = lrkmsave
      rkhsave(is:ie,:) = lrkhsave  
      tkesave(is:ie,:) = ltkesave
      epssave(is:ie,:) = lepssave
      wth_flux(is:ie,:) = lwth_flux 
      wq_flux(is:ie,:) = lwq_flux 
      uw_flux(is:ie,:) = luw_flux 
      vw_flux(is:ie,:) = lvw_flux 
      mfsave(is:ie,:) = lmfsave 
      buoyproduction(is:ie,:) = lbuoyproduction
      shearproduction(is:ie,:) = lshearproduction
      totaltransport(is:ie,:) = ltotaltransport
#endif

    end do ! tile = 1,ntiles
!!$acc end parallel
!$omp end do nowait

  case default  
      ! JLM's local Ri scheme
    
!$omp do schedule(static) private(is,ie,k),               &
!$omp private(lt,lqg,lqfg,lqlg),                          &
!$omp private(lcfrac,lstratcloud,lxtg,lu,lv,lsavu,lsavv), &
!$omp private(lat,lct,idjd_t,mydiag_t)
    do tile = 1,ntiles
      is = (tile-1)*imax + 1
      ie = tile*imax

      idjd_t = mod(idjd-1,imax)+1
      mydiag_t = ((idjd-1)/imax==tile-1).and.mydiag
  
      lt = t(is:ie,:)
      lqg = qg(is:ie,:)
      lqfg = qfg(is:ie,:)
      lqlg = qlg(is:ie,:)
      lcfrac = cfrac(is:ie,:)
      lstratcloud = stratcloud(is:ie,:)
      if ( abs(iaero)>=2 ) then
        lxtg = xtg(is:ie,:,:)
      end if
      ! Adjustment for moving ocean surface
      do k = 1,kl
        lu(:,k) = u(is:ie,k) - uadj(is:ie)
        lv(:,k) = v(is:ie,k) - vadj(is:ie)
        lsavu(:,k) = savu(is:ie,k) - uadj(is:ie)
        lsavv(:,k) = savv(is:ie,k) - vadj(is:ie)
      end do  
    
      call vertmix_work(lt,tss(is:ie),eg(is:ie),fg(is:ie),kbsav(is:ie),ktsav(is:ie),convpsav(is:ie),               &
                        ps(is:ie),lqg,lqfg,lqlg,lstratcloud,                                                       &
                        condc(is:ie),lcfrac,lxtg,cduv(is:ie),lu,lv,pblh(is:ie),lsavu,lsavv,land(is:ie),            &
                        tscrn(is:ie),qgscrn(is:ie),ustar(is:ie),f(is:ie),condx(is:ie),zs(is:ie),                   &
                        lat,lct,                                                                                   &
#ifdef scm
                        lwth_flux,lwq_flux,luw_flux,lvw_flux,lmfsave,lrkmsave,lrkhsave,                            &
                        lbuoyproduction,lshearproduction,ltotaltransport,                                          &
#endif
                        idjd_t,mydiag_t)
                        
      t(is:ie,:)          = lt
      qg(is:ie,:)         = lqg
      qfg(is:ie,:)        = lqfg
      qlg(is:ie,:)        = lqlg
      cfrac(is:ie,:)      = lcfrac
      stratcloud(is:ie,:) = lstratcloud
      at_save(is:ie,:)    = lat
      ct_save(is:ie,:)    = lct
      if ( abs(iaero)>=2 ) then
        xtg(is:ie,:,:) = lxtg
      end if
      do k = 1,kl  
        u(is:ie,k) = lu(:,k) + uadj(is:ie)
        v(is:ie,k) = lv(:,k) + vadj(is:ie)
      end do  
#ifdef scm
#ifdef _OPENMP
      write(6,*) "ERROR: scm requires OMP is disabled"
      stop
#endif
      rkmsave(is:ie,:) = lrkmsave
      rkhsave(is:ie,:) = lrkhsave  
      wth_flux(is:ie,:) = lwth_flux 
      wq_flux(is:ie,:) = lwq_flux 
      uw_flux(is:ie,:) = luw_flux 
      vw_flux(is:ie,:) = lvw_flux 
      mfsave(is:ie,:) = lmfsave 
      buoyproduction(is:ie,:) = lbuoyproduction
      shearproduction(is:ie,:) = lshearproduction
      totaltransport(is:ie,:) = ltotaltransport
#endif

    end do ! tile = 1,ntiles
!$omp end do nowait

end select
    

#ifndef scm
if ( ngas>0 ) then 
!$omp do schedule(static) private(is,ie,k),  &
!$omp private(lt,lat,lct,idjd_t,mydiag_t),   &
!$omp private(ltr,lco2em,loh,lstrloss,ljmcf)
  do tile = 1,ntiles
    is = (tile-1)*imax + 1
    ie = tile*imax
    idjd_t = mod(idjd-1,imax)+1
    mydiag_t = ((idjd-1)/imax==tile-1).and.mydiag
    ltr      = tr(is:ie,:,:)
    lco2em   = co2em(is:ie,:)
    loh      = oh(is:ie,:)
    lstrloss = strloss(is:ie,:)
    ljmcf    = jmcf(is:ie,:)
    lat      = at_save(is:ie,:)
    lct      = ct_save(is:ie,:)
    lt       = t(is:ie,:)
    ! Tracers
    call tracervmix(lat,lct,lt,ps(is:ie),cdtq(is:ie),ltr,fnee(is:ie),fpn(is:ie),             &
                    frp(is:ie),frs(is:ie),lco2em,loh,lstrloss,ljmcf,mcfdep(is:ie),tile,imax)
    tr(is:ie,:,:) = ltr
  end do ! tile = 1,ntiles
!$omp end do nowait
end if 
#endif
   
return
end subroutine vertmix

!--------------------------------------------------------------
! Control subroutine for vertical mixing
subroutine vertmix_work(t,tss,eg,fg,kbsav,ktsav,convpsav,ps,qg,qfg,qlg,stratcloud,condc,cfrac, &
                        xtg,cduv,u,v,pblh,savu,savv,land,tscrn,qgscrn,ustar,f,condx,zs,        &
                        at,ct,                                                                 &
#ifdef scm
                        wth_flux,wq_flux,uw_flux,vw_flux,mfsave,rkmsave,rkhsave,               &
                        buoyproduction,shearproduction,totaltransport,                         &
#endif
                        idjd,mydiag)

use aerosolldr, only : naero        ! LDR prognostic aerosols
use cc_mpi, only : comm_world,       &
    ccmpi_barrier,ccmpi_abort       ! CC MPI routines
use cc_omp                          ! CC OpenMP routines
use const_phys                      ! Physical constants
use diag_m                          ! Diagnostic routines
use estab, only : establ            ! Liquid saturation function
use newmpar_m                       ! Grid parameters
use parm_m, only : diag,ktau,        &
    nvmix,dt,nlv,ia,ib,ja,jb,nmaxpr, &
    iaero,nlocal,av_vmod,            &
    amxlsq,dvmodmin                 ! Model configuration
use sigs_m                          ! Atmosphere sigma levels
use soil_m, only : zmin             ! Soil and surface data

implicit none
      
include 'kuocom.h'                  ! Convection parameters

integer, parameter :: ndvmod=0    ! 0 default, 1+ for dvmod tests
integer, intent(in) :: idjd
integer, dimension(imax), intent(in) :: kbsav, ktsav
integer, dimension(imax) :: kbase,ktop
integer, parameter :: ntest = 0
integer k, nt, iq
real, parameter :: lambda=0.45               ! coefficients for Louis scheme
real, parameter :: vkar4=0.4                 ! coefficients for Louis scheme
real, parameter :: bprmj=5.                  ! coefficients for Louis scheme
real, parameter :: cmj=5.                    ! coefficients for Louis scheme
real, parameter :: chj=2.6                   ! coefficients for Louis scheme
real, dimension(imax,kl,naero), intent(inout) :: xtg
real, dimension(imax,kl), intent(inout) :: t, qg, qfg, qlg
real, dimension(imax,kl), intent(inout) :: stratcloud, u, v, cfrac
real, dimension(imax,kl), intent(out) :: at, ct
real, dimension(imax,kl), intent(in) :: savu, savv
real, dimension(imax), intent(inout) :: pblh, ustar
real, dimension(imax), intent(in) :: tss, eg, fg, convpsav, ps, condc
real, dimension(imax), intent(in) :: cduv, tscrn, qgscrn, f, condx, zs
real, dimension(imax,kl) :: zh
real, dimension(imax,kl) :: rhs, guv, gt
real, dimension(imax,kl) :: au, cu
real, dimension(imax,kl) :: uav, vav
real, dimension(imax,kl) :: rkm, rkh
real, dimension(imax,kl) :: qs, betatt, betaqt, delthet, ri, rk_shal, thee
real, dimension(imax,kl) :: thebas
real, dimension(imax,kl-1) :: tmnht
real, dimension(imax) :: dz, dzr
real, dimension(imax) :: zhv, dvmod, dqtot, x, csq, sqmxl, fm, fh, theeb
real, dimension(imax) :: sigsp
real, dimension(kl) :: sighkap,sigkap,delons,delh
real, dimension(kl) :: prcpv
real rong, rlogs1, rlogs2, rlogh1, rlog12
real delsig, conflux, condrag
real w1, w2, pk, delta, es, dqsdt, betat, betaq, betac, qc, fice, al
real denma, denha, esp, epart, tsp, qbas
logical, intent(in) :: mydiag
logical, dimension(imax), intent(in) :: land

#ifdef scm
real, dimension(imax,kl), intent(inout) :: wth_flux, wq_flux, uw_flux
real, dimension(imax,kl), intent(inout) :: vw_flux
real, dimension(imax,kl), intent(out) :: rkmsave, rkhsave
real, dimension(imax,kl), intent(inout) :: buoyproduction, shearproduction
real, dimension(imax,kl), intent(inout) :: totaltransport
real, dimension(imax,kl-1), intent(inout) :: mfsave
real, dimension(imax,kl) :: mfout
#endif


#ifdef scm
! Initial flux to be added up below.
wth_flux(:,:) = 0.
wq_flux(:,:) = 0.
uw_flux(:,:) = 0.
vw_flux(:,:) = 0.
mfsave(:,:) = 0.
#endif

w1=0.
pk=0.
rkm = 0.
rkh = 0.

! Set-up potential temperature transforms
rong = rdry/grav
do k = 1,kl-1
  sighkap(k) = sigmh(k+1)**(-roncp)
  delons(k)  = rong*((sig(k+1)-sig(k))/sigmh(k+1))
end do      ! k loop
do k = 1,kl
  delh(k)   = -rong*dsig(k)/sig(k)  ! sign of delh defined so always +ve
  sigkap(k) = sig(k)**(-roncp)
end do      ! k loop

if ( diag .or. ntest>=1 .and. ntiles==1 ) then
  call maxmin(u,'%u',ktau,1.,kl)
  call maxmin(v,'%v',ktau,1.,kl)
  call maxmin(t,'%t',ktau,1.,kl)
  call maxmin(qg,'qg',ktau,1.e3,kl)
  call ccmpi_barrier(comm_world)
  if ( mydiag ) then
    write(6,*) 'sig ',sig
    write(6,*) 'dsig ',dsig
    write(6,*) 'delh ',delh
    write(6,*) 'in vertmix'
    write(6,"('uin ',9f8.3/4x,9f8.3)") u(idjd,:) 
    write(6,"('vin ',9f8.3/4x,9f8.3)") v(idjd,:) 
  end if
end if

! Calculate half level heights and temperatures
rlogs1=log(sig(1))
rlogs2=log(sig(2))
rlogh1=log(sigmh(2))
rlog12=1./(rlogs1-rlogs2)
tmnht(:,1)=(t(:,2)*rlogs1-t(:,1)*rlogs2+(t(:,1)-t(:,2))*rlogh1)*rlog12
! n.b. an approximate zh (in m) is quite adequate for this routine
zh(:,1) = t(:,1)*delh(1)
do k = 2,kl-1
  zh(:,k)    = zh(:,k-1) + t(:,k)*delh(k)
  tmnht(:,k) = ratha(k)*t(:,k+1) + rathb(k)*t(:,k)
end do      !  k loop
zh(:,kl) = zh(:,kl-1) + t(1:imax,kl)*delh(kl)

! Calculate theta
do k = 1,kl
  rhs(:,k) = t(:,k)*sigkap(k)  ! rhs is theta here
enddo      !  k loop
    
prcpv(1:kl) = sig(1:kl)**(-roncp)

if ( nmaxpr==1 .and. mydiag ) then
  write (6,"('thet_in',9f8.3/7x,9f8.3)") rhs(idjd,:)
end if
  
! Pre-calculate the buoyancy parameters if using qcloud scheme.
! Follow Smith's (1990) notation; gam() is HBG's notation for (L/cp)dqsdt.
! The factor of (1/sigkap)=T/theta in betatt differs from Smith's formulation
! because we use theta derivative rather than (dry static energy)/cp.
if ( (nvmix>0.and.nvmix<4) .or. nvmix==7 ) then
  delta=1./epsil-1.  ! i.e. 1/.622 -1., i.e. .6077
  do k=1,kl
    if ( sig(k)>.8 ) then ! change made 17/1/06
      do iq=1,imax
        es=establ(t(iq,k))
        pk=ps(iq)*sig(k)
        !qs(iq,k)=.622*es/max(1.,pk-es)  
        qs(iq,k)=.622*es/(pk-es)
        !dqsdt=qs(iq,k)*pk*(hl/rvap)/(t(iq,k)**2*max(1.,pk-es))
        dqsdt=qs(iq,k)*pk*(hl/rvap)/(t(iq,k)**2*(pk-es))
        rhs(iq,k)=rhs(iq,k)-(hlcp*qlg(iq,k)+hlscp*qfg(iq,k))*sigkap(k)   !Convert to thetal - used only to calc Ri
        betat=1./t(iq,k)
        qc=qlg(iq,k)+qfg(iq,k)
        fice=qfg(iq,k)/max(qc,1.e-12)
        betaq=delta/(1.+delta*qg(iq,k)-qc)
        al=1./(1.+hlcp*dqsdt)
        betac=cfrac(iq,k)*al * ((hlcp+fice*hlfcp)*betat - betaq/(1.-epsil) )
        betatt(iq,k)=(betat-dqsdt*betac)/sigkap(k)  !Beta_t_tilde
        betaqt(iq,k)=betaq+betac                    !Beta_q_tilde
      enddo   ! iq loop
    else  ! i.e. (sig(k)<.8)
      do iq=1,imax
        es=establ(t(iq,k))
        pk=ps(iq)*sig(k)
        qs(iq,k)=.622*es/max(1.,pk-es)  ! still need qs(); max for k=kl
        betat=1./t(iq,k)
        ! qc=qlg(iq,k)+qfg(iq,k)
        ! betaq=delta/(1.+delta*qg(iq,k)-qc)
        betaq=delta/(1.+delta*qg(iq,k))
        betatt(iq,k)=betat/sigkap(k)    !Beta_t_tilde
        betaqt(iq,k)=betaq              !Beta_q_tilde
      enddo   ! iq loop
    endif  ! (sig(k)>.8)
    if(diag.and.mydiag)then
      iq=idjd
      dqsdt=qs(iq,k)*pk*(hl/rvap)/(t(iq,k)**2*max(pk-es,1.))
      betat=1./t(iq,k)
      qc=qlg(iq,k)+qfg(iq,k)
      fice=qfg(iq,k)/max(qc,1.e-12)
      betaq=delta/(1+delta*qg(iq,k)-qc)
      ! al=1./(1.+gam(iq,k))
      al=1./(1.+hlcp*dqsdt)
      betac=cfrac(iq,k)*al*((hlcp+fice*hlfcp)*betat - betaq/(1.-epsil) )
      write(6,*)'k,qg,qs,qlg,qfg,qc ',k,qg(iq,k),qs(iq,k),qlg(iq,k),qfg(iq,k),qc
      write(6,*)'t,rhs,cfrac,fice ',t(iq,k),rhs(iq,k),cfrac(iq,k),fice
      write(6,*)'al,delta,betaq,betac ',al,delta,betaq,betac
      write(6,*)'betat,betatt,betaqt ',betat,betatt(iq,k),betaqt(iq,k)
    endif 
  enddo    !  k loop
else       ! other nvmix values (0 or 4+) still need qs()
  do k=1,kl
    do iq=1,imax
      es=establ(t(iq,k))
      qs(iq,k)=.622*es/max(1.,ps(iq)*sig(k)-es)  ! max for k=kl
    enddo   ! iq loop
  enddo    !  k loop
endif      ! (nvmix>0.and.nvmix<4).or.nvmix==7

do k=1,kl-1
  delthet(:,k)=rhs(:,k+1)-rhs(:,k)  ! rhs is theta or thetal here
enddo      !  k loop

!     ****** section for Geleyn shallow convection; others moved lower****
select case(ksc)
  case (-99)
    do k=kscbase,ksctop    ! new usage of ksc thu  02-17-2000
      do iq=1,imax
        delthet(iq,k)=delthet(iq,k)-hlcp*max(0.,qs(iq,k+1)-qg(iq,k+1)-qs(iq,k)+qg(iq,k) )
      enddo  ! iq loop
    enddo   !  k loop
  case (-98) ! modified Geleyn Jan '08
    do k=kscbase,ksctop    ! new usage of ksc thu  02-17-2000
      do iq=1,imax
        if(qg(iq,k)>tied_rh*qs(iq,k).or.qg(iq,k+1)>tied_rh*qs(iq,k+1))then
          delthet(iq,k)=delthet(iq,k)-hlcp*max(0.,qs(iq,k+1)-qg(iq,k+1)-qs(iq,k)+qg(iq,k) )
        endif
      enddo  ! iq loop
    enddo   !  k loop
  case (-97) ! modified Geleyn Jan '08
    do k=kscbase,ksctop    ! new usage of ksc thu  02-17-2000
      do iq=1,imax
        if(qg(iq,k)>tied_rh*qs(iq,k))then
          delthet(iq,k)=delthet(iq,k)-hlcp*max(0.,qs(iq,k+1)-qg(iq,k+1)-qs(iq,k)+qg(iq,k) )
        endif
      enddo  ! iq loop
    enddo   !  k loop
  case (-96) ! combination of Geleyn and jlm 83 (Feb 08)
    do k=1,ksctop    
      do iq=1,imax
        if(k<ktsav(iq).and.k>=kbsav(iq).and.condc(iq)<1.e-20)then  
          delthet(iq,k)=delthet(iq,k)-hlcp*max(0.,qs(iq,k+1)-qg(iq,k+1)-qs(iq,k)+qg(iq,k) )
        endif 
      enddo  ! iq loop
      do iq=1,imax
        if(k<ktsav(iq).and.k>=kbsav(iq).and.condc(iq)<1.e-20)then  
          write(6,*)'-96 iq,k,kbsav,ktsav ',iq,k,kbsav(iq),ktsav(iq),hlcp*(qs(iq,k+1)-qg(iq,k+1)-qs(iq,k)+qg(iq,k)),delthet(iq,k)
        endif 
      enddo  ! iq loop
    enddo   !  k loop
  case (-95) ! same as -99 but has tied_rh (e.g. .75) factor
    ! kshal(:)=0
    do k=kscbase,ksctop     
      do iq=1,imax
        delthet(iq,k)=delthet(iq,k)-tied_rh*hlcp*max(0.,(qs(iq,k+1)-qg(iq,k+1)-qs(iq,k)+qg(iq,k)) )
      enddo  ! iq loop
    enddo   !  k loop
  case (-94) ! combination of Geleyn and jlm 83 (Feb 08)
    do k=1,ksctop    
      do iq=1,imax
        if(k<ktsav(iq).and.k>=kbsav(iq).and.condc(iq)<1.e-20)then  
          delthet(iq,k)=0.
          write(6,*)'-94 iq,k,kbsav,ktsav ',iq,k,kbsav(iq),ktsav(iq),hlcp*(qs(iq,k+1)-qg(iq,k+1)-qs(iq,k)+qg(iq,k)),delthet(iq,k)
        endif 
      enddo  ! iq loop
    enddo   !  k loop
  case (-93) ! single-layer (pblh) version of -95
    do k=kscbase,kl/2
      do iq=1,imax  
        if(zh(iq,k)<pblh(iq).and.zh(iq,k+1)>pblh(iq))then
          delthet(iq,k)=delthet(iq,k)-tied_rh*hlcp*max(0.,(qs(iq,k+1)-qg(iq,k+1)-qs(iq,k)+qg(iq,k)) )
        endif
      enddo   ! iq loop
    enddo  ! k loop
  case (-92) ! capped-by-pblh version of -95
    do k=kscbase,kl/2
      do iq=1,imax  
        if(zh(iq,k)<pblh(iq))then
          delthet(iq,k)=delthet(iq,k)-tied_rh*hlcp*max(0.,(qs(iq,k+1)-qg(iq,k+1)-qs(iq,k)+qg(iq,k)) )
        endif
      enddo   ! iq loop  
    enddo  ! k loop 
  case (-91) ! capped-by-pblh (anywhere in layer) version of -95
    do k=2,kl/2
      do iq=1,imax  
        if(zh(iq,k-1)<pblh(iq))then
          delthet(iq,k)=delthet(iq,k)-tied_rh*hlcp*max(0.,(qs(iq,k+1)-qg(iq,k+1)-qs(iq,k)+qg(iq,k)) )
        endif
      enddo   ! iq loop
    enddo  ! k loop       
end select
!     ********* end of Geleyn shallow convection section ****************

! following now defined in vertmix (don't need to pass from sflux)
uav(:,:) = av_vmod*u(:,:) + (1.-av_vmod)*savu(:,:) 
vav(:,:) = av_vmod*v(:,:) + (1.-av_vmod)*savv(:,:) 
do k = 1,kl-1
  do iq = 1,imax
    dz(iq) =-tmnht(iq,k)*delons(k)  ! this is z(k+1)-z(k)
    dzr(iq)=1./dz(iq)
    zhv(iq)=1./zh(iq,k)
  end do      
  if(ndvmod==0)then
    do iq = 1,imax  
      dvmod(iq)=sqrt( (uav(iq,k+1)-uav(iq,k))**2+(vav(iq,k+1)-vav(iq,k))**2 )
    end do  
  else
    do iq = 1,imax  
      dvmod(iq)=ndvmod  ! just for tests
    end do  
  endif    ! (ndvmod==0)

  ! x is bulk ri *(dvmod **2); used to calc ri(k), rkh(k) etc
  select case (nvmix)
    case (1)
      ! usually nvmix=3       
      if ( sig(k)>.8 ) then ! change made 17/1/06
        dqtot(:)=qg(:,k+1)+qlg(:,k+1)+qfg(:,k+1)-(qg(:,k)+qlg(:,k)+qfg(:,k))
      else
        dqtot(:)=0.
      end if
      w1=dsig(k+1)/(dsig(k)+dsig(k+1)) 
      w2=1.-w1             !weight for upper level
      do iq=1,imax
        x(iq)=grav*dz(iq)*((w1*betatt(iq,k)+w2*betatt(iq,k+1))*delthet(iq,k) + (w1*betaqt(iq,k)+w2*betaqt(iq,k+1))*dqtot(iq) )
      enddo ! iq loop	 
      rhs(:,k)=t(:,k)*sigkap(k)   !need to re-set theta for nvmix=1-3
    case (2) !  jlm May '05
      ! usually nvmix=3       
      if ( sig(k)>.8 ) then ! change made 17/1/06
        dqtot(:)=qg(:,k+1)+qlg(:,k+1)+qfg(:,k+1)-(qg(:,k)+qlg(:,k)+qfg(:,k))
      else
        dqtot(:)=0.
      end if
      do iq=1,imax
        x(iq)=grav*dz(iq)*((min(betatt(iq,k),betatt(iq,k+1)))*delthet(iq,k) + (max(betaqt(iq,k),betaqt(iq,k+1)))*dqtot(iq) )
      enddo ! iq loop	
      rhs(:,k)=t(:,k)*sigkap(k)   !need to re-set theta for nvmix=1-3
    case (3,7)
      ! usually nvmix=3       
      if ( sig(k)>.8 ) then ! change made 17/1/06
        dqtot(:)=qg(:,k+1)+qlg(:,k+1)+qfg(:,k+1)-(qg(:,k)+qlg(:,k)+qfg(:,k))
      else
        dqtot(:)=0.
      end if
      w1 = 1.    !weight for lower level  usual  
      w2=1.-w1   !weight for upper level
      do iq=1,imax
        x(iq)=grav*dz(iq)*((w1*betatt(iq,k)+w2*betatt(iq,k+1))*delthet(iq,k) + (w1*betaqt(iq,k)+w2*betaqt(iq,k+1))*dqtot(iq) )
      enddo ! iq loop
      rhs(:,k)=t(:,k)*sigkap(k)   !need to re-set theta for nvmix=1-3
    case (5) ! non-cloudy x with qfg, qlg
      x(:)=grav*dz(:)*(delthet(:,k)/(tmnht(:,k)*sighkap(k))+.61*(qg(1:imax,k+1)-qg(1:imax,k)) &
       -qfg(1:imax,k+1)-qlg(1:imax,k+1)+qfg(1:imax,k)+qlg(1:imax,k) )
    case default ! original non-cloudy x, nvmix=4
      x(:)=grav*dz(:)*(delthet(:,k)/(tmnht(:,k)*sighkap(k))+.61*(qg(1:imax,k+1)-qg(1:imax,k)))  
  end select 

  ! fm and fh denote f(Louis style)*dvmod
  ! nb. an error exists in the Louis expression for c; this is corrected
  ! in the current code
  csq(:) = zhv(:)*(((1.+dz(:)*zhv(:))**(1./3.)-1.)*dzr(:))**3

  ! newest code, stable same as csiro9 here (originally for nvmix=4)
  do iq=1,imax
    sqmxl(iq)=(vkar4*zh(iq,k)/(1.+vkar4*zh(iq,k)/amxlsq))**2
    dvmod(iq)=max( dvmod(iq) , dvmodmin )
    ri(iq,k)=x(iq)/dvmod(iq)**2
    if(ri(iq,k)< 0.)then  ! unstable case
      ! first do momentum
      denma=dvmod(iq)+cmj*( 2.*bprmj*sqmxl(iq)*sqrt(-x(iq)*csq(iq)) )
      fm(iq)=dvmod(iq)-(2.*bprmj *x(iq))/denma
      ! now heat
      denha=dvmod(iq)+chj*( 2.*bprmj*sqmxl(iq)*sqrt(-x(iq)*csq(iq)) )
      fh(iq)=dvmod(iq)-(2.*bprmj *x(iq))/denha
    else                     ! stable case
      ! the following is the original Louis stable formula
      fm(iq)=dvmod(iq)/(1.+4.7*ri(iq,k))**2
      fh(iq)=fm(iq)
    endif
  enddo   ! iq loop

  if(nvmix<0)then   ! use quasi-neutral mixing for test runs only
    do iq=1,imax
      fm(iq)=1.
      fh(iq)=1.
    enddo   ! iq loop
    if(k<kl/2)then
      do iq=1,imax
        if(land(iq))fh(iq)=abs(nvmix)
      enddo   ! iq loop
    endif
  endif     !  (nvmix<0)

  ! calculate k's, guv and gt
  ! nonlocal scheme usually applied to temperature and moisture, not momentum
  ! (i.e. local scheme is applied to momentum for nlocal=0,1)
  if ( nvmix==7 ) then
    !added by Jing Huang on 4 Feb 2016
    rkm(:,k)=fm(:)*sqmxl(:)*dzr(:)
    rkh(:,k)=fh(:)*sqmxl(:)*dzr(:)
    where ( ri(:,k)>0. )
      sqmxl(:)=(vkar4*zh(:,k)/(1.+vkar4*zh(:,k)/amxlsq+vkar4*zh(:,k)*ri(:,k)/lambda))**2
      rkm(:,k)=dvmod(:)*dzr(:)*sqmxl(:)
      rkh(:,k)=rkm(:,k)/0.85
    end where
  else if ( nvmix/=0 ) then    ! use nvmix=0 only for special test runs
    rkm(:,k)=fm(:)*sqmxl(:)*dzr(:)
    rkh(:,k)=fh(:)*sqmxl(:)*dzr(:)
  else
    rkm(:,k)=0.
    rkh(:,k)=0.
  endif     ! (nvmix/=0)

  if((diag.or.ntest>=1).and.mydiag)then
    iq=idjd
    write(6,*)'k,dt,sqmxl,dzr ',k,dt,sqmxl(idjd),dzr(idjd)
    write(6,*)'k,t,t+,ps ',k,t(idjd,k),t(idjd,k+1),ps(idjd)
    write(6,*)'k,qg,qg+ ',k,qg(idjd,k),qg(idjd,k+1)
    write(6,*)'k,qs,qs+ ',k,qs(idjd,k),qs(idjd,k+1)
    es=establ(t(idjd,k))
    esp=establ(t(idjd,k+1))
    write(6,*)'k,es,es+,delthet ',k,es,esp,delthet(idjd,k)
    write(6,*)'k,fm,dvmod,ri,csq ',k,fm(idjd),dvmod(idjd),ri(idjd,k),csq(idjd)
    write(6,*)'qfg,qfg+ ',qfg(idjd,k),qfg(idjd,k+1)
    write(6,*)'qlg,qlg+,dqtot ',qlg(idjd,k),qlg(idjd,k+1),dqtot(iq)
    write(6,*)'cfrac,cfrac+ ',cfrac(idjd,k),cfrac(idjd,k+1)
    write(6,*)'betatt,betatt+,betatt*delthet ',betatt(iq,k),betatt(iq,k+1),betatt(iq,k)*delthet(iq,k)
    write(6,*)'betaqt,betaqt+,betaqt*dqtot   ',betaqt(iq,k),betaqt(iq,k+1),betaqt(iq,k)*dqtot(iq)
    write(6,*)'x,zh,tmnht,dz,sighkap ',x(idjd),zh(idjd,k),tmnht(idjd,k),dz(idjd),sighkap(k)
  endif     ! (diag.or.ntest>=1)
enddo      ! end of k loop

if( (diag.or.ntest>=1) .and. mydiag )then
  write(6,*)'before possible call to pbldif in vertmix'
  write (6,"('uav ',9f8.3/4x,9f8.3)") uav(idjd,:) 
  write (6,"('vav ',9f8.3/4x,9f8.3)") vav(idjd,:)
  write (6,"('t   ',9f8.3/4x,9f8.3)") t(idjd,:)
  write (6,"('qg  ',9f8.3/4x,9f8.3)") qg(idjd,:)
  write (6,"('qs  ',9f8.3/4x,9f8.3)") qs(idjd,:)
  write (6,"('thee',9f8.3/4x,9f8.3)") (prcpv(k)*t(idjd,k)*(t(idjd,k) + .5*hlcp*qs(idjd,k)) &
                                      /(t(idjd,k) - .5*hlcp*qs(idjd,k)),k=1,kl)
endif
if(nmaxpr==1.and.mydiag)then
  write (6,"('rino_v',9f9.3/6x,9f9.3)") ri(idjd,1:kl-1)
  write (6,"('rkh0 ',9f9.3/5x,9f9.3)") rkh(idjd,1:kl-2)
  write (6,"('rkm0 ',9f9.3/5x,9f9.3)") rkm(idjd,1:kl-2)
endif

if (nlocal/=0) then
  call pbldif(rhs,rkm,rkh,uav,vav,                   &
              t,pblh,ustar,f,ps,fg,eg,qg,land,cfrac  &
#ifdef scm
              ,wth_flux,wq_flux                      &
#endif
              )  ! rhs is theta or thetal
  ! n.b. *** pbldif partially updates qg and theta (t done during trim)	 
  ! and updates rkh and rkm arrays
  if(nmaxpr==1.and.mydiag)then
    write (6,"('pblh ',f8.2)") pblh(idjd)
    write (6,"('rkh1 ',9f9.3/5x,9f9.3)") rkh(idjd,1:kl-2)
    write (6,"('rkm1 ',9f9.3/5x,9f9.3)") rkm(idjd,1:kl-2)
  endif
  if(nmaxpr==1.and.mydiag) write (6,"('thet_pbl',9f8.3/8x,9f8.3)") rhs(idjd,:)
  if( (diag.or.ntest>=1) .and. mydiag ) write (6,"('qg ',9f8.3/4x,9f8.3)") qg(idjd,:)
  if(diag.and.ntiles==1)then
    call printa('rkh ',rkh,ktau,nlv,ia,ib,ja,jb,0.,1.)
    call printa('cond',condx,ktau,1,ia,ib,ja,jb,0.,1.)
    call printa('zs  ',zs,ktau,1,ia,ib,ja,jb,0.,1.)
  endif
else
  pblh(:) = zmin
endif      ! (nlocal>0)

rk_shal(:,:)=0.
!     ***** ***** section for jlm shallow convection v4 *****************
select case (ksc)
  case (81) 
    do k=1,ksctop-1   ! or ksctop?  
      do iq=1,imax
        if(sig(ktsav(iq))>sig_ct.and.k<ktsav(iq).and.condc(iq)<1.e-20)then  
          rk_shal(iq,k)=tied_con
          rk_shal(iq,k+1)=tied_over
        endif ! (sig(ktsav(iq))>sig_ct.and.k<ktsav(iq).and....)
      enddo  ! iq loop
    enddo   !  k loop
  case (82)  
    do k=1,ksctop-1    
      do iq=1,imax
        if(sig(ktsav(iq))>sig_ct.and.k<ktsav(iq).and.k>=kbsav(iq).and.condc(iq)<1.e-20)then  
          rk_shal(iq,k)=tied_con
          rk_shal(iq,k+1)=tied_over
        endif ! (sig(ktsav(iq))>sig_ct.and.k<ktsav(iq).and....)
      enddo  ! iq loop
    enddo   !  k loop
  case (83)  
    do k=1,ksctop-1    
      do iq=1,imax
        if(sig(k)>sig_ct.and.k<ktsav(iq).and.k>=kbsav(iq).and.condc(iq)<1.e-20)then  
          rk_shal(iq,k)=tied_con
          rk_shal(iq,k+1)=tied_over
        endif ! (sig(ktsav(iq))>sig_ct.and.k<ktsav(iq).and....)
      enddo  ! iq loop
    enddo   !  k loop
  case (91)  
    do k=1,ksctop ! May 08
      do iq=1,imax
        if(ktsav(iq)<kl-1.and.k<ktsav(iq))then  ! April 04
          rk_shal(iq,k)=tied_con
          rk_shal(iq,k+1)=tied_over
        endif ! (ktsav(iq)<0.and.k<abs(ktsav(iq)))
      enddo  ! iq loop
    enddo   !  k loop
  case (92)  
    do k=1,ksctop ! May 08
      do iq=1,imax
        if(k>=kbsav(iq).and.k<ktsav(iq).and.condc(iq)<1.e-20)then  ! May 08
          rk_shal(iq,k)=tied_con
          rk_shal(iq,k+1)=tied_over
        endif ! (ktsav(iq)<0.and.k<abs(ktsav(iq)))
      enddo  ! iq loop
    enddo   !  k loop
end select
!     *********** end of jlm shallow convection section *****************

!     ************ section for Tiedtke shallow convection *******************
if(ksc==99)then
  do iq=1,imax
    theeb(iq)=prcpv(kscbase)*t(iq,kscbase)*(t(iq,kscbase)+.5*hlcp*qs(iq,kscbase))/(t(iq,kscbase)-.5*hlcp*qs(iq,kscbase))
  enddo    ! iq loop
  if(kscsea==1)then  ! Tiedtke done only over sea
    do k=kscbase+1,ksctop
      do iq=1,imax
        if(.not.land(iq))then 
          thee(iq,k)=prcpv(k)*t(iq,k)*(t(iq,k) + .5*hlcp*qs(iq,k))/(t(iq,k) - .5*hlcp*qs(iq,k))
          if(qg(iq,kscbase)>tied_rh*qs(iq,kscbase).and.thee(iq,k)<theeb(iq))then         !  default tied_rh=.75
            rk_shal(iq,k-1)=tied_con     !  m**2/sec  6., originally 10.
            rk_shal(iq,k)=tied_over      !  m**2/sec
          endif ! (qg(iq,kscbase)>rhscon*qs(iq,kscbase).....
        endif  ! (.not.land(iq)) 
      enddo   ! iq loop
    enddo    ! end of k=kscbase+1,ksctop loop
  else       !  i.e. Tiedtke original scheme over land and sea
    do k=kscbase+1,ksctop  ! typically kscbase=3 & ksctop=6
      do iq=1,imax
        thee(iq,k)=prcpv(k)*t(iq,k)*(t(iq,k) + .5*hlcp*qs(iq,k))/(t(iq,k) - .5*hlcp*qs(iq,k))
        if(qg(iq,kscbase)>tied_rh*qs(iq,kscbase).and.thee(iq,k)<theeb(iq))then !  default tied_rh=.75
          rk_shal(iq,k-1)=tied_con     !  m**2/sec  6., originally 10.
          rk_shal(iq,k)=tied_over      !  m**2/sec
          if(ntest==3.and.k==ksctop)then
            write(6,*) 'ktau,iq,theeb,thee,delthee ',ktau,iq,theeb(iq),thee(iq,k),theeb(iq)-thee(iq,k)
          endif
        endif ! (qg(iq,kscbase)>rhscon*qs(iq,kscbase).....
      enddo  ! iq loop
    enddo   ! end of k=kscbase+1,ksctop loop
  endif     ! (kscsea==1)  .. else ..
endif       ! (ksc==99)
!     *********** end of Tiedtke shallow convection section **************

!     *********** Tiedtke_ldr-style shallow convection 93 ************
if(ksc>=93.and.ksc<=96)then   
  ! Calculate LCL for near surface air
  ! Assume qstar=qtg(iq,1) but allow some sub-grid variability
  ! The formula for tsp is eqn (21) of Bolton (1980), MWR 108, 1046-1053.
  if(sigkscb>0.)then   ! uses qg1
    do iq=1,imax
      ! Leon used factor 1.01 over sea; here div by tied_rh (e.g. .99)	 
      epart=qg(iq,1)*.01*ps(iq)/(tied_rh*epsil) !in hPa 
      tsp=2840./(3.5*log(tscrn(iq))-log(epart)-4.805) + 55.   
      sigsp(iq)=(tsp/tscrn(iq))**(1./roncp)
    enddo
  else                    ! uses qgscrn
    do iq=1,imax
      ! Leon used factor 1.01 over sea; here div by tied_rh (e.g. .99)	 
      epart=qgscrn(iq)*.01*ps(iq)/(tied_rh*epsil) !in hPa   
      tsp=2840./(3.5*log(tscrn(iq))-log(epart)-4.805) + 55.   
      sigsp(iq)=(tsp/tscrn(iq))**(1./roncp)
    enddo
  endif  ! (sigkscb>0.) .. else ..
  ! Look for the lowest layer s.t. the top of the layer is above the LCL,
  ! and call that layer the cloud base. 
  kbase(:)=kl
  do k=ksctop-1,1,-1
    do iq=1,imax
      if(sigsp(iq)>sigmh(k+1))kbase(iq)=k
    enddo
  enddo
  if(nlocal/=0)then  ! like old ksc=95, but ensures LCL within PBL
    do iq=1,imax
      if(zh(iq,kbase(iq))>pblh(iq))then 
        kbase(iq)=kl
      endif
    enddo  ! iq loop
  endif  ! (nlocal/=0)
! following has some jlm shortcuts	 
  do iq=1,imax
    k=kbase(iq)
    qbas=qs(iq,k)
    theeb(iq)=prcpv(k)*t(iq,k)*(t(iq,k) + .5*hlcp*qbas)/(t(iq,k) - .5*hlcp*qbas)
  enddo  ! iq loop
  ktop(:)=kbase(:)
  do k=2,ksctop+1
    do iq=1,imax
      thee(iq,k)=prcpv(k)*t(iq,k)*(t(iq,k) + .5*hlcp*qs(iq,k))/(t(iq,k) - .5*hlcp*qs(iq,k))
      if(theeb(iq)>thee(iq,k).and.ktop(iq)==k-1)then
        ktop(iq)=k       ! also checking if contiguous
      endif
    enddo  ! iq loop
  enddo   ! k loop	 
  if(ksc==94)then  ! from April same as 93* 
    do k=2,ksctop
      do iq=1,imax
!       if(ktop(iq)>kbase(iq).and.k<=ktop(iq).and.k>kbase(iq))then
        if(k>kbase(iq).and.k<=ktop(iq))then
          rk_shal(iq,k-1)=tied_con     !  m**2/sec  6., originally 10.
          rk_shal(iq,k)=tied_over      !  m**2/sec
        endif
      enddo  ! iq loop
    enddo   ! k loop
  endif     ! (ksc==94)
  if(ksc==93)then
    do k=2,ksctop
      do iq=1,imax
        if(ktop(iq)>kbase(iq).and.k<=ktop(iq).and.k>kbase(iq).and.ktop(iq)<=ksctop)then
          rk_shal(iq,k-1)=tied_con     !  m**2/sec  6., originally 10.
          rk_shal(iq,k)=tied_over      !  m**2/sec
        endif
      enddo  ! iq loop
    enddo   ! k loop
  endif     ! (ksc==93)
  if(ksc==95)then   ! new from 7 April
    do k=2,ksctop
      do iq=1,imax
        if(k>kbase(iq).and.k<=ktop(iq).and.condc(iq)<1.e-20)then
          rk_shal(iq,k-1)=tied_con     !  m**2/sec  6., originally 10.
          rk_shal(iq,k)=tied_over      !  m**2/sec
        endif
      enddo  ! iq loop
    enddo   ! k loop
  endif     ! (ksc==95)
  if(ksc==96)then   ! applied from sfce up
    do k=2,ksctop
      do iq=1,imax
        if(ktop(iq)>kbase(iq).and.k<=ktop(iq).and.condc(iq)<1.e-20)then  
          rk_shal(iq,k-1)=tied_con     !  m**2/sec  6., originally 10.
          rk_shal(iq,k)=tied_over      !  m**2/sec
        endif
      enddo  ! iq loop
    enddo   ! k loop
  endif     ! (ksc==96)
endif       ! (ksc>=93.and.ksc<=96)
!     *********** end of Tiedtke_ldr shallow convection 93-96 *********

!     *************** Tiedtke_jlm shallow convection 97 ***************
if(ksc==97)then
! this one is similar to ksc=99, but uses enhanced qbas <= qs
! and also works up the column with other possible thebas values
  do k=kscbase,ksctop
    do iq=1,imax
      qbas=min(qg(iq,k)/tied_rh,qs(iq,k))
!     thebas(iq,k)=t(iq,k)+hlcp*qbas  + hght
      thebas(iq,k)=prcpv(k)*t(iq,k)*(t(iq,k) + .5*hlcp*qbas)/(t(iq,k) - .5*hlcp*qbas)
    enddo  ! iq loop
  enddo   ! k loop
  theeb(:)=thebas(:,kscbase)
  do k=kscbase+1,ksctop+1
    do iq=1,imax
      thee(iq,k)=prcpv(k)*t(iq,k)*(t(iq,k) + .5*hlcp*qs(iq,k))/(t(iq,k) - .5*hlcp*qs(iq,k))
      if(theeb(iq)>thee(iq,k))then
        rk_shal(iq,k-1)=tied_con       !  m**2/sec  6., originally 10.
        rk_shal(iq,k)=tied_over        !  m**2/sec
      else
        theeb(iq)=thebas(iq,k) ! ready for next k in k-loop
      endif
    enddo  ! iq loop
  enddo   ! k loop
  if(ntest==3)then
    do iq=1,imax
      if(rk_shal(iq,ksctop-1)>.9*tied_con)then 
        write(6,*) 'iq,land,rk_shal ',iq,land(iq),rk_shal(iq,ksctop-1)
      endif
    enddo
  endif   ! (ntest==3)
endif     ! (ksc==97)
!     *********** end of Tiedtke_jlm shallow convection 97 *************

! add in effects of shallow convection
do k=1,kl-1  
  rkh(:,k)=rkh(:,k)+rk_shal(:,k)
enddo   !  k loop
if(kscmom==1)then
  do k=1,kl-1  
    rkm(:,k)=rkm(:,k)+rk_shal(:,k)
  enddo   !  k loop
endif     ! (kscmom==1)
      
if(ksc/=0.and.(ntest/=0.or.diag).and.nproc==1)then
  do iq=1,imax
    if(rk_shal(iq,1)>rk_shal(iq,2))then
      write(6,*) 'iq,rk_shal1,rk_shal2',iq,rk_shal(iq,1),rk_shal(iq,2)
    endif
  enddo
  if (mydiag) then 
    iq=idjd
    write(6,*)'for shallow convection in vertmix '
    write(6,*)'ktsav,ksctop ',ktsav(idjd),ksctop
    write(6,*)'kbase,ktop ',kbase(idjd),ktop(idjd)
    write(6,*)'kbsav,ktsav,theeb: ',kbsav(iq),ktsav(iq),theeb(iq)
    write (6,"('rk_shal ',9f7.2/(9x,9f7.2))") (rk_shal(idjd,k),k=1,kl)
    write (6,"('rh   ',9f7.2/(5x,9f7.2))") (100.*qg(idjd,k)/qs(idjd,k),k=1,kl)
    write (6,"('qs   ',9f7.3/(5x,9f7.3))") (1000.*qs(idjd,k),k=1,kl)
    write (6,"('qg   ',9f7.3/(5x,9f7.3))") (1000.*qg(idjd,k),k=1,kl)
    write (6,"('qbas ',9f7.3/(5x,9f7.3))") (1000.*qg(idjd,k)/tied_rh,k=1,kl)
    write (6,"('t    ',9f7.2/(5x,9f7.2))") (t(idjd,k),k=1,kl)
    write (6,"('thebas',9f7.2/(6x,9f7.2))") (thebas(iq,k),k=1,kl)
!   write (6,"('hs  ',9f7.2/(4x,9f7.2))") (t(idjd,k)+hlcp*qs(idjd,k),k=1,kl)
    write (6,"('thee',9f7.2/(4x,9f7.2))") (thee(idjd,k),k=1,kl)
  endif ! mydiag
endif


do k = 1,kl-1
  delsig   = sig(k+1) - sig(k)
  dz(:)    = -tmnht(:,k)*delons(k)  ! this is z(k+1)-z(k)
  dzr(:)   = 1./dz(:)
  guv(:,k) = rkm(:,k)*dt*delsig*dzr(:)**2
  gt(:,k)  = rkh(:,k)*dt*delsig*dzr(:)**2
end do      ! k loop
guv(:,kl) = 0.
gt(:,kl)  = 0.

if ( diag .and. ntiles==1 ) then
  call maxmin(rkh,'rk',ktau,.01,kl-1)
  if ( mydiag ) then
    write(6,*)'vertmix guv ',(guv(idjd,k),k=1,kl)
    write(6,*)'vertmix gt ',(gt(idjd,k),k=1,kl)
  end if
  call printa('t   ',t,ktau,nlv,ia,ib,ja,jb,200.,1.)
  call printa('tss ',tss,ktau,nlv,ia,ib,ja,jb,200.,1.)
  call printa('eg  ',eg,ktau,nlv,ia,ib,ja,jb,0.,1.)
  call printa('fg  ',fg,ktau,nlv,ia,ib,ja,jb,0.,1.)
  call printa('thet',rhs,ktau,nlv,ia,ib,ja,jb,200.,1.)
end if

if ( ncvmix>0 ) then  ! cumulus mixing of momentum - jlm version
  do k = kuocb+1-ncvmix,kl-3
    ! for no-conv-points, doing k=1-ncvmix,1 with convpsav=0.
    !  1 below for ncvmix=2
    where ( kbsav>0 .and. k>=kbsav+1-ncvmix .and. k<ktsav )
      guv(:,k) = guv(:,k) - convpsav(:)*.5 ! with factor of .5
    end where
  enddo    ! k loop
  if ( diag .and. mydiag ) then
    write(6,*)'vertmix after conv; kb,kt,convp',kbsav(idjd),ktsav(idjd),convpsav(idjd)
    write(6,*)'new guv',(guv(idjd,k),k=1,kl)
  end if
end if      !   (ncvmix>0)

conflux = grav*dt/dsig(1)
condrag = grav*dt/(dsig(1)*rdry)
! first do theta (then convert back to t)
at(:,1)  = 0.
ct(:,kl) = 0.

do k = 2,kl
  at(:,k) = -gt(:,k-1)/dsig(k)
end do
do k = 1,kl-1
  ct(:,k) = -gt(:,k)/dsig(k)
enddo    !  k loop
if ( ( diag .or. ntest==2 ) .and. mydiag ) then
  write(6,*)'ktau,fg,tss,ps ',ktau,fg(idjd),tss(idjd),ps(idjd)
  write(6,*)'at ',(at(idjd,k),k=1,kl)
  write(6,*)'ct ',(ct(idjd,k),k=1,kl)
  write(6,*)'rhs ',(rhs(idjd,k),k=1,kl)
end if      ! (ntest==2)

!--------------------------------------------------------------
! Temperature
if ( nmaxpr==1 .and. mydiag ) write (6,"('thet_inx',9f8.3/8x,9f8.3)") rhs(idjd,:)
rhs(:,1) = rhs(:,1) - (conflux/cp)*fg(:)/ps(1:imax)
call trim(at,ct,rhs,imax,kl)   ! for t
if ( nmaxpr==1 .and. mydiag ) write (6,"('thet_out',9f8.3/8x,9f8.3)") rhs(idjd,:)
do k = 1,kl
  t(1:imax,k) = rhs(:,k)/sigkap(k)
enddo    !  k loop
if ( diag ) then
  if ( mydiag ) then
    write(6,*)'vertmix eg,fg ',eg(idjd),fg(idjd)
  end if
  if ( ntiles==1 ) then
    call printa('thet',rhs,ktau,nlv,ia,ib,ja,jb,200.,1.)
  end if  
end if
 
#ifdef scm
! Save Km and Kh for output
rkmsave(:,:) = rkm(:,:)
rkhsave(:,:) = rkh(:,:)
! counter-gradied included in pbldif.f90
wth_flux(:,1) = fg(:)*rdry*t(1:imax,1)/(ps(1:imax)*cp)
do k = 1,kl-1
  wth_flux(:,k+1) = rkh(:,k)*(rhs(1:imax,k+1)-rhs(1:imax,k))*(grav/rdry)*sig(k)/(t(1:imax,k)*dsig(k))
end do    
#endif

!--------------------------------------------------------------
! Moisture
rhs = qg(1:imax,:)
rhs(:,1) = rhs(:,1) - (conflux/hl)*eg/ps(1:imax)
! could add extra sfce moisture flux term for crank-nicholson
call trim(at,ct,rhs,imax,kl)    ! for qg
qg(1:imax,:) = rhs
if ( diag .and. mydiag ) then
  write(6,*)'vertmix rhs & qg after trim ',(rhs(idjd,k),k=1,kl)
  write (6,"('qg ',9f7.3/(8x,9f7.3))") (1000.*qg(idjd,k),k=1,kl)
end if

#ifdef scm  
! counter-gradied included in pbldif.f90
wq_flux(:,1) = eg(:)*rdry*t(1:imax,1)/(ps(1:imax)*hl)
do k = 1,kl-1
  wq_flux(:,k+1) = rkh(:,k)*(qg(1:imax,k+1)-qg(1:imax,k))*(grav/rdry)*sig(k)/(t(1:imax,k)*dsig(k))
end do  
#endif
  
!--------------------------------------------------------------
! Cloud microphysics terms
if ( ldr/=0 ) then
  ! now do qfg
  call trim(at,ct,qfg,imax,kl)
  ! now do qlg
  call trim(at,ct,qlg,imax,kl)
  ! now do cfrac
  call trim(at,ct,cfrac,imax,kl)
  ! now do stratcloud
  call trim(at,ct,stratcloud,imax,kl)
end if    ! (ldr/=0)
  
!--------------------------------------------------------------
! Aerosols
if ( abs(iaero)>=2 ) then
  do nt = 1,naero
    call trim(at,ct,xtg(:,:,nt),imax,kl)
    xtg(:,:,nt) = max( xtg(:,:,nt), 0. )
  end do
end if ! (abs(iaero)>=2)

!--------------------------------------------------------------
! Momentum terms
au(:,1) = cduv(:)*condrag/tss(:)
do k = 2,kl
  au(:,k) = -guv(:,k-1)/dsig(k)
enddo    !  k loop
do k = 1,kl-1
  cu(:,k) = -guv(:,k)/dsig(k)
enddo    !  k loop
cu(:,kl) = 0.
if ( ( diag .or. ntest==2 ) .and. mydiag ) then
  write(6,*)'au ',(au(idjd,k),k=1,kl)
  write(6,*)'cu ',(cu(idjd,k),k=1,kl)
end if      ! (ntest==2)

! first do u
call trim(au,cu,u,imax,kl)
  
#ifdef scm
uw_flux(:,1) = -cduv(:)*u(1:imax,1)
do k = 1,kl-1
  uw_flux(:,k+1) = rkm(:,k)*(u(1:imax,k+1)-u(1:imax,k))*(grav/rdry)*sig(k)/(t(1:imax,k)*dsig(k))
end do
#endif
  
! now do v; with properly unstaggered au,cu
call trim(au,cu,v,imax,kl)    ! note now that au, cu unstaggered globpea

#ifdef scm
vw_flux(:,1) = -cduv(:)*v(1:imax,1)
do k = 1,kl-1
  vw_flux(:,k+1) = rkm(:,k)*(v(1:imax,k+1)-v(1:imax,k))*(grav/rdry)*sig(k)/(t(1:imax,k)*dsig(k))
end do
#endif
  
if ( ( diag .or. ntest>=1 ) .and. mydiag ) then
  write(6,*)'after trim in vertmix '
  write (6,"('thet',9f7.2/(8x,9f7.2))") (sigkap(k)*t(idjd,k),k=1,kl) 
  write (6,"('t   ',9f7.2/(8x,9f7.2))") (t(idjd,k),k=1,kl) 
  write (6,"('u   ',9f7.2/(8x,9f7.2))") (u(idjd,k),k=1,kl) 
  write (6,"('v   ',9f7.2/(8x,9f7.2))") (v(idjd,k),k=1,kl) 
  write(6,*)'cduv,cduv+1,cduvj+1,tss ',cduv(idjd),cduv(idjd+1),cduv(idjd+il),tss(idjd)
  write (6,"('au ',9f7.3/(8x,9f7.3))") (au(idjd,k),k=1,kl) 
  write (6,"('cu ',9f7.3/(8x,9f7.3))") (cu(idjd,k),k=1,kl) 
end if
if ( nmaxpr==1 .and. mydiag ) then
  write (6,"('qg_vm ',9f7.3/(6x,9f7.3))") (1000.*qg(idjd,k),k=1,kl)
end if
           
return
end subroutine vertmix_work

subroutine pbldif(theta,rkm,rkh,uav,vav,                &
                  t,pblh,ustar,f,ps,fg,eg,qg,land,cfrac &
#ifdef scm
                  ,wth_flux,wq_flux                     &
#endif
                  )   

use cc_mpi, only : mydiag, myid
use cc_omp
use const_phys
use newmpar_m
use parm_m
use sigs_m     !sig,sigmh

implicit none

include 'kuocom.h'

integer, parameter :: ntest=0
integer, parameter :: nrkmin=1   ! 1 original (& from 0510); 2 new; 3 newer
integer, parameter :: npblmin=4  ! 1 original (best for Oz); 2 new ; 3,4 newer
integer kmax,iq
integer k                 ! level index
integer, dimension(imax) :: iflag

!------------------------------------------------------------------------
! 
! Atmospheric boundary layer computation.
!
! Nonlocal scheme that determines eddy diffusivities based on a
! diagnosed boundary layer height and a turbulent velocity scale;
! also, countergradient effects for heat and moisture, and constituents
! are included, along with temperature and humidity perturbations which 
! measure the strength of convective thermals in the lower part of the 
! atmospheric boundary layer.
!
! For more information, see Holtslag, A.A.M., and B.A. Boville, 1993:
! Local versus Nonlocal Boundary-Layer Diffusion in a Global Climate
! Model. J. Clim., vol. 6., p. 1825--1842.
!
! Updated by Holtslag and Hack to exclude the surface layer from the
! definition of the boundary layer Richardson number. Ri is now defined
! across the outer layer of the pbl (between the top of the surface
! layer and the pbl top) instead of the full pbl (between the surface and
! the pbl top). For simplicity, the surface layer is assumed to be the
! region below the first model level (otherwise the boundary layer depth 
! determination would require iteration).
!
! NOTE that all calculation in this module is at temperature points (DARLAM)
!------------------------------Code history--------------------------------
!
! Original version:  B. Boville
! Standardized:      J. Rosinski, June 1992
! Reviewed:          B. Boville, P. Rasch, August 1992
! Reviewed:          B. Boville, P. Rasch, April 1996
!
! Modified for boundary layer height diagnosis: Bert Holtslag, june 1994
! >>>>>>>>>  (Use ricr = 0.3 in this formulation)
!
!------------------------------Arguments--------------------------------
!
! Input arguments:u,v,fg,eg,theta,ustar,uav,vav
      
! Input & Output arguments
real, dimension(imax,kl) :: rkm            ! eddy diffusivity for momentum [m2/s]
real, dimension(imax,kl) :: rkh            ! eddy diffusivity for heat [m2/s]
real, dimension(imax,kl) :: theta          ! potential temperature [K]
!     also qg                              ! mixing ratio [kg/kg}

real, dimension(imax,kl) :: cgh            ! counter-gradient term for heat [K/m]
real, dimension(imax,kl) :: cgq            ! counter-gradient term for constituents
real, dimension(imax,kl) :: zg
real ztodtgor,delsig,tmp1,sigotbk,sigotbkm1
real cgs                     ! counter-gradient star (cg/flux)
real, dimension(imax,kl), intent(in) :: t
real, dimension(imax), intent(inout) :: pblh
real, dimension(imax), intent(inout) :: ustar
real, dimension(imax), intent(in) :: f
real, dimension(imax), intent(in) :: ps
real, dimension(imax), intent(in) :: fg
real, dimension(imax), intent(in) :: eg
real, dimension(imax,kl), intent(inout) :: qg
logical, dimension(imax), intent(in) :: land
real, dimension(imax,kl), intent(in) :: cfrac
#ifdef scm
real, dimension(imax,kl), intent(inout) :: wth_flux
real, dimension(imax,kl), intent(inout) :: wq_flux
#endif
!
!---------------------------Local parameters----------------------------
!
real, parameter :: tiny=1.e-36             ! lower bound for wind magnitude
!
!---------------------------Local workspace-----------------------------
!

real, dimension(imax) :: heatv          ! surface virtual heat flux
real, dimension(imax) :: thvref         ! reference level virtual temperature
real, dimension(imax) :: phiminv        ! inverse phi function for momentum
real, dimension(imax) :: phihinv        ! inverse phi function for heat 
real, dimension(imax) :: wm             ! turbulent velocity scale for momentum
real, dimension(imax) :: rkhfs          ! surface kinematic heat flux [mK/s]
real, dimension(imax) :: rkqfs          ! sfc kinematic constituent flux [m/s]
real, dimension(imax,kl) :: rino        ! bulk Richardson no. from level to ref lev
real, dimension(imax) :: tlv            ! ref. level pot tmp + tmp excess
real, dimension(imax) :: wstr           ! w*, convective velocity scale
real, dimension(imax) :: obklen         ! Obukhov length
real tkv                                ! model level potential temperature
real therm                              ! thermal virtual temperature excess
real pmid                               ! midpoint pressures
real vvk                                ! velocity magnitude squared
real zmzp                               ! level height halfway between zm and zp
real fak1                               ! k*ustar*pblh
real fak2                               ! k*wm*pblh
real fak3                               ! fakn*wstr/wm 
real pblk                               ! level eddy diffusivity for momentum
real pr                                 ! Prandtl number for eddy diffusivities
real zm                                 ! current level height
real zp                                 ! current level height + one level up
real zl                                 ! zmzp / Obukhov length
real zh                                 ! zmzp / pblh      at half levels
real zzh                                ! (1-(zmzp/pblh))**2      at half levels
real rrho                               ! 1./bottom level density (temporary)
real term                               ! intermediate calculation
real fac                                ! interpolation factor

!------------------------------Commons----------------------------------
real, dimension(imax,kl) :: uav,vav

real, parameter :: betam  = 15.0  ! Constant in wind gradient expression
real, parameter :: betas  = 5.0   ! Constant in surface layer gradient expression
real, parameter :: betah  = 15.0  ! Constant in temperature gradient expression
real, parameter :: fak    = 8.5   ! Constant in surface temperature excess
real, parameter :: fakn   = 7.2   ! Constant in turbulent prandtl number
real, parameter :: ricr   = 0.25  ! Critical richardson number
real, parameter :: sffrac = 0.1   ! Surface layer fraction of boundary layer
real, parameter :: vk     = 0.4   ! Von Karman's constant
real ccon    ! fak * sffrac * vk
real binm    ! betam * sffrac
real binh    ! betah * sffrac
real rkmin ! minimum eddy coeffs based on Hourdin et. al. (2001)
      
kmax=kl-1
if (nlocal==6) then
  fac=10.
else
  fac=100.
end if
binh   = betah*sffrac
binm   = betam*sffrac
ccon   = fak*sffrac*vk
!****************************************************************
zg(1:imax,1) = bet(1)*t(1:imax,1)/grav
do k = 2,kl
  zg(1:imax,k) = zg(1:imax,k-1)+(bet(k)*t(1:imax,k)+betm(k)*t(1:imax,k-1))/grav
enddo         ! k  loop
cgh(:,:) = 0.   ! 3D
cgq(:,:) = 0.   ! 3D
if ( ktau==1 .and. myid==0 .and. ntiles==1 ) then
  write(6,*) 'in pbldif nrkmin,npblmin: ',nrkmin,npblmin 
end if
      
! Compute kinematic surface fluxes
do iq=1,imax
  pmid=ps(iq)*sigmh(1) 
  rrho = rdry*t(iq,1)/pmid
  ustar(iq) = max(ustar(iq),0.01)
  rkhfs(iq) = fg(iq)*rrho/cp           !khfs=w'theta'
  rkqfs(iq) = eg(iq)*rrho/hl           !kqfs=w'q'

  ! Compute various arrays for use later:

  thvref(iq) = theta(iq,1)*(1.0 + 0.61*qg(iq,1))
  heatv(iq)  = rkhfs(iq) + 0.61*theta(iq,1)*rkqfs(iq)
  wm(iq)     = 0.
  ! obklen at t point
  obklen(iq) = -thvref(iq)*ustar(iq)**3/(grav*vk*(heatv(iq) + sign(1.e-10,heatv(iq))))
         
  ! >>>> Define first a new factor fac=100 for use in Richarson number
  !      Calculate virtual potential temperature first level
  !      and initialize pbl height to z1 i.e  1st full level

  pblh(iq) = zg(iq,1)    
  rino(iq,1) = 0.
enddo

! PBL height calculation:
! Search for level of pbl. Scan upward until the Richardson number between
! the first level and the current level exceeds the "critical" value.
! Richardson no. is computed using eq. (4.d.18) NCAR technical report, CCM3)

iflag(:)=0
do k=2,kmax
  do iq=1,imax
    vvk = (uav(iq,k) - uav(iq,1))**2 + (vav(iq,k) - vav(iq,1))**2 + fac*ustar(iq)**2
    tkv = theta(iq,k)*(1. + 0.61*qg(iq,k))
    rino(iq,k) = grav*(tkv - thvref(iq))*(zg(iq,k)-zg(iq,1))/max(thvref(iq)*vvk,tiny)
    if(rino(iq,k)>=ricr.and.iflag(iq)==0)then
      pblh(iq) = zg(iq,k-1) + (ricr - rino(iq,k-1))/(rino(iq,k) - rino(iq,k-1))*(zg(iq,k) - zg(iq,k-1))
      iflag(iq)=1
    endif  ! (rino(iq,k)>=ricr.and.iflag(iq)==0)
  enddo  ! iq loop
enddo   ! k loop
if(nmaxpr==1.and.mydiag.and.ntiles==1)then
  write (6,"('zg',9f8.1)") zg(idjd,1:kmax)
  write (6,"('rino_pa',9f8.3)") rino(idjd,1:kmax)
endif

! Set pbl height to maximum value where computation exceeds number of
! layers allowed
 
! Improve estimate of pbl height for the unstable points.
! Find unstable points (virtual heat flux is positive):
!
phiminv=0. ! MJT bug fix for IBM compiler
tlv=thvref

do iq=1,imax
  if(heatv(iq)>0.)then  ! unstable case
    phiminv(iq) = (1. - binm*pblh(iq)/obklen(iq))**(1./3.)
    wm(iq)= ustar(iq)*phiminv(iq)
    ! therm: 2nd term in eq. (4.d.19):
    ! temperature excess due to convective thermal
    therm = heatv(iq)*fak/wm(iq)
    ! eq. (4.d.19) : tlv then used in eq. (4.d.18) to improve pblh
    tlv(iq) = thvref(iq) + therm
  end if  
end do

! Improve pblh estimate for unstable conditions using the
! convective temperature excess:

do k=2,kmax
  do iq=1,imax
    vvk = (uav(iq,k) - uav(iq,1))**2 + (vav(iq,k) - vav(iq,1))**2 + fac*ustar(iq)**2
    vvk = max(vvk,tiny)
    tkv = theta(iq,k)*(1. + 0.61*qg(iq,k))
    rino(iq,k) = grav*(tkv - tlv(iq))*(zg(iq,k)-zg(iq,1))/max(thvref(iq)*vvk,tiny)     ! (see (4.d.18)
  enddo  !  i loop
enddo   !  k loop

iflag(:)=0
do k=2,kmax
  do iq=1,imax
    if(heatv(iq)>0..and.iflag(iq)==0)then  ! unstable case
      pblh(iq) = zg(iq,kl)    ! large default for unstable case
      if(rino(iq,k)>=ricr)then
        pblh(iq) = zg(iq,k-1) + (ricr - rino(iq,k-1))/(rino(iq,k) - rino(iq,k-1))*(zg(iq,k) - zg(iq,k-1))
        iflag(iq)=1  ! i.e. found it
      endif  ! (rino(iq,k)>=ricr)
    endif    ! (heatv(iq)>0..and.iflag(iq)==0)
  enddo     ! i loop
enddo      ! k loop

! Points for which pblh exceeds number of pbl layers allowed;
! set to maximum
 
! PBL height must be greater than some minimum mechanical mixing depth
! Several investigators have proposed minimum mechanical mixing depth
! relationships as a function of the local friction velocity, u*.  We 
! make use of a linear relationship of the form h = c u* where c=700.
! The scaling arguments that give rise to this relationship most often 
! represent the coefficient c as some constant over the local coriolis
! parameter.  Here we make use of the experimental results of Koracin 
! and Berkowicz (1988) [BLM, Vol 43] for which they recommend 0.07/f
! where f was evaluated at 39.5 N and 52 N.  Thus we use a typical mid
! latitude value for f so that c = 0.07/f = 700.
 
if(npblmin==1)pblh(:) = max(pblh(:),min(200.,700.*ustar(:)))
if(npblmin==2)pblh(:) = max(pblh(:),.07*ustar(:)/max(.5e-4,abs(f(:))))
if(npblmin==3)pblh(:) = max(pblh(:),.07*ustar(:)/max(1.e-4,abs(f(:)))) ! to ~agree 39.5N
if(npblmin==4)pblh(:) = max(pblh(1:imax),50.)

! pblh is now available; do preparation for diffusivity calculation:

! Do additional preparation for unstable cases only, set temperature
! and moisture perturbations depending on stability.

do iq=1,imax
  if(heatv(iq)>0.)then  ! unstable case
    phiminv(iq) =     (1. - binm*pblh(iq)/obklen(iq))**(1./3.)
    phihinv(iq) = sqrt(1. - binh*pblh(iq)/obklen(iq))
    wm(iq)      = ustar(iq)*phiminv(iq)
    wstr(iq)    = (heatv(iq)*grav*pblh(iq)/thvref(iq))**(1./3.)
  end if
end do

! Main level loop to compute the diffusivities and 
! counter-gradient terms:

if(nlocal==3)then
  ! suppress nonlocal scheme over the sea   jlm
  do iq=1,imax
    if(.not.land(iq))pblh(iq)=0.
  enddo
endif  !  (nlocal==3)

if(nlocal==4)then
  ! suppress nonlocal scheme for column if cloudy layer in pbl   jlm
  do k=1,kl/2
    do iq=1,imax
      if(zg(iq,k)<pblh(iq).and.cfrac(iq,k)>0.)pblh(iq)=0.
    enddo
  enddo
endif  !  (nlocal==4)

if(nlocal==5)then
  ! suppress nonlocal scheme for column if cloudy layer in pbl   jlm
  ! restores pblh at the bottom to have it available in convjlm/vertmix
  do k=1,kl/2
    do iq=1,imax
      if(zg(iq,k)<pblh(iq).and.cfrac(iq,k)>0.) pblh(iq)=-pblh(iq)  
    enddo
  enddo
endif  !  (nlocal==5)

if(nlocal==2)then
  do k=1,kmax-1        
    ! suppress nonlocal scheme if cloudy layers in pbl   jlm
    ! note this allows layers below to be done as nonlocal
    do iq=1,imax
      if(zg(iq,k)<pblh(iq).and.cfrac(iq,k)>0.)pblh(iq)=0.
    enddo
  end do
endif  !  (nlocal==2)

! Find levels within boundary layer:
! This is where Kh is at half model levels 
! zmzp = 0.5*(zm + zp)

do k=1,kmax-1
  do iq=1,imax
    fak1 = ustar(iq)*pblh(iq)*vk
    zm = zg(iq,k)
    zp = zg(iq,k+1)
    if (zm < pblh(iq)) then
      zmzp = 0.5*(zm + zp)
      zh = zmzp/pblh(iq)
      zl = zmzp/obklen(iq)
      zzh= 0.
      if (zh<=1.0) zzh = (1. - zh)**2

! stblev for points zm < plbh and stable and neutral
! unslev for points zm < plbh and unstable
      if(heatv(iq)>0.)then  ! unstable case
        fak2   = wm(iq)*pblh(iq)*vk
! unssrf, unstable within surface layer of pbl
! Unstable for surface layer; counter-gradient terms zero
        if (zh<sffrac) then
          term =     (1. - betam*zl)**(1./3.)
          pblk = fak1*zh*zzh*term
          pr = term/sqrt(1. - betah*zl)
        else
! unsout, unstable within outer   layer of pbl
! Unstable for outer layer; counter-gradient terms non-zero:
          pblk = fak2*zh*zzh
          fak3 = fakn*wstr(iq)/wm(iq)
          cgs     = fak3/(pblh(iq)*wm(iq))
          cgh(iq,k) = rkhfs(iq)*cgs                 !eq. (4.d.17)
          cgq(iq,k) = rkqfs(iq)*cgs                 !eq. (4.d.17)
          pr = phiminv(iq)/phihinv(iq) + ccon*fak3/fak
        end if
        rkm(iq,k) = max(pblk,rkm(iq,k))
        rkh(iq,k) = max(pblk/pr,rkh(iq,k))
      elseif(nlocal>0)then    ! following are stable or neutral
! Stable and neutral points; set diffusivities; counter-gradient
! terms zero for stable case:
! term: pblk is Kc in eq. (4.d.16)
! but reverts to Louis stable treatment for nlocal=-1
        if (zl<=1.) then   ! 0 < z/L < 1.
          pblk = fak1*zh*zzh/(1. + betas*zl)
        else
          pblk = fak1*zh*zzh/(betas + zl)
        endif
        if(nrkmin==2)rkmin=vk*ustar(iq)*zmzp*zzh
        if(nrkmin==3)rkmin=max(rkh(iq,k),vk*ustar(iq)*zmzp*zzh)
        if(nrkmin==1.or.nlocal==6)rkmin=rkh(iq,k)
        if(ntest==1.and.mydiag.and.ntiles==1)then
          if(iq==idjd)then
            write(6,*) 'in pbldif k,ustar,zmzp,zh,zl,zzh ',k,ustar(iq),zmzp,zh,zl,zzh
            write(6,*) 'rkh_L,rkmin,pblk,fak1,pblh ',rkh(iq,k),rkmin,pblk,fak1,pblh(iq)
          endif  ! (iq==idjd)
        endif    ! (ntest==1)
        rkm(iq,k) = max(pblk,rkmin)        
        rkh(iq,k) = rkm(iq,k)
      endif      ! (heatv(iq)>0.)    unstbl(i)
    endif        ! zm < pblh(iq)
  enddo         ! iq=1,imax
enddo             !end of k loop
if(diag.and.mydiag.and.ntiles==1)then
  if(heatv(idjd)>0.) write (6,"('rino_pb',9f8.3)") rino(idjd,1:kmax) ! not meaningful or used otherwise
  write (6,*) 'ricr,obklen,heatv,pblh ',ricr,obklen(idjd),heatv(idjd),pblh(idjd)
  write (6,"('rkh_p',9f9.3/5x,9f9.3)") rkh(idjd,1:kl-2)
endif


ztodtgor = dtin*grav/rdry
!     update theta and qtg due to counter gradient
do k=2,kmax-1
  do iq=1,imax
    delsig = sigmh(k+1)-sigmh(k)
    tmp1 = ztodtgor/delsig
    sigotbk=sigmh(k+1)/(0.5*(t(iq,k+1) + t(iq,k)))
    sigotbkm1=sigmh(k)/(0.5*(t(iq,k-1) + t(iq,k)))
    theta(iq,k) = theta(iq,k) + tmp1*(sigotbk*rkh(iq,k)*cgh(iq,k) - sigotbkm1*rkh(iq,k-1)*cgh(iq,k-1))
    qg(iq,k) = qg(iq,k) + tmp1*(sigotbk*rkh(iq,k)*cgq(iq,k) - sigotbkm1*rkh(iq,k-1)*cgq(iq,k-1))
#ifdef scm
    wth_flux(iq,k+1) = wth_flux(iq,k+1) + tmp1*sigotbk*rkh(iq,k)*cgh(iq,k)/dtin
    wq_flux(iq,k+1) = wq_flux(iq,k+1) + tmp1*sigotbk*rkh(iq,k)*cgq(iq,k)/dtin
#endif
  end do
end do
k=1
do iq=1,imax
  delsig = sigmh(k+1)-sigmh(k)
  tmp1 = ztodtgor/delsig
  sigotbk=sigmh(k+1)/(0.5*(t(iq,k+1) + t(iq,k)))
  theta(iq,k) = theta(iq,k) + tmp1*sigotbk*rkh(iq,k)*cgh(iq,k)
  qg(iq,k) = qg(iq,k) + tmp1*sigotbk*rkh(iq,k)*cgq(iq,k)
#ifdef scm
  wth_flux(iq,k+1) = wth_flux(iq,k+1) + tmp1*sigotbk*rkh(iq,k)*cgh(iq,k)/dtin
  wq_flux(iq,k+1) = wq_flux(iq,k+1) + tmp1*sigotbk*rkh(iq,k)*cgq(iq,k)/dtin
#endif
end do

if (ntest>0.and.mydiag.and.ntiles==1) then
  write(6,*) 'pbldif'
  write(6,*) 'rkh= ',(rkh(idjd,k),k=1,kl)
  write(6,*) 'theta= ',(theta(idjd,k),k=1,kl)
  write(6,*) 'qg= ',(qg(idjd,k),k=1,kl)
  write(6,*) 'cgh= ',(cgh(idjd,k),k=1,kl)
  write(6,*) 'cgq= ',(cgq(idjd,k),k=1,kl)
endif
!
!    Check for neg qtg's and put the original vertical
!    profile back if a neg value is found. A neg value implies that the
!    quasi-equilibrium conditions assumed for the countergradient term are
!    strongly violated.
!    Original code rewritten by Rosinski 7/8/91 to vectorize in longitude.

if(nlocal==5)then
  ! restoring pblh to have it available in convjlm/vertmix jlm
  pblh(:)=abs(pblh(:))
endif  !  (nlocal==5)

return
end subroutine pbldif

subroutine trim(a,c,rhs,imax,kl)
!!$acc routine vector

implicit none

integer, intent(in) :: imax, kl
integer k, iq
real, dimension(imax,kl), intent(in) :: a, c
real, dimension(imax,kl), intent(inout) :: rhs
real, dimension(imax,kl) :: e, g
real temp, b

! this routine solves the system
!   a(k)*u(k-1)+b(k)*u(k)+c(k)*u(k+1)=rhs(k)    for k=2,kl-1
!   with  b(k)*u(k)+c(k)*u(k+1)=rhs(k)          for k=1
!   and   a(k)*u(k-1)+b(k)*u(k)=rhs(k)          for k=kl

! the Thomas algorithm is used
do concurrent (iq = 1:imax)
  b=1.-a(iq,1)-c(iq,1)
  e(iq,1)=c(iq,1)/b
  g(iq,1)=rhs(iq,1)/b
end do
do k = 2,kl-1
  do concurrent (iq = 1:imax)
    b=1.-a(iq,k)-c(iq,k)
    temp= 1./(b-a(iq,k)*e(iq,k-1))
    e(iq,k)=c(iq,k)*temp
    g(iq,k)=(rhs(iq,k)-a(iq,k)*g(iq,k-1))*temp
  end do
end do

! do back substitution to give answer now
do concurrent (iq = 1:imax)
  b=1.-a(iq,kl)-c(iq,kl)
  rhs(iq,kl)=(rhs(iq,kl)-a(iq,kl)*g(iq,kl-1))/(b-a(iq,kl)*e(iq,kl-1))
end do
do k = kl-1,1,-1
  do concurrent (iq = 1:imax)
    rhs(iq,k) = g(iq,k)-e(iq,k)*rhs(iq,k+1)
  end do
end do

return
end subroutine trim

subroutine tkeeps_work(t,em,tss,eg,fg,ps,qg,qfg,qlg,stratcloud,                                 &
                       xtg,cduv,u,v,pblh,ustar,tke,eps,shear,at,ct,                             &
#ifdef scm
                       wth_flux,wq_flux,uw_flux,vw_flux,mfsave,tkesave,epssave,rkmsave,rkhsave, &
                       buoyproduction,shearproduction,totaltransport,                           &
#endif
                       imax,kl,naero)
!!$acc routine vector

use const_phys                   ! Physical constants
use parm_m, only : ds, nlocal, iaero, dt, qgmin, cqmix
                                 ! Model configuration
use sigs_m                       ! Atmosphere sigma levels
use tkeeps, only : tkemix        ! TKE-EPS boundary layer

implicit none

integer, intent(in) :: imax, kl, naero
integer k, nt, iq
real, dimension(imax,kl,naero), intent(inout) :: xtg
real, dimension(imax,kl), intent(inout) :: t, qg, qfg, qlg
real, dimension(imax,kl), intent(inout) :: stratcloud, u, v
real, dimension(imax,kl), intent(inout) :: tke, eps
real, dimension(imax,kl), intent(out) :: at, ct
real, dimension(imax,kl), intent(in) :: shear
real, dimension(imax,kl) :: zh
real, dimension(imax,kl) :: rhs, zg
real, dimension(imax,kl) :: rkm
real tmnht
real, dimension(imax), intent(inout) :: pblh, ustar
real, dimension(imax), intent(in) :: em, tss, eg, fg, ps
real, dimension(imax), intent(in) :: cduv
real dz, gt
real, dimension(imax) :: rhos, dx
real, dimension(kl) :: sigkap, delh
real rong, rlogs1, rlogs2, rlogh1, rlog12
#ifdef scm
real, dimension(imax,kl), intent(inout) :: wth_flux, wq_flux, uw_flux
real, dimension(imax,kl), intent(inout) :: vw_flux, tkesave, epssave
real, dimension(imax,kl), intent(out) :: rkmsave, rkhsave
real, dimension(imax,kl), intent(inout) :: buoyproduction, shearproduction
real, dimension(imax,kl), intent(inout) :: totaltransport
real, dimension(imax,kl-1), intent(inout) :: mfsave
#endif


! estimate grid spacing for scale aware MF
dx(:) = ds/em

! Set-up potential temperature transforms
rong = rdry/grav
do k = 1,kl
  delh(k)   = -rong*dsig(k)/sig(k)  ! sign of delh defined so always +ve
  sigkap(k) = sig(k)**(-rdry/cp)
end do      ! k loop
      
! note ksc/=0 options are clobbered when nvmix=6
! However, nvmix=6 with nlocal=7 supports its own shallow
! convection options

! calculate height on full levels
zg(:,1) = bet(1)*t(:,1)/grav
zh(:,1) = t(:,1)*delh(1)
do k = 2,kl-1
  zg(:,k) = zg(:,k-1) + (bet(k)*t(:,k)+betm(k)*t(:,k-1))/grav
  zh(:,k) = zh(:,k-1) + t(:,k)*delh(k)
end do ! k  loop
zg(:,kl) = zg(:,kl-1) + (bet(kl)*t(:,kl)+betm(kl)*t(:,kl-1))/grav
zh(:,kl) = zh(:,kl-1) + t(:,kl)*delh(kl)
       
! near surface air density (see sflux.f and cable_ccam2.f90)
rhos(:) = ps(:)/(rdry*tss(:))
    
! transform to ocean reference frame and temp to theta
do k = 1,kl
  rhs(:,k) = t(:,k)*sigkap(k) ! theta
end do

#ifdef scm
! Initial flux for SCM to be added up below.
do k = 1,kl
  wth_flux(:,k) = 0.
  wq_flux(:,k) = 0.
  uw_flux(:,k) = 0.
  vw_flux(:,k) = 0.
  mfsave(:,k) = 0.
  tkesave(:,k) = -1. ! missing value
end do
! Evaluate EDMF scheme
select case(nlocal)
  case(0) ! no counter gradient
    call tkemix(rkm,rhs,qg,qlg,qfg,stratcloud,u,v,pblh,fg,eg,cduv,ps,zg,zh,sig,rhos, &
                ustar,dt,qgmin,1,tke,eps,shear,dx,                                   &
                wth_flux,wq_flux,uw_flux,vw_flux,mfsave,buoyproduction,              &
                shearproduction,totaltransport,imax,kl)
  case(7) ! mass-flux counter gradient
    call tkemix(rkm,rhs,qg,qlg,qfg,stratcloud,u,v,pblh,fg,eg,cduv,ps,zg,zh,sig,rhos, &
                ustar,dt,qgmin,0,tke,eps,shear,dx,                                   &
                wth_flux,wq_flux,uw_flux,vw_flux,mfsave,buoyproduction,              &
                shearproduction,totaltransport,imax,kl)
    
end select
do k = 1,kl
  ! save Km and Kh for output
  rkmsave(:,k) = rkm(:,k)
  rkhsave(:,k) = rkm(:,k)
  tkesave(:,k) = tke(:,k)
  epssave(:,k) = eps(:,k)
end do
#else
! Evaluate EDMF scheme
select case(nlocal)
  case(0) ! no counter gradient
    call tkemix(rkm,rhs,qg,qlg,qfg,stratcloud,u,v,pblh,fg,eg,cduv,ps,zg,zh,sig,rhos, &
                ustar,dt,qgmin,1,tke,eps,shear,dx,imax,kl) 
  case(7) ! mass-flux counter gradient
    call tkemix(rkm,rhs,qg,qlg,qfg,stratcloud,u,v,pblh,fg,eg,cduv,ps,zg,zh,sig,rhos, &
                ustar,dt,qgmin,0,tke,eps,shear,dx,imax,kl)     
end select
#endif

! replace counter gradient term
do k = 1,kl
  rkm(:,k) = rkm(:,k)*cqmix
end do

! tracers
at(:,1) = 0.
ct(:,kl) = 0.
rlogs1=log(sig(1))
rlogs2=log(sig(2))
rlogh1=log(sigmh(2))
rlog12=1./(rlogs1-rlogs2)
do concurrent (iq = 1:imax)
  tmnht=(t(iq,2)*rlogs1-t(iq,1)*rlogs2+(t(iq,1)-t(iq,2))*rlogh1)*rlog12  
  dz = -tmnht*rong*((sig(2)-sig(1))/sigmh(2))  ! this is z(k+1)-z(k)
  gt = rkm(iq,1)*dt*(sig(2)-sig(1))/(dz**2)
  at(iq,2) =-gt/dsig(2)  
  ct(iq,1) = -gt/dsig(1)
end do
do concurrent (k = 2:kl-1)
  do concurrent (iq = 1:imax)
    ! Calculate half level heights and temperatures
    ! n.b. an approximate zh (in m) is quite adequate for this routine
    tmnht = ratha(k)*t(iq,k+1) + rathb(k)*t(iq,k)
    dz = -tmnht*rong*((sig(k+1)-sig(k))/sigmh(k+1))  ! this is z(k+1)-z(k)
    gt = rkm(iq,k)*dt*(sig(k+1)-sig(k))/(dz**2)
    at(iq,k+1) =-gt/dsig(k+1)  
    ct(iq,k) = -gt/dsig(k)
  end do
end do
  
! Aerosols
if ( abs(iaero)>=2 ) then
  do nt = 1,naero
    call trim(at,ct,xtg(:,:,nt),imax,kl)
  end do
end if ! (abs(iaero)>=2)  

! transform winds back to Earth reference frame and theta to temp
do k = 1,kl
  t(:,k) = rhs(:,k)/sigkap(k)
enddo    !  k loop
      
return
end subroutine tkeeps_work


end module vertmix_m
