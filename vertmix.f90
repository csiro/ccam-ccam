
! Conformal Cubic Atmospheric Model
    
! Copyright 2015-2016 Commonwealth Scientific Industrial Research Organisation (CSIRO)
    
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
module  vertmix_m
implicit none

private

public vertmix,vertmix_init

integer, save :: nb,imax
real, dimension(:), allocatable, save :: prcpv
integer, save :: kscbase,ksctop

contains

subroutine vertmix_init(ifull,kl,nbin)
use cc_mpi, only : myid
use const_phys, only : roncp        ! Physical constants
use sigs_m, only : sig              ! Atmosphere sigma levels
use vertmixdata_m
use tracers_m, only : ntracmax
use tracermodule, only : numtracer
use aerosolldr, only : naero

implicit none
include 'kuocom.h'                  ! Convection parameters
integer, intent(in) :: ifull,kl,nbin
integer :: k

nb=nbin
imax=ifull/nb

allocate(prcpv(kl))

if(ksc/=0)then
  ! set ksctop for shallow convection
  ksctop=1    ! ksctop will be first level below sigkcst
  do while(sig(ksctop+1)>sigksct)  !  e.g. sigksct=.75
    ksctop=ksctop+1
  enddo
  kscbase=1  ! kscbase will be first level above sigkcsb
  do while(sig(kscbase)>sigkscb.and.sigkscb>0.) ! e.g. sigkscb=.99
    kscbase=kscbase+1
  enddo
  if ( myid == 0 ) then
    write(6,*)'For shallow convection:'
    write(6,*)'ksc,kscbase,ksctop,kscsea ', ksc,kscbase,ksctop,kscsea
    write (6,"(' sigkscb,sigksct,tied_con,tied_over,tied_rh:',5f8.3)")sigkscb,sigksct,tied_con,tied_over,tied_rh
  end if
  do k=1,kl
    prcpv(k)=sig(k)**(-roncp)
  enddo  ! k loop
endif    ! (ktau==1.and.ksc/=0)

call vertmixdata_init(nb,imax,kl,ntracmax,numtracer,naero)

end subroutine vertmix_init

subroutine vertmix
use cc_mpi, only : start_log,end_log,vertmx_begin,vertmx_end,myid
use aerosolldr, only : naero,xtg                                                                                              ! LDR prognostic aerosols
use arrays_m, only : t,u,v,qg,ps                                                                                              ! Atmosphere dyamics prognostic arrays
use cfrac_m, only : cfrac                                                                                                     ! Cloud fraction
use cloudmod, only : stratcloud,combinecloudfrac,convectivecloudfrac                                                          ! Prognostic strat cloud
use const_phys, only : rdry,grav,roncp,cp,hl                                                                                  ! Physical constants
use diag_m, only : maxmin,printa                                                                                              ! Diagnostic routines
use extraout_m, only : ustar                      ! Additional diagnostics
use kuocomb_m, only : kbsav,ktsav,convpsav                                                                                    ! JLM convection
use liqwpar_m, only : qfg,qlg                                                                                                 ! Cloud water mixing ratios
use map_m, only : em                                                                                                          ! Grid map arrays
use mlo, only : mloexport,mloexpice                                                                                           ! Ocean physics and prognostic arrays
use morepbl_m, only : eg,fg,pblh,condc                                                                                        ! Additional boundary layer diagnostics
use newmpar_m, only : il,kl                                                                                             ! Grid parameters
use nharrs_m, only : phi_nh                                                                                                   ! Non-hydrostatic atmosphere arrays
use parm_m, only : cgmap_offset,ds,cgmap_scale,diag,ktau,idjd,nmlo,nvmix,dt,nlv,ia,ib,ja,jb,nmaxpr,iaero,nlocal,qgmin,av_vmod ! Model configuration
use pbl_m, only : tss,cdtq,cduv                                                                                               ! Boundary layer arrays
use savuvt_m, only : savu,savv                                                                                                ! Saved dynamic arrays
use sigs_m, only : bet,betm,sigmh,sig,dsig,ratha,rathb                                                                        ! Atmosphere sigma levels
use soilsnow_m, only : fracice                                                                                                ! Soil, snow and surface data
use tkeeps, only : cq,tkemix,tke,eps,shear                                                                                                  ! TKE-EPS boundary layer
#ifndef scm
use tracers_m, only : ngas,tr                                                                                                    ! Tracer data
use trvmix, only : tracervmix                                                                                                 ! Tracer mixing routines
#endif
use work2_m, only : zo                                                                                                        ! Diagnostic arrays
use soil_m, only : land
use kuocomb_m, only : kbsav,ktsav
use vertmixdata_m
use carbpools_m, only : fnee,fpn,frp,frs
use nsibd_m, only : ivegt
use xyzinfo_m, only : wts
use tracermodule, only : mcfdep,co2em,oh,strloss,jmcf
use map_m, only : f
use screen_m, only : tscrn,qgscrn

implicit none
integer :: i, is, ie

do i=1,nb
  is=(i-1)*imax+1
  ie=i*imax

  if (allocated(land)) b_land(i)%data=land(is:ie)

  if (allocated(kbsav)) b_kbsav(i)%data=kbsav(is:ie)
  if (allocated(ktsav)) b_ktsav(i)%data=ktsav(is:ie)

  if (allocated(ps)) b_ps(i)%data=ps(is:ie)
  if (allocated(cdtq)) b_cdtq(i)%data=cdtq(is:ie)
  if (allocated(fnee)) b_fnee(i)%data=fnee(is:ie)
  if (allocated(fpn)) b_fpn(i)%data=fpn(is:ie)
  if (allocated(frp)) b_frp(i)%data=frp(is:ie)
  if (allocated(frs)) b_frs(i)%data=frs(is:ie)
  if (allocated(ivegt)) b_ivegt(i)%data=ivegt(is:ie)
  if (allocated(wts)) b_wts(i)%data=wts(is:ie)
  if (allocated(mcfdep)) b_mcfdep(i)%data=mcfdep(is:ie)
  if (allocated(em)) b_em(i)%data=em(is:ie)
  if (allocated(fracice)) b_fracice(i)%data=fracice(is:ie)
  if (allocated(tss)) b_tss(i)%data=tss(is:ie)
  if (allocated(eg)) b_eg(i)%data=eg(is:ie)
  if (allocated(fg)) b_fg(i)%data=fg(is:ie)
  if (allocated(convpsav)) b_convpsav(i)%data=convpsav(is:ie)
  if (allocated(ustar)) b_ustar(i)%data=ustar(is:ie)
  if (allocated(pblh)) b_pblh(i)%data=pblh(is:ie)
  if (allocated(f)) b_f(i)%data=f(is:ie)
  if (allocated(condc)) b_condc(i)%data=condc(is:ie)
  if (allocated(cduv)) b_cduv(i)%data=cduv(is:ie)
  if (allocated(zo)) b_zo(i)%data=zo(is:ie)
  if (allocated(tscrn)) b_tscrn(i)%data=tscrn(is:ie)
  if (allocated(qgscrn)) b_qgscrn(i)%data=qgscrn(is:ie)

  if (allocated(stratcloud)) b_stratcloud(i)%data=stratcloud(is:ie,:)
  if (allocated(tke)) b_tke(i)%data=tke(is:ie,:)
  if (allocated(eps)) b_eps(i)%data=eps(is:ie,:)
  if (allocated(shear)) b_shear(i)%data=shear(is:ie,:)
  if (allocated(phi_nh)) b_phi_nh(i)%data=phi_nh(is:ie,:)
  if (allocated(t)) b_t(i)%data=t(is:ie,:)
  if (allocated(co2em)) b_co2em(i)%data=co2em(is:ie,:)
  if (allocated(oh)) b_oh(i)%data=oh(is:ie,:)
  if (allocated(strloss)) b_strloss(i)%data=strloss(is:ie,:)
  if (allocated(jmcf)) b_jmcf(i)%data=jmcf(is:ie,:)
  if (allocated(qg)) b_qg(i)%data=qg(is:ie,:)
  if (allocated(qfg)) b_qfg(i)%data=qfg(is:ie,:)
  if (allocated(qlg)) b_qlg(i)%data=qlg(is:ie,:)
  if (allocated(cfrac)) b_cfrac(i)%data=cfrac(is:ie,:)
  if (allocated(u)) b_u(i)%data=u(is:ie,:)
  if (allocated(v)) b_v(i)%data=v(is:ie,:)
  if (allocated(savu)) b_savu(i)%data=savu(is:ie,:)
  if (allocated(savv)) b_savv(i)%data=savv(is:ie,:)

  if (allocated(tr)) b_tr(i)%data=tr(is:ie,:,:)
  if (allocated(xtg)) b_xtg(i)%data=xtg(is:ie,:,:)

end do

call start_log(vertmx_begin)
!$omp parallel do
do i=1,nb
  call vertmix_work(i,imax)
end do
call end_log(vertmx_end)

do i=1,nb
  is=(i-1)*imax+1
  ie=i*imax

  if (allocated(land)) land(is:ie)=b_land(i)%data

  if (allocated(kbsav)) kbsav(is:ie)=b_kbsav(i)%data
  if (allocated(ktsav)) ktsav(is:ie)=b_ktsav(i)%data

  if (allocated(ps)) ps(is:ie)=b_ps(i)%data
  if (allocated(cdtq))cdtq(is:ie)=b_cdtq(i)%data
  if (allocated(fnee)) fnee(is:ie)=b_fnee(i)%data
  if (allocated(fpn)) fpn(is:ie)=b_fpn(i)%data
  if (allocated(frp)) frp(is:ie)=b_frp(i)%data
  if (allocated(frs)) frs(is:ie)=b_frs(i)%data
  if (allocated(ivegt)) ivegt(is:ie)=b_ivegt(i)%data
  if (allocated(wts)) wts(is:ie)=b_wts(i)%data
  if (allocated(mcfdep)) mcfdep(is:ie)=b_mcfdep(i)%data
  if (allocated(em)) em(is:ie)=b_em(i)%data
  if (allocated(fracice)) fracice(is:ie)=b_fracice(i)%data
  if (allocated(tss)) tss(is:ie)=b_tss(i)%data
  if (allocated(eg)) eg(is:ie)=b_eg(i)%data
  if (allocated(fg)) fg(is:ie)=b_fg(i)%data
  if (allocated(convpsav)) convpsav(is:ie)=b_convpsav(i)%data
  if (allocated(ustar)) ustar(is:ie)=b_ustar(i)%data
  if (allocated(pblh)) pblh(is:ie)=b_pblh(i)%data
  if (allocated(f)) f(is:ie)=b_f(i)%data
  if (allocated(condc)) condc(is:ie)=b_condc(i)%data
  if (allocated(cduv)) cduv(is:ie)=b_cduv(i)%data
  if (allocated(zo)) zo(is:ie)=b_zo(i)%data
  if (allocated(tscrn)) tscrn(is:ie)=b_tscrn(i)%data
  if (allocated(qgscrn)) qgscrn(is:ie)=b_qgscrn(i)%data

  if (allocated(stratcloud)) stratcloud(is:ie,:)=b_stratcloud(i)%data
  if (allocated(tke)) tke(is:ie,:)=b_tke(i)%data
  if (allocated(eps)) eps(is:ie,:)=b_eps(i)%data
  if (allocated(shear)) shear(is:ie,:)=b_shear(i)%data
  if (allocated(phi_nh)) phi_nh(is:ie,:)=b_phi_nh(i)%data
  if (allocated(t)) t(is:ie,:)=b_t(i)%data
  if (allocated(co2em)) co2em(is:ie,:)=b_co2em(i)%data
  if (allocated(oh)) oh(is:ie,:)=b_oh(i)%data
  if (allocated(strloss)) strloss(is:ie,:)=b_strloss(i)%data
  if (allocated(jmcf)) jmcf(is:ie,:)=b_jmcf(i)%data
  if (allocated(qg)) qg(is:ie,:)=b_qg(i)%data
  if (allocated(qfg)) qfg(is:ie,:)=b_qfg(i)%data
  if (allocated(qlg)) qlg(is:ie,:)=b_qlg(i)%data
  if (allocated(cfrac)) cfrac(is:ie,:)=b_cfrac(i)%data
  if (allocated(u)) u(is:ie,:)=b_u(i)%data
  if (allocated(v)) v(is:ie,:)=b_v(i)%data
  if (allocated(savu)) savu(is:ie,:)=b_savu(i)%data
  if (allocated(savv)) savv(is:ie,:)=b_savv(i)%data

  if (allocated(tr)) tr(is:ie,:,:)=b_tr(i)%data
  if (allocated(xtg)) xtg(is:ie,:,:)=b_xtg(i)%data

end do

end subroutine vertmix

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

!--------------------------------------------------------------
! Control subroutine for vertical mixing
subroutine vertmix_work(tile,imax)

use aerosolldr, only : naero,xtg                                                                                              ! LDR prognostic aerosols
use arrays_m, only : t,u,v,qg,ps                                                                                              ! Atmosphere dyamics prognostic arrays
use cc_mpi, only : ccmpi_barrier,comm_world,mydiag,ccmpi_abort                                                                ! CC MPI routines
use cfrac_m, only : cfrac                                                                                                     ! Cloud fraction
use cloudmod, only : stratcloud,combinecloudfrac,convectivecloudfrac                                                          ! Prognostic strat cloud
use const_phys, only : rdry,grav,roncp,cp,hl                                                                                  ! Physical constants
use diag_m, only : maxmin,printa                                                                                              ! Diagnostic routines
!use extraout_m                      ! Additional diagnostics
use kuocomb_m, only : kbsav,ktsav,convpsav                                                                                    ! JLM convection
use liqwpar_m, only : qfg,qlg                                                                                                 ! Cloud water mixing ratios
use map_m, only : em                                                                                                          ! Grid map arrays
use mlo, only : mloexport,mloexpice                                                                                           ! Ocean physics and prognostic arrays
use morepbl_m, only : eg,fg,pblh,condc                                                                                        ! Additional boundary layer diagnostics
use newmpar_m, only : il,kl                                                                                             ! Grid parameters
use nharrs_m, only : phi_nh                                                                                                   ! Non-hydrostatic atmosphere arrays
use parm_m, only : cgmap_offset,ds,cgmap_scale,diag,ktau,idjd,nmlo,nvmix,dt,nlv,ia,ib,ja,jb,nmaxpr,iaero,nlocal,qgmin,av_vmod ! Model configuration
use pbl_m, only : tss,cdtq,cduv                                                                                               ! Boundary layer arrays
use savuvt_m, only : savu,savv                                                                                                ! Saved dynamic arrays
use sigs_m, only : bet,betm,sigmh,sig,dsig,ratha,rathb                                                                        ! Atmosphere sigma levels
use soilsnow_m, only : fracice                                                                                                ! Soil, snow and surface data
use tkeeps, only : cq,tkemix                                                                                                  ! TKE-EPS boundary layer
#ifndef scm
use tracers_m, only : ngas                                                                                                    ! Tracer data
use trvmix, only : tracervmix                                                                                                 ! Tracer mixing routines
#endif
use work2_m, only : zo                                                                                                        ! Diagnostic arrays
use vertmixdata_m
      
implicit none
      
include 'kuocom.h'                  ! Convection parameters

integer, intent(in) :: tile,imax
integer, parameter :: ntest = 0
integer k, tnaero, nt
real rong, rlogs1, rlogs2, rlogh1, rlog12
real delsig, conflux, condrag
real, dimension(imax,kl) :: cnhs_fl, zh, clcon
real, dimension(imax,kl) :: rhs, guv, gt
real, dimension(imax,kl) :: at, ct, au, cu, zg, cldtmp
real, dimension(imax,kl) :: uav, vav
real, dimension(imax,kl) :: rkm, rkh
real, dimension(imax,kl-1) :: tmnht, cnhs_hl
real, dimension(imax) :: ou, ov, iu, iv
real, dimension(imax) :: dz, dzr
real, dimension(imax) :: cgmap, tnhs_fl
real, dimension(imax) :: rhos
real, dimension(kl) :: sighkap,sigkap,delons,delh
#ifdef scm
real, dimension(imax,kl) :: mfout
#endif
integer :: is, ie

is=(tile-1)*imax+1
ie=tile*imax

! Non-hydrostatic terms
tnhs_fl(1:imax) = b_phi_nh(tile)%data(:,1)/bet(1)
cnhs_fl(1:imax,1) = max( 1. + tnhs_fl(:)/b_t(tile)%data(:,1), 0.001 )
do k = 2,kl
  ! representing non-hydrostatic term as a correction to air temperature
  tnhs_fl(1:imax) = (b_phi_nh(tile)%data(:,k)-b_phi_nh(tile)%data(:,k-1)-betm(k)*tnhs_fl(:))/bet(k)
  cnhs_fl(1:imax,k) = max( 1. + tnhs_fl(:)/b_t(tile)%data(:,k), 0.001 )
end do

! Weight as a function of grid spacing for turning off CG term
!cgmap = 0.982, 0.5, 0.018 for 1000m, 600m, 200m when cgmap_offset=600 and cgmap_scale=200.
if ( cgmap_offset>0. ) then
  cgmap(1:imax) = 0.5*tanh((ds/b_em(tile)%data(:)-cgmap_offset)/cgmap_scale) + 0.5 ! MJT suggestion
else
  cgmap(1:imax) = 1.
end if

#ifdef scm
! Initial flux to be added up below.
wth_flux(tile)%data(:,:) = 0.
wq_flux(tile)%data(:,:) = 0.
uw_flux(tile)%data(:,:) = 0.
vw_flux(tile)%data(:,:) = 0.
mfsave(tile)%data(:,:) = 0.
tkesave(tile)%data(:,:) = -1. ! missing value
#endif

! Set-up potential temperature transforms
rong = rdry/grav
do k = 1,kl-1
  sighkap(k)=sigmh(k+1)**(-roncp)
  delons(k) =rong *((sig(k+1)-sig(k))/sigmh(k+1))
end do      ! k loop
do k = 1,kl
  delh(k)  =-rong *dsig(k)/sig(k)  ! sign of delh defined so always +ve
  sigkap(k)=sig(k)**(-roncp)
end do      ! k loop

if ( diag .or. ntest>=1 ) then
  call maxmin(u,'%u',ktau,1.,kl)
  call maxmin(v,'%v',ktau,1.,kl)
  call maxmin(t,'%t',ktau,1.,kl)
  call maxmin(qg,'qg',ktau,1.e3,kl)
  call ccmpi_barrier(comm_world)  ! stop others going past    
  if ( mydiag ) then
    write(6,*)'sig ',sig
    write(6,*)'dsig ',dsig
    write(6,*)'delh ',delh
    write(6,*)'in vertmix'
    write (6,"('uin ',9f8.3/4x,9f8.3)") u(idjd,:) 
    write (6,"('vin ',9f8.3/4x,9f8.3)") v(idjd,:) 
  end if
end if

! Adjustment for moving ocean surface
ou=0.
ov=0.
if ( nmlo/=0 ) then
  iu=0.
  iv=0.
  call mloexport(2,ou,1,0,tile,imax)
  call mloexport(3,ov,1,0,tile,imax)
  call mloexpice(iu, 9,0,tile,imax)
  call mloexpice(iv,10,0,tile,imax)
  ou = (1.-b_fracice(tile)%data(:))*ou + b_fracice(tile)%data(:)*iu
  ov = (1.-b_fracice(tile)%data(:))*ov + b_fracice(tile)%data(:)*iv
end if
      
if ( nvmix/=6 ) then

  !--------------------------------------------------------------
  ! JLM's local Ri scheme

  ! Calculate half level heights, temperatures and NHS correction
  rlogs1=log(sig(1))
  rlogs2=log(sig(2))
  rlogh1=log(sigmh(2))
  rlog12=1./(rlogs1-rlogs2)
  tmnht(:,1)=(b_t(tile)%data(:,2)*rlogs1-b_t(tile)%data(:,1)*rlogs2+(b_t(tile)%data(:,1)-b_t(tile)%data(:,2))*rlogh1)*rlog12
  cnhs_hl(:,1)=(cnhs_fl(1:imax,2)*rlogs1-cnhs_fl(1:imax,1)*rlogs2+(cnhs_fl(1:imax,1)-cnhs_fl(1:imax,2))*rlogh1)*rlog12
  ! n.b. an approximate zh (in m) is quite adequate for this routine
  zh(:,1)   =b_t(tile)%data(:,1)*cnhs_fl(:,1)*delh(1)
  do k = 2,kl-1
    zh(:,k)   =zh(:,k-1)+b_t(tile)%data(:,k)*cnhs_fl(:,k)*delh(k)
    tmnht(:,k)=ratha(k)*b_t(tile)%data(:,k+1)+rathb(k)*b_t(tile)%data(:,k)
    ! non-hydrostatic temperature correction at half level height
    cnhs_hl(:,k) =ratha(k)*cnhs_fl(1:imax,k+1)+rathb(k)*cnhs_fl(1:imax,k)
  end do      !  k loop
  zh(:,kl)=zh(:,kl-1)+b_t(tile)%data(:,kl)*cnhs_fl(:,kl)*delh(kl)

  ! Calculate theta
  do k = 1,kl
    rhs(:,k)=b_t(tile)%data(:,k)*sigkap(k)  ! rhs is theta here
  enddo      !  k loop
    
  call vertjlm(rkm,rkh,rhs,sigkap,sighkap,delons,zh,tmnht,cnhs_hl,ntest,cgmap,tile,imax)

  do k = 1,kl-1
    delsig  =(sig(k+1)-sig(k))
    dz(:)   =-tmnht(:,k)*delons(k)  ! this is z(k+1)-z(k)
    dzr(:)  =1./dz(:)
    guv(:,k)=rkm(:,k)*dt*delsig*dzr(:)**2/cnhs_hl(:,k)
    gt(:,k) =rkh(:,k)*dt*delsig*dzr(:)**2/cnhs_hl(:,k)
  end do      ! k loop
  guv(:,kl)=0.
  gt(:,kl) =0.

  if ( diag ) then
    call maxmin(rkh,'rk',ktau,.01,kl-1)
    if ( mydiag ) then
      write(6,*)'vertmix guv ',(guv(idjd,k),k=1,kl)
      write(6,*)'vertmix gt ',(gt(idjd,k),k=1,kl)
    end if
    call printa('t   ',b_t(tile)%data(:,:),ktau,nlv,ia,ib,ja,jb,200.,1.)
    call printa('tss ',b_tss(tile)%data(:),ktau,nlv,ia,ib,ja,jb,200.,1.)
    call printa('eg  ',b_eg(tile)%data(:),ktau,nlv,ia,ib,ja,jb,0.,1.)
    call printa('fg  ',b_fg(tile)%data(:),ktau,nlv,ia,ib,ja,jb,0.,1.)
    call printa('thet',rhs,ktau,nlv,ia,ib,ja,jb,200.,1.)
  end if

  if ( ncvmix>0 ) then  ! cumulus mixing of momentum - jlm version
    do k = kuocb+1-ncvmix,kl-3
      ! for no-conv-points, doing k=1-ncvmix,1 with convpsav=0.
      !  1 below for ncvmix=2
      where ( b_kbsav(tile)%data(:)>0 .and. k>=b_kbsav(tile)%data(:)+1-ncvmix .and. k<b_ktsav(tile)%data(:) )
        guv(:,k)=guv(:,k)-b_convpsav(tile)%data(:)*.5 ! with factor of .5
      end where
    enddo    ! k loop
    if ( diag .and. mydiag ) then
      write(6,*)'vertmix after conv; kb,kt,convp',kbsav(idjd),ktsav(idjd),convpsav(idjd)
      write(6,*)'new guv',(guv(idjd,k),k=1,kl)
    end if
  end if      !   (ncvmix>0)

  conflux=grav*dt/dsig(1)
  condrag=grav*dt/(dsig(1)*rdry)
  ! first do theta (then convert back to t)
  at(:,1) =0.
  ct(:,kl)=0.

  do k = 2,kl
    at(:,k)=-gt(:,k-1)/dsig(k)/cnhs_fl(:,k)
  end do
  do k = 1,kl-1
    ct(:,k)=-gt(:,k)/dsig(k)/cnhs_fl(:,k)
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
  rhs(:,1) = rhs(:,1) - (conflux/cp)*b_fg(tile)%data(:)/(b_ps(tile)%data(:)*cnhs_fl(:,1))
  call trim(at,ct,rhs)   ! for t
  if ( nmaxpr==1 .and. mydiag ) write (6,"('thet_out',9f8.3/8x,9f8.3)") rhs(idjd,:)
  do k = 1,kl
    b_t(tile)%data(:,k) = rhs(:,k)/sigkap(k)
  enddo    !  k loop
  if ( diag ) then
    if ( mydiag ) then
      write(6,*)'vertmix eg,fg,cdtq ',eg(idjd),fg(idjd),cdtq(idjd)
    end if
    call printa('thet',rhs,ktau,nlv,ia,ib,ja,jb,200.,1.)
  end if
  
#ifdef scm
  rkmsave(tile)%data(:,:) = rkm(:,:)
  rkhsave(tile)%data(:,:) = rkh(:,:)

  ! counter-gradied included in pbldif.f90
  wth_flux(tile)%data(:,1) = fg(:)*rdry*b_t(tile)%data(:,1)/(b_ps(tile)%data(:)*cp)
  do k = 1,kl-1
    wth_flux(tile)%data(:,k+1) = rkh(:,k)*(rhs(1:imax,k+1)-rhs(1:imax,k))*(grav/rdry)*sig(k)/(b_t(tile)%data(:,k)*dsig(k))
  end do    
#endif

  !--------------------------------------------------------------
  ! Moisture
  rhs=b_qg(tile)%data(:,:)
  rhs(:,1)=rhs(:,1)-(conflux/hl)*b_eg(tile)%data(:)/(b_ps(tile)%data(:)*cnhs_fl(:,1))
  ! could add extra sfce moisture flux term for crank-nicholson
  call trim(at,ct,rhs)    ! for qg
  b_qg(tile)%data(:,:)=rhs
  if ( diag .and. mydiag ) then
    write(6,*)'vertmix rhs & qg after trim ',(rhs(idjd,k),k=1,kl)
    write (6,"('qg ',9f7.3/(8x,9f7.3))") (1000.*qg(idjd,k),k=1,kl)
  end if

#ifdef scm  
  ! counter-gradied included in pbldif.f90
  wq_flux(tile)%data(:,1) = eg(:)*rdry*b_t(tile)%data(:,1)/(b_ps(tile)%data(:)*hl)
  do k = 1,kl-1
    wq_flux(tile)%data(:,k+1) = rkh(:,k)*(b_qg(tile)%data(:,k+1)-b_qg(tile)%data(:,k))*(grav/rdry)*sig(k)/(b_t(tile)%data(:,k)*dsig(k))
  end do  
#endif
  
  !--------------------------------------------------------------
  ! Cloud microphysics terms
  if ( ldr/=0 ) then
    ! now do qfg
    rhs=b_qfg(tile)%data(:,:)
    call trim(at,ct,rhs)       ! for qfg
    b_qfg(tile)%data(:,:)=rhs
    ! now do qlg
    rhs=b_qlg(tile)%data(:,:)
    call trim(at,ct,rhs)       ! for qlg
    b_qlg(tile)%data(:,:)=rhs
    if ( ncloud>=4 ) then
      ! now do cldfrac
      rhs=b_stratcloud(tile)%data(:,:)
      call trim(at,ct,rhs)    ! for cldfrac
      b_stratcloud(tile)%data(:,:)=rhs
      call combinecloudfrac(b_cfrac(tile)%data(:,:),b_condc(tile)%data(:),b_kbsav(tile)%data(:),b_ktsav(tile)%data(:),tile,imax)
    else
      ! now do cfrac
      call convectivecloudfrac(clcon,b_condc(tile)%data(:),b_kbsav(tile)%data(:),b_ktsav(tile)%data(:),imax)
      rhs = (b_cfrac(tile)%data(:,:)-clcon(:,:))/(1.-clcon(:,:))
      call trim(at,ct,rhs)    ! for cfrac
      b_cfrac(tile)%data(:,:)=min(max(rhs+clcon-rhs*clcon,0.),1.)
    end if  ! (ncloud>=4)
  end if        ! (ldr/=0)
  
  !--------------------------------------------------------------
  ! Aerosols
  if ( abs(iaero)>=2 ) then
    do nt = 1,naero
      rhs(:,:) = b_xtg(tile)%data(:,:,nt) ! Total grid-box
      call trim(at,ct,rhs)
      b_xtg(tile)%data(:,:,nt) = rhs(:,:)
    end do
  end if ! (abs(iaero)>=2)

  !--------------------------------------------------------------
  ! Momentum terms
  au(:,1) = b_cduv(tile)%data(:)*condrag/(b_tss(tile)%data(:)*cnhs_fl(:,1))
  cu(:,kl) = 0.
  do k = 2,kl
    au(:,k) = -guv(:,k-1)/(dsig(k)*cnhs_fl(:,k))
  enddo    !  k loop
  do k = 1,kl-1
    cu(:,k) = -guv(:,k)/(dsig(k)*cnhs_fl(:,k))
  enddo    !  k loop
  if ( ( diag .or. ntest==2 ) .and. mydiag ) then
    write(6,*)'au ',(au(idjd,k),k=1,kl)
    write(6,*)'cu ',(cu(idjd,k),k=1,kl)
  end if      ! (ntest==2)

  ! first do u
  do k = 1,kl
    rhs(:,k) = b_u(tile)%data(:,k) - ou(:)
  end do
  call trim(au,cu,rhs)
  do k = 1,kl
    b_u(tile)%data(:,k) = rhs(:,k) + ou(:)
  end do
  if ( diag .and. mydiag ) then
    write(6,*)'vertmix au ',(au(idjd,k),k=1,kl)
  end if
  
#ifdef scm
  uw_flux(tile)%data(:,1) = -cduv(:)*(b_u(tile)%data(:,1)-ou(:))
  do k = 1,kl-1
    uw_flux(tile)%data(:,k+1) = rkm(:,k)*(b_u(tile)%data(:,k+1)-b_u(tile)%data(:,k))*(grav/rdry)*sig(k)/(b_t(tile)%data(:,k)*dsig(k))
  end do
#endif
  
  ! now do v; with properly unstaggered au,cu
  do k = 1,kl
    rhs(:,k) = b_v(tile)%data(:,k) - ov(:)
  end do
  call trim(au,cu,rhs)    ! note now that au, cu unstaggered globpea
  do k = 1,kl
    b_v(tile)%data(:,k) = rhs(:,k) + ov(:)
  end do

#ifdef scm
  vw_flux(tile)%data(:,1) = -cduv(:)*(b_v(tile)%data(:,1)-ov(:))
  do k = 1,kl-1
    vw_flux(tile)%data(:,k+1) = rkm(:,k)*(b_v(tile)%data(:,k+1)-b_v(tile)%data(:,k))*(grav/rdry)*sig(k)/(b_t(tile)%data(:,k)*dsig(k))
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

#ifndef scm
  !--------------------------------------------------------------
  ! Tracers
  if ( ngas>0 ) then
    call tracervmix(at,ct,tile,imax)
  end if ! (ngas>0)
#endif
      
else
      
  !-------------------------------------------------------------
  ! k-e + MF closure scheme
      
  ! note ksc/=0 options are clobbered when nvmix=6
  ! However, nvmix=6 with nlocal=7 supports its own shallow
  ! convection options

  ! calculate height on full levels
  zg(:,1) = bet(1)*b_t(tile)%data(:,1)/grav + b_phi_nh(tile)%data(:,1)/grav
  zg(:,1) = max( 1., zg(:,1) )
  zh(:,1) = b_t(tile)%data(:,1)*cnhs_fl(:,1)*delh(1)
  zh(:,1) = max( zg(:,1)+1., zh(:,1) )
  do k = 2,kl-1
    zg(:,k) = zg(:,k-1) + (bet(k)*b_t(tile)%data(:,k)+betm(k)*b_t(tile)%data(:,k-1))/grav + b_phi_nh(tile)%data(:,k)/grav
    zg(:,k) = max( zh(:,k-1)+1., zg(:,k) )
    zh(:,k) = zh(:,k-1) + b_t(tile)%data(:,k)*cnhs_fl(:,k)*delh(k)
    zh(:,k) = max( zg(:,k)+1., zh(:,k) )
  end do ! k  loop
  zg(:,kl) = zg(:,kl-1) + (bet(kl)*b_t(tile)%data(:,kl)+betm(kl)*b_t(tile)%data(:,kl-1))/grav + b_phi_nh(tile)%data(:,kl)/grav
  zg(:,kl) = max( zh(:,kl-1)+1., zg(:,kl) )
  zh(:,kl) = zh(:,kl-1) + b_t(tile)%data(:,kl)*cnhs_fl(:,kl)*delh(kl)
  zh(:,kl) = max( zg(:,kl)+1., zh(:,kl) )
       
  ! near surface air density (see sflux.f and cable_ccam2.f90)
  rhos(1:imax) = b_ps(tile)%data(:)/(rdry*b_tss(tile)%data(:))
    
  ! Use counter gradient for aerosol tracers
  if ( abs(iaero)>=2 ) then 
    tnaero = naero
  else
    tnaero = 0
  end if
  
  ! Special treatment for prognostic cloud fraction
  if ( ncloud>=4 ) then
    cldtmp = b_stratcloud(tile)%data(:,:)
  else
    call convectivecloudfrac(clcon,b_condc(tile)%data(:),b_kbsav(tile)%data(:),b_ktsav(tile)%data(:),imax)
    cldtmp(:,:) = (b_cfrac(tile)%data(:,:)-clcon(:,:))/(1.-clcon(:,:))
  end if
       
  ! transform to ocean reference frame and temp to theta
  do k = 1,kl
    b_u(tile)%data(:,k) = b_u(tile)%data(:,k) - ou
    b_v(tile)%data(:,k) = b_v(tile)%data(:,k) - ov
    rhs(:,k) = b_t(tile)%data(:,k)*sigkap(k) ! theta
  end do

#ifdef scm
  ! Evaluate EDMF scheme
  select case(nlocal)
    case(0) ! no counter gradient
      call tkemix(rkm,rhs,b_qg(tile)%data(:,:),b_qlg(tile)%data(:,:),b_qfg(tile)%data(:,:),cldtmp,b_u(tile)%data(:,:),b_v(tile)%data(:,:),b_pblh(tile)%data(:),b_fg(tile)%data(:),b_eg(tile)%data(:),b_ps(tile)%data(:),b_zo(tile)%data(:),zg,zh,sig,rhos,        &
                  dt,qgmin,1,0,tnaero,b_xtg(tile)%data(:,:,:),cgmap,wth_flux(tile)%data(:,:),wq_flux(tile)%data(:,:),uw_flux(tile)%data(:,:),vw_flux(tile)%data(:,:),mfout,tile,imax)
      rkh = rkm
    case(1,2,3,4,5,6) ! KCN counter gradient method
      call tkemix(rkm,rhs,b_qg(tile)%data(:,:),b_qlg(tile)%data(:,:),b_qfg(tile)%data(:,:),cldtmp,b_u(tile)%data(:,:),b_v(tile)%data(:,:),b_pblh(tile)%data(:),b_fg(tile)%data(:),b_eg(tile)%data(:),b_ps(tile)%data(:),b_zo(tile)%data(:),zg,zh,sig,rhos,        &
                  dt,qgmin,1,0,tnaero,b_xtg(tile)%data(:,:,:),cgmap,wth_flux(tile)%data(:,:),wq_flux(tile)%data(:,:),uw_flux(tile)%data(:,:),vw_flux(tile)%data(:,:),mfout,tile,imax)
      rkh = rkm
      do k = 1,kl
        uav(1:imax,k) = av_vmod*b_u(tile)%data(:,k) + (1.-av_vmod)*(b_savu(tile)%data(:,k)-ou)
        vav(1:imax,k) = av_vmod*b_v(tile)%data(:,k) + (1.-av_vmod)*(b_savv(tile)%data(:,k)-ov)
      end do
      call pbldif(rkm,rkh,rhs,uav,vav,cgmap,tile,imax)
    case(7) ! mass-flux counter gradient
      call tkemix(rkm,rhs,b_qg(tile)%data(:,:),b_qlg(tile)%data(:,:),b_qfg(tile)%data(:,:),cldtmp,b_u(tile)%data(:,:),b_v(tile)%data(:,:),b_pblh(tile)%data(:),b_fg(tile)%data(:),b_eg(tile)%data(:),b_ps(tile)%data(:),b_zo(tile)%data(:),zg,zh,sig,rhos,        &
                  dt,qgmin,0,0,tnaero,b_xtg(tile)%data(:,:,:),cgmap,wth_flux(tile)%data(:,:),wq_flux(tile)%data(:,:),uw_flux(tile)%data(:,:),vw_flux(tile)%data(:,:),mfout,tile,imax)
      rkh = rkm
    case DEFAULT
      write(6,*) "ERROR: Unknown nlocal option for nvmix=6"
      call ccmpi_abort(-1)
    end select
#else
  ! Evaluate EDMF scheme
  select case(nlocal)
    case(0) ! no counter gradient
      call tkemix(rkm,rhs,b_qg(tile)%data(:,:),b_qlg(tile)%data(:,:),b_qfg(tile)%data(:,:),cldtmp,b_u(tile)%data(:,:),b_v(tile)%data(:,:),b_pblh(tile)%data(:),b_fg(tile)%data(:),b_eg(tile)%data(:),b_ps(tile)%data(:),b_zo(tile)%data(:),zg,zh,sig,rhos, &
                  dt,qgmin,1,0,tnaero,b_xtg(tile)%data(:,:,:),cgmap,tile,imax) 
      rkh = rkm
    case(1,2,3,4,5,6) ! KCN counter gradient method
      call tkemix(rkm,rhs,b_qg(tile)%data(:,:),b_qlg(tile)%data(:,:),b_qfg(tile)%data(:,:),cldtmp,b_u(tile)%data(:,:),b_v(tile)%data(:,:),b_pblh(tile)%data(:),b_fg(tile)%data(:),b_eg(tile)%data(:),b_ps(tile)%data(:),b_zo(tile)%data(:),zg,zh,sig,rhos, &
                  dt,qgmin,1,0,tnaero,b_xtg(tile)%data(:,:,:),cgmap,tile,imax) 
      rkh = rkm
      do k = 1,kl
        uav(1:imax,k) = av_vmod*b_u(tile)%data(:,k) + (1.-av_vmod)*(b_savu(tile)%data(:,k)-ou)
        vav(1:imax,k) = av_vmod*b_v(tile)%data(:,k) + (1.-av_vmod)*(b_savv(tile)%data(:,k)-ov)
      end do
      call pbldif(rkm,rkh,rhs,uav,vav,cgmap,tile,imax)
    case(7) ! mass-flux counter gradient
      call tkemix(rkm,rhs,b_qg(tile)%data(:,:),b_qlg(tile)%data(:,:),b_qfg(tile)%data(:,:),cldtmp,b_u(tile)%data(:,:),b_v(tile)%data(:,:),b_pblh(tile)%data(:),b_fg(tile)%data(:),b_eg(tile)%data(:),b_ps(tile)%data(:),b_zo(tile)%data(:),zg,zh,sig,rhos, &
                  dt,qgmin,0,0,tnaero,b_xtg(tile)%data(:,:,:),cgmap,tile,imax) 
      rkh = rkm
    case DEFAULT
      write(6,*) "ERROR: Unknown nlocal option for nvmix=6"
      call ccmpi_abort(-1)
    end select
#endif
  
  ! special treatment for prognostic cloud fraction  
  if ( ncloud>=4 ) then
    b_stratcloud(tile)%data(:,:) = cldtmp
    call combinecloudfrac(b_cfrac(tile)%data(:,:),b_condc(tile)%data(:),b_kbsav(tile)%data(:),b_ktsav(tile)%data(:),tile,imax)
  else
    b_cfrac(tile)%data(:,:) = min(max(cldtmp+clcon-cldtmp*clcon,0.),1.)
  endif
  
  ! transform winds back to Earth reference frame and theta to temp
  do k = 1,kl
    b_u(tile)%data(:,k) = b_u(tile)%data(:,k) + ou
    b_v(tile)%data(:,k) = b_v(tile)%data(:,k) + ov
    b_t(tile)%data(:,k) = rhs(1:imax,k)/sigkap(k)
  enddo    !  k loop
  
#ifdef scm
  rkmsave(tile)%data(:,:) = rkm(:,:)
  rkhsave(tile)%data(:,:) = rkh(:,:)
  tkesave(tile)%data(:,1:kl) = tke(tile)%data(:,1:kl)
  mfsave(tile)%data(:,:) = mfout(:,:)
#else
  ! tracers
  if ( ngas>0 ) then
    do k = 1,kl-1
      delsig       =sig(k+1)-sig(k)
      dz(1:imax)  =-tmnht(1:imax,k)*delons(k)*cnhs_hl(1:imax,k)  ! this is z(k+1)-z(k)
      dzr(1:imax) =1./dz(1:imax)
      gt(1:imax,k)=rkh(1:imax,k)*dt*delsig*dzr(1:imax)**2
    end do      ! k loop
    gt(:,kl)=0.
    at(:,1) =0.
    ct(:,kl)=0.
    do k = 2,kl
      at(1:imax,k)=-gt(1:imax,k-1)/dsig(k)
    end do
    do k = 1,kl-1
      ct(1:imax,k)=-gt(1:imax,k)/dsig(k)
    end do
    ! increase mixing to replace counter gradient term      
    if ( nlocal>0 ) then
      at=cq*at
      ct=cq*ct
    end if
    call tracervmix(at,ct,tile,imax)
  end if
#endif
       
end if ! nvmix/=6 ..else..
      
return
end subroutine vertmix_work

subroutine vertjlm(rkm,rkh,rhs,sigkap,sighkap,delons,zh,tmnht,cnhs_hl,ntest,cgmap,tile,imax)

use arrays_m, only : ps,qg,t,u,v,zs                                                          ! Atmosphere dyamics prognostic arrays
use cc_mpi, only : mydiag                                                                    ! CC MPI routines
use cfrac_m, only : cfrac                                                                    ! Cloud fraction
use const_phys, only : epsil,grav,hl,hlcp,hlfcp,hlscp,roncp,rvap                             ! Physical constants
use diag_m, only : printa                                                                    ! Diagnostic routines
use estab, only : establ                                                                     ! Liquid saturation function
use kuocomb_m, only : kbsav,ktsav                                                            ! JLM convection
use liqwpar_m, only : qfg,qlg                                                                ! Cloud water mixing ratios
use morepbl_m, only : condc,condx,pblh                                                       ! Additional boundary layer diagnostics
use newmpar_m, only : kl,nproc                                                               ! Grid parameters
use parm_m, only : amxlsq,av_vmod,diag,dt,ia,ib,idjd,ja,jb,ktau,nlocal,nlv,nmaxpr,ntau,nvmix ! Model configuration
use savuvt_m, only : savu,savv                                                               ! Saved dynamic arrays
use screen_m, only : qgscrn,tscrn                                                            ! Screen level diagnostics
use sigs_m, only : dsig,sig,sigmh                                                            ! Atmosphere sigma levels
use soil_m, only : land,zmin                                                                 ! Soil and surface data
use vertmixdata_m

implicit none

include 'kuocom.h'                  ! Convection parameters

integer, intent(in) :: tile,imax
integer, parameter :: ndvmod=0    ! 0 default, 1+ for dvmod tests
integer, intent(in) :: ntest
integer iq,k,iqmax
integer, dimension(imax) :: kbase,ktop
real, parameter :: lambda=0.45               ! coefficients for Louis scheme
real, parameter :: vkar4=0.4,bprmj=5.,cmj=5. ! coefficients for Louis scheme
real, parameter :: chj=2.6                   ! coefficients for Louis scheme
real delta,es,pk,dqsdt,betat,betaq,betac,al,qc,fice
real w1,w2,diffmax,rhsk,rhskp,delthet_old,xold,diff
real denma,denha,esp,epart,tsp,qbas
real, dimension(imax,kl), intent(inout) :: rhs, rkm, rkh
real, dimension(imax,kl), intent(in) :: zh
real, dimension(imax,kl-1), intent(in) :: tmnht, cnhs_hl
real, dimension(imax,kl) :: qs,betatt,betaqt,delthet
real, dimension(imax,kl) :: uav,vav,ri,rk_shal
real, dimension(imax,kl) :: thee,thebas
real, dimension(imax), intent(in) :: cgmap
real, dimension(imax) :: dz,dzr,dvmod,dqtot,x,zhv
real, dimension(imax) :: csq,sqmxl,fm,fh,sigsp
real, dimension(imax) :: theeb
real, dimension(kl), intent(in) :: sigkap,sighkap,delons
integer :: is, ie

is=(tile-1)*imax+1
ie=tile*imax

w1=0.
pk=0.

if ( nmaxpr==1 .and. mydiag ) write (6,"('thet_in',9f8.3/7x,9f8.3)") rhs(idjd,:)
  
! Pre-calculate the buoyancy parameters if using qcloud scheme.
! Follow Smith's (1990) notation; gam() is HBG's notation for (L/cp)dqsdt.
! The factor of (1/sigkap)=T/theta in betatt differs from Smith's formulation
! because we use theta derivative rather than (dry static energy)/cp.
if ( (nvmix>0.and.nvmix<4) .or. nvmix==7 ) then
  delta=1./epsil-1.  ! i.e. 1/.622 -1., i.e. .6077
  do k=1,kl
    if ( sig(k)>.8 ) then ! change made 17/1/06
      do iq=1,imax
        es=establ(b_t(tile)%data(iq,k))
        pk=b_ps(tile)%data(iq)*sig(k)
        qs(iq,k)=.622*es/max(1.,pk-es)  
        dqsdt=qs(iq,k)*pk*(hl/rvap)/(b_t(tile)%data(iq,k)**2*max(1.,pk-es))
        rhs(iq,k)=rhs(iq,k)-(hlcp*b_qlg(tile)%data(iq,k)+hlscp*b_qfg(tile)%data(iq,k))*sigkap(k)   !Convert to thetal - used only to calc Ri
        betat=1./b_t(tile)%data(iq,k)
        qc=b_qlg(tile)%data(iq,k)+b_qfg(tile)%data(iq,k)
        fice=b_qfg(tile)%data(iq,k)/max(qc,1.e-12)
        betaq=delta/(1.+delta*b_qg(tile)%data(iq,k)-qc)
        al=1./(1.+hlcp*dqsdt)
        betac=b_cfrac(tile)%data(iq,k)*al * ((hlcp+fice*hlfcp)*betat - betaq/(1.-epsil) )
        betatt(iq,k)=(betat-dqsdt*betac)/sigkap(k)  !Beta_t_tilde
        betaqt(iq,k)=betaq+betac                    !Beta_q_tilde
      enddo   ! iq loop
    else  ! i.e. (sig(k)<.8)
      do iq=1,imax
        es=establ(b_t(tile)%data(iq,k))
        pk=b_ps(tile)%data(iq)*sig(k)
        qs(iq,k)=.622*es/max(1.,pk-es)  ! still need qs(); max for k=kl
        betat=1./b_t(tile)%data(iq,k)
        ! qc=b_qlg(tile)%data(iq,k)+b_qfg(tile)%data(iq,k)
        ! betaq=delta/(1.+delta*b_qg(tile)%data(iq,k)-qc)
        betaq=delta/(1.+delta*b_qg(tile)%data(iq,k))
        betatt(iq,k)=betat/sigkap(k)    !Beta_t_tilde
        betaqt(iq,k)=betaq              !Beta_q_tilde
      enddo   ! iq loop
    endif  ! (sig(k)>.8)
    if(diag.and.mydiag)then
      iq=idjd
      dqsdt=qs(iq,k)*pk*(hl/rvap)/(b_t(tile)%data(iq,k)**2*max(pk-es,1.))
      betat=1./b_t(tile)%data(iq,k)
      qc=b_qlg(tile)%data(iq,k)+b_qfg(tile)%data(iq,k)
      fice=b_qfg(tile)%data(iq,k)/max(qc,1.e-12)
      betaq=delta/(1+delta*b_qg(tile)%data(iq,k)-qc)
      ! al=1./(1.+gam(iq,k))
      al=1./(1.+hlcp*dqsdt)
      betac=b_cfrac(tile)%data(iq,k)*al*((hlcp+fice*hlfcp)*betat - betaq/(1.-epsil) )
      write(6,*)'k,qg,qs,qlg,qfg,qc ',k,b_qg(tile)%data(iq,k),qs(iq,k),b_qlg(tile)%data(iq,k),b_qfg(tile)%data(iq,k),qc
      write(6,*)'t,rhs,cfrac,fice ',b_t(tile)%data(iq,k),rhs(iq,k),b_cfrac(tile)%data(iq,k),fice
      write(6,*)'al,delta,betaq,betac ',al,delta,betaq,betac
      write(6,*)'betat,betatt,betaqt ',betat,betatt(iq,k),betaqt(iq,k)
    endif 
  enddo    !  k loop
else       ! other nvmix values (0 or 4+) still need qs()
  do k=1,kl
    do iq=1,imax
      es=establ(b_t(tile)%data(iq,k))
      qs(iq,k)=.622*es/max(1.,b_ps(tile)%data(iq)*sig(k)-es)  ! max for k=kl
    enddo   ! iq loop
  enddo    !  k loop
endif      ! (nvmix>0.and.nvmix<4).or.nvmix==7

do k=1,kl-1
  delthet(:,k)=rhs(:,k+1)-rhs(:,k)  ! rhs is theta or thetal here
enddo      !  k loop

!     ****** section for Geleyn shallow convection; others moved lower****
if(ksc==-99)then
  do k=kscbase,ksctop    ! new usage of ksc thu  02-17-2000
    do iq=1,imax
      delthet(iq,k)=delthet(iq,k)-hlcp*max(0.,qs(iq,k+1)-b_qg(tile)%data(iq,k+1)-qs(iq,k)+b_qg(tile)%data(iq,k) )
    enddo  ! iq loop
  enddo   !  k loop
endif     ! (ksc==-99)
if(ksc==-98)then    ! modified Geleyn Jan '08
  do k=kscbase,ksctop    ! new usage of ksc thu  02-17-2000
    do iq=1,imax
      if(b_qg(tile)%data(iq,k)>tied_rh*qs(iq,k).or.b_qg(tile)%data(iq,k+1)>tied_rh*qs(iq,k+1))then
        delthet(iq,k)=delthet(iq,k)-hlcp*max(0.,qs(iq,k+1)-b_qg(tile)%data(iq,k+1)-qs(iq,k)+b_qg(tile)%data(iq,k) )
      endif
    enddo  ! iq loop
  enddo   !  k loop
endif     ! (ksc==-98)
if(ksc==-97)then    ! modified Geleyn Jan '08
  do k=kscbase,ksctop    ! new usage of ksc thu  02-17-2000
    do iq=1,imax
      if(b_qg(tile)%data(iq,k)>tied_rh*qs(iq,k))then
        delthet(iq,k)=delthet(iq,k)-hlcp*max(0.,qs(iq,k+1)-b_qg(tile)%data(iq,k+1)-qs(iq,k)+b_qg(tile)%data(iq,k) )
      endif
    enddo  ! iq loop
  enddo   !  k loop
endif     ! (ksc==-97)
if(ksc==-96)then   ! combination of Geleyn and jlm 83 (Feb 08)
  do k=1,ksctop    
    do iq=1,imax
      if(k<b_ktsav(tile)%data(iq).and.k>=b_kbsav(tile)%data(iq).and.b_condc(tile)%data(iq)<1.e-20)then  
        delthet(iq,k)=delthet(iq,k)-hlcp*max(0.,qs(iq,k+1)-b_qg(tile)%data(iq,k+1)-qs(iq,k)+b_qg(tile)%data(iq,k) )
        write(6,*)'-96 iq,k,kbsav,ktsav ',iq,k,b_kbsav(tile)%data(iq),b_ktsav(tile)%data(iq),hlcp*(qs(iq,k+1)-b_qg(tile)%data(iq,k+1)-qs(iq,k)+b_qg(tile)%data(iq,k)),delthet(iq,k)
      endif 
    enddo  ! iq loop
  enddo   !  k loop
endif     ! (ksc==-96)
if(ksc==-95)then ! same as -99 but has tied_rh (e.g. .75) factor
  ! kshal(:)=0
  do k=kscbase,ksctop     
    do iq=1,imax
      delthet(iq,k)=delthet(iq,k)-tied_rh*hlcp*max(0.,(qs(iq,k+1)-b_qg(tile)%data(iq,k+1)-qs(iq,k)+b_qg(tile)%data(iq,k)) )
    enddo  ! iq loop
  enddo   !  k loop
endif     ! (ksc==-95)
if(ksc==-94)then   ! combination of Geleyn and jlm 83 (Feb 08)
  do k=1,ksctop    
    do iq=1,imax
      if(k<b_ktsav(tile)%data(iq).and.k>=b_kbsav(tile)%data(iq).and.b_condc(tile)%data(iq)<1.e-20)then  
        delthet(iq,k)=0.
        write(6,*)'-94 iq,k,kbsav,ktsav ',iq,k,b_kbsav(tile)%data(iq),b_ktsav(tile)%data(iq),hlcp*(qs(iq,k+1)-b_qg(tile)%data(iq,k+1)-qs(iq,k)+b_qg(tile)%data(iq,k)),delthet(iq,k)
      endif 
    enddo  ! iq loop
  enddo   !  k loop
endif     ! (ksc==-94)
if(ksc==-93)then ! single-layer (pblh) version of -95
  do iq=1,imax
    do k=kscbase,kl/2
      if(zh(iq,k)<b_pblh(tile)%data(iq).and.zh(iq,k+1)>b_pblh(tile)%data(iq))then
        delthet(iq,k)=delthet(iq,k)-tied_rh*hlcp*max(0.,(qs(iq,k+1)-b_qg(tile)%data(iq,k+1)-qs(iq,k)+b_qg(tile)%data(iq,k)) )
      endif
    enddo  ! k loop 
  enddo   ! iq loop
endif     ! (ksc==-93)
if(ksc==-92)then ! capped-by-pblh version of -95
  do iq=1,imax
    do k=kscbase,kl/2
      if(zh(iq,k)<b_pblh(tile)%data(iq))then
        delthet(iq,k)=delthet(iq,k)-tied_rh*hlcp*max(0.,(qs(iq,k+1)-b_qg(tile)%data(iq,k+1)-qs(iq,k)+b_qg(tile)%data(iq,k)) )
      endif
    enddo  ! k loop 
  enddo   ! iq loop
endif     ! (ksc==-92)
if(ksc==-91)then ! capped-by-pblh (anywhere in layer) version of -95
  do iq=1,imax
    do k=2,kl/2
      if(zh(iq,k-1)<b_pblh(tile)%data(iq))then
        delthet(iq,k)=delthet(iq,k)-tied_rh*hlcp*max(0.,(qs(iq,k+1)-b_qg(tile)%data(iq,k+1)-qs(iq,k)+b_qg(tile)%data(iq,k)) )
      endif
    enddo  ! k loop 
  enddo   ! iq loop
endif     ! (ksc==-91)
!     ********* end of Geleyn shallow convection section ****************

! following now defined in vertmix (don't need to pass from sflux)
uav(1:imax,:)=av_vmod*b_u(tile)%data(:,:)+(1.-av_vmod)*b_savu(tile)%data(:,:) 
vav(1:imax,:)=av_vmod*b_v(tile)%data(:,:)+(1.-av_vmod)*b_savv(tile)%data(:,:) 
do k=1,kl-1
  do iq=1,imax
    dz(iq) =-tmnht(iq,k)*delons(k)*cnhs_hl(iq,k)  ! this is z(k+1)-z(k)
    dzr(iq)=1./dz(iq)
    zhv(iq)=1./zh(iq,k)
    if(ndvmod==0)then
      dvmod(iq)=sqrt( (uav(iq,k+1)-uav(iq,k))**2+(vav(iq,k+1)-vav(iq,k))**2 )
    else
      dvmod(iq)=ndvmod  ! just for tests
    endif    ! (ndvmod==0)
  enddo ! iq loop

  ! x is bulk ri *(dvmod **2); used to calc ri(k), rkh(k) etc
  if((nvmix>0.and.nvmix<4).or.nvmix==7)then  ! new one allowing for cloudy air
    ! usually nvmix=3       
    if(sig(k)>.8)then ! change made 17/1/06
      dqtot(:)=b_qg(tile)%data(:,k+1)+b_qlg(tile)%data(:,k+1)+b_qfg(tile)%data(:,k+1)-(b_qg(tile)%data(:,k)  +b_qlg(tile)%data(:,k)  +b_qfg(tile)%data(:,k))
    else
      dqtot(:)=0.
    endif
    if(nvmix==2)then !  jlm May '05
      do iq=1,imax
        x(iq)=grav*dz(iq)*((min(betatt(iq,k),betatt(iq,k+1)))*delthet(iq,k) + (max(betaqt(iq,k),betaqt(iq,k+1)))*dqtot(iq) )
      enddo ! iq loop	
    else    ! i.e. nvmix=1 or 3
      if(nvmix==1)w1=dsig(k+1)/(dsig(k)+dsig(k+1)) 
      if(nvmix==3)w1=1.    !weight for lower level  usual  
      if(nvmix==7)w1=1.
      w2=1.-w1             !weight for upper level
      do iq=1,imax
        x(iq)=grav*dz(iq)*((w1*betatt(iq,k)+w2*betatt(iq,k+1))*delthet(iq,k) + (w1*betaqt(iq,k)+w2*betaqt(iq,k+1))*dqtot(iq) )
      enddo ! iq loop	 
    endif  !  (nvmix==2) .. else ..
    if(ntest==4.and.k<=9.and.ktau==ntau)then
      diffmax=0.
      do iq=1,imax
        rhsk=b_t(tile)%data(iq,k)*sigkap(k)
        rhskp=b_t(tile)%data(iq,k+1)*sigkap(k+1)
        delthet_old=rhs(iq,k+1)-rhs(iq,k)
        xold=grav*dz(iq)*(delthet_old/(tmnht(iq,k)*sighkap(k))+.61*(b_qg(tile)%data(iq,k+1)-b_qg(tile)%data(iq,k)))
        diff=abs(xold-x(iq))
        if(diff>diffmax)then
          diffmax=diff
          iqmax=iq
        endif
        write(47,'(3g13.4,i7,i4)') xold,x(iq),diff,iq,k
      enddo
      write(6,*)'k,iqmax,diffmax ',k,iqmax,diffmax
    endif   ! (ntest==4.and.k<=9.and.ktau==ntau)
    rhs(:,k)=b_t(tile)%data(:,k)*sigkap(k)   !need to re-set theta for nvmix=1-3
  elseif(nvmix==5)then        ! non-cloudy x with qfg, qlg
    x(:)=grav*dz(:)*(delthet(:,k)/(tmnht(:,k)*sighkap(k))+.61*(b_qg(tile)%data(:,k+1)-b_qg(tile)%data(:,k)) &
         -b_qfg(tile)%data(:,k+1)-b_qlg(tile)%data(:,k+1)+b_qfg(tile)%data(:,k)+b_qlg(tile)%data(:,k) )
  else                 ! original non-cloudy x, nvmix=4
    x(:)=grav*dz(:)*(delthet(:,k)/(tmnht(:,k)*sighkap(k))+.61*(b_qg(tile)%data(:,k+1)-b_qg(tile)%data(:,k)))
  endif     ! (nvmix>0.and.nvmix<4) .. else ..

  ! fm and fh denote f(Louis style)*dvmod
  ! nb. an error exists in the Louis expression for c; this is corrected
  ! in the current code
  csq(:) = zhv(:)*(((1.+dz(:)*zhv(:))**(1./3.)-1.)*dzr(:))**3

  ! newest code, stable same as csiro9 here (originally for nvmix=4)
  do iq=1,imax
    sqmxl(iq)=(vkar4*zh(iq,k)/(1.+vkar4*zh(iq,k)/amxlsq))**2
    dvmod(iq)=max( dvmod(iq) , 1. )
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
        if(b_land(tile)%data(iq))fh(iq)=abs(nvmix)
      enddo   ! iq loop
    endif
  endif     !  (nvmix<0)

  ! calculate k's, guv and gt
  ! nonlocal scheme usually applied to temperature and moisture, not momentum
  ! (i.e. local scheme is applied to momentum for nlocal=0,1)
  if (nvmix/=0) then    ! use nvmix=0 only for special test runs
    rkm(:,k)=fm(:)*sqmxl(:)*dzr(:)
    rkh(:,k)=fh(:)*sqmxl(:)*dzr(:)
    !added by Jing Huang on 4 Feb 2016
    if (nvmix==7) then
      where ( ri(:,k)>0. )
        sqmxl(:)=(vkar4*zh(:,k)/(1.+vkar4*zh(:,k)/amxlsq+vkar4*zh(:,k)*ri(:,k)/lambda))**2
        rkm(:,k)=dvmod(:)*dzr(:)*sqmxl(:)
        rkh(:,k)=rkh(:,k)
      end where
    end if
  else
    rkm(:,:)=0.
    rkh(:,:)=0.
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
  do k=1,kl
    prcpv(k)=sig(k)**(-roncp)
  enddo       
  write (6,"('thee',9f8.3/4x,9f8.3)") (prcpv(k)*t(idjd,k)*(t(idjd,k) + .5*hlcp*qs(idjd,k)) &
                                      /(t(idjd,k) - .5*hlcp*qs(idjd,k)),k=1,kl)
endif
if(nmaxpr==1.and.mydiag)then
  write (6,"('rino_v',9f9.3/6x,9f9.3)") ri(idjd,1:kl-1)
  write (6,"('rkh0 ',9f9.3/5x,9f9.3)") rkh(idjd,1:kl-2)
  write (6,"('rkm0 ',9f9.3/5x,9f9.3)") rkm(idjd,1:kl-2)
endif

if(nlocal/=0)then
  call pbldif(rkm,rkh,rhs,uav,vav,cgmap,tile,imax)  ! rhs is theta or thetal
  ! n.b. *** pbldif partially updates qg and theta (t done during trim)	 
  ! and updates rkh and rkm arrays
  if(nmaxpr==1.and.mydiag)then
    write (6,"('pblh ',f8.2)") pblh(idjd)
    write (6,"('rkh1 ',9f9.3/5x,9f9.3)") rkh(idjd,1:kl-2)
    write (6,"('rkm1 ',9f9.3/5x,9f9.3)") rkm(idjd,1:kl-2)
  endif
  if(nmaxpr==1.and.mydiag) write (6,"('thet_pbl',9f8.3/8x,9f8.3)") rhs(idjd,:)
  if( (diag.or.ntest>=1) .and. mydiag ) write (6,"('qg ',9f8.3/4x,9f8.3)") qg(idjd,:)
  if(diag)then
    call printa('rkh ',rkh,ktau,nlv,ia,ib,ja,jb,0.,1.)
    call printa('cond',condx,ktau,1,ia,ib,ja,jb,0.,1.)
    call printa('zs  ',zs,ktau,1,ia,ib,ja,jb,0.,1.)
  endif
else
  b_pblh(tile)%data(:) = zmin
endif      ! (nlocal>0)

rk_shal(:,:)=0.
!     ***** ***** section for jlm shallow convection v4 *****************
if(ksc==81)then
  do k=1,ksctop-1   ! or ksctop?  
    do iq=1,imax
      if(sig(b_ktsav(tile)%data(iq))>sig_ct.and.k<b_ktsav(tile)%data(iq).and.b_condc(tile)%data(iq)<1.e-20)then  
        rk_shal(iq,k)=tied_con
        rk_shal(iq,k+1)=tied_over
      endif ! (sig(b_ktsav(tile)%data(iq))>sig_ct.and.k<b_ktsav(tile)%data(iq).and....)
    enddo  ! iq loop
  enddo   !  k loop
endif     ! (ksc==81)
if(ksc==82)then
  do k=1,ksctop-1    
    do iq=1,imax
      if(sig(b_ktsav(tile)%data(iq))>sig_ct.and.k<b_ktsav(tile)%data(iq).and.k>=b_kbsav(tile)%data(iq).and.b_condc(tile)%data(iq)<1.e-20)then  
        rk_shal(iq,k)=tied_con
        rk_shal(iq,k+1)=tied_over
      endif ! (sig(b_ktsav(tile)%data(iq))>sig_ct.and.k<b_ktsav(tile)%data(iq).and....)
    enddo  ! iq loop
  enddo   !  k loop
endif     ! (ksc==82)
if(ksc==83)then
  do k=1,ksctop-1    
    do iq=1,imax
      if(sig(k)>sig_ct.and.k<b_ktsav(tile)%data(iq).and.k>=b_kbsav(tile)%data(iq).and.b_condc(tile)%data(iq)<1.e-20)then  
        rk_shal(iq,k)=tied_con
        rk_shal(iq,k+1)=tied_over
      endif ! (sig(b_ktsav(tile)%data(iq))>sig_ct.and.k<b_ktsav(tile)%data(iq).and....)
    enddo  ! iq loop
  enddo   !  k loop
endif     ! (ksc==83)
if(ksc==91)then
  do k=1,ksctop ! May 08
    do iq=1,imax
      if(b_ktsav(tile)%data(iq)<kl-1.and.k<b_ktsav(tile)%data(iq))then  ! April 04
        rk_shal(iq,k)=tied_con
        rk_shal(iq,k+1)=tied_over
      endif ! (b_ktsav(tile)%data(iq)<0.and.k<abs(b_ktsav(tile)%data(iq)))
    enddo  ! iq loop
  enddo   !  k loop
endif     ! (ksc==91)
if(ksc==92)then
  do k=1,ksctop ! May 08
    do iq=1,imax
      if(k>=b_kbsav(tile)%data(iq).and.k<b_ktsav(tile)%data(iq).and.b_condc(tile)%data(iq)<1.e-20)then  ! May 08
        rk_shal(iq,k)=tied_con
        rk_shal(iq,k+1)=tied_over
      endif ! (b_ktsav(tile)%data(iq)<0.and.k<abs(b_ktsav(tile)%data(iq)))
    enddo  ! iq loop
  enddo   !  k loop
endif     ! (ksc==92)
!     *********** end of jlm shallow convection section *****************

!     ************ section for Tiedtke shallow convection *******************
if(ksc==99)then
  do iq=1,imax
    theeb(iq)=prcpv(kscbase)*b_t(tile)%data(iq,kscbase)*(b_t(tile)%data(iq,kscbase)+.5*hlcp*qs(iq,kscbase))/(b_t(tile)%data(iq,kscbase)-.5*hlcp*qs(iq,kscbase))
  enddo    ! iq loop
  if(kscsea==1)then  ! Tiedtke done only over sea
    do k=kscbase+1,ksctop
      do iq=1,imax
        if(.not.b_land(tile)%data(iq))then 
          thee(iq,k)=prcpv(k)*b_t(tile)%data(iq,k)*(b_t(tile)%data(iq,k) + .5*hlcp*qs(iq,k))/(b_t(tile)%data(iq,k) - .5*hlcp*qs(iq,k))
          if(b_qg(tile)%data(iq,kscbase)>tied_rh*qs(iq,kscbase).and.thee(iq,k)<theeb(iq))then         !  default tied_rh=.75
            rk_shal(iq,k-1)=tied_con     !  m**2/sec  6., originally 10.
            rk_shal(iq,k)=tied_over      !  m**2/sec
          endif ! (b_qg(tile)%data(iq,kscbase)>rhscon*qs(iq,kscbase).....
        endif  ! (.not.b_land(tile)%data(iq)) 
      enddo   ! iq loop
    enddo    ! end of k=kscbase+1,ksctop loop
  else       !  i.e. Tiedtke original scheme over land and sea
    do k=kscbase+1,ksctop  ! typically kscbase=3 & ksctop=6
      do iq=1,imax
        thee(iq,k)=prcpv(k)*b_t(tile)%data(iq,k)*(b_t(tile)%data(iq,k) + .5*hlcp*qs(iq,k))/(b_t(tile)%data(iq,k) - .5*hlcp*qs(iq,k))
        if(b_qg(tile)%data(iq,kscbase)>tied_rh*qs(iq,kscbase).and.thee(iq,k)<theeb(iq))then !  default tied_rh=.75
          rk_shal(iq,k-1)=tied_con     !  m**2/sec  6., originally 10.
          rk_shal(iq,k)=tied_over      !  m**2/sec
          if(ntest==3.and.k==ksctop)then
            write(6,*) 'ktau,iq,theeb,thee,delthee ',ktau,iq,theeb(iq),thee(iq,k),theeb(iq)-thee(iq,k)
          endif
        endif ! (b_qg(tile)%data(iq,kscbase)>rhscon*qs(iq,kscbase).....
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
      epart=b_qg(tile)%data(iq,1)*.01*b_ps(tile)%data(iq)/(tied_rh*epsil) !in hPa 
      tsp=2840./(3.5*log(b_tscrn(tile)%data(iq))-log(epart)-4.805) + 55.   
      sigsp(iq)=(tsp/b_tscrn(tile)%data(iq))**(1./roncp)
    enddo
  else                    ! uses qgscrn
    do iq=1,imax
      ! Leon used factor 1.01 over sea; here div by tied_rh (e.g. .99)	 
      epart=b_qgscrn(tile)%data(iq)*.01*b_ps(tile)%data(iq)/(tied_rh*epsil) !in hPa   
      tsp=2840./(3.5*log(b_tscrn(tile)%data(iq))-log(epart)-4.805) + 55.   
      sigsp(iq)=(tsp/b_tscrn(tile)%data(iq))**(1./roncp)
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
      if(zh(iq,kbase(iq))>b_pblh(tile)%data(iq))then 
        kbase(iq)=kl
      endif
    enddo  ! iq loop
  endif  ! (nlocal/=0)
! following has some jlm shortcuts	 
  do iq=1,imax
    k=kbase(iq)
    qbas=qs(iq,k)
    theeb(iq)=prcpv(k)*b_t(tile)%data(iq,k)*(b_t(tile)%data(iq,k) + .5*hlcp*qbas)/(b_t(tile)%data(iq,k) - .5*hlcp*qbas)
  enddo  ! iq loop
  ktop(:)=kbase(:)
  do k=2,ksctop+1
    do iq=1,imax
      thee(iq,k)=prcpv(k)*b_t(tile)%data(iq,k)*(b_t(tile)%data(iq,k) + .5*hlcp*qs(iq,k))/(b_t(tile)%data(iq,k) - .5*hlcp*qs(iq,k))
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
        if(k>kbase(iq).and.k<=ktop(iq).and.b_condc(tile)%data(iq)<1.e-20)then
          rk_shal(iq,k-1)=tied_con     !  m**2/sec  6., originally 10.
          rk_shal(iq,k)=tied_over      !  m**2/sec
        endif
      enddo  ! iq loop
    enddo   ! k loop
  endif     ! (ksc==95)
  if(ksc==96)then   ! applied from sfce up
    do k=2,ksctop
      do iq=1,imax
        if(ktop(iq)>kbase(iq).and.k<=ktop(iq).and.b_condc(tile)%data(iq)<1.e-20)then  
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
      qbas=min(b_qg(tile)%data(iq,k)/tied_rh,qs(iq,k))
!     thebas(iq,k)=b_t(tile)%data(iq,k)+hlcp*qbas  + hght
      thebas(iq,k)=prcpv(k)*b_t(tile)%data(iq,k)*(b_t(tile)%data(iq,k) + .5*hlcp*qbas)/(b_t(tile)%data(iq,k) - .5*hlcp*qbas)
    enddo  ! iq loop
  enddo   ! k loop
  theeb(:)=thebas(:,kscbase)
  do k=kscbase+1,ksctop+1
    do iq=1,imax
      thee(iq,k)=prcpv(k)*b_t(tile)%data(iq,k)*(b_t(tile)%data(iq,k) + .5*hlcp*qs(iq,k))/(b_t(tile)%data(iq,k) - .5*hlcp*qs(iq,k))
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
        write(6,*) 'iq,land,rk_shal ',iq,b_land(tile)%data(iq),rk_shal(iq,ksctop-1)
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
    write(6,*)'kbsav,ktsav,theeb: ',b_kbsav(tile)%data(iq),b_ktsav(tile)%data(iq),theeb(iq)
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

return
end subroutine vertjlm

subroutine trim(a,c,rhs)

use newmpar_m

implicit none

integer k
real, dimension(imax,kl), intent(in) :: a, c
real, dimension(imax,kl), intent(inout) :: rhs
real, dimension(imax,kl) :: e, g
real, dimension(imax) :: b, temp

! this routine solves the system
!   a(k)*u(k-1)+b(k)*u(k)+c(k)*u(k+1)=rhs(k)    for k=2,kl-1
!   with  b(k)*u(k)+c(k)*u(k+1)=rhs(k)          for k=1
!   and   a(k)*u(k-1)+b(k)*u(k)=rhs(k)          for k=kl

! the Thomas algorithm is used
b(:)=1.-a(:,1)-c(:,1)
e(:,1)=c(:,1)/b(:)
g(:,1)=rhs(:,1)/b(:)
do k = 2,kl-1
  b(:)=1.-a(:,k)-c(:,k)
  temp(:)= 1./(b(:)-a(:,k)*e(:,k-1))
  e(:,k)=c(:,k)*temp(:)
  g(:,k)=(rhs(:,k)-a(:,k)*g(:,k-1))*temp(:)
end do

! do back substitution to give answer now
b(:)=1.-a(:,kl)-c(:,kl)
rhs(:,kl)=(rhs(:,kl)-a(:,kl)*g(:,kl-1))/(b(:)-a(:,kl)*e(:,kl-1))
do k = kl-1,1,-1
  rhs(:,k)=g(:,k)-e(:,k)*rhs(:,k+1)
end do

return
end subroutine trim

end module vertmix_m
  
subroutine trim(a,c,rhs)

use newmpar_m

implicit none

integer k
real, dimension(ifull,kl), intent(in) :: a, c
real, dimension(ifull,kl), intent(inout) :: rhs
real, dimension(ifull,kl) :: e, g
real, dimension(ifull) :: b, temp

! this routine solves the system
!   a(k)*u(k-1)+b(k)*u(k)+c(k)*u(k+1)=rhs(k)    for k=2,kl-1
!   with  b(k)*u(k)+c(k)*u(k+1)=rhs(k)          for k=1
!   and   a(k)*u(k-1)+b(k)*u(k)=rhs(k)          for k=kl

! the Thomas algorithm is used
b(:)=1.-a(:,1)-c(:,1)
e(:,1)=c(:,1)/b(:)
g(:,1)=rhs(:,1)/b(:)
do k = 2,kl-1
  b(:)=1.-a(:,k)-c(:,k)
  temp(:)= 1./(b(:)-a(:,k)*e(:,k-1))
  e(:,k)=c(:,k)*temp(:)
  g(:,k)=(rhs(:,k)-a(:,k)*g(:,k-1))*temp(:)
end do

! do back substitution to give answer now
b(:)=1.-a(:,kl)-c(:,kl)
rhs(:,kl)=(rhs(:,kl)-a(:,kl)*g(:,kl-1))/(b(:)-a(:,kl)*e(:,kl-1))
do k = kl-1,1,-1
  rhs(:,k)=g(:,k)-e(:,k)*rhs(:,k+1)
end do

return
end subroutine trim
