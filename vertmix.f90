! Conformal Cubic Atmospheric Model
    
! Copyright 2015 Commonwealth Scientific Industrial Research Organisation (CSIRO)
    
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

! nlocal=0 No counter gradient term
! nlocal=6 Holtslag and Boville non-local term
! nlocal=7 Mass flux based counter gradient (requires nvmix=6)
      
! kscmom   0 off, 1 turns on shal_conv momentum (usual) (requires nvmix<6)

!--------------------------------------------------------------
! Control subroutine for vertical mixing
subroutine vertmix

use aerosolldr                      ! LDR prognostic aerosols
use arrays_m                        ! Atmosphere dyamics prognostic arrays
use cc_mpi                          ! CC MPI routines
use cfrac_m                         ! Cloud fraction
use cloudmod                        ! Prognostic strat cloud
use diag_m                          ! Diagnostic routines
use extraout_m                      ! Additional diagnostics
use kuocomb_m                       ! JLM convection
use liqwpar_m                       ! Cloud water mixing ratios
use map_m                           ! Grid map arrays
use mlo                             ! Ocean physics and prognostic arrays
use morepbl_m                       ! Additional boundary layer diagnostics
use nharrs_m                        ! Non-hydrostatic atmosphere arrays
use pbl_m                           ! Boundary layer arrays
use savuvt_m                        ! Saved dynamic arrays
use sigs_m                          ! Atmosphere sigma levels
use soilsnow_m, only : fracice      ! Soil, snow and surface data
use tkeeps                          ! TKE-EPS boundary layer
use tracers_m, only : ngas          ! Tracer data
use trvmix, only : tracervmix       ! Tracer mixing routines
use work2_m                         ! Diagnostic arrays
      
implicit none
      
include 'newmpar.h'                 ! Grid parameters
include 'const_phys.h'              ! Physical constants
include 'kuocom.h'                  ! Convection parameters
include 'parm.h'                    ! Model configuration

integer, parameter :: ntest = 0
integer k, tnaero, nt
real rong, rlogs1, rlogs2, rlogh1, rlog12
real delsig, conflux, condrag
real, dimension(ifull,kl) :: tnhs, tv, zh
real, dimension(ifull,kl) :: rhs, guv, gt, rkm, rkh
real, dimension(ifull,kl) :: at, ct, au, cu, zg, cldtmp
real, dimension(ifull,kl) :: uav, vav
real, dimension(ifull,kl-1) :: tmnht, cnhs
real, dimension(ifull) :: ou, ov, iu, iv, rhos
real, dimension(ifull) :: dz, dzr
real, dimension(ifull) :: cgmap
real, dimension(kl) :: sighkap,sigkap,delons,delh

! Non-hydrostatic terms
tnhs(1:ifull,1) = phi_nh(:,1)/bet(1)
do k = 2,kl
  ! representing non-hydrostatic term as a correction to air temperature
  tnhs(1:ifull,k) = (phi_nh(:,k)-phi_nh(:,k-1)-betm(k)*tnhs(:,k-1))/bet(k)
end do
tv(1:ifull,1:kl) = t(1:ifull,:)*(1.+0.61*qg(1:ifull,:)-qlg(1:ifull,:)-qfg(1:ifull,:) &
                                -qrg(1:ifull,:)-qsng(1:ifull,:)-qgrg(1:ifull,:))

! Weight as a function of grid spacing for turning off CG term
!cgmap = 0.982, 0.5, 0.018 for 1000m, 600m, 200m when cgmap_offset=600 and cgmap_scale=200.
if ( cgmap_offset>0. ) then
  cgmap(1:ifull) = 0.5*tanh((ds/em(1:ifull)-cgmap_offset)/cgmap_scale) + 0.5 ! MJT suggestion
else
  cgmap(1:ifull) = 1.
end if

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

! Calculate half level heights, temperatures and NHS correction
rlogs1=log(sig(1))
rlogs2=log(sig(2))
rlogh1=log(sigmh(2))
rlog12=1./(rlogs1-rlogs2)
tmnht(:,1)=(tv(:,2)*rlogs1-tv(:,1)*rlogs2+(tv(:,1)-tv(:,2))*rlogh1)*rlog12
cnhs(:,1) =1.+(tnhs(:,2)*rlogs1-tnhs(:,1)*rlogs2+(tnhs(:,1)-tnhs(:,2))*rlogh1)*rlog12/tmnht(:,1)
! n.b. an approximate zh (in m) is quite adequate for this routine
zh(:,1)   =(tv(:,1)+tnhs(:,1))*delh(1)
do k = 2,kl-1
  zh(:,k)   =zh(:,k-1)+(tv(:,k)+tnhs(:,k))*delh(k)
  tmnht(:,k)=ratha(k)*tv(:,k+1)+rathb(k)*tv(:,k)                      ! MJT suggestion
  ! non-hydrostatic temperature correction at half level height
  cnhs(:,k) =1.+(ratha(k)*tnhs(:,k+1)+rathb(k)*tnhs(:,k))/tmnht(:,k)  ! MJT suggestion
end do      !  k loop
zh(:,kl)=zh(:,kl-1)+(tv(:,kl)+tnhs(:,kl))*delh(kl)                    ! MJT suggestion

! Adjustment for moving ocean surface
ou=0.
ov=0.
if ( nmlo/=0 ) then
  iu=0.
  iv=0.
  call mloexport(2,ou,1,0)
  call mloexport(3,ov,1,0)
  call mloexpice(iu, 9,0)
  call mloexpice(iv,10,0)
  ou=(1.-fracice)*ou+fracice*iu
  ov=(1.-fracice)*ov+fracice*iv
end if
      
if ( nvmix/=6 ) then

  !--------------------------------------------------------------
  ! JLM's local Ri scheme

  do k = 1,kl
    rhs(:,k)=t(1:ifull,k)*sigkap(k)  ! rhs is theta here
  enddo      !  k loop
    
  call vertjlm(rkm,rkh,rhs,sigkap,sighkap,delons,zh,tmnht,cnhs,ntest,cgmap)

  do k = 1,kl-1
    delsig  =(sig(k+1)-sig(k))
    dz(:)   =-tmnht(:,k)*delons(k)*cnhs(:,k)  ! this is z(k+1)-z(k)
    dzr(:)  =1./dz(:)
    guv(:,k)=rkm(:,k)*dt*delsig*dzr(:)**2
    gt(:,k) =rkh(:,k)*dt*delsig*dzr(:)**2
  end do      ! k loop
  guv(:,kl)=0.
  gt(:,kl) =0.

  if ( diag ) then
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
        guv(:,k)=guv(:,k)-convpsav(:)*.5 ! with factor of .5
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
    at(:,k)=-gt(:,k-1)/dsig(k)
  end do
  do k = 1,kl-1
    ct(:,k)=-gt(:,k)/dsig(k)
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
  rhs(:,1) = rhs(:,1) - (conflux/cp)*fg(:)/ps(1:ifull)
  call trim(at,ct,rhs)   ! for t
  if ( nmaxpr==1 .and. mydiag ) write (6,"('thet_out',9f8.3/8x,9f8.3)") rhs(idjd,:)
  do k = 1,kl
    t(1:ifull,k) = rhs(:,k)/sigkap(k)
  enddo    !  k loop
  if ( diag ) then
    if ( mydiag ) then
      write(6,*)'vertmix eg,fg,cdtq ',eg(idjd),fg(idjd),cdtq(idjd)
    end if
    call printa('thet',rhs,ktau,nlv,ia,ib,ja,jb,200.,1.)
  end if

  !--------------------------------------------------------------
  ! Moisture
  rhs=qg(1:ifull,:)
  rhs(:,1)=rhs(:,1)-(conflux/hl)*eg/ps(1:ifull)
  ! could add extra sfce moisture flux term for crank-nicholson
  call trim(at,ct,rhs)    ! for qg
  qg(1:ifull,:)=rhs
  if ( diag .and. mydiag ) then
    write(6,*)'vertmix rhs & qg after trim ',(rhs(idjd,k),k=1,kl)
    write (6,"('qg ',9f7.3/(8x,9f7.3))") (1000.*qg(idjd,k),k=1,kl)
  end if

  !--------------------------------------------------------------
  ! Cloud microphysics terms
  if ( ldr/=0 ) then
    ! now do qfg
    rhs=qfg(1:ifull,:)
    call trim(at,ct,rhs)       ! for qfg
    qfg(1:ifull,:)=rhs
    ! now do qlg
    rhs=qlg(1:ifull,:)
    call trim(at,ct,rhs)       ! for qlg
    qlg(1:ifull,:)=rhs
    if ( ncloud>=2 ) then
      ! now do qrg
      rhs=qrg(1:ifull,:)
      call trim(at,ct,rhs)      ! for qrg
      qrg(1:ifull,:)=rhs
      if ( ncloud>=3 ) then
        ! now do qsng
        rhs=qsng(1:ifull,:)
        call trim(at,ct,rhs)      ! for qsng
        qsng(1:ifull,:)=rhs
        ! now do qgrg
        rhs=qgrg(1:ifull,:)
        call trim(at,ct,rhs)      ! for qgrg
        qgrg(1:ifull,:)=rhs
        ! now do rfrac
        rhs=rfrac(1:ifull,:)
        call trim(at,ct,rhs)      ! for rfrac
        rfrac(1:ifull,:)=min(max(rhs,0.),1.)
        ! now do sfrac
        rhs=sfrac(1:ifull,:)
        call trim(at,ct,rhs)      ! for sfrac
        sfrac(1:ifull,:)=min(max(rhs,0.),1.)
        ! now do gfrac
        rhs=gfrac(1:ifull,:)
        call trim(at,ct,rhs)      ! for gfrac
        gfrac(1:ifull,:)=min(max(rhs,0.),1.)
        if ( ncloud>=4 ) then
          ! now do cldfrac
          rhs=stratcloud(1:ifull,:)
          call trim(at,ct,rhs)    ! for cldfrac
          stratcloud(1:ifull,:)=rhs
          call combinecloudfrac
        else
          ! now do cfrac
          rhs=cfrac(1:ifull,:)
          call trim(at,ct,rhs)    ! for cfrac
          cfrac(1:ifull,:)=rhs
        end if  ! (ncloud>=4)
      end if    ! (ncloud>=3)
    end if      ! (ncloud>=2)
  end if        ! (ldr/=0)
      
  !--------------------------------------------------------------
  ! Aerosols
  if ( abs(iaero)==2 ) then
    do nt = 1,naero
      rhs(:,:) = xtg(1:ifull,:,nt) ! Total grid-box
      call trim(at,ct,rhs)
      xtg(1:ifull,:,nt) = rhs(:,:)
    end do
  end if ! (abs(iaero)==2)

  !--------------------------------------------------------------
  ! Momentum terms
  au(:,1) = cduv(:)*condrag/tss(:)
  cu(:,kl) = 0.
  do k = 2,kl
    au(:,k) = -guv(:,k-1)/dsig(k)
  enddo    !  k loop
  do k = 1,kl-1
    cu(:,k) = -guv(:,k)/dsig(k)
  enddo    !  k loop
  if ( ( diag .or. ntest==2 ) .and. mydiag ) then
    write(6,*)'au ',(au(idjd,k),k=1,kl)
    write(6,*)'cu ',(cu(idjd,k),k=1,kl)
  end if      ! (ntest==2)

  ! first do u
  do k = 1,kl
    rhs(:,k) = u(1:ifull,k) - ou(:)
  end do
  call trim(au,cu,rhs)
  do k = 1,kl
    u(1:ifull,k) = rhs(:,k) + ou(:)
  end do
  if ( diag .and. mydiag ) then
    write(6,*)'vertmix au ',(au(idjd,k),k=1,kl)
  end if

  ! now do v; with properly unstaggered au,cu
  do k = 1,kl
    rhs(:,k) = v(1:ifull,k) - ov(:)
  end do
  call trim(au,cu,rhs)    ! note now that au, cu unstaggered globpea
  do k = 1,kl
    v(1:ifull,k) = rhs(:,k) + ov(:)
  end do

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

  !--------------------------------------------------------------
  ! Tracers
  if ( ngas>0 ) then
    call tracervmix(at,ct)
  end if ! (ngas>0)
      
else
      
  !-------------------------------------------------------------
  ! k-e + MF closure scheme
      
  ! note ksc/=0 options are clobbered when nvmix=6
  ! However, nvmix=6 with nlocal=7 supports its own shallow
  ! convection options
       
  ! calculate height on full levels
  zg(:,1)=bet(1)*tv(:,1)/grav
  do k = 2,kl
    zg(:,k)=zg(:,k-1)+(bet(k)*tv(:,k)+betm(k)*tv(:,k-1))/grav
  end do ! k  loop
  zg=zg+phi_nh/grav ! add non-hydrostatic component
       
  ! near surface air density (see sflux.f and cable_ccam2.f90)
  rhos=sig(1)*ps(1:ifull)/(rdry*t(1:ifull,1))
  
  ! Use counter gradient for aerosol tracers
  if ( abs(iaero)==2 ) then 
    tnaero=naero
  else
    tnaero=0
  end if
  
  ! Special treatment for prognostic cloud fraction
  if ( ncloud>=4 ) then
    cldtmp=stratcloud(1:ifull,:)
  else
    cldtmp=cfrac(1:ifull,:)
  end if
       
  ! transform to ocean reference frame and temp to theta
  do k = 1,kl
    u(1:ifull,k)=u(1:ifull,k)-ou
    v(1:ifull,k)=v(1:ifull,k)-ov
    rhs(:,k)=t(1:ifull,k)*sigkap(k) ! theta
  end do

  ! Evaluate EDMF scheme
  select case(nlocal)
    case(0) ! no counter gradient
      call tkemix(rkm,rhs,qg,qlg,qfg,qrg,qsng,qgrg,cldtmp,rfrac,sfrac,gfrac,u,v, &
                  pblh,fg,eg,ps,zo,zg,zh,sig,rhos,dt,qgmin,1,0,tnaero,xtg,cgmap)
      rkh=rkm
    case(1,2,3,4,5,6) ! KCN counter gradient method
      call tkemix(rkm,rhs,qg,qlg,qfg,qrg,qsng,qgrg,cldtmp,rfrac,sfrac,gfrac,u,v, &
                  pblh,fg,eg,ps,zo,zg,zh,sig,rhos,dt,qgmin,1,0,tnaero,xtg,cgmap)
      rkh=rkm
      do k=1,kl
        uav(1:ifull,k)=av_vmod*u(1:ifull,k)+(1.-av_vmod)*(savu(1:ifull,k)-ou)
        vav(1:ifull,k)=av_vmod*v(1:ifull,k)+(1.-av_vmod)*(savv(1:ifull,k)-ov)
      end do
      call pbldif(rhs,rkh,rkm,uav,vav,cgmap)
    case(7) ! mass-flux counter gradient
      call tkemix(rkm,rhs,qg,qlg,qfg,qrg,qsng,qgrg,cldtmp,rfrac,sfrac,gfrac,u,v, &
                  pblh,fg,eg,ps,zo,zg,zh,sig,rhos,dt,qgmin,0,0,tnaero,xtg,cgmap)
      rkh=rkm
    case DEFAULT
      write(6,*) "ERROR: Unknown nlocal option for nvmix=6"
      call ccmpi_abort(-1)
  end select
  
  ! special treatment for prognostic cloud fraction  
  if ( ncloud>=4 ) then
    stratcloud(1:ifull,:)=cldtmp
    call combinecloudfrac
  else
    cfrac(1:ifull,:)=cldtmp
  endif
  
  ! transform winds back to Earth reference frame and theta to temp
  do k = 1,kl
    u(1:ifull,k)=u(1:ifull,k)+ou
    v(1:ifull,k)=v(1:ifull,k)+ov
    t(1:ifull,k)=rhs(1:ifull,k)/sigkap(k)
  enddo    !  k loop

  ! tracers
  if ( ngas>0 ) then
    do k = 1,kl-1
      delsig       =sig(k+1)-sig(k)
      dz(1:ifull)  =-tmnht(1:ifull,k)*delons(k)*cnhs(1:ifull,k)  ! this is z(k+1)-z(k)
      dzr(1:ifull) =1./dz(1:ifull)
      gt(1:ifull,k)=rkh(1:ifull,k)*dt*delsig*dzr(1:ifull)**2
    end do      ! k loop
    gt(:,kl)=0.
    at(:,1) =0.
    ct(:,kl)=0.
    do k = 2,kl
      at(1:ifull,k)=-gt(1:ifull,k-1)/dsig(k)
    end do
    do k = 1,kl-1
      ct(1:ifull,k)=-gt(1:ifull,k)/dsig(k)
    end do
    ! increase mixing to replace counter gradient term      
    if ( nlocal>0 ) then
      at=cq*at
      ct=cq*ct
    end if
    call tracervmix(at,ct)
  end if
       
end if ! nvmix/=6 ..else..
      
return
end subroutine vertmix

subroutine vertjlm(rkm,rkh,rhs,sigkap,sighkap,delons,zh,tmnht,cnhs,ntest,cgmap)

use arrays_m                        ! Atmosphere dyamics prognostic arrays
use cc_mpi                          ! CC MPI routines
use cfrac_m                         ! Cloud fraction
use diag_m                          ! Diagnostic routines
use estab                           ! Liquid saturation function
use kuocomb_m                       ! JLM convection
use liqwpar_m                       ! Cloud water mixing ratios
use morepbl_m                       ! Additional boundary layer diagnostics
use savuvt_m                        ! Saved dynamic arrays
use screen_m                        ! Screen level diagnostics
use sigs_m                          ! Atmosphere sigma levels
use soil_m, only : land             ! Soil and surface data

implicit none

include 'newmpar.h'                 ! Grid parameters
include 'const_phys.h'              ! Physical constants
include 'kuocom.h'                  ! Convection parameters
include 'parm.h'                    ! Model configuration

integer, parameter :: ndvmod=0    ! 0 default, 1+ for dvmod tests
integer, intent(in) :: ntest
integer iq,k,iqmax
integer, save :: kscbase,ksctop
integer, dimension(ifull) :: kbase,ktop
real, parameter :: bprm=4.7,cm=7.4,ch=5.3,amxlsq=100.   ! coefficients for Louis scheme
real, parameter :: vkar3=0.35,vkar4=0.4,bprmj=5.,cmj=5. ! coefficients for Louis scheme
real, parameter :: chj=2.6                              ! coefficients for Louis scheme
real delta,es,pk,dqsdt,betat,betaq,betac,al,qc,fice
real w1,w2,diffmax,rhsk,rhskp,delthet_old,xold,diff
real denma,denha,esp,epart,tsp,qbas
real, dimension(ifull,kl), intent(inout) :: rkm,rkh,rhs
real, dimension(ifull,kl), intent(in) :: zh
real, dimension(ifull,kl-1), intent(in) :: tmnht,cnhs
real, dimension(ifull,kl) :: qs,betatt,betaqt,delthet
real, dimension(ifull,kl) :: uav,vav,ri,rk_shal
real, dimension(ifull,kl) :: thee,thebas
real, dimension(ifull), intent(in) :: cgmap
real, dimension(ifull) :: dz,dzr,dvmod,dqtot,x,zhv
real, dimension(ifull) :: csq,sqmxl,fm,fh,sigsp
real, dimension(ifull) :: theeb
real, dimension(kl), intent(in) :: sigkap,sighkap,delons
real, dimension(:), allocatable, save :: prcpv

w1=0.
pk=0.

if ( .not.allocated(prcpv) ) allocate(prcpv(kl))

if ( nmaxpr==1 .and. mydiag ) write (6,"('thet_in',9f8.3/7x,9f8.3)") rhs(idjd,:)
  
! Pre-calculate the buoyancy parameters if using qcloud scheme.
! Follow Smith's (1990) notation; gam() is HBG's notation for (L/cp)dqsdt.
! The factor of (1/sigkap)=T/theta in betatt differs from Smith's formulation
! because we use theta derivative rather than (dry static energy)/cp.
if ( nvmix>0 .and. nvmix<4 ) then
  delta=1./epsil-1.  ! i.e. 1/.622 -1., i.e. .6077
  do k=1,kl
    if ( sig(k)>.8 ) then ! change made 17/1/06
      do iq=1,ifull
        es=establ(t(iq,k))
        pk=ps(iq)*sig(k)
        qs(iq,k)=.622*es/(pk-es)  
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
      do iq=1,ifull
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
    do iq=1,ifull
      es=establ(t(iq,k))
      qs(iq,k)=.622*es/max(1.,ps(iq)*sig(k)-es)  ! max for k=kl
    enddo   ! iq loop
  enddo    !  k loop
endif      ! (nvmix>0.and.nvmix<4)

do k=1,kl-1
  delthet(:,k)=rhs(:,k+1)-rhs(:,k)  ! rhs is theta or thetal here
enddo      !  k loop

if(ktau==1.and.ksc/=0)then
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

!     ****** section for Geleyn shallow convection; others moved lower****
if(ksc==-99)then
  do k=kscbase,ksctop    ! new usage of ksc thu  02-17-2000
    do iq=1,ifull
      delthet(iq,k)=delthet(iq,k)-hlcp*max(0.,qs(iq,k+1)-qg(iq,k+1)-qs(iq,k)+qg(iq,k) )
    enddo  ! iq loop
  enddo   !  k loop
endif     ! (ksc==-99)
if(ksc==-98)then    ! modified Geleyn Jan '08
  do k=kscbase,ksctop    ! new usage of ksc thu  02-17-2000
    do iq=1,ifull
      if(qg(iq,k)>tied_rh*qs(iq,k).or.qg(iq,k+1)>tied_rh*qs(iq,k+1))then
        delthet(iq,k)=delthet(iq,k)-hlcp*max(0.,qs(iq,k+1)-qg(iq,k+1)-qs(iq,k)+qg(iq,k) )
      endif
    enddo  ! iq loop
  enddo   !  k loop
endif     ! (ksc==-98)
if(ksc==-97)then    ! modified Geleyn Jan '08
  do k=kscbase,ksctop    ! new usage of ksc thu  02-17-2000
    do iq=1,ifull
      if(qg(iq,k)>tied_rh*qs(iq,k))then
        delthet(iq,k)=delthet(iq,k)-hlcp*max(0.,qs(iq,k+1)-qg(iq,k+1)-qs(iq,k)+qg(iq,k) )
      endif
    enddo  ! iq loop
  enddo   !  k loop
endif     ! (ksc==-97)
if(ksc==-96)then   ! combination of Geleyn and jlm 83 (Feb 08)
  do k=1,ksctop    
    do iq=1,ifull
      if(k<ktsav(iq).and.k>=kbsav(iq).and.condc(iq)==0.)then  
        delthet(iq,k)=delthet(iq,k)-hlcp*max(0.,qs(iq,k+1)-qg(iq,k+1)-qs(iq,k)+qg(iq,k) )
        write(6,*)'-96 iq,k,kbsav,ktsav ',iq,k,kbsav(iq),ktsav(iq),hlcp*(qs(iq,k+1)-qg(iq,k+1)-qs(iq,k)+qg(iq,k)),delthet(iq,k)
      endif 
    enddo  ! iq loop
  enddo   !  k loop
endif     ! (ksc==-96)
if(ksc==-95)then ! same as -99 but has tied_rh (e.g. .75) factor
  ! kshal(:)=0
  do k=kscbase,ksctop     
    do iq=1,ifull
      delthet(iq,k)=delthet(iq,k)-tied_rh*hlcp*max(0.,(qs(iq,k+1)-qg(iq,k+1)-qs(iq,k)+qg(iq,k)) )
    enddo  ! iq loop
  enddo   !  k loop
endif     ! (ksc==-95)
if(ksc==-94)then   ! combination of Geleyn and jlm 83 (Feb 08)
  do k=1,ksctop    
    do iq=1,ifull
      if(k<ktsav(iq).and.k>=kbsav(iq).and.condc(iq)==0.)then  
        delthet(iq,k)=0.
        write(6,*)'-94 iq,k,kbsav,ktsav ',iq,k,kbsav(iq),ktsav(iq),hlcp*(qs(iq,k+1)-qg(iq,k+1)-qs(iq,k)+qg(iq,k)),delthet(iq,k)
      endif 
    enddo  ! iq loop
  enddo   !  k loop
endif     ! (ksc==-94)
if(ksc==-93)then ! single-layer (pblh) version of -95
  do iq=1,ifull
    do k=kscbase,kl/2
      if(zh(iq,k)<pblh(iq).and.zh(iq,k+1)>pblh(iq))then
        delthet(iq,k)=delthet(iq,k)-tied_rh*hlcp*max(0.,(qs(iq,k+1)-qg(iq,k+1)-qs(iq,k)+qg(iq,k)) )
      endif
    enddo  ! k loop 
  enddo   ! iq loop
endif     ! (ksc==-93)
if(ksc==-92)then ! capped-by-pblh version of -95
  do iq=1,ifull
    do k=kscbase,kl/2
      if(zh(iq,k)<pblh(iq))then
        delthet(iq,k)=delthet(iq,k)-tied_rh*hlcp*max(0.,(qs(iq,k+1)-qg(iq,k+1)-qs(iq,k)+qg(iq,k)) )
      endif
    enddo  ! k loop 
  enddo   ! iq loop
endif     ! (ksc==-92)
if(ksc==-91)then ! capped-by-pblh (anywhere in layer) version of -95
  do iq=1,ifull
    do k=2,kl/2
      if(zh(iq,k-1)<pblh(iq))then
        delthet(iq,k)=delthet(iq,k)-tied_rh*hlcp*max(0.,(qs(iq,k+1)-qg(iq,k+1)-qs(iq,k)+qg(iq,k)) )
      endif
    enddo  ! k loop 
  enddo   ! iq loop
endif     ! (ksc==-91)
!     ********* end of Geleyn shallow convection section ****************

! following now defined in vertmix (don't need to pass from sflux)
uav(1:ifull,:)=av_vmod*u(1:ifull,:)+(1.-av_vmod)*savu(1:ifull,:) 
vav(1:ifull,:)=av_vmod*v(1:ifull,:)+(1.-av_vmod)*savv(1:ifull,:) 
do k=1,kl-1
  do iq=1,ifull
    dz(iq) =-tmnht(iq,k)*delons(k)*cnhs(iq,k)  ! this is z(k+1)-z(k)
    dzr(iq)=1./dz(iq)
    zhv(iq)=1./zh(iq,k)
    if(ndvmod==0)then
      dvmod(iq)=sqrt( (uav(iq,k+1)-uav(iq,k))**2+(vav(iq,k+1)-vav(iq,k))**2 )
    else
      dvmod(iq)=ndvmod  ! just for tests
    endif    ! (ndvmod==0)
  enddo ! iq loop

  ! x is bulk ri *(dvmod **2); used to calc ri(k), rkh(k) etc
  if(nvmix>0.and.nvmix<4)then  ! new one allowing for cloudy air
    ! usually nvmix=3       
    if(sig(k)>.8)then ! change made 17/1/06
      dqtot(:)=qg(1:ifull,k+1)+qlg(1:ifull,k+1)+qfg(1:ifull,k+1)-(qg(1:ifull,k)  +qlg(1:ifull,k)  +qfg(1:ifull,k))
    else
      dqtot(:)=0.
    endif
    if(nvmix==2)then !  jlm May '05
      do iq=1,ifull
        x(iq)=grav*dz(iq)*((min(betatt(iq,k),betatt(iq,k+1)))*delthet(iq,k) + (max(betaqt(iq,k),betaqt(iq,k+1)))*dqtot(iq) )
      enddo ! iq loop	
    else    ! i.e. nvmix=1 or 3
      if(nvmix==1)w1=dsig(k+1)/(dsig(k)+dsig(k+1)) 
      if(nvmix==3)w1=1.    !weight for lower level  usual           
      w2=1.-w1             !weight for upper level
      do iq=1,ifull
        x(iq)=grav*dz(iq)*((w1*betatt(iq,k)+w2*betatt(iq,k+1))*delthet(iq,k) + (w1*betaqt(iq,k)+w2*betaqt(iq,k+1))*dqtot(iq) )
      enddo ! iq loop	 
    endif  !  (nvmix==2) .. else ..
    if(ntest==4.and.k<=9.and.ktau==ntau)then
      diffmax=0.
      do iq=1,ifull
        rhsk=t(iq,k)*sigkap(k)
        rhskp=t(iq,k+1)*sigkap(k+1)
        delthet_old=rhs(iq,k+1)-rhs(iq,k)
        xold=grav*dz(iq)*(delthet_old/(tmnht(iq,k)*sighkap(k))+.61*(qg(iq,k+1)-qg(iq,k)))
        diff=abs(xold-x(iq))
        if(diff>diffmax)then
          diffmax=diff
          iqmax=iq
        endif
        write(47,'(3g13.4,i7,i4)') xold,x(iq),diff,iq,k
      enddo
      write(6,*)'k,iqmax,diffmax ',k,iqmax,diffmax
    endif   ! (ntest==4.and.k<=9.and.ktau==ntau)
    rhs(:,k)=t(1:ifull,k)*sigkap(k)   !need to re-set theta for nvmix=1-3
  elseif(nvmix==5)then        ! non-cloudy x with qfg, qlg
    x(:)=grav*dz(:)*(delthet(:,k)/(tmnht(:,k)*sighkap(k))+.61*(qg(1:ifull,k+1)-qg(1:ifull,k)) &
         -qfg(1:ifull,k+1)-qlg(1:ifull,k+1)+qfg(1:ifull,k)+qlg(1:ifull,k) )
  else                 ! original non-cloudy x, nvmix=4
    x(:)=grav*dz(:)*(delthet(:,k)/(tmnht(:,k)*sighkap(k))+.61*(qg(1:ifull,k+1)-qg(1:ifull,k)))
  endif     ! (nvmix>0.and.nvmix<4) .. else ..

  ! fm and fh denote f(Louis style)*dvmod
  ! nb. an error exists in the Louis expression for c; this is corrected
  ! in the current code
  csq(:) = zhv(:)*(((1.+dz(:)*zhv(:))**(1./3.)-1.)*dzr(:))**3

  ! newest code, stable same as csiro9 here (originally for nvmix=4)
  do iq=1,ifull
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
    do iq=1,ifull
      fm(iq)=1.
      fh(iq)=1.
    enddo   ! iq loop
    if(k<kl/2)then
      do iq=1,ifull
        if(land(iq))fh(iq)=abs(nvmix)
      enddo   ! iq loop
    endif
  endif     !  (nvmix<0)

  ! calculate k's, guv and gt
  ! nonlocal scheme usually applied to temperature and moisture, not momentum
  ! (i.e. local scheme is applied to momentum for nlocal=0,1)
  if(nvmix/=0)then    ! use nvmix=0 only for special test runs
    do iq=1,ifull
      rkm(iq,k)=fm(iq)*sqmxl(iq)*dzr(iq)
      rkh(iq,k)=fh(iq)*sqmxl(iq)*dzr(iq)
    enddo   ! iq loop
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
  call pbldif(rhs,rkh,rkm,uav,vav,cgmap)  ! rhs is theta or thetal
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
endif      ! (nlocal>0)

rk_shal(:,:)=0.
!     ***** ***** section for jlm shallow convection v4 *****************
if(ksc==81)then
  do k=1,ksctop-1   ! or ksctop?  
    do iq=1,ifull
      if(sig(ktsav(iq))>sig_ct.and.k<ktsav(iq).and.condc(iq)==0.)then  
        rk_shal(iq,k)=tied_con
        rk_shal(iq,k+1)=tied_over
      endif ! (sig(ktsav(iq))>sig_ct.and.k<ktsav(iq).and....)
    enddo  ! iq loop
  enddo   !  k loop
endif     ! (ksc==81)
if(ksc==82)then
  do k=1,ksctop-1    
    do iq=1,ifull
      if(sig(ktsav(iq))>sig_ct.and.k<ktsav(iq).and.k>=kbsav(iq).and.condc(iq)==0.)then  
        rk_shal(iq,k)=tied_con
        rk_shal(iq,k+1)=tied_over
      endif ! (sig(ktsav(iq))>sig_ct.and.k<ktsav(iq).and....)
    enddo  ! iq loop
  enddo   !  k loop
endif     ! (ksc==82)
if(ksc==83)then
  do k=1,ksctop-1    
    do iq=1,ifull
      if(sig(k)>sig_ct.and.k<ktsav(iq).and.k>=kbsav(iq).and.condc(iq)==0.)then  
        rk_shal(iq,k)=tied_con
        rk_shal(iq,k+1)=tied_over
      endif ! (sig(ktsav(iq))>sig_ct.and.k<ktsav(iq).and....)
    enddo  ! iq loop
  enddo   !  k loop
endif     ! (ksc==83)
if(ksc==91)then
  do k=1,ksctop ! May 08
    do iq=1,ifull
      if(ktsav(iq)<kl-1.and.k<ktsav(iq))then  ! April 04
        rk_shal(iq,k)=tied_con
        rk_shal(iq,k+1)=tied_over
      endif ! (ktsav(iq)<0.and.k<abs(ktsav(iq)))
    enddo  ! iq loop
  enddo   !  k loop
endif     ! (ksc==91)
if(ksc==92)then
  do k=1,ksctop ! May 08
    do iq=1,ifull
      if(k>=kbsav(iq).and.k<ktsav(iq).and.condc(iq)==0.)then  ! May 08
        rk_shal(iq,k)=tied_con
        rk_shal(iq,k+1)=tied_over
      endif ! (ktsav(iq)<0.and.k<abs(ktsav(iq)))
    enddo  ! iq loop
  enddo   !  k loop
endif     ! (ksc==92)
!     *********** end of jlm shallow convection section *****************

!     ************ section for Tiedtke shallow convection *******************
if(ksc==99)then
  do iq=1,ifull
    theeb(iq)=prcpv(kscbase)*t(iq,kscbase)*(t(iq,kscbase)+.5*hlcp*qs(iq,kscbase))/(t(iq,kscbase)-.5*hlcp*qs(iq,kscbase))
  enddo    ! iq loop
  if(kscsea==1)then  ! Tiedtke done only over sea
    do k=kscbase+1,ksctop
      do iq=1,ifull
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
      do iq=1,ifull
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
    do iq=1,ifull
      ! Leon used factor 1.01 over sea; here div by tied_rh (e.g. .99)	 
      epart=qg(iq,1)*.01*ps(iq)/(tied_rh*epsil) !in hPa 
      tsp=2840./(3.5*log(tscrn(iq))-log(epart)-4.805) + 55.   
      sigsp(iq)=(tsp/tscrn(iq))**(1./roncp)
    enddo
  else                    ! uses qgscrn
    do iq=1,ifull
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
    do iq=1,ifull
      if(sigsp(iq)>sigmh(k+1))kbase(iq)=k
    enddo
  enddo
  if(nlocal/=0)then  ! like old ksc=95, but ensures LCL within PBL
    do iq=1,ifull
      if(zh(iq,kbase(iq))>pblh(iq))then 
        kbase(iq)=kl
      endif
    enddo  ! iq loop
  endif  ! (nlocal/=0)
! following has some jlm shortcuts	 
  do iq=1,ifull
    k=kbase(iq)
    qbas=qs(iq,k)
    theeb(iq)=prcpv(k)*t(iq,k)*(t(iq,k) + .5*hlcp*qbas)/(t(iq,k) - .5*hlcp*qbas)
  enddo  ! iq loop
  ktop(:)=kbase(:)
  do k=2,ksctop+1
    do iq=1,ifull
      thee(iq,k)=prcpv(k)*t(iq,k)*(t(iq,k) + .5*hlcp*qs(iq,k))/(t(iq,k) - .5*hlcp*qs(iq,k))
      if(theeb(iq)>thee(iq,k).and.ktop(iq)==k-1)then
        ktop(iq)=k       ! also checking if contiguous
      endif
    enddo  ! iq loop
  enddo   ! k loop	 
  if(ksc==94)then  ! from April same as 93* 
    do k=2,ksctop
      do iq=1,ifull
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
      do iq=1,ifull
        if(ktop(iq)>kbase(iq).and.k<=ktop(iq).and.k>kbase(iq).and.ktop(iq)<=ksctop)then
          rk_shal(iq,k-1)=tied_con     !  m**2/sec  6., originally 10.
          rk_shal(iq,k)=tied_over      !  m**2/sec
        endif
      enddo  ! iq loop
    enddo   ! k loop
  endif     ! (ksc==93)
  if(ksc==95)then   ! new from 7 April
    do k=2,ksctop
      do iq=1,ifull
        if(k>kbase(iq).and.k<=ktop(iq).and.condc(iq)==0.)then
          rk_shal(iq,k-1)=tied_con     !  m**2/sec  6., originally 10.
          rk_shal(iq,k)=tied_over      !  m**2/sec
        endif
      enddo  ! iq loop
    enddo   ! k loop
  endif     ! (ksc==95)
  if(ksc==96)then   ! applied from sfce up
    do k=2,ksctop
      do iq=1,ifull
        if(ktop(iq)>kbase(iq).and.k<=ktop(iq).and.condc(iq)==0.)then  
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
    do iq=1,ifull
      qbas=min(qg(iq,k)/tied_rh,qs(iq,k))
!     thebas(iq,k)=t(iq,k)+hlcp*qbas  + hght
      thebas(iq,k)=prcpv(k)*t(iq,k)*(t(iq,k) + .5*hlcp*qbas)/(t(iq,k) - .5*hlcp*qbas)
    enddo  ! iq loop
  enddo   ! k loop
  theeb(:)=thebas(:,kscbase)
  do k=kscbase+1,ksctop+1
    do iq=1,ifull
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
    do iq=1,ifull
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
      
if(ksc/=0.and.(ntest/=0.or.diag).and.nproc==1.)then
  do iq=1,ifull
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

return
end subroutine vertjlm
  
subroutine trim(a,c,rhs)

implicit none

include 'newmpar.h'

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