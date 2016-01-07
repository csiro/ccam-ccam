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
    
! CCAM nudging/assimilation routines
      
! These routines preturb the regional model with the large scale circulation of the host model.
! Currently, relaxation, far-field and scale-selective filter options are supported for both
! the atmosphere and ocean.
      
! We support both 1D and 2D versions of the scale-selective filter.  2D is exact, but expensive.
! Current tests suggest the 1D is a good approximation of the 2D filter.
      
! Now use global sparse arrays with 1D filter (see globalpack in cc_mpi.f90).  This considerably
! reduces the memory used by the filter.  Nevertheless, MPI_Win_Allocate_Shared in MPI-3 can
! still be used to make shared copies of some arrays such as x_g, y_g, z_g and em_g.

! nbd/=0       Far-field or relaxation nudging
! mbd/=0       Spectral filter (1D and 2D versions, see nud_uv)
! nud_uv =1    Nudge winds (=9 for 2D filter)
! nud_t  =1    Nudge air temperature
! nud_qg =1    Nudge mixing ratio
! nud_p  =1    Nudge surface pressure
! nud_sst=1    Nudge water temperature (numbers greater than mbd control strength)
! nud_sss=1    Nudge salinity
! nud_ouv=1    Nudge water currents
! nud_sfh=1    Nudge water free surface height
! kbotdav      Lowest atmospheric level to nudge
! ktopdav      Highest atmospheric level to nudge
! ktopmlo      Deepest water level to nudge
! kbotmlo      Shallowest water level to nudge
! mloalpha     Weight of water nudging strength

module nesting

private
public nestin, nestinb, mlofilterhub, mlonudge, specinit

contains

!--------------------------------------------------------------
! FAR-FIELD NUDGING AND RELAXATION ROUTINES
! Called for nbd/=0
subroutine nestin
      
use aerosolldr                   ! LDR prognostic aerosols
use arrays_m                     ! Atmosphere dyamics prognostic arrays
use cc_mpi                       ! CC MPI routines
use daviesnudge                  ! Far-field nudging
use diag_m                       ! Diagnostic routines
use indices_m                    ! Grid index arrays
use latlong_m                    ! Lat/lon coordinates
use mlo                          ! Ocean physics and prognostic arrays
use onthefly_m                   ! Input interpolation routines
use pbl_m                        ! Boundary layer arrays
use soil_m                       ! Soil and surface data
use soilsnow_m                   ! Soil, snow and surface data
      
implicit none
      
include 'newmpar.h'              ! Grid parameters
include 'dates.h'                ! Date data
include 'parm.h'                 ! Model configuration
include 'stime.h'                ! File date data

integer, dimension(ifull) :: dumm
integer, save :: mtimea = -1
integer, save :: mtimeb = -1
integer, save :: wl = -1
integer i
integer kdate_r, ktime_r, kdhour, kdmin
real timerm, cona, conb
real, dimension(2) :: depthcheck
real, dimension(:,:), allocatable, save :: ta, ua, va, qa
real, dimension(:,:), allocatable, save :: tb, ub, vb, qb, ocndep
real, dimension(:), allocatable, save :: psla, pslb, tssa, tssb
real, dimension(:), allocatable, save :: sicedepb, fraciceb
real, dimension(:,:,:), allocatable, save :: sssa, sssb
real, dimension(:,:,:), allocatable, save :: xtghosta, xtghostb
real, dimension(ifull) :: zsb, timelt
real, dimension(ifull,3) :: duma
real, dimension(ifull,wlev,4) :: dumaa
real, dimension(ifull,ms,3) :: dumg
real, dimension(ifull,kl,5) :: dumv
real, dimension(ifull,3,3) :: dums
      
!     mtimer, mtimeb are in minutes
if ( ktau<100 .and. myid==0 ) then
  write(6,*) 'in nestin ktau,mtimer,mtimea,mtimeb ',ktau,mtimer,mtimea,mtimeb
  write(6,*) 'with kdate_s,ktime_s >= ',kdate_s,ktime_s
end if

! Load next host model timestep for nudging
if( mtimer>mtimeb ) then  ! allows for dt<1 minute
      
  ! Intialise nudging
  if ( .not.allocated(ta) ) then
    ! Allocate host data arrays
    allocate( ta(ifull,kl), ua(ifull,kl), va(ifull,kl), qa(ifull,kl) )
    allocate( tb(ifull,kl), ub(ifull,kl), vb(ifull,kl), qb(ifull,kl) )
    allocate( psla(ifull), pslb(ifull), tssa(ifull), tssb(ifull) )
    allocate( sicedepb(ifull), fraciceb(ifull) )
    allocate( sssa(ifull,wlev,4), sssb(ifull,wlev,4), ocndep(ifull,2) )
    allocate( xtghosta(ifull,kl,naero) )
    allocate( xtghostb(ifull,kl,naero) )

    ! Save host atmospheric data
    if ( myid==0 ) write(6,*) 'set nesting fields to those already read in via indata'
    pslb(:) = psl(:)
    tssb(:) = tss(:)
    sicedepb(:) = sicedep(:)
    fraciceb(:) = fracice(:)
    tb(:,:) = t(1:ifull,:)
    qb(:,:) = qg(1:ifull,:)
    ub(:,:) = u(1:ifull,:)
    vb(:,:) = v(1:ifull,:)

    ! Save host ocean data
    if ( nmlo/=0 ) then
      ocndep(:,:) = 0.
      sssb(:,:,1) = 293.16
      sssb(:,:,2) = 34.72
      sssb(:,:,3) = 0.
      sssb(:,:,4) = 0.
      do i = 1,4
        call mloexport3d(i-1,sssb(:,:,i),0)
      end do
    end if
          
    ! Save host aerosol data
    if ( abs(iaero)>=2 .and. nud_aero/=0 ) then
      xtghostb(:,:,:) = xtg(1:ifull,:,:)
    end if
        
    call setdavvertwgt
    
    ! record time of saved data
    mtimeb = mtimer
  endif       ! (.not.allocated(ta))
      
! transfer mtimeb fields to mtimea and update sice variables
  mtimea = mtimeb
  psla(:) = pslb(:)
  tssa(:) = tssb(:)
  ta(1:ifull,:) = tb(1:ifull,:)
  qa(1:ifull,:) = qb(1:ifull,:)
  ua(1:ifull,:) = ub(1:ifull,:)
  va(1:ifull,:) = vb(1:ifull,:)
  if ( nmlo/=0 ) then
    sssa(:,:,:) = sssb(:,:,:)
  end if
  if ( abs(iaero)>=2 .and. nud_aero/=0 ) then
    xtghosta(:,:,:) = xtghostb(:,:,:)
  end if

  ! Read sea-ice data from host when not using
  ! AMIP SSTs or Mixed-Layer-Ocean sea-ice      
  if ( namip==0 .and. nmlo==0 ) then
    ! check whether present ice points should change to/from sice points
    where ( fraciceb(:)>0. .and. fracice(:)==0 .and. .not.land(1:ifull) )
      ! N.B. if already a sice point, keep present tice (in tggsn)
      tggsn(:,1) = min( 271.2, tssb(:) )
    end where
    sicedep(:) = sicedepb(:)
    fracice(:) = fraciceb(:)
    ! ensure that sice is only over sea
    where ( fracice(:)<0.02 )
      fracice(:) = 0.
    end where
    where ( land(1:ifull) )
      sicedep(:) = 0.
      fracice(:) = 0.
    elsewhere ( fracice(:)>0. .and. sicedep(:)==0. .and. rlatt(:)>0. )
      ! assign to 2m in NH and 1m in SH (according to spo)
      sicedep(:) = 2.
    elsewhere ( fracice(:)>0. .and. sicedep(:)==0. )
      ! assign to 2m in NH and 1m in SH (according to spo)
      sicedep(:) = 1.
    elsewhere ( fracice(:)==0. .and. sicedep(:)>0. )
        fracice(:) = 1.
    end where
  endif ! (namip==0.and.nmlo==0)

  ! Read host atmospheric and ocean data for nudging      
  if ( abs(io_in)==1 ) then
    call onthefly(1,kdate_r,ktime_r,                            &
                  pslb,zsb,tssb,sicedepb,fraciceb,tb,ub,vb,qb,  &
                  dumg(:,:,1),dumg(:,:,2),dumg(:,:,3),          &
                  duma(:,1),dumv(:,:,1),dumv(:,:,2),            &
                  dumv(:,:,3),dumv(:,:,4),dumv(:,:,5),          &
                  dums(:,:,1),dums(:,:,2),dums(:,:,3),          &
                  duma(:,2),duma(:,3),dumm,sssb,ocndep,xtghostb)
  else
    write(6,*) 'ERROR: Nudging requires abs(io_in)=1'
    call ccmpi_abort(-1)
  endif   ! (io_in==1)
  tssb(:) = abs(tssb(:))
#ifdef debug
  if ( mydiag ) then
    write (6,"('zsb# nestin  ',9f7.1)") diagvals(zsb)
    write (6,"('tssb# nestin ',9f7.1)") diagvals(tssb) 
  end if
#endif
   
  ! determine time corrosponding to new host nudging data
  kdhour = ktime_r/100 - ktime/100
  kdmin = (ktime_r-100*(ktime_r/100)) - (ktime-100*(ktime/100))
  mtimeb = 60*24*(iabsdate(kdate_r,kdate)-iabsdate(kdate,kdate)) + 60*kdhour+kdmin

#ifdef debug
  if ( myid==0 ) then
    write(6,*) 'nesting file has: kdate_r,ktime_r,kdhour,kdmin ',kdate_r,ktime_r,kdhour,kdmin
    write(6,*) 'kdate_r,iabsdate ',kdate_r,iabsdate(kdate_r,kdate)
    write(6,*) 'giving mtimeb = ',mtimeb
!   print additional information
    write(6,*) ' kdate ',kdate,' ktime ',ktime
    write(6,*) 'timeg,mtimer,mtimea,mtimeb: ',timeg,mtimer,mtimea,mtimeb
    write(6,*) 'ds ',ds
  end if
#endif

! ensure qb big enough, but not too big in top levels (from Sept '04)
  qb(1:ifull,:) = max(qb(1:ifull,:), 0.)

#ifdef debug
! following is useful if troublesome data is read in
  if ( mod(ktau,nmaxpr)==0 .or. ktau==2 .or. diag ) then
    if ( myid==0 ) then
      write(6,*) 'following max/min values printed from nestin'
    end if
    call maxmin(ub,'ub',ktau,1.,kl)
    call maxmin(vb,'vb',ktau,1.,kl)
    call maxmin(tb,'tb',ktau,1.,kl)
    call maxmin(qb,'qb',ktau,1.e3,kl)
    if ( myid==0 ) then
      write(6,*) 'following are really psl not ps'
    end if
    call maxmin(pslb,'ps',ktau,100.,1)
  endif

! in these cases redefine pslb, tb and (effectively) zsb using zs
! this keeps fine-mesh land mask & zs
! presently simplest to do whole pslb, tb (& qb) arrays
  if ( nmaxpr==1 .and. mydiag ) then
    write(6,*) 'zs (idjd) :',zs(idjd)
    write(6,*) 'zsb (idjd) :',zsb(idjd)
    write (6,"('100*psl.wesn ',2p5f8.3)") psl(idjd),psl(iw(idjd)),psl(ie(idjd)),psl(is(idjd)),psl(in(idjd))
    write (6,"('ps.wesn ',-2p5f9.3)") ps(idjd),ps(iw(idjd)),ps(ie(idjd)),ps(is(idjd)),ps(in(idjd))
    write(6,*) 'pslb in(idjd) :',pslb(idjd)
    write(6,*) 'now call retopo from nestin'
  end if
#endif
  call retopo(pslb,zsb,zs(1:ifull),tb,qb)
#ifdef debug
  if ( nmaxpr==1. and. mydiag ) then
    write (6,"('100*pslb.wesn ',2p5f8.3)") pslb(idjd),pslb(iw(idjd)),pslb(ie(idjd)),pslb(is(idjd)),pslb(in(idjd))
    write(6,*) 'pslb out(idjd) :',pslb(idjd)
    write(6,*) 'after pslb print; num= ',num
  end if

  ! display diagnostics      
  if ( num==0 ) then
    num = 1
    call printa('zs  ',zs        ,ktau,0  ,ia,ib,ja,jb,0.,.01)
    call printa('zsb ',zsb       ,ktau,0  ,ia,ib,ja,jb,0.,.01)
    call printa('psl ',psl       ,ktau,0  ,ia,ib,ja,jb,0.,1.e2)
    call printa('pslb',pslb      ,ktau,0  ,ia,ib,ja,jb,0.,1.e2)
    call printa('t   ',t,ktau,nlv,ia,ib,ja,jb,200.,1.)
    call printa('tb  ',tb,ktau,nlv,ia,ib,ja,jb,200.,1.)
    call printa('u   ',u,ktau,nlv,ia,ib,ja,jb,0.,1.)
    call printa('ub  ',ub,ktau,nlv,ia,ib,ja,jb,0.,1.)
    call printa('v   ',v,ktau,nlv,ia,ib,ja,jb,0.,1.)
    call printa('vb  ',vb,ktau,nlv,ia,ib,ja,jb,0.,1.)
    return
  end if   !  num==0
#endif
      
end if ! (mtimer>mtimeb)

! now use tt, uu, vv arrays for time interpolated values
timerm = ktau*dt/60.   ! real value in minutes (in case dt < 60 seconds)
cona = (mtimeb-timerm)/real(mtimeb-mtimea)
conb = (timerm-mtimea)/real(mtimeb-mtimea)
psls(:) = cona*psla(:) + conb*pslb(:)
tt (:,:) = cona*ta(:,:) + conb*tb(:,:)
qgg(:,:) = cona*qa(:,:) + conb*qb(:,:)
uu (:,:) = cona*ua(:,:) + conb*ub(:,:)
vv (:,:) = cona*va(:,:) + conb*vb(:,:)

! calculate time interpolated tss 
if ( namip==0 ) then     ! namip SSTs/sea-ice take precedence
  if ( nmlo==0 ) then
    ! SSTs read from host model
    where ( .not.land )
      tss = cona*tssa + conb*tssb
      tgg(:,1) = tss
    end where  ! (.not.land)
  else
    if ( nud_sst/=0 .or. nud_sss/=0 .or. nud_ouv/=0 .or. nud_sfh/=0 ) then
      ! nudge mlo
      dumaa = cona*sssa + conb*sssb
      if ( wl<1 ) then
        ! determine if multiple levels of ocean data exist in host
        depthcheck(1) = maxval(ocndep(:,1)) ! check if 3D data exists
        call ccmpi_allreduce(depthcheck(1:1),depthcheck(2:2),"max",comm_world)
        if ( depthcheck(2)<0.5 ) then
          wl = 1
        else
          wl = wlev
        end if
      end if
      if ( wl==1 ) then ! switch to 2D if 3D data is missing
        call mloexpmelt(timelt)
        timelt = min( timelt, 271.2 )
        dumaa(:,1,1) = (cona*tssa+conb*tssb)*(1.-fraciceb) + timelt*fraciceb
        dumaa(:,1,1) = dumaa(:,1,1) - wrtemp
      end if
      call mlonudge(dumaa(:,:,1),dumaa(:,:,2),dumaa(:,:,3:4),ocndep(:,2),wl)
    end if  ! nud_sst/=0.or.nud_sss/=0.or.nud_ouv/=0.or.nud_sfh/=0
  endif     ! nmlo==0 ..else..
endif       ! namip==0
     
if ( abs(iaero)>=2 .and. nud_aero/=0 ) then
  xtgdav(:,:,:) = cona*xtghosta(:,:,:) + conb*xtghostb(:,:,:)
end if
     
return
end subroutine nestin


!--------------------------------------------------------------
! SCALE SELECTIVE FILTER ASSIMILATION
! Called for mbd/=0
subroutine nestinb

use aerosolldr                   ! LDR prognostic aerosols
use arrays_m                     ! Atmosphere dyamics prognostic arrays
use cc_mpi                       ! CC MPI routines
use daviesnudge                  ! Far-field nudging
use diag_m                       ! Diagnostic routines
use indices_m                    ! Grid index arrays
use latlong_m                    ! Lat/lon coordinates
use mlo                          ! Ocean physics and prognostic arrays
use onthefly_m                   ! Input interpolation routines
use pbl_m                        ! Boundary layer arrays
use savuvt_m                     ! Saved dynamic arrays
use soil_m                       ! Soil and surface data
use soilsnow_m                   ! Soil, snow and surface data
 
implicit none
 
include 'newmpar.h'              ! Grid parameters
include 'dates.h'                ! Date data
include 'parm.h'                 ! Model configuration
include 'stime.h'                ! File date data
 
integer, dimension(ifull) :: dumm
integer, save :: mtimeb = -1
integer, save :: wl = -1
integer kdate_r, ktime_r
integer kdhour, kdmin
real, dimension(2) :: depthcheck
real, dimension(:,:), allocatable, save :: tb, ub, vb, qb, ocndep
real, dimension(:), allocatable, save :: pslb, tssb, fraciceb
real, dimension(:), allocatable, save :: sicedepb
real, dimension(:,:,:), allocatable, save :: sssb
real, dimension(:,:,:), allocatable, save :: xtghostb
real, dimension(ifull) :: zsb, timelt
real, dimension(ifull,3) :: duma
real, dimension(ifull,ms,3) :: dumg
real, dimension(ifull,kl,5) :: dumv
real, dimension(ifull,3,3) :: dums
 
! mtimer, mtimeb are in minutes
if ( ktau<100 .and. myid==0 ) then
  write(6,*) 'in nestinb ktau,mtimer,mtimeb ',ktau,mtimer,mtimeb
  write(6,*) 'with kdate_s,ktime_s >= ',kdate_s,ktime_s
end if

! Load new host data to be ready for next call to filter
if ( mtimer>mtimeb ) then

  ! allocate arrays on first call     
  if ( .not.allocated(tb) ) then
    allocate( tb(ifull,kl), ub(ifull,kl), vb(ifull,kl), qb(ifull,kl) )
    allocate( pslb(ifull), tssb(ifull), fraciceb(ifull) )
    allocate( sicedepb(ifull), ocndep(ifull,2) )
    allocate( sssb(ifull,wlev,4) )
    allocate( xtghostb(ifull,kl,naero) )
    if ( nud_uv/=9 ) then
      ! initialise arrays for 1D filter
      call specinit
    end if
    ! define vertical weights
    call setdavvertwgt
  end if
          
! following (till end of subr) reads in next bunch of data in readiness
! read tb etc  - for globpea, straight into tb etc
  if ( abs(io_in)==1 ) then
    call onthefly(1,kdate_r,ktime_r,                            &
                  pslb,zsb,tssb,sicedepb,fraciceb,tb,ub,vb,qb,  &
                  dumg(:,:,1),dumg(:,:,2),dumg(:,:,3),          &
                  duma(:,1),dumv(:,:,1),dumv(:,:,2),            &
                  dumv(:,:,3),dumv(:,:,4),dumv(:,:,5),          &
                  dums(:,:,1),dums(:,:,2),dums(:,:,3),          &
                  duma(:,2),duma(:,3),dumm,sssb,ocndep,xtghostb)
  else
    write(6,*) 'ERROR: Scale-selective filter requires abs(io_in)=1'
    call ccmpi_abort(-1)
  endif   ! (abs(io_in)==1)
  tssb(:) = abs(tssb(:))  ! moved here Mar '03

  ! calculate time for next filter call   
  kdhour = ktime_r/100 - ktime/100   ! integer hour diff from Oct '05
  kdmin = (ktime_r-100*(ktime_r/100)) - (ktime-100*(ktime/100))
  mtimeb = 60*24*(iabsdate(kdate_r,kdate)-iabsdate(kdate,kdate)) + 60*kdhour + kdmin
  ! adjust input data for change in orography
  call retopo(pslb,zsb,zs(1:ifull),tb,qb)

end if ! ((mtimer>mtimeb).or.firstcall)

! Apply filter to model data using previously loaded host data
if ( mtimer==mtimeb .and. mod(nint(ktau*dt),60)==0 ) then

  ! atmospheric nudging if required
  if ( nud_p/=0 .or. nud_t/=0 .or. nud_uv/=0 .or. nud_q/=0 .or. nud_aero/=0 ) then
    pslb(:) = pslb(:) - psl(1:ifull)
    ub(:,:) = ub(:,:) - u(1:ifull,:)
    vb(:,:) = vb(:,:) - v(1:ifull,:)
    tb(:,:) = tb(:,:) - t(1:ifull,:)
    qb(:,:) = qb(:,:) - qg(1:ifull,:)
    if ( abs(iaero)>=2 .and. nud_aero/=0 ) then
      xtghostb(:,:,:) = xtghostb(:,:,:) - xtg(1:ifull,:,:)
    end if
    call getspecdata(pslb,ub,vb,tb,qb,xtghostb)
  end if

  ! specify sea-ice if not AMIP or Mixed-Layer-Ocean
  if ( namip==0 ) then  ! namip SSTs/sea-ice take precedence
    if ( nmlo==0 ) then
      ! following sice updating code copied from nestin June '08      
      ! check whether present ice points should change to/from sice points
      where ( fraciceb(:)>0. .and. fracice(:)==0 .and. .not.land(1:ifull) )
        ! N.B. if already a sice point, keep present tice (in tggsn)
        tggsn(:,1) = min( 271.2, tssb(:) )
      end where
      sicedep(:) = sicedepb(:)
      fracice(:) = fraciceb(:)
      ! because of new zs etc, ensure that sice is only over sea
      where ( fracice(:)<0.02 )
        fracice(:) = 0.
      end where
      where ( land(1:ifull) )
        sicedep(:) = 0.
        fracice(:) = 0.
      elsewhere ( fracice(:)>0. .and. sicedep(:)==0. .and. rlatt(:)>0. )
        ! assign to 2m in NH and 1m in SH (according to spo)
        sicedep(:) = 2.
      elsewhere ( fracice(:)>0. .and. sicedep(:)==0. )
        ! assign to 2m in NH and 1m in SH (according to spo)
        sicedep(:) = 1.
      elsewhere ( fracice(:)==0. .and. sicedep(:)>0. )
        fracice(:) = 1.
      end where
      ! update tss 
      where ( .not.land(1:ifull) )
        tss(:) = tssb(:)
        tgg(:,1) = tss(:)
      end where  ! (.not.land(iq))
    else
      ! nudge Mixed-Layer-Ocean
      if ( nud_sst/=0 .or. nud_sss/=0 .or. nud_ouv/=0 .or. nud_sfh/=0 ) then
        ! check host for 2D or 3D data
        if ( wl<1 ) then
          depthcheck(1) = maxval(ocndep(:,1)) ! check for 3D ocean data in host model
          call ccmpi_allreduce(depthcheck(1:1),depthcheck(2:2),"max",comm_world)
          if ( depthcheck(2)<0.5 ) then
            wl = 1
          else
            wl = wlev
          end if
        end if
        if ( wl==1 ) then ! switch to 2D data if 3D is missing
          call mloexpmelt(timelt)
          timelt(:) = min( timelt(:), 271.2 )
          sssb(:,1,1) = tssb(:)*(1.-fraciceb(:)) + timelt(:)*fraciceb(:)
          sssb(:,1,1) = sssb(:,1,1) - wrtemp
        end if
        call mlofilterhub(sssb(:,:,1),sssb(:,:,2),sssb(:,:,3:4),ocndep(:,2),wl)
      end if
    end if ! (nmlo==0)
  end if   ! (namip==0)
end if     ! (mod(nint(ktau*dt),60)==0)

return
end subroutine nestinb

!--------------------------------------------------------------
! This subroutine gathers and distributes data for the
! scale-selective filter
subroutine getspecdata(pslb,ub,vb,tb,qb,xtgb)

use aerosolldr                   ! Aerosol interface
use arrays_m                     ! Atmosphere dyamics prognostic arrays
use cc_mpi                       ! CC MPI routines
use daviesnudge                  ! Far-field nudging
use liqwpar_m                    ! Cloud water mixing ratios
use nharrs_m                     ! Non-hydrostatic atmosphere arrays
use savuvt_m                     ! Saved dynamic arrays
use savuv1_m                     ! Saved dynamic arrays
use sigs_m                       ! Atmosphere sigma levels
use vecsuv_m                     ! Map to cartesian coordinates
use work3sav_m                   ! Water and tracer saved arrays
use xyzinfo_m, only : x,y,z      ! Grid coordinate arrays
      
implicit none

include 'newmpar.h'              ! Grid parameters
include 'const_phys.h'           ! Physical constants
include 'kuocom.h'               ! Convection parameters
include 'parm.h'                 ! Model configuration
include 'parmdyn.h'              ! Dynamics parameters
include 'parmgeom.h'             ! Coordinate data

integer iq, k, ntr, kb, kln, klx, klt
real, dimension(ifull), intent(inout) :: pslb
real, dimension(ifull) :: costh,sinth
real, dimension(ifull,kl), intent(inout) :: ub, vb, tb, qb
real, dimension(ifull,kl,naero), intent(inout) :: xtgb
real, dimension(ifull,kl) :: tv
real, dimension(ifull) :: dum
real den, polenx, poleny, polenz, zonx, zony, zonz
logical lblock

! nud_uv=0 (no preturbing of winds)
! nud_uv=1 (1D scale-selective filter)
! nud_uv=3 (JLM preturb zonal winds with 1D filter)
! nud_uv=9 (2D scale-selective filter)

! zonal wind option
if ( nud_uv==3 ) then
  polenx = -cos(rlat0*pi/180.)
  poleny = 0.
  polenz = sin(rlat0*pi/180.)
  do iq = 1,ifull
   zonx = real(            -polenz*y(iq))
   zony = real(polenz*x(iq)-polenx*z(iq))
   zonz = real(polenx*y(iq)             )
   den = sqrt( max( zonx**2 + zony**2 + zonz**2, 1.e-7 ) ) 
   costh(iq) =  (zonx*ax(iq)+zony*ay(iq)+zonz*az(iq))/den
   sinth(iq) = -(zonx*bx(iq)+zony*by(iq)+zonz*bz(iq))/den
  enddo
  do k = kbotdav,ktopdav
    dum(:) = costh(:)*ub(:,k) - sinth(:)*vb(:,k)
    ub(:,k) = dum(:)
  end do
end if
    
! Loop over maximum block size
! kblock can be reduced to save memory
do kb = kbotdav,ktopdav,kblock
  if ( myid==0 ) then     
    write(6,*) "Gather data for spectral filter      ",kb
  end if      
  kln = kb                          ! lower limit of block
  klx = min( kb+kblock-1, ktopdav ) ! upper limit of block
  klt = klx - kln + 1               ! number of levels in block
  lblock = (kb==kbotdav)            ! flag for first loop (include psl)
        
  !-----------------------------------------------------------------------
  ! select nudging option
  if ( nud_uv==9 ) then 
    if ( myid==0 ) then
      write(6,*) "Two dimensional spectral filter      ",kb
    end if
    call slowspecmpi(.1*real(mbd)/(pi*schmidt),pslb,ub,vb,tb,qb,xtgb,lblock,klt,kln,klx)
  else
    if ( myid==0 ) then
#ifdef uniform_decomp
      write(6,*) "Separable 1D filter                  ",kb
#else
      write(6,*) "Separable 1D filter (optimised)      ",kb
#endif
    end if
    call specfastmpi(.1*real(mbd)/(pi*schmidt),pslb,ub,vb,tb,qb,xtgb,lblock,klt,kln,klx)
  endif  ! (nud_uv==9) .. else ..
  !-----------------------------------------------------------------------

  if ( myid==0 ) then
    write(6,*) "Distribute data from spectral filter ",kb
  end if
        
end do

      
if ( nud_p>0 ) then
  psl(1:ifull) = psl(1:ifull) + pslb(:)
  ps(1:ifull) = 1.e5*exp(psl(1:ifull))
end if
if ( nud_uv/=0 ) then
  if ( nud_uv==3 ) then
    do k = kbotdav,ktopdav
      dum(:) = ub(:,k)
      ub(1:ifull,k) = costh(:)*dum(:)
      vb(1:ifull,k) = -sinth(:)*dum(:)
    end do
  end if
  do k = kbotdav,ktopdav
    ub(:,k) = ub(:,k)*vertwgt(k)
    vb(:,k) = vb(:,k)*vertwgt(k)
  end do
  u(1:ifull,kbotdav:ktopdav) = u(1:ifull,kbotdav:ktopdav) + ub(:,kbotdav:ktopdav)
  v(1:ifull,kbotdav:ktopdav) = v(1:ifull,kbotdav:ktopdav) + vb(:,kbotdav:ktopdav)
  savu(1:ifull,kbotdav:ktopdav) = savu(1:ifull,kbotdav:ktopdav) + ub(:,kbotdav:ktopdav)
  savu1(1:ifull,kbotdav:ktopdav) = savu1(1:ifull,kbotdav:ktopdav) + ub(:,kbotdav:ktopdav)
  savu2(1:ifull,kbotdav:ktopdav) = savu2(1:ifull,kbotdav:ktopdav) + ub(:,kbotdav:ktopdav)
  savv(1:ifull,kbotdav:ktopdav) = savv(1:ifull,kbotdav:ktopdav) + vb(:,kbotdav:ktopdav)
  savv1(1:ifull,kbotdav:ktopdav) = savv1(1:ifull,kbotdav:ktopdav) + vb(:,kbotdav:ktopdav)
  savv2(1:ifull,kbotdav:ktopdav) = savv2(1:ifull,kbotdav:ktopdav) + vb(:,kbotdav:ktopdav)
end if
if ( nud_t>0 ) then
  do k = kbotdav,ktopdav
    tb(:,k) = tb(:,k)*vertwgt(k)
  end do
  t(1:ifull,kbotdav:ktopdav) = t(1:ifull,kbotdav:ktopdav) + tb(:,kbotdav:ktopdav)
end if
if ( nud_q>0 ) then
  do k = kbotdav,ktopdav
    qb(:,k) = qb(:,k)*vertwgt(k)
  end do
  qg(1:ifull,kbotdav:ktopdav) = max(qg(1:ifull,kbotdav:ktopdav)+qb(:,kbotdav:ktopdav), 0.)
end if
if ( nud_t>0 .or. nud_q>0 ) then
  tv(:,:) = t(1:ifull,:)*(1.+0.61*qg(1:ifull,:)-qlg(1:ifull,:)-qfg(1:ifull,:) &
                         -qsng(1:ifull,:)-qgrg(1:ifull,:))
  phi(:,1) = bet(1)*tv(1:ifull,1)
  do k = 2,kl
    phi(:,k) = phi(:,k-1) + bet(k)*tv(1:ifull,k) + betm(k)*tv(1:ifull,k-1)
  end do
  phi(:,:) = phi(:,:) + phi_nh(:,:)
end if
if ( abs(iaero)>=2 .and. nud_aero>0 ) then
  do ntr = 1,size(xtg,3)
    do k = kbotdav,ktopdav
      xtgb(:,k,ntr) = xtgb(:,k,ntr)*vertwgt(k)
    end do
  end do
  xtg(1:ifull,kbotdav:ktopdav,:) = max(xtg(1:ifull,kbotdav:ktopdav,:)+xtgb(:,kbotdav:ktopdav,:), 0.)
end if

return
end subroutine getspecdata

!---------------------------------------------------------------------------------
! Slow 2D spectral downscaling - MPI version
! This option is an exact treatment of the filter
subroutine slowspecmpi(cin,pslb,ub,vb,tb,qb,xtgb,lblock,klt,kln,klx)

use aerosolldr        ! Aerosol interface
use cc_mpi            ! CC MPI routines
use vecsuv_m          ! Map to cartesian coordinates
      
implicit none
      
include 'newmpar.h'   ! Grid parameters
include 'parm.h'      ! Model configuration

integer, intent(in) :: klt, kln, klx
integer k, n
real, intent(in) :: cin
real, dimension(ifull), intent(inout) :: pslb
real, dimension(ifull,kl), intent(inout) :: ub, vb
real, dimension(ifull,kl), intent(inout) :: tb, qb
real, dimension(ifull,kl,naero), intent(inout) :: xtgb
real, dimension(ifull,kln:klx) :: wb
real, dimension(ifull_g,klt) :: tt ! large common array
real, dimension(ifull) :: da, db
real cq
logical, intent(in) :: lblock

cq = sqrt(4.5)*cin

if ( nud_p>0 .and. lblock ) then
  ! Create global copy of host data on each processor.  This can require a lot of memory
  call ccmpi_gatherall(pslb(:), tt(:,1))
  ! Apply 2D filter
  call slowspecmpi_work(cq,tt(:,1),pslb,1)
end if
if ( nud_uv==3 ) then
  call ccmpi_gatherall(ub(:,kln:klx),tt)
  call slowspecmpi_work(cq,tt,ub(:,kln:klx),klt)
else if ( nud_uv>0 ) then
  ! vectors are processed as Cartesian coordinates (X,Y,Z),
  ! avoiding complications along panel boundaries
  do k = kln,klx
    da(:) = ub(:,k)
    db(:) = vb(:,k)
    ub(:,k) = ax(1:ifull)*da(:) + bx(1:ifull)*db(:)
    vb(:,k) = ay(1:ifull)*da(:) + by(1:ifull)*db(:)
    wb(:,k) = az(1:ifull)*da(:) + bz(1:ifull)*db(:)
  end do
  call ccmpi_gatherall(ub(:,kln:klx),tt)
  call slowspecmpi_work(cq,tt,ub(:,kln:klx),klt)
  call ccmpi_gatherall(vb(:,kln:klx),tt)
  call slowspecmpi_work(cq,tt,vb(:,kln:klx),klt)
  call ccmpi_gatherall(wb(:,kln:klx),tt)
  call slowspecmpi_work(cq,tt,wb(:,kln:klx),klt)
  ! Convert Cartesian vectors back to Conformal Cubic vectors
  do k = kln,klx
    da(:) = ax(1:ifull)*ub(:,k) + ay(1:ifull)*vb(:,k) + az(1:ifull)*wb(:,k)
    db(:) = bx(1:ifull)*ub(:,k) + by(1:ifull)*vb(:,k) + bz(1:ifull)*wb(:,k)
    ub(:,k) = da(:)
    vb(:,k) = db(:)
  end do
end if
if ( nud_t>0 ) then
  call ccmpi_gatherall(tb(:,kln:klx),tt)
  call slowspecmpi_work(cq,tt,tb(:,kln:klx),klt)
end if
if ( nud_q>0 ) then
  call ccmpi_gatherall(qb(:,kln:klx),tt)
  call slowspecmpi_work(cq,tt,qb(:,kln:klx),klt)
end if
if ( abs(iaero)>=2 .and. nud_aero>0 ) then
  do n = 1,naero
    call ccmpi_gatherall(xtgb(:,kln:klx,n),tt)
    call slowspecmpi_work(cq,tt,xtgb(:,kln:klx,n),klt)
  end do
end if      

return
end subroutine slowspecmpi

subroutine slowspecmpi_work(cq,tt,tb,klt)

use cc_mpi            ! CC MPI routines
use map_m             ! Grid map arrays
use xyzinfo_m         ! Grid coordinate arrays
      
implicit none
      
include 'newmpar.h'   ! Grid parameters

integer, intent(in) :: klt
integer i, j, n, iq, iqg, k
real, intent(in) :: cq
real, dimension(ifull,klt), intent(out) :: tb
real, dimension(ifull_g,klt), intent(in) :: tt
real, dimension(ifull_g) :: r ! large working array
real psum

#ifdef debug
if ( myid==0 .and. nmaxpr==1 ) then
  write(6,*) "Start 2D filter"
end if
#endif

! evaluate the 2D convolution
call START_LOG(nestcalc_begin)
do n = 1,npan
  do j = 1,jpan
    do i = 1,ipan
      iqg = i + ioff + (j+joff-1)*il_g + (n-noff)*il_g*il_g
      iq  = i + (j-1)*ipan + (n-1)*ipan*jpan
      ! calculate distance between targer grid point and all other grid points
      r(:) = real(x_g(iqg)*x_g(:)+y_g(iqg)*y_g(:)+z_g(iqg)*z_g(:))
      r(:) = acos(max( min( r(:), 1. ), -1. ))
      ! evaluate Gaussian weights as a function of distance
      r(:) = exp(-(cq*r(:))**2)/(em_g(:)*em_g(:))
      ! discrete normalisation factor
      psum = sum(r(:))
      ! apply low band pass filter
      do k = 1,klt
        tb(iq,k) = sum(r(:)*tt(:,k))/psum
      end do
    end do
  end do
end do
call END_LOG(nestcalc_end)
 
#ifdef debug
if ( myid==0 .and. nmaxpr==1 ) write(6,*) "End 2D filter"
#endif

return
end subroutine slowspecmpi_work
!---------------------------------------------------------------------------------

!---------------------------------------------------------------------------------
! Four pass spectral downscaling
subroutine specfastmpi(cin,psls,uu,vv,tt,qgg,xtgg,lblock,klt,kln,klx)
      
use aerosolldr         ! Aerosol interface
use cc_mpi             ! CC MPI routines
      
implicit none
      
include 'newmpar.h'    ! Grid parameters
include 'parm.h'       ! Model configuration
      
integer, intent(in) :: klt, kln, klx
real, intent(in) :: cin
real, dimension(ifull), intent(inout) :: psls
real, dimension(ifull,kl), intent(inout) :: uu, vv
real, dimension(ifull,kl), intent(inout) :: tt, qgg
real, dimension(ifull,kl,naero), intent(inout) :: xtgg
logical, intent(in) :: lblock
      
if ( npta==1 ) then
  ! face version (nproc>=6)
  call spechost_n(cin,psls,uu,vv,tt,qgg,xtgg,lblock,klt,kln,klx)
else
  ! normal version
  call spechost(cin,psls,uu,vv,tt,qgg,xtgg,lblock,klt,kln,klx)
end if

return
end subroutine specfastmpi
!---------------------------------------------------------------------------------


!---------------------------------------------------------------------------------
! This is the main routine for the scale-selective filter
! (see spechost_n for a reduced memory version)
subroutine spechost(cin,pslb,ub,vb,tb,qb,xtgb,lblock,klt,kln,klx)

use aerosolldr        ! Aerosol interface
use cc_mpi            ! CC MPI routines
use vecsuv_m          ! Map to cartesian coordinates
      
implicit none
      
include 'newmpar.h'   ! Grid parameters
include 'parm.h'      ! Model configuration
      
integer, intent(in) :: klt, kln, klx
integer i, k, n, ppass
real, intent(in) :: cin
real, dimension(ifull), intent(inout) :: pslb
real, dimension(ifull,kl), intent(inout) :: ub, vb
real, dimension(ifull,kl), intent(inout) :: tb, qb
real, dimension(ifull,kl,naero), intent(inout) :: xtgb
real, dimension(ifull,kln:klx) :: wb
real, dimension(ipan*jpan,klt) :: qt
real, dimension(ifull) :: da, db
logical, intent(in) :: lblock

if ( nud_p>0 .and. lblock ) then
  call ccmpi_gathermap(pslb,1)            ! gather data onto global sparse array (1)
  do ppass = pprocn,pprocx
    call copyglobalpack(1,0,1)            ! copy sparse array data (1) to (0)
    call fastspecmpi_work(cin,qt,1,ppass) ! filter sparse array (0)
    do n = 1,ipan*jpan
      pslb(n+ipan*jpan*(ppass-pprocn)) = qt(n,1)
    end do
  end do
end if
if ( nud_t>0 ) then
  call ccmpi_gathermap(tb(:,kln:klx),klt)   ! gather data onto global sparse array (1)
  do ppass = pprocn,pprocx
    call copyglobalpack(klt,0,klt)          ! copy sparse array data (1) to (0)
    call fastspecmpi_work(cin,qt,klt,ppass) ! filter sparse array (0)
    do k = 1,klt
      do n = 1,ipan*jpan
        tb(n+ipan*jpan*(ppass-pprocn),k+kln-1) = qt(n,k)
      end do
    end do
  end do
end if
if ( nud_q>0 ) then
  call ccmpi_gathermap(qb(:,kln:klx),klt)   ! gather data onto global sparse array (1)
  do ppass = pprocn,pprocx
    call copyglobalpack(klt,0,klt)          ! copy sparse array data (1) to (0)
    call fastspecmpi_work(cin,qt,klt,ppass) ! filter sparse array (0)
    do k = 1,klt
      do n = 1,ipan*jpan
        qb(n+ipan*jpan*(ppass-pprocn),k+kln-1) = qt(n,k)
      end do
    end do
  end do
end if
if ( nud_uv==3 ) then
  call ccmpi_gathermap(ub(:,kln:klx),klt)        ! gather data onto global sparse array (1)
  do ppass = pprocn,pprocx
    call copyglobalpack(klt,0,klt)               ! copy sparse array data (1) to (0)
    call fastspecmpi_work(cin,qt(:,1),klt,ppass) ! filter sparse array (0)
    do k = 1,klt
      do n = 1,ipan*jpan
        ub(n+ipan*jpan*(ppass-pprocn),k+kln-1) = qt(n,k)
      end do
    end do
  end do
else if ( nud_uv>0 ) then
  do k = kln,klx
    da(:) = ub(:,k)
    db(:) = vb(:,k)
    ub(:,k) = ax(1:ifull)*da(:) + bx(1:ifull)*db(:)
    vb(:,k) = ay(1:ifull)*da(:) + by(1:ifull)*db(:)
    wb(:,k) = az(1:ifull)*da(:) + bz(1:ifull)*db(:)
  end do
  call ccmpi_gathermap(ub(:,kln:klx),klt)   ! gather data onto global sparse array (1)
  do ppass = pprocn,pprocx
    call copyglobalpack(klt,0,klt)          ! copy sparse array data (1) to (0)
    call fastspecmpi_work(cin,qt,klt,ppass) ! filter sparse array (0)
    do k = 1,klt
      do n = 1,ipan*jpan
        ub(n+ipan*jpan*(ppass-pprocn),k+kln-1) = qt(n,k)
      end do
    end do
  end do
  call ccmpi_gathermap(vb(:,kln:klx),klt)   ! gather data onto global sparse array (1)
  do ppass = pprocn,pprocx
    call copyglobalpack(klt,0,klt)          ! copy sparse array data (1) to (0)
    call fastspecmpi_work(cin,qt,klt,ppass) ! filter sparse array (0)
    do k = 1,klt
      do n = 1,ipan*jpan
        vb(n+ipan*jpan*(ppass-pprocn),k+kln-1) = qt(n,k)
      end do
    end do
  end do
  call ccmpi_gathermap(wb(:,kln:klx),klt)   ! gather data onto global sparse array (1)
  do ppass = pprocn,pprocx
    call copyglobalpack(klt,0,klt)          ! copy sparse array data (1) to (0)
    call fastspecmpi_work(cin,qt,klt,ppass) ! filter sparse array (0)
    do k = 1,klt
      do n = 1,ipan*jpan
        wb(n+ipan*jpan*(ppass-pprocn),k+kln-1) = qt(n,k)
      end do
    end do
  end do
  do k = kln,klx
    da(:) = ax(1:ifull)*ub(:,k) + ay(1:ifull)*vb(:,k) + az(1:ifull)*wb(:,k)
    db(:) = bx(1:ifull)*ub(:,k) + by(1:ifull)*vb(:,k) + bz(1:ifull)*wb(:,k)
    ub(:,k) = da(:)
    vb(:,k) = db(:)
  end do
endif
if ( abs(iaero)>=2 .and. nud_aero>0 ) then
  do i = 1,naero
    call ccmpi_gathermap(xtgb(:,kln:klx,i),klt) ! gather data onto global sparse array (1)
    do ppass = pprocn,pprocx
      call copyglobalpack(klt,0,klt)            ! copy sparse array data (1) to (0)
      call fastspecmpi_work(cin,qt,klt,ppass)   ! filter sparse array (0)
      do k = 1,klt
        do n = 1,ipan*jpan
          xtgb(n+ipan*jpan*(ppass-pprocn),k+kln-1,i) = qt(n,k)
        end do
      end do
    end do
  end do
end if      

return
end subroutine spechost
!---------------------------------------------------------------------------------

! This version of spechost is for one panel per processor (reduced memory)
subroutine spechost_n(cin,pslb,ub,vb,tb,qb,xtgb,lblock,klt,kln,klx)

use aerosolldr        ! Aerosol interface
use cc_mpi            ! CC MPI routines
use vecsuv_m          ! Map to cartesian coordinates
      
implicit none
      
include 'newmpar.h'   ! Grid parameters
include 'parm.h'      ! Model configuration
      
integer, intent(in) :: klt,kln,klx
integer k,n
real, intent(in) :: cin
real, dimension(ifull), intent(inout) :: pslb
real, dimension(ifull,kl), intent(inout) :: ub, vb
real, dimension(ifull,kl), intent(inout) :: tb, qb
real, dimension(ifull,kl,naero), intent(inout) :: xtgb
real, dimension(ifull,kln:klx) :: wb
real, dimension(ipan*jpan,klt) :: qt
real, dimension(ifull) :: da, db
logical, intent(in) :: lblock
      
if ( nud_p>0 .and. lblock ) then
  call ccmpi_gathermap(pslb,0)                ! gather data onto global sparse array (0)
  call fastspecmpi_work(cin,qt(:,1),1,pprocn) ! filter sparse array (0)
  pslb(:) = qt(:,1)
end if
if ( nud_uv==3 ) then
  call ccmpi_gathermap(ub(:,kln:klx),0)    ! gather data onto global sparse array (0)
  call fastspecmpi_work(cin,qt,klt,pprocn) ! filter sparse array (0)
  ub(:,kln:klx) = qt(:,:)
else if ( nud_uv>0 ) then
  do k = kln,klx
    da(:) = ub(:,k)
    db(:) = vb(:,k)
    ub(:,k) = ax(1:ifull)*da(:) + bx(1:ifull)*db(:)
    vb(:,k) = ay(1:ifull)*da(:) + by(1:ifull)*db(:)
    wb(:,k) = az(1:ifull)*da(:) + bz(1:ifull)*db(:)
  end do
  call ccmpi_gathermap(ub(:,kln:klx),0)    ! gather data onto global sparse array (0)
  call fastspecmpi_work(cin,qt,klt,pprocn) ! filter sparse array (0)
  ub(:,kln:klx) = qt(:,:)
  call ccmpi_gathermap(vb(:,kln:klx),0)    ! gather data onto global sparse array (0)
  call fastspecmpi_work(cin,qt,klt,pprocn) ! filter sparse array (0)
  vb(:,kln:klx) = qt(:,:)
  call ccmpi_gathermap(wb(:,kln:klx),0)    ! gather data onto global sparse array (0)
  call fastspecmpi_work(cin,qt,klt,pprocn) ! filter sparse array (0)
  wb(:,kln:klx) = qt(:,:)
  do k = kln,klx
    da(:) = ax(1:ifull)*ub(:,k) + ay(1:ifull)*vb(:,k) + az(1:ifull)*wb(:,k)
    db(:) = bx(1:ifull)*ub(:,k) + by(1:ifull)*vb(:,k) + bz(1:ifull)*wb(:,k)
    ub(:,k) = da(:)
    vb(:,k) = db(:)
  end do
endif
if ( nud_t>0 ) then
  call ccmpi_gathermap(tb(:,kln:klx),0)    ! gather data onto global sparse array (0)
  call fastspecmpi_work(cin,qt,klt,pprocn) ! filter sparse array (0)
  tb(:,kln:klx) = qt(:,:)
end if
if ( nud_q>0 ) then
  call ccmpi_gathermap(qb(:,kln:klx),0)    ! gather data onto global sparse array (0)
  call fastspecmpi_work(cin,qt,klt,pprocn) ! filter sparse array (0)
  qb(:,kln:klx) = qt(:,:)
end if
if ( abs(iaero)>=2 .and. nud_aero>0 ) then
  do n = 1,naero
    call ccmpi_gathermap(xtgb(:,kln:klx,n),0) ! gather data onto global sparse array (0)
    call fastspecmpi_work(cin,qt,klt,pprocn)  ! filter sparse array (0)
    xtgb(:,kln:klx,n) = qt(:,:)
  end do
end if      

return
end subroutine spechost_n

! determine if panel is 'left' or 'right'.  This
! basically reorders the convolution so that the
! final convolution leaves the local processors
! data in the correct orientation.

subroutine fastspecmpi_work(cin,tt,klt,ppass)

use cc_mpi            ! CC MPI routines
      
implicit none
      
include 'newmpar.h'   ! Grid parameters

integer, intent(in) :: klt,ppass
integer xpan
real, intent(in) :: cin
real, dimension(ipan*jpan,klt), intent(inout) :: tt
real cq

cq = sqrt(4.5)*cin ! filter length scale
xpan = max(ipan, jpan)

! computations for the local processor group
select case(ppass)
  case(1, 2, 3)
    call speclocal_left(cq,ppass,tt,klt,xpan)
  case(0, 4, 5)
    call speclocal_right(cq,ppass,tt,klt,xpan)
end select
     
return
end subroutine fastspecmpi_work
      
!---------------------------------------------------------------------------------
! This code runs between the local processors
      
subroutine speclocal_left(cq,ppass,qt,klt,xpan)

use cc_mpi            ! CC MPI routines
use map_m             ! Grid map arrays
use xyzinfo_m         ! Grid coordinate arrays

implicit none
      
include 'newmpar.h'   ! Grid parameters
include 'parm.h'      ! Model configuration
      
integer, intent(in) :: ppass, klt, xpan
integer j, k, n, ipass
integer jpoff, ibase
integer me, ns, ne, os, oe
integer til, a, b, c, sn, sy, jj, nn
integer ibeg, iend
integer, dimension(0:3) :: astr, bstr, cstr
integer, dimension(0:3) :: maps
real, intent(in) :: cq
real, dimension(ipan*jpan,klt), intent(out) :: qt
real, dimension(4*il_g,klt) :: at             ! subset of sparse array
real, dimension(4*il_g) :: asum               ! subset of sparse array
real, dimension(4*il_g) :: ra                 ! large working array
real, dimension(xpan,xpan,klt) :: pt
real, dimension(xpan,xpan) :: psum
real, dimension(il_g*xpan*(klt+1)) :: dd      ! subset of sparse array
real, dimension(xpan*xpan*(klt+1)) :: ff
real(kind=8), dimension(4*il_g) :: xa, ya, za ! subset of shared array
      
! matched for panels 1,2 and 3
      
maps = (/ il_g, il_g, 4*il_g, 3*il_g /)
til = il_g*il_g
astr = 0
bstr = 0
cstr = 0

ns = joff + 1
ne = joff + jpan
os = ioff + 1
oe = ioff + ipan
      
do ipass = 0,2
  me = maps(ipass)
  call getiqa(astr,bstr,cstr,me,ipass,ppass,il_g)

#ifdef debug        
  if ( myid==0 ) write(6,*) "Start convolution"
#endif
  
  do j = 1,jpan
  
    ! pack data from sparse arrays
    call START_LOG(nestpack_begin)
    jj = j + ns - 1
    do sn = 1,me,il_g
      sy = (sn-1)/il_g
      a = astr(sy)
      b = bstr(sy)
      c = cstr(sy)
      ibeg = a*sn + b*jj + c
      iend = a*(sn+il_g-1) + b*jj + c
      xa(sn:sn+il_g-1) = x_g(ibeg:iend:a)
      ya(sn:sn+il_g-1) = y_g(ibeg:iend:a)
      za(sn:sn+il_g-1) = z_g(ibeg:iend:a)
      asum(sn:sn+il_g-1) = 1./em_g(ibeg:iend:a)**2
      do k = 1,klt
        call getglobalpack(at(:,k),sn,ibeg,iend,k)
        at(sn:sn+il_g-1,k) = at(sn:sn+il_g-1,k)*asum(sn:sn+il_g-1)
      end do
    end do
    call END_LOG(nestpack_end)
    
    ! start convolution
    call START_LOG(nestcalc_begin)
    do n = 1,ipan
      nn = n + os - 1
      ra(1:me) = real(xa(nn)*xa(1:me)+ya(nn)*ya(1:me)+za(nn)*za(1:me))
      ra(1:me) = acos(max(min(ra(1:me), 1.), -1.))
      ra(1:me) = exp(-(cq*ra(1:me))**2)
      ! can also use the lines below which integrate the gaussian
      ! analytically over the length element (but slower)
      !ra(1) = 2.*erf(cq*0.5*(ds/rearth)
      !ra(2:me) = erf(cq*(ra(2:me)+0.5*(ds/rearth)))-erf(cq*(ra(2:me)-0.5*(ds/rearth)))
      psum(n,j) = sum(ra(1:me)*asum(1:me))
      do k = 1,klt
        pt(n,j,k) = sum(ra(1:me)*at(1:me,k))
      end do
    end do
    call END_LOG(nestcalc_end)
    
  end do

#ifdef debug
  if ( myid==0 ) then
    write(6,*) "End convolution"
    write(6,*) "Send arrays to local host"
  end if
#endif

  ! unpacking grid
  a = astr(0)
  b = bstr(0)
  c = cstr(0)

  ! gather data for final pass
  call START_LOG(nestcomm_begin)
  ff(1:ipan*jpan*klt) = reshape( pt(1:ipan,1:jpan,1:klt), (/ ipan*jpan*klt /) )
  ff(ipan*jpan*klt+1:ipan*jpan*klt+ipan*jpan) = reshape( psum(1:ipan,1:jpan), (/ ipan*jpan /) )  
  call ccmpi_allgatherx(dd(1:il_g*ipan*(klt+1)),ff(1:ipan*jpan*(klt+1)),comm_cols)
  call END_LOG(nestcomm_end)
  
  ! unpack to sparse arrays
  call START_LOG(nestunpack_begin)
  do jpoff = 0,il_g-1,jpan
    sy = jpoff/jpan
    do k = 1,klt
      do j = jpoff+1,jpoff+jpan
        ibase = 1 + ipan*(j-jpoff-1) + ipan*jpan*(k-1) + ipan*jpan*(klt+1)*sy
        ibeg = a*os + b*j + c
        iend = a*oe + b*j + c
        call setglobalpack(dd(:),ibase,ibeg,iend,k,trans=.false.)
      end do
    end do
    do j = jpoff+1,jpoff+jpan
      ibase = 1 + ipan*(j-jpoff-1) + ipan*jpan*klt + ipan*jpan*(klt+1)*sy
      ibeg = a*os + b*j + c
      iend = a*oe + b*j + c
      call setglobalpack(dd(:),ibase,ibeg,iend,0,trans=.false.)
    end do
  end do
  call END_LOG(nestunpack_end)

end do

ns = ioff + 1
ne = ioff + ipan
os = joff + 1
oe = joff + jpan

ipass = 3
me = maps(ipass)
call getiqa(astr,bstr,cstr,me,ipass,ppass,il_g)

#ifdef debug
if ( myid==0 ) write(6,*) "Start convolution"
#endif

do j = 1,ipan
    
  ! pack from sparse arrays
  call START_LOG(nestpack_begin)
  jj = j + ns - 1
  do sn = 1,me,il_g
    sy = (sn-1)/il_g
    a = astr(sy)
    b = bstr(sy)
    c = cstr(sy)
    ibeg = a*sn + b*jj + c
    iend = a*(sn+il_g-1) + b*jj + c
    xa(sn:sn+il_g-1) = x_g(ibeg:iend:a)
    ya(sn:sn+il_g-1) = y_g(ibeg:iend:a)
    za(sn:sn+il_g-1) = z_g(ibeg:iend:a)
    call getglobalpack(asum(:),sn,ibeg,iend,0)
    do k = 1,klt
      call getglobalpack(at(:,k),sn,ibeg,iend,k)
    end do
  end do
  call END_LOG(nestpack_end)
  
  ! start convolution
  call START_LOG(nestcalc_begin)
  do n = 1,jpan
    nn = n + os - 1
    ra(1:me) = real(xa(nn)*xa(1:me)+ya(nn)*ya(1:me)+za(nn)*za(1:me))
    ra(1:me) = acos(max(min(ra(1:me), 1.), -1.))
    ra(1:me) = exp(-(cq*ra(1:me))**2)
    ! can also use the lines below which integrate the gaussian
    ! analytically over the length element (but slower)
    !ra(1) = 2.*erf(cq*0.5*(ds/rearth)
    !ra(2:me) = erf(cq*(ra(2:me)+0.5*(ds/rearth)))-erf(cq*(ra(2:me)-0.5*(ds/rearth)))
    psum(n,j) = sum(ra(1:me)*asum(1:me))
    do k = 1,klt
      pt(n,j,k) = sum(ra(1:me)*at(1:me,k))
    end do
  end do
  call END_LOG(nestcalc_end)
  
end do

#ifdef debug
if ( myid==0 ) then
  write(6,*) "End convolution"
  write(6,*) "Send arrays to local host"
end if
#endif

! unpack data to local array
call START_LOG(nestunpack_begin)
do k = 1,klt
  do j = 1,ipan
    do n = 1,jpan
      qt(j+ipan*(n-1),k) = pt(n,j,k)/psum(n,j)
    end do
  end do
end do
call END_LOG(nestunpack_end)

return  
end subroutine speclocal_left

subroutine speclocal_right(cq,ppass,qt,klt,xpan)

use cc_mpi            ! CC MPI routines
use map_m             ! Grid map arrays
use xyzinfo_m         ! Grid coordinate arrays

implicit none
      
include 'newmpar.h'   ! Grid parameters
include 'parm.h'      ! Model configuration
      
integer, intent(in) :: ppass, klt, xpan
integer j, k, n, ipass
integer jpoff, ibase
integer me, ns, ne, os, oe
integer til, a, b, c, sn, sy, jj, nn
integer ibeg, iend
integer, dimension(0:3) :: astr, bstr, cstr
integer, dimension(0:3) :: maps
real, intent(in) :: cq
real, dimension(ipan*jpan,klt), intent(out) :: qt
real, dimension(4*il_g,klt) :: at
real, dimension(4*il_g) :: asum, ra
real, dimension(xpan,xpan,klt) :: pt
real, dimension(xpan,xpan) :: psum
real, dimension(il_g*xpan*(klt+1)) :: dd
real, dimension(xpan*xpan*(klt+1)) :: ff
real(kind=8), dimension(4*il_g) :: xa, ya, za
      
! matched for panels 0, 4 and 5
      
maps = (/ il_g, il_g, 4*il_g, 3*il_g /)
til = il_g*il_g
astr = 0
bstr = 0
cstr = 0

ns = ioff + 1
ne = ioff + ipan
os = joff + 1
oe = joff + jpan
      
do ipass = 0,2
  me = maps(ipass)
  call getiqa(astr,bstr,cstr,me,ipass,ppass,il_g)

#ifdef debug
  if ( myid==0 ) write(6,*) "Start convolution"
#endif

  do j = 1,ipan
      
    ! pack data from sparse arrays
    call START_LOG(nestpack_begin)
    jj = j + ns - 1
    do sn = 1,me,il_g
      sy = (sn-1)/il_g
      a = astr(sy)
      b = bstr(sy)
      c = cstr(sy)
      ibeg = a*sn + b*jj + c
      iend = a*(sn+il_g-1) + b*jj + c
      xa(sn:sn+il_g-1) = x_g(ibeg:iend:a)
      ya(sn:sn+il_g-1) = y_g(ibeg:iend:a)
      za(sn:sn+il_g-1) = z_g(ibeg:iend:a)
      asum(sn:sn+il_g-1) = 1./em_g(ibeg:iend:a)**2
      do k = 1,klt
        call getglobalpack(at(:,k),sn,ibeg,iend,k)
        at(sn:sn+il_g-1,k) = at(sn:sn+il_g-1,k)*asum(sn:sn+il_g-1)
      end do
    end do
    call END_LOG(nestpack_end)
    
    ! start convolution
    call START_LOG(nestcalc_begin)
    do n = 1,jpan
      nn = n + os - 1
      ra(1:me) = real(xa(nn)*xa(1:me)+ya(nn)*ya(1:me)+za(nn)*za(1:me))
      ra(1:me) = acos(max( min(ra(1:me), 1.), -1.))
      ra(1:me) = exp(-(cq*ra(1:me))**2)
      ! can also use the lines below which integrate the gaussian
      ! analytically over the length element (but slower)
      !ra(1) = 2.*erf(cq*0.5*(ds/rearth)
      !ra(2:me) = erf(cq*(ra(2:me)+0.5*(ds/rearth)))-erf(cq*(ra(2:me)-0.5*(ds/rearth)))
      psum(n,j) = sum(ra(1:me)*asum(1:me))
      do k = 1,klt
        pt(n,j,k) = sum(ra(1:me)*at(1:me,k))
      end do
    end do
    call END_LOG(nestcalc_end)
    
  end do

#ifdef debug
  if ( myid==0 ) then
    write(6,*) "End convolution"
    write(6,*) "Send arrays to local host"
  end if
#endif

  ! unpacking grid
  a = astr(0)
  b = bstr(0)
  c = cstr(0)

  call START_LOG(nestcomm_begin)
  ff(1:ipan*jpan*klt) = reshape( pt(1:jpan,1:ipan,1:klt), (/ ipan*jpan*klt /) )
  ff(ipan*jpan*klt+1:ipan*jpan*klt+ipan*jpan) = reshape( psum(1:jpan,1:ipan), (/ ipan*jpan /) )
  call ccmpi_allgatherx(dd(1:il_g*jpan*(klt+1)),ff(1:ipan*jpan*(klt+1)),comm_rows)
  call END_LOG(nestcomm_end)
  
  ! unpack data to sparse arrays
  call START_LOG(nestunpack_begin)
  do jpoff = 0,il_g-1,ipan
    sy = jpoff/ipan
    do k = 1,klt
      do j = jpoff+1,jpoff+ipan
        ibase = 1 + jpan*(j-jpoff-1) + ipan*jpan*(k-1) + ipan*jpan*(klt+1)*sy
        ibeg = a*os + b*j + c
        iend = a*oe + b*j + c
        call setglobalpack(dd(:),ibase,ibeg,iend,k,trans=.true.)
      end do
    end do
    do j = jpoff+1,jpoff+ipan
      ibase = 1 + jpan*(j-jpoff-1) + ipan*jpan*klt + ipan*jpan*(klt+1)*sy
      ibeg = a*os + b*j + c
      iend = a*oe + b*j + c
      call setglobalpack(dd(:),ibase,ibeg,iend,0,trans=.true.)
    end do
  end do
  call END_LOG(nestunpack_end)

end do

ns = joff + 1
ne = joff + jpan
os = ioff + 1
oe = ioff + ipan

ipass = 3
me = maps(ipass)
call getiqa(astr,bstr,cstr,me,ipass,ppass,il_g)

#ifdef debug
if ( myid == 0 ) write(6,*) "Start convolution"
#endif

do j = 1,jpan
    
  ! pack data from sparse arrays
  call START_LOG(nestpack_begin)  
  jj = j + ns - 1
  do sn = 1,me,il_g
    sy = (sn-1)/il_g
    a = astr(sy)
    b = bstr(sy)
    c = cstr(sy)
    ibeg = a*sn + b*jj + c
    iend = a*(sn+il_g-1) + b*jj + c
    xa(sn:sn+il_g-1) = x_g(ibeg:iend:a)
    ya(sn:sn+il_g-1) = y_g(ibeg:iend:a)
    za(sn:sn+il_g-1) = z_g(ibeg:iend:a)
    call getglobalpack(asum(:),sn,ibeg,iend,0)
    do k = 1,klt
      call getglobalpack(at(:,k),sn,ibeg,iend,k)
    end do
  end do
  call END_LOG(nestpack_end)
  
  ! start convolution
  call START_LOG(nestcalc_begin)
  do n = 1,ipan
    nn = n + os - 1
    ra(1:me) = real(xa(nn)*xa(1:me)+ya(nn)*ya(1:me)+za(nn)*za(1:me))
    ra(1:me) = acos(max(min( ra(1:me), 1.), -1.))
    ra(1:me) = exp(-(cq*ra(1:me))**2)
    ! can also use the lines below which integrate the gaussian
    ! analytically over the length element (but slower)
    !ra(1) = 2.*erf(cq*0.5*(ds/rearth)
    !ra(2:me) = erf(cq*(ra(2:me)+0.5*(ds/rearth)))-erf(cq*(ra(2:me)-0.5*(ds/rearth)))
    psum(n,j) = sum(ra(1:me)*asum(1:me))
    do k = 1,klt
      pt(n,j,k) = sum(ra(1:me)*at(1:me,k))
    end do
  end do
  call END_LOG(nestcalc_end)

end do

#ifdef debug
if ( myid==0 ) then
  write(6,*) "End convolution"
  write(6,*) "Send arrays to local host"
end if
#endif

! unpack data
call START_LOG(nestunpack_begin)
do k = 1,klt
  do j = 1,jpan
    do n = 1,ipan
      qt(n+ipan*(j-1),k) = pt(n,j,k)/psum(n,j)
    end do
  end do
end do
call END_LOG(nestunpack_end)

return  
end subroutine speclocal_right

!---------------------------------------------------------------------------------
! Map from 1D convolution to global index
subroutine getiqa(a,b,c,ne,ipass,ppass,il_g)
      
implicit none
      
integer, intent(in) :: ne,ipass,ppass,il_g
integer, dimension(0:3), intent(out) :: a,b,c
integer sn,sy
      
do sn = 1,ne,il_g
  sy = (sn-1)/il_g

  select case(ppass*100+ipass*10+sy)
    case(0)                                ! panel 5   - x pass
      a(sy) = il_g
      b(sy) = 1
      c(sy) = il_g*(5*il_g-1)
    case(10)                               ! panel 2   - x pass
      a(sy) = -1
      b(sy) = il_g
      c(sy) = 2*il_g*il_g + 1
    case(20,21,321)                        ! panel 0,1 - y pass
      a(sy) = il_g
      b(sy) = 1
      c(sy) = -il_g
    case(22)                               ! panel 3   - y pass
      a(sy) = 1
      b(sy) = -il_g
      c(sy) = il_g*(4*il_g-2)
    case(23,323)                           ! panel 4   - y pass
      a(sy) = 1
      b(sy) = -il_g
      c(sy) = il_g*(5*il_g-3)
    case(30,100)                           ! panel 0   - z pass
      a(sy) = 1
      b(sy) = il_g
      c(sy) = -il_g
    case(31,223,523,532)                   ! panel 2   - z pass ! panel 4   - x pass ! panel 3   - z pass
      a(sy) = il_g
      b(sy) = -1
      c(sy) = il_g*il_g + 1
    case(32)                               ! panel 5   - z pass
      a(sy) = 1
      b(sy) = il_g
      c(sy) = il_g*(5*il_g-3)
    case(110)                              ! panel 3   - z pass
      a(sy) = -il_g
      b(sy) = 1
      c(sy) = 4*il_g*il_g
    case(120)                              ! panel 1   - x pass
      a(sy) = 1
      b(sy) = il_g
      c(sy) = il_g*(il_g-1)
    case(121,421)                          ! panel 2   - x pass
      a(sy) = 1
      b(sy) = il_g
      c(sy) = 2*il_g*(il_g-1)
    case(122,123,423)                      ! panel 4,5 - x pass ! panel 2   - z pass
      a(sy) = il_g
      b(sy) = -1
      c(sy) = 2*il_g*il_g+1
    case(130)                              ! panel 1   - y pass
      a(sy) = il_g
      b(sy) = 1
      c(sy) = il_g*(il_g-1)
    case(131)                              ! panel 3   - y pass
      a(sy) = 1
      b(sy) = -il_g
      c(sy) = il_g*(4*il_g-1)
    case(132,322)                          ! panel 0,1 - y pass
      a(sy) = il_g
      b(sy) = 1
      c(sy) = -il_g*(2*il_g+1)
    case(200)                              ! panel 0   - y pass
      a(sy) = -il_g
      b(sy) = 1
      c(sy) = il_g*il_g
    case(210)                              ! panel 3   - y pass
      a(sy) = 1
      b(sy) = il_g
      c(sy) = il_g*(3*il_g-1)
    case(220)                              ! panel 2   - x pass
      a(sy) = 1
      b(sy) = il_g
      c(sy) = il_g*(2*il_g-1)
    case(221,521)                          ! panel 1   - x pass
      a(sy) = 1
      b(sy) = il_g
      c(sy) = il_g*(il_g-2)
    case(222,410)                          ! panel 3   - z pass ! panel 5   - x pass
      a(sy) = il_g
      b(sy) = -1
      c(sy) = 3*il_g*il_g+1
    case(230)                              ! panel 2   - z pass
      a(sy) = il_g
      b(sy) = 1
      c(sy) = il_g*(2*il_g-1)
    case(231)                              ! panel 0   - z pass
      a(sy) = 1
      b(sy) = -il_g
      c(sy) = il_g*(il_g-1)
    case(232)                              ! panel 3   - z pass
      a(sy) = il_g
      b(sy) = 1
      c(sy) = il_g*(il_g-1)
    case(300)                              ! panel 5   - x pass
      a(sy) = -il_g
      b(sy) = -1
      c(sy) = il_g*(6*il_g+1) + 1
    case(310)                              ! panel 2   - x pass
      a(sy) = 1
      b(sy) = -il_g
      c(sy) = 3*il_g*il_g
    case(320)                              ! panel 3   - y pass
      a(sy) = 1
      b(sy) = -il_g
      c(sy) = 4*il_g*il_g
    case(330)                              ! panel 3   - z pass
      a(sy) = il_g
      b(sy) = 1
      c(sy) = il_g*(3*il_g-1)
    case(331)                              ! panel 2   - z pass
      a(sy) = il_g
      b(sy) = 1
      c(sy) = il_g*(il_g-1)
    case(332)                              ! panel 5   - z pass
      a(sy) = 1
      b(sy) = -il_g
      c(sy) = il_g*(6*il_g-2)
    case(400)                              ! panel 0   - z pass
      a(sy) = -1
      b(sy) = -il_g
      c(sy) = il_g*(il_g+1)+1
    case(420)                              ! panel 4   - x pass
      a(sy) = il_g
      b(sy) = -1
      c(sy) = 4*il_g*il_g+1
    case(422)                              ! panel 1   - x pass
      a(sy) = 1
      b(sy) = il_g
      c(sy) = il_g*(il_g-3)
    case(430)                              ! panel 4   - y pass
      a(sy) = 1
      b(sy) = il_g
      c(sy) = il_g*(4*il_g-1)
    case(431)                              ! panel 3   - y pass
      a(sy) = 1
      b(sy) = il_g
      c(sy) = il_g*(3*il_g-2)
    case(432)                              ! panel 0   - y pass
      a(sy) = il_g
      b(sy) = -1
      c(sy) = -2*il_g*il_g+1
    case(500)                              ! panel 0   - y pass
      a(sy) = il_g
      b(sy) = -1
      c(sy) = 1
    case(510)                              ! panel 3   - y pass
      a(sy) = -1
      b(sy) = -il_g
      c(sy) = il_g*(4*il_g+1) + 1
    case(520)                              ! panel 5   - x pass
      a(sy) = il_g
      b(sy) = -1
      c(sy) = 5*il_g*il_g+1
    case(522)                              ! panel 2   - x pass
      a(sy) = 1
      b(sy) = il_g
      c(sy) = il_g*(2*il_g-3)
    case(530)                              ! panel 5   - z pass
      a(sy) = 1
      b(sy) = il_g
      c(sy) = il_g*(5*il_g-1)            
    case(531)                              ! panel 0   - z pass
      a(sy) = 1
      b(sy) = il_g
      c(sy) = -2*il_g
    case DEFAULT
      write(6,*) "Invalid index ",ppass,ipass,sn,ppass*100+ipass*10+sy
      stop
  end select

end do

return
end subroutine getiqa

!---------------------------------------------------------------------------------
! This subroutine gathers and distributes data for the
! MLO scale-selective filter
subroutine mlofilterhub(sstb,sssb,suvb,sfh,wl)

use cc_mpi                                          ! CC MPI routines
use mlo, only : mloimport,mloexport,mloexpdep, &    ! Ocean physics and prognostic arrays
                wlev,wrtemp
use mlodynamics                                     ! Ocean dynamics routines
use soil_m                                          ! Soil and surface data
use vecsuv_m                                        ! Map to cartesian coordinates

implicit none

include 'newmpar.h'                                 ! Grid parameters      
include 'parm.h'                                    ! Model configuration

integer, intent(in) :: wl
integer k,ka,kb,kc,kln,klx,klt,kbb
real, dimension(ifull), intent(in) :: sfh
real, dimension(ifull,wlev), intent(in) :: sstb,sssb
real, dimension(ifull,wlev,2), intent(in) :: suvb
real, dimension(ifull,1) :: diffh_l
real, dimension(ifull,kblock) :: diff_l,diffs_l
real, dimension(ifull,kblock) :: diffu_l,diffv_l
real, dimension(ifull) :: old
logical lblock
real, parameter :: miss = 999999.
 
ka=0
kb=0
kc=min(kbotmlo,ktopmlo+wl-1)

if (nud_sfh/=0) then
  old=sfh
  call mloexport(4,old,0,0)
  where (.not.land)
    diffh_l(:,1)=sfh-old
  elsewhere
    diffh_l(:,1)=miss
  end where
end if
      
do kbb=ktopmlo,kc,kblock
      
  if (myid==0) then
    write(6,*) "Gather data for MLO filter     ",kbb
  end if
            
  kln=kbb
  klx=min(kbb+kblock-1,kc)
  klt=klx-kln+1
  lblock=(kbb==ktopmlo)
      
  if (nud_sst/=0) then
    do k=kln,klx
      kb=k-kln+1
      old=sstb(:,k)
      call mloexport(0,old,k,0)
      where (.not.land)
        diff_l(:,kb)=sstb(:,k)-old
      elsewhere
        diff_l(:,kb)=miss
      end where
    end do
  end if

  if (nud_sss/=0) then
    do k=kln,klx
      kb=k-kln+1
      old=sssb(:,k)
      call mloexport(1,old,k,0)
      where (.not.land)
        diffs_l(:,kb)=sssb(:,k)-old
      elsewhere
        diffs_l(:,kb)=miss
      end where
    end do
  end if

  if (nud_ouv/=0) then
    do k=kln,klx
      kb=k-kln+1
      old=suvb(:,k,1)
      call mloexport(2,old,k,0)
      where (.not.land)
        diffu_l(:,kb)=suvb(:,k,1)-old
      elsewhere
        diffu_l(:,kb)=miss
      end where
      old=suvb(:,k,2)
      call mloexport(3,old,k,0)
      where (.not.land)
        diffv_l(:,kb)=suvb(:,k,2)-old
      elsewhere
        diffv_l(:,kb)=miss
      end where
    end do
  end if

  if ((nud_uv/=9.and.abs(nmlo)/=1).or.namip/=0) then
    call mlofilterfast(diff_l(:,1:klt),diffs_l(:,1:klt),diffu_l(:,1:klt),diffv_l(:,1:klt),diffh_l(:,1),miss,lblock,klt)
  else
    call mlofilter(diff_l(:,1:klt),diffs_l(:,1:klt),diffu_l(:,1:klt),diffv_l(:,1:klt),diffh_l(:,1),miss,lblock,klt)
  end if

  if (myid==0) then
    write(6,*) "Distribute data for MLO filter ",kbb
  end if

  if (nud_sst/=0) then
    do k=kln,klx
      ka=min(wl,k)
      kb=k-kln+1
      old=sstb(:,ka)
      call mloexport(0,old,k,0)
      old=old+diff_l(:,kb)*10./real(mloalpha)
      call mloimport(0,old,k,0)
    end do
    if (klx==kc) then
      do k=kc+1,kbotmlo
        old=sstb(:,ka)
        call mloexport(0,old,k,0)
        old=old+diff_l(:,kb)*10./real(mloalpha) ! kb saved from above loop
        call mloimport(0,old,k,0)
      end do
    end if
  end if

  if (nud_sss/=0) then
    do k=kln,klx
      ka=min(wl,k)
      kb=k-kln+1
      old=sssb(:,ka)
      call mloexport(1,old,k,0)
      old=old+diffs_l(:,kb)*10./real(mloalpha)
      old=max(old,0.)
      call mloimport(1,old,k,0)
    end do
    if (klx==kc) then
      do k=kc+1,kbotmlo
        old=sssb(:,ka)
        call mloexport(1,old,k,0)
        old=old+diffs_l(:,kb)*10./real(mloalpha) ! kb saved from above loop
        old=max(old,0.)
        call mloimport(1,old,k,0)
      end do
    end if
  end if

  if (nud_ouv/=0) then
    do k=kln,klx
      ka=min(wl,k)
      kb=k-kln+1
      old=suvb(:,ka,1)
      call mloexport(2,old,k,0)
      old=old+diffu_l(:,kb)*10./real(mloalpha)
      call mloimport(2,old,k,0)
      if (allocated(oldu1)) then
        oldu1(:,k)=oldu1(:,k)+diffu_l(:,kb)*10./real(mloalpha)
        oldu2(:,k)=oldu2(:,k)+diffu_l(:,kb)*10./real(mloalpha)
      end if
    end do
    if (klx==kc) then
      do k=kc+1,kbotmlo
        old=suvb(:,ka,1)
        call mloexport(2,old,k,0)
        old=old+diffu_l(:,kb)*10./real(mloalpha) ! kb saved from above loop
        call mloimport(2,old,k,0)
        if (allocated(oldu1)) then
          oldu1(:,k)=oldu1(:,k)+diffu_l(:,kb)*10./real(mloalpha)
          oldu2(:,k)=oldu2(:,k)+diffu_l(:,kb)*10./real(mloalpha)
        end if
      end do
    end if
    do k=kln,klx
      ka=min(wl,k)
      kb=k-kln+1
      old=suvb(:,ka,2)
      call mloexport(3,old,k,0)
      old=old+diffv_l(:,kb)*10./real(mloalpha)
      call mloimport(3,old,k,0)
      if (allocated(oldv1)) then
        oldv1(:,k)=oldv1(:,k)+diffv_l(:,kb)*10./real(mloalpha)
        oldv2(:,k)=oldv2(:,k)+diffv_l(:,kb)*10./real(mloalpha)
      end if
    end do
    if (klx==kc) then
      do k=kc+1,kbotmlo
        old=suvb(:,ka,2)
        call mloexport(3,old,k,0)
        old=old+diffv_l(:,kb)*10./real(mloalpha)
        call mloimport(3,old,k,0)
        if (allocated(oldv1)) then
          oldv1(:,k)=oldv1(:,k)+diffv_l(:,kb)*10./real(mloalpha)
          oldv2(:,k)=oldv2(:,k)+diffv_l(:,kb)*10./real(mloalpha)
        end if
      end do
    end if
  end if
      
end do
     
if (nud_sfh/=0) then
  old=sfh
  call mloexport(4,old,0,0)
  old=old+diffh_l(:,1)*10./real(mloalpha)
  call mloimport(4,old,0,0)
end if

return
end subroutine mlofilterhub
      
!---------------------------------------------------------------------------------
! 2D Filter for MLO 
subroutine mlofilter(diff_l,diffs_l,diffu_l,diffv_l,diffh_l,miss,lblock,kd)

use cc_mpi                  ! CC MPI routines
use vecsuv_m                ! Map to cartesian coordinates

implicit none

include 'newmpar.h'         ! Grid parameters
include 'parm.h'            ! Model configuration

integer, intent(in) :: kd
integer k
real, intent(in) :: miss
real, dimension(ifull,1), intent(inout) :: diffh_l
real, dimension(ifull,kd), intent(inout) :: diff_l, diffs_l
real, dimension(ifull,kd), intent(inout) :: diffu_l, diffv_l
real, dimension(ifull,kd) :: diffw_l
real, dimension(ifull_g,kd) :: diff_g ! large common array
real, dimension(ifull) :: xa_l, xb_l
logical, intent(in) :: lblock
logical, dimension(ifull_g) :: landg  ! large common array

if (myid==0) then
  write(6,*) "MLO 2D scale-selective filter"
  if (kd==1) then
    write(6,*) "Single level filter"
  else
    write(6,*) "Multiple level filter"
  end if
end if

if (nud_sst/=0) then
  call ccmpi_gatherall(diff_l(:,1:kd),diff_g(:,1:kd))
  landg=abs(diff_g(:,1)-miss)<0.1
  call mlofilterhost(diff_g,diff_l,kd,landg)
end if
if (nud_sss/=0) then
  call ccmpi_gatherall(diffs_l(:,1:kd),diff_g(:,1:kd))
  landg=abs(diff_g(:,1)-miss)<0.1
  call mlofilterhost(diff_g,diffs_l,kd,landg)
end if
if (nud_ouv/=0) then
  do k=1,kd
    xa_l=diffu_l(:,k)
    xb_l=diffv_l(:,k)
    diffu_l(:,k)=ax(1:ifull)*xa_l+bx(1:ifull)*xb_l
    diffv_l(:,k)=ay(1:ifull)*xa_l+by(1:ifull)*xb_l
    diffw_l(:,k)=az(1:ifull)*xa_l+bz(1:ifull)*xb_l
    where (abs(xa_l-miss)<0.1)
      diffu_l(:,k)=miss
      diffv_l(:,k)=miss
      diffw_l(:,k)=miss
    end where
  end do        
  call ccmpi_gatherall(diffu_l(:,1:kd),diff_g(:,1:kd))
  landg=abs(diff_g(:,1)-miss)<0.1
  call mlofilterhost(diff_g,diffu_l,kd,landg)
  call ccmpi_gatherall(diffv_l(:,1:kd),diff_g(:,1:kd))
  call mlofilterhost(diff_g,diffv_l,kd,landg)
  call ccmpi_gatherall(diffw_l(:,1:kd),diff_g(:,1:kd))
  call mlofilterhost(diff_g,diffw_l,kd,landg)
end if
if (nud_sfh/=0.and.lblock) then
  call ccmpi_gatherall(diffh_l(:,1),diff_g(:,1))
  landg=abs(diff_g(:,1)-miss)<0.1
  call mlofilterhost(diff_g,diffh_l,1,landg)
end if

#ifdef debug
if (myid==0.and.nmaxpr==1) then
  write(6,*) "MLO end 2D filter"
end if
#endif
      
return
end subroutine mlofilter

subroutine mlofilterhost(diff_g,dd,kd,landg)

use cc_mpi             ! CC MPI routines
use map_m              ! Grid map arrays
use vecsuv_m           ! Map to cartesian coordinates
use xyzinfo_m          ! Grid coordinate arrays

implicit none

include 'newmpar.h'    ! Grid parameters
include 'const_phys.h' ! Physical constants
include 'parm.h'       ! Model configuration
include 'parmgeom.h'   ! Coordinate data

integer, intent(in) :: kd
integer i, j, n, iqq, iqqg, k
real nsum, cq
real, dimension(ifull_g,kd), intent(inout) :: diff_g ! large common array
real, dimension(ifull_g) :: rr
real, dimension(ifull_g) :: mm, nn                   ! large common array
real, dimension(ifull,kd), intent(out) :: dd
logical, dimension(ifull_g), intent(in) :: landg     ! large common array

! eventually will be replaced with mbd once full ocean coupling is complete
cq = sqrt(4.5)*.1*real(max( nud_sst, nud_sss, nud_ouv, nud_sfh, mbd ))/(pi*schmidt)

call START_LOG(nestpack_begin)
dd(:,:) = 0.
mm(:) = 1./(em_g(:)*em_g(:))
where( .not.landg(:) )
  nn(:) = mm(:)
elsewhere
  nn(:) = 0.
end where
do k = 1,kd
  diff_g(:,k) = diff_g(:,k)*nn(:)
end do
call END_LOG(nestpack_end)

call START_LOG(nestcalc_begin)
do n = 1,npan
  do j = 1,jpan
    do i = 1,ipan
      iqqg = i + ioff + (j+joff-1)*il_g + (n-noff)*il_g*il_g
      iqq = i + (j-1)*ipan + (n-1)*ipan*jpan
      if ( .not.landg(iqqg) ) then
        rr(:) = real(x_g(iqqg)*x_g(:)+y_g(iqqg)*y_g(:)+z_g(iqqg)*z_g(:))
        rr(:) = acos(max( min( rr(:), 1. ), -1. ))
        rr(:) = exp(-(cq*rr(:))**2)
        nsum = sum(rr(:)*mm(:))
        do k = 1,kd
          dd(iqq,k) = sum(rr(:)*diff_g(:,k))/nsum
        end do
      end if
    end do
  end do
end do
call END_LOG(nestcalc_end)

return
end subroutine mlofilterhost

! 1D filer for mlo
subroutine mlofilterfast(diff_l,diffs_l,diffu_l,diffv_l,diffh_l,miss,lblock,kd)

use cc_mpi                  ! CC MPI routines

implicit none

include 'newmpar.h'         ! Grid parameters
include 'const_phys.h'      ! Physical constants
include 'parm.h'            ! Model configuration
include 'parmgeom.h'        ! Coordinate data

integer, intent(in) :: kd
real, intent(in) :: miss      
real, dimension(ifull,1), intent(inout) :: diffh_l
real, dimension(ifull,kd), intent(inout) :: diff_l, diffs_l
real, dimension(ifull,kd), intent(inout) :: diffu_l, diffv_l
real cq
logical, intent(in) :: lblock
      
! eventually will be replaced with mbd once full ocean coupling is complete
cq = sqrt(4.5)*.1*real(max( nud_sst, nud_sss, nud_ouv, nud_sfh, mbd ))/(pi*schmidt)
      
#ifdef uniform_decomp
if ( myid==0 ) then
  write(6,*) "MLO 1D scale-selective filter"
  if ( kd==1 ) then
    write(6,*) "Single level filter"
  else
    write(6,*) "Multiple level filter"
  end if
end if        
#else
if ( myid==0 ) then
  write(6,*) "MLO 1D scale-selective filter (optimised)"
  if ( kd==1 ) then
    write(6,*) "Single level filter"
  else
    write(6,*) "Multiple level filter"
  end if
end if
#endif

if ( pprocn==pprocx ) then
  call mlospechost_n(cq,diff_l,diffs_l,diffu_l,diffv_l,diffh_l,miss,lblock,kd)
else
  call mlospechost(cq,diff_l,diffs_l,diffu_l,diffv_l,diffh_l,miss,lblock,kd)
end if

#ifdef debug
if ( myid==0 .and. nmaxpr==1 ) then
  write(6,*) "MLO end 1D filter"
end if
#endif

return
end subroutine mlofilterfast

subroutine mlospechost(cq,diff_l,diffs_l,diffu_l,diffv_l,diffh_l,miss,lblock,kd)

use cc_mpi             ! CC MPI routines
use vecsuv_m           ! Map to cartesian coordinates
      
implicit none
      
include 'newmpar.h'    ! Grid parameters
include 'parm.h'       ! Model configuration
      
integer, intent(in) :: kd
integer n,k,ppass
real, intent(in) :: cq,miss
real, dimension(ifull,1), intent(inout) :: diffh_l
real, dimension(ifull,kd), intent(inout) :: diff_l,diffs_l
real, dimension(ifull,kd), intent(inout) :: diffu_l,diffv_l
real, dimension(ifull,kd) :: diffw_l
real, dimension(ipan*jpan,kd) :: qp
real, dimension(ifull) :: xa_l,xb_l
logical, intent(in) :: lblock

if (nud_sst/=0) then
  call ccmpi_gathermap(diff_l(:,:),kd)         ! gather data onto global sparse array (1)
  do ppass=pprocn,pprocx
    call copyglobalpack(kd,0,kd)               ! copy sparse array (1) to (0)
    call mlofastspec_work(cq,qp,kd,ppass,miss) ! filter sparse array (0)
    do k=1,kd
      do n=1,ipan*jpan
        diff_l(n+ipan*jpan*(ppass-pprocn),k)=qp(n,k)
      end do
    end do
  end do
end if
if (nud_sss/=0) then
  call ccmpi_gathermap(diffs_l(:,:),kd)        ! gather data onto global sparse array (1)
  do ppass=pprocn,pprocx
    call copyglobalpack(kd,0,kd)               ! copy sparse array (1) to (0)
    call mlofastspec_work(cq,qp,kd,ppass,miss) ! filter sparse array (0)
    do k=1,kd
      do n=1,ipan*jpan
        diffs_l(n+ipan*jpan*(ppass-pprocn),k)=qp(n,k)
      end do
    end do
  end do
end if
if (nud_ouv/=0) then
  do k=1,kd
    xa_l=diffu_l(:,k)
    xb_l=diffv_l(:,k)
    diffu_l(:,k)=ax(1:ifull)*xa_l+bx(1:ifull)*xb_l
    diffv_l(:,k)=ay(1:ifull)*xa_l+by(1:ifull)*xb_l
    diffw_l(:,k)=az(1:ifull)*xa_l+bz(1:ifull)*xb_l
    where (abs(xa_l-miss)<0.1)
      diffu_l(:,k)=miss
      diffv_l(:,k)=miss
      diffw_l(:,k)=miss
    end where
  end do
  call ccmpi_gathermap(diffu_l(:,:),kd)        ! gather data onto global sparse array (1)
  do ppass=pprocn,pprocx
    call copyglobalpack(kd,0,kd)               ! copy sparse array (1) to (0)
    call mlofastspec_work(cq,qp,kd,ppass,miss) ! filter sparse array (0)
    do k=1,kd
      do n=1,ipan*jpan
        diffu_l(n+ipan*jpan*(ppass-pprocn),k)=qp(n,k)
      end do
    end do
  end do
  call ccmpi_gathermap(diffv_l(:,:),kd)        ! gather data onto global sparse array (1)
  do ppass=pprocn,pprocx
    call copyglobalpack(kd,0,kd)               ! copy sparse array (1) to (0)
    call mlofastspec_work(cq,qp,kd,ppass,miss) ! filter sparse array (0)
    do k=1,kd
      do n=1,ipan*jpan
        diffv_l(n+ipan*jpan*(ppass-pprocn),k)=qp(n,k)
      end do
    end do
  end do
  call ccmpi_gathermap(diffw_l(:,:),kd)        ! gather data onto global sparse array (1)
  do ppass=pprocn,pprocx
    call copyglobalpack(kd,0,kd)               ! copy sparse array (1) to (0)
    call mlofastspec_work(cq,qp,kd,ppass,miss) ! filter sparse array (0)
    do k=1,kd
      do n=1,ipan*jpan
        diffw_l(n+ipan*jpan*(ppass-pprocn),k)=qp(n,k)
      end do
    end do
  end do
  do k=1,kd
    xa_l=ax(1:ifull)*diffu_l(:,k)+ay(1:ifull)*diffv_l(:,k)+az(1:ifull)*diffw_l(:,k)
    xb_l=bx(1:ifull)*diffu_l(:,k)+by(1:ifull)*diffv_l(:,k)+bz(1:ifull)*diffw_l(:,k)
    diffu_l(:,k)=xa_l
    diffv_l(:,k)=xb_l
  end do
endif
if (nud_sfh/=0.and.lblock) then
  call ccmpi_gathermap(diffh_l(:,1),1)        ! gather data onto global sparse array (1)
  do ppass=pprocn,pprocx
    call copyglobalpack(1,0,1)                ! copy sparse array (1) to (0)
    call mlofastspec_work(cq,qp,1,ppass,miss) ! filter sparse array (0)
    do n=1,ipan*jpan
      diffh_l(n+ipan*jpan*(ppass-pprocn),1)=qp(n,1)
    end do
  end do
end if

return
end subroutine mlospechost

!---------------------------------------------------------------------------------
! memory reduced version of mlospechost
subroutine mlospechost_n(cq,diff_l,diffs_l,diffu_l,diffv_l,diffh_l,miss,lblock,kd)

use cc_mpi             ! CC MPI routines
use vecsuv_m           ! Map to cartesian coordinates
      
implicit none
      
include 'newmpar.h'    ! Grid parameters
include 'parm.h'       ! Model configuration
      
integer, intent(in) :: kd
integer k
real, intent(in) :: cq,miss
real, dimension(ifull,1), intent(inout) :: diffh_l
real, dimension(ifull,kd), intent(inout) :: diff_l,diffs_l
real, dimension(ifull,kd), intent(inout) :: diffu_l,diffv_l
real, dimension(ifull,kd) :: diffw_l
real, dimension(ifull) :: xa_l,xb_l
logical, intent(in) :: lblock

if (nud_sst/=0) then
  call ccmpi_gathermap(diff_l(:,:),0)             ! gather data onto global sparse array (0)
  call mlofastspec_work(cq,diff_l,kd,pprocn,miss) ! filter sparse array (0)
end if
if (nud_sss/=0) then
  call ccmpi_gathermap(diffs_l(:,:),0)             ! gather data onto global sparse array (0)
  call mlofastspec_work(cq,diffs_l,kd,pprocn,miss) ! filter sparse array (0)
end if
if (nud_ouv/=0) then
  do k=1,kd
    xa_l=diffu_l(:,k)
    xb_l=diffv_l(:,k)
    diffu_l(:,k)=ax(1:ifull)*xa_l+bx(1:ifull)*xb_l
    diffv_l(:,k)=ay(1:ifull)*xa_l+by(1:ifull)*xb_l
    diffw_l(:,k)=az(1:ifull)*xa_l+bz(1:ifull)*xb_l
    where (abs(xa_l-miss)<0.1)
      diffu_l(:,k)=miss
      diffv_l(:,k)=miss
      diffw_l(:,k)=miss
    end where
  end do
  call ccmpi_gathermap(diffu_l(:,:),0)             ! gather data onto global sparse array (0)
  call mlofastspec_work(cq,diffu_l,kd,pprocn,miss) ! filter sparse array (0)
  call ccmpi_gathermap(diffv_l(:,:),0)             ! gather data onto global sparse array (0)
  call mlofastspec_work(cq,diffv_l,kd,pprocn,miss) ! filter sparse array (0)
  call ccmpi_gathermap(diffw_l(:,:),0)             ! gather data onto global sparse array (0)
  call mlofastspec_work(cq,diffw_l,kd,pprocn,miss) ! filter sparse array (0)
  do k=1,kd
    xa_l=ax(1:ifull)*diffu_l(:,k)+ay(1:ifull)*diffv_l(:,k)+az(1:ifull)*diffw_l(:,k)
    xb_l=bx(1:ifull)*diffu_l(:,k)+by(1:ifull)*diffv_l(:,k)+bz(1:ifull)*diffw_l(:,k)
    diffu_l(:,k)=xa_l
    diffv_l(:,k)=xb_l
  end do
endif
if (nud_sfh/=0.and.lblock) then
  call ccmpi_gathermap(diffh_l(:,1),0)           ! gather data onto global sparse array (0)
  call mlofastspec_work(cq,diffh_l(:,1),1,pprocn,miss) ! filter sparse array (0)
end if

return
end subroutine mlospechost_n

subroutine mlofastspec_work(cq,diff_g,kd,ppass,miss)

use cc_mpi

implicit none
      
include 'newmpar.h'

integer, intent(in) :: kd,ppass
integer xpan
real, intent(in) :: cq,miss
real, dimension(ipan*jpan,kd), intent(inout) :: diff_g

! computations for the local processor group
xpan = max( ipan, jpan )
select case(ppass)
  case(1,2,3)
    call mlospeclocal_left(cq,ppass,diff_g,kd,xpan,miss)
  case(0,4,5)
    call mlospeclocal_right(cq,ppass,diff_g,kd,xpan,miss)
end select

return
end subroutine mlofastspec_work
      
!---------------------------------------------------------------------------------
! This version is for asymmetric decomposition
subroutine mlospeclocal_left(cq,ppass,qp,kd,xpan,miss)

use cc_mpi             ! CC MPI routines
use map_m              ! Grid map arrays
use xyzinfo_m          ! Grid coordinate arrays
     
implicit none
      
include 'newmpar.h'    ! Grid parameters
include 'parm.h'       ! Model configuration
      
integer, intent(in) :: ppass, kd, xpan
integer j, n, ipass, ns, ne, os, oe
integer jpoff, ibase
integer me, k, til, sn, sy, a, b, c, jj, nn
integer ibeg, iend
integer, dimension(0:3) :: astr, bstr, cstr
integer, dimension(0:3) :: maps
real, intent(in) :: cq, miss
real, dimension(ipan*jpan,kd), intent(out) :: qp
real, dimension(4*il_g,kd) :: ap      
real, dimension(4*il_g) :: rr, asum
real, dimension(xpan,xpan,kd) :: pp      
real, dimension(xpan,xpan) :: psum
real, dimension(il_g*xpan*(kd+1)) :: zz
real, dimension(xpan*xpan*(kd+1)) :: yy
real(kind=8), dimension(4*il_g) :: xa, ya, za
      
maps = (/ il_g, il_g, 4*il_g, 3*il_g /)
til = il_g*il_g
astr = 0
bstr = 0
cstr = 0

ns = joff + 1
ne = joff + jpan
os = ioff + 1
oe = ioff + ipan
      
do ipass = 0,2
  me = maps(ipass)
  call getiqa(astr,bstr,cstr,me,ipass,ppass,il_g)

#ifdef debug
  if ( myid == 0 ) write(6,*) "MLO start convolution"
#endif

  do j = 1,jpan
      
    ! pack data from sparse arrays
    call START_LOG(nestpack_begin)  
    jj = j + ns - 1
    do sn = 1,me,il_g
      sy = (sn-1)/il_g
      a = astr(sy)
      b = bstr(sy)
      c = cstr(sy)
      ibeg = a*sn + b*jj + c
      iend = a*(sn+il_g-1) + b*jj + c
      xa(sn:sn+il_g-1) = x_g(ibeg:iend:a)
      ya(sn:sn+il_g-1) = y_g(ibeg:iend:a)
      za(sn:sn+il_g-1) = z_g(ibeg:iend:a)
      asum(sn:sn+il_g-1) = 1./em_g(ibeg:iend:a)**2
      do k = 1,kd
        call getglobalpack(ap(:,k),sn,ibeg,iend,k)
        where ( abs(ap(sn:sn+il_g-1,k)-miss)<0.1 ) ! land
          ap(sn:sn+il_g-1,k) = 0.
        elsewhere
          ap(sn:sn+il_g-1,k) = ap(sn:sn+il_g-1,k)*asum(sn:sn+il_g-1)
        end where
      end do
    end do
    call END_LOG(nestpack_end)
    
    ! start convolution
    call START_LOG(nestcalc_begin)
    do n = 1,ipan
      nn = n + os - 1
      rr(1:me) = real(xa(nn)*xa(1:me)+ya(nn)*ya(1:me)+za(nn)*za(1:me))
      rr(1:me) = acos(max( min( rr(1:me), 1. ), -1. ))
      rr(1:me) = exp(-(cq*rr(1:me))**2)
      psum(n,j) = sum(rr(1:me)*asum(1:me))
      do k = 1,kd
        pp(n,j,k) = sum(rr(1:me)*ap(1:me,k))
      end do
    end do
    call END_LOG(nestcalc_end)
    
  end do

#ifdef debug
  if ( myid == 0 ) then
    write(6,*) "MLO end conv"
    write(6,*) "MLO Send arrays to local host"
  end if
#endif

  ! unpack grid
  a = astr(0)
  b = bstr(0)
  c = cstr(0)

  ! gather data on host processors
  call START_LOG(nestcomm_begin)
  yy(1:ipan*jpan*kd) = reshape( pp(1:ipan,1:jpan,1:kd), (/ ipan*jpan*kd /) )
  yy(ipan*jpan*kd+1:ipan*jpan*kd+ipan*jpan) = reshape( psum(1:ipan,1:jpan), (/ ipan*jpan /) )
  call ccmpi_allgatherx(zz(1:il_g*ipan*(kd+1)),yy(1:ipan*jpan*(kd+1)),comm_cols)
  call END_LOG(nestcomm_end)
  
  ! unpack data to sparse arrays
  call START_LOG(nestunpack_begin)
  do jpoff = 0,il_g-1,jpan
    sy = jpoff/jpan
    do k = 1,kd
      do j = jpoff+1,jpoff+jpan
        ibase = 1 + ipan*(j-jpoff-1) + ipan*jpan*(k-1) + ipan*jpan*(kd+1)*sy
        ibeg = a*os + b*j + c
        iend = a*oe + b*j + c
        call setglobalpack(zz(:),ibase,ibeg,iend,k,trans=.false.)
      end do
    end do
    do j = jpoff+1,jpoff+jpan
      ibase = 1 + ipan*(j-jpoff-1) + ipan*jpan*kd + ipan*jpan*(kd+1)*sy
      ibeg = a*os + b*j + c
      iend = a*oe + b*j + c
      call setglobalpack(zz(:),ibase,ibeg,iend,0,trans=.false.)
    end do
  end do
  call END_LOG(nestunpack_end)
          
end do

ns = ioff + 1
ne = ioff + ipan
os = joff + 1
oe = joff + jpan

ipass = 3
me = maps(ipass)
call getiqa(astr,bstr,cstr,me,ipass,ppass,il_g)

#ifdef debug
if ( myid==0 ) write(6,*) "MLO start convolution"
#endif

do j = 1,ipan
    
  ! pack data from sparse arrays
  call START_LOG(nestpack_begin)
  jj = j + ns - 1
  do sn = 1,me,il_g
    sy = (sn-1)/il_g
    a = astr(sy)
    b = bstr(sy)
    c = cstr(sy)
    ibeg = a*sn + b*jj + c
    iend = a*(sn+il_g-1) + b*jj + c
    xa(sn:sn+il_g-1) = x_g(ibeg:iend:a)
    ya(sn:sn+il_g-1) = y_g(ibeg:iend:a)
    za(sn:sn+il_g-1) = z_g(ibeg:iend:a)
    call getglobalpack(asum(:),sn,ibeg,iend,0)
    do k = 1,kd
      call getglobalpack(ap(:,k),sn,ibeg,iend,k)
    end do
  end do
  call END_LOG(nestpack_end)
  
  ! start convolution
  call START_LOG(nestcalc_begin)
  do n = 1,jpan
    nn = n + os - 1
    rr(1:me) = real(xa(nn)*xa(1:me)+ya(nn)*ya(1:me)+za(nn)*za(1:me))
    rr(1:me) = acos(max( min( rr(1:me), 1. ), -1. ))
    rr(1:me) = exp(-(cq*rr(1:me))**2)
    psum(n,j) = sum(rr(1:me)*asum(1:me))
    do k = 1,kd
      pp(n,j,k) = sum(rr(1:me)*ap(1:me,k))
    end do
  end do
  call END_LOG(nestcalc_end)
  
end do

#ifdef debug
if ( myid == 0 ) then
  write(6,*) "MLO end conv"
  write(6,*) "MLO Send arrays to local host"
end if
#endif

! unpack data
call START_LOG(nestunpack_begin)
do k = 1,kd
  do j = 1,ipan
    do n = 1,jpan
      if ( psum(n,j)>1.E-8 ) then
        qp(j+ipan*(n-1),k) = pp(n,j,k)/psum(n,j)
      else
        qp(j+ipan*(n-1),k) = 0.
      end if
    end do
  end do
end do
call END_LOG(nestunpack_end)
      
return  
end subroutine mlospeclocal_left

 subroutine mlospeclocal_right(cq,ppass,qp,kd,xpan,miss)

use cc_mpi             ! CC MPI routines
use map_m              ! Grid map arrays
use xyzinfo_m          ! Grid coordinate arrays
     
implicit none
      
include 'newmpar.h'    ! Grid parameters
include 'parm.h'       ! Model configuration
      
integer, intent(in) :: ppass, kd, xpan
integer j, n, ipass, ns, ne, os, oe
integer jpoff, ibase
integer me, k, til, sn, sy, a, b, c, jj, nn
integer ibeg, iend
integer, dimension(0:3) :: astr, bstr, cstr
integer, dimension(0:3) :: maps
real, intent(in) :: cq, miss
real, dimension(ipan*jpan,kd), intent(out) :: qp
real, dimension(4*il_g,kd) :: ap      
real, dimension(4*il_g) :: rr, asum
real, dimension(xpan,xpan,kd) :: pp      
real, dimension(xpan,xpan) :: psum
real, dimension(il_g*xpan*(kd+1)) :: zz
real, dimension(xpan*xpan*(kd+1)) :: yy
real(kind=8), dimension(4*il_g) :: xa, ya, za
      
maps = (/ il_g, il_g, 4*il_g, 3*il_g /)
til = il_g*il_g
astr = 0
bstr = 0
cstr = 0

ns = ioff + 1
ne = ioff + ipan
os = joff + 1
oe  =joff + jpan
      
do ipass = 0,2
  me = maps(ipass)
  call getiqa(astr,bstr,cstr,me,ipass,ppass,il_g)

#ifdef debug
  if ( myid == 0 ) write(6,*) "MLO start convolution"
#endif

  do j = 1,ipan
      
    ! pack data from sparse arrays
    call START_LOG(nestpack_begin)
    jj = j + ns - 1
    do sn = 1,me,il_g
      sy = (sn-1)/il_g
      a = astr(sy)
      b = bstr(sy)
      c = cstr(sy)
      ibeg = a*sn + b*jj + c
      iend = a*(sn+il_g-1) + b*jj + c
      xa(sn:sn+il_g-1) = x_g(ibeg:iend:a)
      ya(sn:sn+il_g-1) = y_g(ibeg:iend:a)
      za(sn:sn+il_g-1) = z_g(ibeg:iend:a)
      asum(sn:sn+il_g-1) = 1./em_g(ibeg:iend:a)**2
      do k = 1,kd
        call getglobalpack(ap(:,k),sn,ibeg,iend,k)
        where ( abs(ap(sn:sn+il_g-1,k)-miss)<0.1 )
          ap(sn:sn+il_g-1,k) = 0.
        elsewhere
          ap(sn:sn+il_g-1,k) = ap(sn:sn+il_g-1,k)*asum(sn:sn+il_g-1)
        end where
      end do
    end do
    call END_LOG(nestpack_end)
    
    ! start convolution
    call START_LOG(nestcalc_begin)
    do n = 1,jpan
      nn = n + os - 1
      rr(1:me) = real(xa(nn)*xa(1:me)+ya(nn)*ya(1:me)+za(nn)*za(1:me))
      rr(1:me) = acos(max( min( rr(1:me), 1. ), -1. ))
      rr(1:me) = exp(-(cq*rr(1:me))**2)
      psum(n,j) = sum(rr(1:me)*asum(1:me))
      do k = 1,kd
        pp(n,j,k) = sum(rr(1:me)*ap(1:me,k))
      end do
    end do
    call END_LOG(nestcalc_end)
    
  end do

#ifdef debug
  if ( myid == 0 ) then
    write(6,*) "MLO end conv"
    write(6,*) "MLO Send arrays to local host"
  end if
#endif

  ! unpack grid
  a = astr(0)
  b = bstr(0)
  c = cstr(0)

  ! gather data on host processors
  call START_LOG(nestcomm_begin)
  yy(1:ipan*jpan*kd) = reshape( pp(1:jpan,1:ipan,1:kd), (/ ipan*jpan*kd /) )
  yy(ipan*jpan*kd+1:ipan*jpan*kd+ipan*jpan) = reshape( psum(1:jpan,1:ipan), (/ ipan*jpan /) )
  call ccmpi_allgatherx(zz(1:il_g*jpan*(kd+1)),yy(1:ipan*jpan*(kd+1)),comm_rows)
  call END_LOG(nestcomm_end)
  
  ! unpack data to sparse arrays
  call START_LOG(nestunpack_begin)
  do jpoff = 0,il_g-1,ipan
    sy = jpoff/ipan
    do k = 1,kd
      do j = jpoff+1,jpoff+ipan
        ibase = 1 + jpan*(j-jpoff-1) + ipan*jpan*(k-1) + ipan*jpan*(kd+1)*sy
        ibeg = a*os + b*j + c
        iend = a*oe + b*j + c
        call setglobalpack(zz(:),ibase,ibeg,iend,k,trans=.true.)
      end do
    end do
    do j = jpoff+1,jpoff+ipan
      ibase = 1 + jpan*(j-jpoff-1) + ipan*jpan*kd + ipan*jpan*(kd+1)*sy
      ibeg = a*os + b*j + c
      iend = a*oe + b*j + c
      call setglobalpack(zz(:),ibase,ibeg,iend,0,trans=.true.)
    end do
  end do
  call END_LOG(nestunpack_end)
          
end do

ns = joff + 1
ne = joff + jpan
os = ioff + 1
oe = ioff + ipan

ipass = 3
me = maps(ipass)
call getiqa(astr,bstr,cstr,me,ipass,ppass,il_g)

#ifdef debug
if ( myid == 0 ) write(6,*) "MLO start convolution"
#endif

do j = 1,jpan
    
  ! pack data from sparse arrays
  call START_LOG(nestpack_begin)  
  jj = j + ns - 1
  do sn = 1,me,il_g
    sy = (sn-1)/il_g
    a = astr(sy)
    b = bstr(sy)
    c = cstr(sy)
    ibeg = a*sn + b*jj + c
    iend = a*(sn+il_g-1) + b*jj + c
    xa(sn:sn+il_g-1) = x_g(ibeg:iend:a)
    ya(sn:sn+il_g-1) = y_g(ibeg:iend:a)
    za(sn:sn+il_g-1) = z_g(ibeg:iend:a)
    call getglobalpack(asum(:),sn,ibeg,iend,0)
    do k = 1,kd
      call getglobalpack(ap(:,k),sn,ibeg,iend,k)
    end do
  end do
  call END_LOG(nestpack_end)
  
  ! start convolution
  call START_LOG(nestcalc_begin)
  do n = 1,ipan
    nn = n + os - 1
    rr(1:me) = real(xa(nn)*xa(1:me)+ya(nn)*ya(1:me)+za(nn)*za(1:me))
    rr(1:me) = acos(max( min( rr(1:me), 1. ), -1. ))
    rr(1:me) = exp(-(cq*rr(1:me))**2)
    psum(n,j) = sum(rr(1:me)*asum(1:me))
    do k = 1,kd
      pp(n,j,k) = sum(rr(1:me)*ap(1:me,k))
    end do
  end do
  call END_LOG(nestcalc_end)
  
end do

#ifdef debug
if ( myid==0 ) then
  write(6,*) "MLO end conv"
  write(6,*) "MLO Send arrays to local host"
end if
#endif

! gather data on host processors
call START_LOG(nestunpack_begin)
do k = 1,kd
  do j = 1,jpan
    do n = 1,ipan
      if ( psum(n,j)>1.E-8 ) then
        qp(n+ipan*(j-1),k) = pp(n,j,k)/psum(n,j)
      else
        qp(n+ipan*(j-1),k) = 0.
      end if
    end do
  end do
end do
call END_LOG(nestunpack_end)
      
return  
end subroutine mlospeclocal_right


! Relaxtion method for mlo
subroutine mlonudge(new,sssb,suvb,sfh,wl)

use mlo, only : mloimport,mloexport,wlev ! Ocean physics and prognostic arrays
      
implicit none

include 'newmpar.h'        ! Grid parameters
include 'parm.h'           ! Model configuration

integer, intent(in) :: wl
integer k,ka,i
real, dimension(ifull), intent(in) :: sfh
real, dimension(ifull,wlev), intent(in) :: new,sssb
real, dimension(ifull,wlev,2), intent(in) :: suvb
real, dimension(ifull) :: old
real wgt
      
wgt=dt/real(nud_hrs*3600)
if (nud_sst/=0) then
  do k=ktopmlo,kbotmlo
    ka=min(k,wl)
    old=new(:,ka)
    call mloexport(0,old,k,0)
    old=old*(1.-wgt)+new(:,ka)*wgt
    call mloimport(0,old,k,0)
  end do
end if
      
if (nud_sss/=0) then
  do k=ktopmlo,kbotmlo
    ka=min(k,wl)
    old=sssb(:,ka)
    call mloexport(1,old,k,0)
    old=old*(1.-wgt)+sssb(:,ka)*wgt
    old=max(old,0.)
    call mloimport(1,old,k,0)
  end do
end if
      
if (nud_ouv/=0) then
  do i=2,3
    do k=ktopmlo,kbotmlo
      ka=min(k,wl)
      old=suvb(:,ka,i-1)
      call mloexport(i,old,k,0)
      old=old*(1.-wgt)+suvb(:,ka,i-1)*wgt
      call mloimport(i,old,k,0)
    end do
  end do
end if

if (nud_sfh/=0) then
  old=sfh
  call mloexport(4,old,0,0)
  old=old*(1.-wgt)+sfh*wgt
  call mloimport(4,old,0,0)
end if
      
return
end subroutine mlonudge

!--------------------------------------------------------------
! Initialise RMA windows used by 1D convolutions to retrieve
! data from other processor ranks
subroutine specinit
      
use cc_mpi
      
implicit none
      
include 'newmpar.h'
include 'parm.h'
      
integer ncount,ipass,ppass,me
integer n,j,jj,sn,sy,kx
integer iqg,ng,jg,ig
integer a,b,c,ns
integer iproc,rproc
integer, dimension(0:3) :: maps
integer, dimension(0:3) :: astr,bstr,cstr
logical, dimension(0:nproc-1) :: lproc
      
! length of the 1D convolution for each 'pass'
maps = (/ il_g, il_g, 4*il_g, 3*il_g /)
! flag for data required from processor rank
lproc(:) = .false.
      
! loop over 1D convolutions and determine rank of the required data
! Note that convolution directions are ordered to minimise message passing
do ppass = pprocn,pprocx
  select case(ppass)
    case(1, 2, 3) ! left
      ns = joff + 1
      do ipass = 0,2
        me = maps(ipass)
        call getiqa(astr,bstr,cstr,me,ipass,ppass,il_g)
        do j = 1,jpan
          jj = j + ns - 1
          do sn = 1,me,il_g
            sy = (sn-1)/il_g
            a = astr(sy)
            b = bstr(sy)
            c = cstr(sy)
            do n = sn,sn+il_g-1
              iqg = a*n + b*jj + c
              ! Global ig, jg, ng
              ng = (iqg - 1)/(il_g*il_g)
              jg = 1 + (iqg - ng*il_g*il_g - 1)/il_g
              ig = iqg - (jg - 1)*il_g - ng*il_g*il_g
              ! fproc converts global indices to processor rank
              lproc(fproc(ig,jg,ng)) = .true.
            end do
          end do
        end do
      end do
    case(0, 4, 5) ! right
      ns = ioff + 1
      do ipass = 0,2
        me = maps(ipass)
        call getiqa(astr,bstr,cstr,me,ipass,ppass,il_g)
        do j = 1,ipan
          jj = j + ns - 1
          do sn = 1,me,il_g
            sy = (sn-1)/il_g
            a = astr(sy)
            b = bstr(sy)
            c = cstr(sy)
            do n = sn,sn+il_g-1
              iqg = a*n + b*jj + c
              ! Global ig, jg, ng
              ng = (iqg - 1)/(il_g*il_g)
              jg = 1 + (iqg - ng*il_g*il_g - 1)/il_g
              ig = iqg - (jg - 1)*il_g - ng*il_g*il_g
              ! fproc converts global indices to processor rank
              lproc(fproc(ig,jg,ng)) = .true.
            end do
          end do
        end do
      end do
  end select
end do
      
! specify required RMA windows from the list of processor ranks in specmap
! cc_mpi employs specmap when calling gathermap
ncount = count(lproc)
allocate(specmap(ncount))
ncount = 0
do iproc = 0,nproc-1
  ! stagger reading of windows - does this make any difference with active RMA?
  rproc = modulo(myid+iproc, nproc)
  if ( lproc(rproc) ) then
    ncount = ncount + 1
    specmap(ncount) = rproc
  end if
end do
      
! Include final filter pass before allocating global sparse arrays
do ppass = pprocn,pprocx
  select case(ppass)
    case(1, 2, 3) ! left
      ns = ioff + 1
      ipass = 3
      me = maps(ipass)
      call getiqa(astr,bstr,cstr,me,ipass,ppass,il_g)
      do j = 1,ipan
        jj = j + ns - 1
        do sn = 1,me,il_g
          sy = (sn-1)/il_g
          a = astr(sy)
          b = bstr(sy)
          c = cstr(sy)
          do n = sn,sn+il_g-1
            iqg = a*n + b*jj + c
            ! Global ig, jg, ng
            ng = (iqg - 1)/(il_g*il_g)
            jg = 1 + (iqg - ng*il_g*il_g - 1)/il_g
            ig = iqg - (jg - 1)*il_g - ng*il_g*il_g
            ! fproc converts global indices to processor rank
            lproc(fproc(ig,jg,ng)) = .true.
          end do
        end do
      end do
    case(0, 4, 5) ! right
      ns = joff + 1
      ipass = 3
      me = maps(ipass)
      call getiqa(astr,bstr,cstr,me,ipass,ppass,il_g)
      do j = 1,jpan
        jj = j + ns - 1
        do sn = 1,me,il_g
          sy = (sn-1)/il_g
          a = astr(sy)
          b = bstr(sy)
          c = cstr(sy)
          do n = sn,sn+il_g-1
            iqg = a*n + b*jj + c
            ! Global ig, jg, ng
            ng = (iqg - 1)/(il_g*il_g)
            jg = 1 + (iqg - ng*il_g*il_g - 1)/il_g
            ig = iqg - (jg - 1)*il_g - ng*il_g*il_g
            ! fproc converts global indices to processor rank
            lproc(fproc(ig,jg,ng)) = .true.
          end do
        end do
      end do
  end select
end do

ncount = count(lproc)
allocate(specmapext(ncount))
ncount = 0
do iproc = 0,nproc-1
   ! stagger reading of windows - does this make any difference with active RMA?
   rproc = modulo(myid+iproc,nproc)
   if ( lproc(rproc) ) then
     ncount = ncount+1
     specmapext(ncount) = rproc
   end if
 end do
   
if ( npta==1 ) then
  ! face version (nproc>=6)
  kx = min(max(kl,ol),kblock)
else
   ! normal version
   kx = 2*min(max(kl,ol),kblock) ! extra memory for copy
 end if
 call allocateglobalpack(kx)
      
return
end subroutine specinit

subroutine setdavvertwgt

use daviesnudge                  ! Far-field nudging
use sigs_m                       ! Atmosphere sigma levels

implicit none

include 'newmpar.h'              ! Grid parameters
include 'parm.h'                 ! Model configuration

integer klow, khigh, k
real siglow, sighigh

siglow = sig(kbotdav) - sigramplow
sighigh = sig(ktopdav) + sigramphigh

do klow = 1,kl-1
  if ( siglow>=sig(klow) ) exit
end do
do khigh = kl,2,-1
  if ( sighigh<=sig(khigh) ) exit
end do

do k = 1,kbotdav-1
  vertwgt(k) = 0.
end do
do k = kbotdav,klow-1
  vertwgt(k) = (sig(kbotdav)-sig(k))/(sig(kbotdav)-siglow)
end do
do k = klow,khigh
  vertwgt(k) = 1.
end do
do k = khigh+1,ktopdav
  vertwgt(k) = (sig(k)-sig(ktopdav))/(sighigh-sig(ktopdav))
end do
do k = ktopdav+1,kl
  vertwgt(k) = 0.
end do

return
end subroutine setdavvertwgt


!     Little function to convert kdate_r in form YYYYMMDD to no. of   !  Y2K
!     days since start of year of kdate. (LDR 3/1992, jlm 15/12/98,15/12/00)
!     Accounts for leap years. N.B. no exception in year 2000

!     N.B. this latest version is designed for nested runs running 1 month
!     at a time, so iyear0=iyear except for the case 
!     when kdate is Dec. and kdate_r is the following Jan.

integer function iabsdate(kdate_r,kdate)
      
implicit none
      
integer leap
common/leap_yr/leap
      
integer, intent(in) :: kdate_r, kdate
integer iyear,iyear0,month,iday
integer months,mon,mnth,nl
integer, dimension(12) :: mdays
data mdays/31,28,31,30,31,30,31,31,30,31,30,31/

iabsdate=0
iyear=kdate_r/10000
iyear0=kdate/10000                ! year of kdate
month=(kdate_r-10000*iyear)/100
iday=(kdate_r-10000*iyear)-100*month

! calculate number of months since start of kdate year
months=(iyear-iyear0)*12+month-1  

! Accumulate days month by month, up to last completed month
do mon=1,months
  mnth=mod(mon-1,12)+1
  iabsdate=iabsdate+mdays(mnth)
  if (leap==1.and.mnth==2) then      ! MJT bug fix
    nl=0                             ! MJT bug fix        
    if (mod(iyear0,4  )==0) nl=1     ! MJT bug fix
    if (mod(iyear0,100)==0) nl=0     ! MJT bug fix
    if (mod(iyear0,400)==0) nl=1     ! MJT bug fix
    iabsdate=iabsdate+nl             ! MJT bug fix
  end if
enddo

! Add days from this current month
iabsdate=iabsdate+iday

end function iabsdate

end module nesting
