! Conformal Cubic Atmospheric Model
    
! Copyright 2015-2024 Commonwealth Scientific Industrial Research Organisation (CSIRO)
    
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

implicit none

private
public nestin, nestinb, mlofilterhub, mlonudge, specinit
public mtimea, mtimeb, nestin_exit

integer, save :: mtimea = -1  ! previous mesonest time (mins)
integer, save :: mtimeb = -1  ! next mesonest time (mins)
integer, save :: mtimec = -1
integer, save :: wl = -1

real, dimension(:), allocatable, save :: pslb, tssb, fraciceb
real, dimension(:), allocatable, save :: psla, tssa
real, dimension(:), allocatable, save :: sicedepb
real, dimension(:,:), allocatable, save :: ta, ua, va, qa
real, dimension(:,:), allocatable, save :: tb, ub, vb, qb, ocndep
real, dimension(:,:,:), allocatable, save :: sssb, xtghostb
real, dimension(:,:,:), allocatable, save :: sssa, xtghosta


contains

!--------------------------------------------------------------
! FAR-FIELD NUDGING AND RELAXATION ROUTINES
! Called for nbd/=0
subroutine nestin
      
use aerosol_arrays               ! Aerosol arrays
use arrays_m                     ! Atmosphere dyamics prognostic arrays
use cc_mpi                       ! CC MPI routines
use dates_m                      ! Date data
use daviesnudge                  ! Far-field nudging
use diag_m                       ! Diagnostic routines
use indices_m                    ! Grid index arrays
use latlong_m                    ! Lat/lon coordinates
use mlo                          ! Ocean physics and prognostic arrays
use newmpar_m                    ! Grid parameters
use onthefly_m                   ! Input interpolation routines
use parm_m                       ! Model configuration
use pbl_m                        ! Boundary layer arrays
use soil_m                       ! Soil and surface data
use soilsnow_m                   ! Soil, snow and surface data
use stime_m                      ! File date data
      
integer, dimension(ifull) :: dumm
integer kdate_r, ktime_r, kdhour, kdmin, kddate
integer khour_r, kmin_r, khour, kmin
integer :: num=0
real, dimension(ifull,kl,8) :: dumv
real, dimension(ifull,wlev,4) :: dumaa
real, dimension(ifull,ms,3) :: dumg
real, dimension(ifull,3,3) :: dums
real, dimension(ifull,3) :: duma
real, dimension(ifull) :: zsb, timelt
real, dimension(2) :: depthcheck
real timerm, cona, conb
      
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
    allocate( sssa(ifull,wlev,4), sssb(ifull,wlev,8), ocndep(ifull,6) )
    sssb = 0.
    allocate( xtghosta(ifull,kl,naero) )
    allocate( xtghostb(ifull,kl,naero) )
    xtghostb = 0.

    if ( abs(io_in)==1 ) then
      call START_LOG(nestotf_begin)
      call onthefly(1,kdate_r,ktime_r,                            &
                    pslb,zsb,tssb,sicedepb,fraciceb,tb,ub,vb,qb,  &
                    dumg(:,:,1),dumg(:,:,2),dumg(:,:,3),          & !unused
                    duma(:,1),dumv(:,:,1),dumv(:,:,2),            & !unused
                    dumv(:,:,3),dumv(:,:,4),dumv(:,:,5),          & !unused
                    dumv(:,:,6),dumv(:,:,7),dumv(:,:,8),          & !unused
                    dums(:,:,1),dums(:,:,2),dums(:,:,3),          & !unused
                    duma(:,2),duma(:,3),dumm,                     & !unused
                    sssb,ocndep,xtghostb)
      call END_LOG(nestotf_end)
      tssb(:) = abs(tssb(:))
      !qb = max(qb,0.)
      call retopo(pslb,zsb,zs(1:ifull),tb,qb)
    else
      write(6,*) 'ERROR: Nudging requires abs(io_in)=1'
      call ccmpi_abort(-1)
    endif   ! (io_in==1)
    call setdavvertwgt
    ! record time of saved data
    !mtimeb = mtimer
    khour_r = ktime_r/100
    khour   = ktime/100
    kdhour = khour_r - khour
    kmin_r = ktime_r - 100*khour_r
    kmin   = ktime   - 100*khour
    kdmin  = kmin_r - kmin
    kddate = iabsdate(kdate_r,kdate) - iabsdate(kdate,kdate)
    mtimeb = 1440*kddate + 60*kdhour + kdmin
    if ( myid==0 ) then
      write(6,*) "Starting mtimea at ",mtimeb  
    end if
  endif       ! (.not.allocated(ta))
      
! transfer mtimeb fields to mtimea and update sice variables
  mtimea = mtimeb
  psla(:) = pslb(:)
  tssa(:) = tssb(:)
  ta(1:ifull,:) = tb(1:ifull,:)
  qa(1:ifull,:) = qb(1:ifull,:)
  ua(1:ifull,:) = ub(1:ifull,:)
  va(1:ifull,:) = vb(1:ifull,:)
  sssa(:,:,1:4) = sssb(:,:,1:4)
  xtghosta(:,:,:) = xtghostb(:,:,:)

  ! Read host atmospheric and ocean data for nudging      
  if ( abs(io_in)==1 ) then
    call START_LOG(nestotf_begin)
    call onthefly(1,kdate_r,ktime_r,                            &
                  pslb,zsb,tssb,sicedepb,fraciceb,tb,ub,vb,qb,  &
                  dumg(:,:,1),dumg(:,:,2),dumg(:,:,3),          & !unused
                  duma(:,1),dumv(:,:,1),dumv(:,:,2),            & !unused
                  dumv(:,:,3),dumv(:,:,4),dumv(:,:,5),          & !unused
                  dumv(:,:,6),dumv(:,:,7),dumv(:,:,8),          & !unused
                  dums(:,:,1),dums(:,:,2),dums(:,:,3),          & !unused
                  duma(:,2),duma(:,3),dumm,                     & !unused
                  sssb,ocndep,xtghostb)
    call END_LOG(nestotf_end)
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
  khour_r = ktime_r/100
  khour   = ktime/100
  kdhour = khour_r - khour
  kmin_r = ktime_r - 100*khour_r
  kmin   = ktime   - 100*khour
  kdmin = kmin_r - kmin
  kddate = iabsdate(kdate_r,kdate) - iabsdate(kdate,kdate)
  mtimeb = 1440*kddate + 60*kdhour + kdmin
  
  if ( myid==0 ) then
    write(6,*) 'nesting file has: kdate_r,ktime_r     ',kdate_r,ktime_r
    write(6,*) '                  kddate,kdhour,kdmin ',kddate,kdhour,kdmin
    write(6,*) 'kdate_r,iabsdate ',kdate_r,iabsdate(kdate_r,kdate)
    write(6,*) 'kdate,iabsdate   ',kdate,iabsdate(kdate,kdate)
    write(6,*) 'giving mtimeb = ',mtimeb
  end if

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
  
  
  call retopo(pslb,zsb,zs(1:ifull),tb,qb)
  
  
  if ( nmaxpr==1 .and. mydiag ) then
    write (6,"('100*pslb.wesn ',2p5f8.3)") pslb(idjd),pslb(iw(idjd)),pslb(ie(idjd)),pslb(is(idjd)),pslb(in(idjd))
    write(6,*) 'pslb out(idjd) :',pslb(idjd)
    write(6,*) 'after pslb print ',num
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

  
end if ! (mtimer>mtimeb)

! now use tt, uu, vv arrays for time interpolated values
timerm = real(ktau)*dt/60.   ! real value in minutes (in case dt < 60 seconds)
cona = (real(mtimeb)-timerm)/real(mtimeb-mtimea)
conb = (timerm-real(mtimea))/real(mtimeb-mtimea)
psls(:) = cona*psla(:) + conb*pslb(:)
tt (:,:) = cona*ta(:,:) + conb*tb(:,:)
qgg(:,:) = cona*qa(:,:) + conb*qb(:,:)
uu (:,:) = cona*ua(:,:) + conb*ub(:,:)
vv (:,:) = cona*va(:,:) + conb*vb(:,:)

! calculate time interpolated tss 
if ( namip==0 ) then     ! namip SSTs/sea-ice take precedence
  if ( nmlo==0 ) then
    where ( fraciceb(:)>0. .and. fracice(:)<1.e-20 .and. .not.land(1:ifull) )
      ! N.B. if already a sice point, keep present tice (in tggsn)
      tggsn(:,1) = min( 271.2, tssb(:) )
    end where
    ! SSTs read from host model
    where ( .not.land(1:ifull) )
      tss = cona*tssa + conb*tssb
      tgg(:,1) = tss
      sicedep(:) = sicedepb(:)
      fracice(:) = fraciceb(:)
    end where  ! (.not.land)
  else
    if ( nud_sst/=0 .or. nud_sss/=0 .or. nud_ouv/=0 .or. nud_sfh/=0 ) then
      ! nudge mlo
      dumaa = cona*sssa(:,:,1:4) + conb*sssb(:,:,1:4)
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
        !timelt = min( timelt, 271.2, tgg(:,1) )
        timelt = min( timelt+0.01, tgg(:,1) )
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

use aerosol_arrays               ! Aerosol arrays
use arrays_m                     ! Atmosphere dyamics prognostic arrays
use cc_mpi                       ! CC MPI routines
use dates_m                      ! Date data
use daviesnudge                  ! Far-field nudging
use diag_m                       ! Diagnostic routines
use indices_m                    ! Grid index arrays
use latlong_m                    ! Lat/lon coordinates
use mlo                          ! Ocean physics and prognostic arrays
use newmpar_m                    ! Grid parameters
use onthefly_m                   ! Input interpolation routines
use parm_m                       ! Model configuration
use pbl_m                        ! Boundary layer arrays
use savuvt_m                     ! Saved dynamic arrays
use soil_m                       ! Soil and surface data
use soilsnow_m                   ! Soil, snow and surface data
use stime_m                      ! File date data
 
integer, dimension(ifull) :: dumm
integer kdate_r, ktime_r, ntr
integer kdhour, kdmin, kddate, khour_r, khour, kmin_r, kmin
real, dimension(ifull,kl,naero) :: xtghostc
real, dimension(ifull,kl,8) :: dumv
real, dimension(ifull,wlev,4) :: sssc
real, dimension(ifull,ms,3) :: dumg
real, dimension(ifull,3,3) :: dums
real, dimension(ifull,kl) :: tc, uc, vc, qc
real, dimension(ifull,3) :: duma
real, dimension(ifull) :: pslc
real, dimension(ifull) :: zsb, timelt
real, dimension(2) :: depthcheck
real cona, timerm
 
! mtimer, mtimeb are in minutes
if ( ktau<100 .and. myid==0 ) then
  write(6,*) 'in nestinb ktau,mtimer,mtimec ',ktau,mtimer,mtimec
  write(6,*) 'with kdate_s,ktime_s >= ',kdate_s,ktime_s
end if

! Load new host data to be ready for next call to filter
if ( mtimer>mtimeb ) then

  ! allocate arrays on first call     
  if ( .not.allocated(tb) ) then
    allocate( tb(ifull,kl), ub(ifull,kl), vb(ifull,kl), qb(ifull,kl) )
    allocate( ta(ifull,kl), ua(ifull,kl), va(ifull,kl), qa(ifull,kl) )
    allocate( pslb(ifull), tssb(ifull), fraciceb(ifull) )
    allocate( psla(ifull), tssa(ifull) )
    allocate( sicedepb(ifull) )
    allocate( ocndep(ifull,6) )
    allocate( sssb(ifull,wlev,8) )
    allocate( sssa(ifull,wlev,4) )
    sssb = 0.
    allocate( xtghostb(ifull,kl,naero) )
    allocate( xtghosta(ifull,kl,naero) )
    xtghostb = 0.
    if ( abs(io_in)==1 ) then
      call START_LOG(nestotf_begin)
      call onthefly(1,kdate_r,ktime_r,                            &
                    pslb,zsb,tssb,sicedepb,fraciceb,tb,ub,vb,qb,  &
                    dumg(:,:,1),dumg(:,:,2),dumg(:,:,3),          & !unused
                    duma(:,1),dumv(:,:,1),dumv(:,:,2),            & !unused
                    dumv(:,:,3),dumv(:,:,4),dumv(:,:,5),          & !unused
                    dumv(:,:,6),dumv(:,:,7),dumv(:,:,8),          & !unused
                    dums(:,:,1),dums(:,:,2),dums(:,:,3),          & !unused
                    duma(:,2),duma(:,3),dumm,                     & !unused
                    sssb,ocndep,xtghostb)
      call END_LOG(nestotf_end)
      tssb(:) = abs(tssb(:))
    else
      write(6,*) 'ERROR: Scale-selective filter requires abs(io_in)=1'
      call ccmpi_abort(-1)
    endif   ! (abs(io_in)==1)
    ! initialise arrays for 1D filter
    if ( nud_uv/=9 ) then
      call specinit
    end if
    ! define vertical weights
    call setdavvertwgt
    khour_r = ktime_r/100
    khour   = ktime/100
    kdhour = khour_r - khour
    kmin_r = ktime_r - 100*khour_r
    kmin   = ktime   - 100*khour
    kdmin  = kmin_r - kmin
    kddate = iabsdate(kdate_r,kdate) - iabsdate(kdate,kdate)
    mtimeb = 1440*kddate + 60*kdhour + kdmin
    if ( myid==0 ) then
      write(6,*) "Starting mtimea at ",mtimeb  
    end if
  end if
  
  mtimea = mtimeb
  psla(:) = pslb(:)
  tssa(:) = tssb(:)
  ta(1:ifull,:) = tb(1:ifull,:)
  qa(1:ifull,:) = qb(1:ifull,:)
  ua(1:ifull,:) = ub(1:ifull,:)
  va(1:ifull,:) = vb(1:ifull,:)
  sssa(:,:,1:4) = sssb(:,:,1:4)
  xtghosta(:,:,:) = xtghostb(:,:,:)
          
  ! following (till end of subr) reads in next bunch of data in readiness
  ! read tb etc  - for globpea, straight into tb etc
  if ( abs(io_in)==1 ) then
    call START_LOG(nestotf_begin)
    call onthefly(1,kdate_r,ktime_r,                            &
                  pslb,zsb,tssb,sicedepb,fraciceb,tb,ub,vb,qb,  &
                  dumg(:,:,1),dumg(:,:,2),dumg(:,:,3),          & !unused
                  duma(:,1),dumv(:,:,1),dumv(:,:,2),            & !unused
                  dumv(:,:,3),dumv(:,:,4),dumv(:,:,5),          & !unused
                  dumv(:,:,6),dumv(:,:,7),dumv(:,:,8),          & !unused
                  dums(:,:,1),dums(:,:,2),dums(:,:,3),          & !unused
                  duma(:,2),duma(:,3),dumm,                     & !unused
                  sssb,ocndep,xtghostb)
    call END_LOG(nestotf_end)
  else
    write(6,*) 'ERROR: Scale-selective filter requires abs(io_in)=1'
    call ccmpi_abort(-1)
  endif   ! (abs(io_in)==1)
  tssb(:) = abs(tssb(:))  ! moved here Mar '03

  ! calculate time for next filter call
  khour_r = ktime_r/100
  khour   = ktime/100
  kdhour = khour_r - khour
  kmin_r = ktime_r - 100*khour_r
  kmin   = ktime   - 100*khour
  kdmin  = kmin_r - kmin
  kddate = iabsdate(kdate_r,kdate) - iabsdate(kdate,kdate)
  mtimeb = 1440*kddate + 60*kdhour + kdmin
  
  if ( myid==0 ) then
    write(6,*) 'nesting file has: kdate_r,ktime_r     ',kdate_r,ktime_r
    write(6,*) '                  kddate,kdhour,kdmin ',kddate,kdhour,kdmin
    write(6,*) 'kdate_r,iabsdate ',kdate_r,iabsdate(kdate_r,kdate)
    write(6,*) 'kdate,iabsdate   ',kdate,iabsdate(kdate,kdate)
    write(6,*) 'giving mtimeb = ',mtimeb
  end if
  
  ! adjust input data for change in orography
  call retopo(pslb,zsb,zs,tb,qb)

  if ( nud_period == -1 ) then
    mtimec = mtimeb
  else
    mtimec = min( mtimea + nud_period, mtimeb )
  end if
  
end if ! (mtimer>mtimeb)

! Apply filter to model data using previously loaded host data
if ( mtimer>=mtimec .and. mod(nint(ktau*dt),60)==0 ) then

  if ( nud_period/=-1 ) then
    mtimec = min( mtimec + nud_period, mtimeb )
  end if

  timerm = real(ktau)*dt/60.   ! real value in minutes (in case dt < 60 seconds)
  cona = (real(mtimeb)-timerm)/real(mtimeb-mtimea)
  
  ! atmospheric nudging if required
  if ( mbd/=0 ) then
    if ( nud_p/=0 .or. nud_t/=0 .or. nud_uv/=0 .or. nud_q/=0 .or. nud_aero/=0 ) then
      pslc(:) = cona*psla(:) + (1.-cona)*pslb(:) - psl(1:ifull)
      uc(:,:) = cona*ua(:,:) + (1.-cona)*ub(:,:) - u(1:ifull,:)
      vc(:,:) = cona*va(:,:) + (1.-cona)*vb(:,:) - v(1:ifull,:)
      tc(:,:) = cona*ta(:,:) + (1.-cona)*tb(:,:) - t(1:ifull,:)
      qc(:,:) = cona*qa(:,:) + (1.-cona)*qb(:,:) - qg(1:ifull,:)
      if ( abs(iaero)>=2 ) then
        do ntr = 1,naero
          xtghostc(:,:,ntr) = cona*xtghosta(:,:,ntr) + (1.-cona)*xtghostb(:,:,ntr) - xtg(1:ifull,:,ntr)        
        end do
      end if
      call getspecdata(pslc,uc,vc,tc,qc,xtghostc)
    end if
  end if  

  ! specify sea-ice if not AMIP or Mixed-Layer-Ocean
  if ( namip==0 ) then  ! namip SSTs/sea-ice take precedence
    if ( nmlo==0 ) then
      ! check whether present ice points should change to/from sice points
      where ( fraciceb(:)>0. .and. fracice(:)<1.e-20 .and. .not.land(1:ifull) )
        ! N.B. if already a sice point, keep present tice (in tggsn)
        tggsn(:,1) = min( 271.2, tssb(:) )
      end where
      ! update tss 
      where ( .not.land(1:ifull) )
        tss(:) = cona*tssa(:) + (1.-cona)*tssb(:)
        tgg(:,1) = tss(:)
        sicedep(:) = sicedepb(:)
        fracice(:) = fraciceb(:)
      end where  ! (.not.land(iq))
    else
      ! nudge Mixed-Layer-Ocean
      if ( mbd_mlo/=0 ) then  
        if ( nud_sst/=0 .or. nud_sss/=0 .or. nud_ouv/=0 .or. nud_sfh/=0 ) then
          sssc(:,:,1:4) = cona*sssa(:,:,1:4) + (1.-cona)*sssb(:,:,1:4)  
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
            !timelt(:) = min( timelt(:), 271.2, tgg(:,1) )
            timelt(:) = min( timelt(:)+0.01, tgg(:,1) )
            sssc(:,1,1) = (cona*tssa(:) + (1.-cona)*tssb(:))*(1.-fraciceb(:)) + timelt*fraciceb(:)
            sssc(:,1,1) = sssc(:,1,1) - wrtemp
          end if
          call mlofilterhub(sssc(:,:,1),sssc(:,:,2),sssc(:,:,3:4),ocndep(:,2),wl)
        end if
      end if  
    end if ! (nmlo==0) ..else..
  end if   ! (namip==0)
  
end if     ! (mtimer==mtimec).and.(mod(nint(ktau*dt),60)==0)

return
end subroutine nestinb

!--------------------------------------------------------------
! This subroutine gathers and distributes data for the
! scale-selective filter
subroutine getspecdata(pslbb,ubb,vbb,tbb,qbb,xtgbb)

use aerosol_arrays               ! Aerosol arrays
use arrays_m                     ! Atmosphere dyamics prognostic arrays
use cc_mpi                       ! CC MPI routines
use const_phys                   ! Physical constants
use daviesnudge                  ! Far-field nudging
use kuocom_m                     ! JLM convection
use liqwpar_m                    ! Cloud water mixing ratios
use newmpar_m                    ! Grid parameters
use nharrs_m                     ! Non-hydrostatic atmosphere arrays
use parm_m                       ! Model configuration
use parmdyn_m                    ! Dynamics parameters
use parmgeom_m                   ! Coordinate data
use savuvt_m                     ! Saved dynamic arrays
use savuv1_m                     ! Saved dynamic arrays
use sigs_m                       ! Atmosphere sigma levels
use vecsuv_m                     ! Map to cartesian coordinates
use work3sav_m                   ! Water and tracer saved arrays
use xyzinfo_m, only : x,y,z      ! Grid coordinate arrays
      
integer iq, k, ntr, kb, kln, klx, klt
real, dimension(ifull), intent(inout) :: pslbb
real, dimension(ifull) :: costh,sinth
real, dimension(ifull,kl), intent(inout) :: ubb, vbb, tbb, qbb
real, dimension(ifull,kl,naero), intent(inout) :: xtgbb
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
    dum(:) = costh(:)*ubb(:,k) - sinth(:)*vbb(:,k)
    ubb(:,k) = dum(:)
  end do
end if
    
! Loop over maximum block size
! kblock can be reduced to save memory
do kb = kbotdav,ktopdav,kblock
     
  kln = kb                          ! lower limit of block
  klx = min( kb+kblock-1, ktopdav ) ! upper limit of block
  klt = klx - kln + 1               ! number of levels in block
  lblock = (kb==kbotdav)            ! flag for first loop (include psl)
        
  !-----------------------------------------------------------------------
  ! select nudging option
  if ( nud_uv==9 ) then 
    if ( myid==0 ) then
      write(6,*) "ATM 2D filter                        ",kb,min(kb+kblock-1,ktopdav)
    end if
    call slowspecmpi(.1*real(mbd)/(pi*schmidt),pslbb,ubb,vbb,tbb,qbb,xtgbb,lblock,klt,kln,klx)
  else
    if ( myid==0 ) then
      write(6,*) "ATM 1D filter                        ",kb,min(kb+kblock-1,ktopdav)
    end if
    call specfastmpi(.1*real(mbd)/(pi*schmidt),pslbb,ubb,vbb,tbb,qbb,xtgbb,lblock,klt,kln,klx)
  endif  ! (nud_uv==9) .. else ..
  !-----------------------------------------------------------------------
        
end do

      
if ( nud_p>0 ) then
  psl(1:ifull) = psl(1:ifull) + pslbb(:)
  ps(1:ifull) = 1.e5*exp(psl(1:ifull))
end if
if ( nud_uv/=0 ) then
  if ( nud_uv==3 ) then
    do k = kbotdav,ktopdav
      dum(:) = ub(:,k)
      ubb(1:ifull,k) = costh(:)*dum(:)
      vbb(1:ifull,k) = -sinth(:)*dum(:)
    end do
  end if
  do k = kbotdav,ktopdav
    u(1:ifull,k) = u(1:ifull,k) + ubb(:,k)*vertwgt(k)
    v(1:ifull,k) = v(1:ifull,k) + vbb(:,k)*vertwgt(k)
    savu(1:ifull,k) = savu(1:ifull,k) + ubb(:,k)*vertwgt(k)
    savu1(1:ifull,k) = savu1(1:ifull,k) + ubb(:,k)*vertwgt(k)
    savu2(1:ifull,k) = savu2(1:ifull,k) + ubb(:,k)*vertwgt(k)
    savv(1:ifull,k) = savv(1:ifull,k) + vbb(:,k)*vertwgt(k)
    savv1(1:ifull,k) = savv1(1:ifull,k) + vbb(:,k)*vertwgt(k)
    savv2(1:ifull,k) = savv2(1:ifull,k) + vbb(:,k)*vertwgt(k)
  end do
end if
if ( nud_t>0 ) then
  do k = kbotdav,ktopdav
    t(1:ifull,k) = t(1:ifull,k) + tbb(:,k)*vertwgt(k)
  end do
end if
if ( nud_q>0 ) then
  do k = kbotdav,ktopdav
    qg(1:ifull,k) = max(qg(1:ifull,k)+qbb(:,k)*vertwgt(k), qgmin)
  end do
end if
if ( nud_t>0 .or. nud_q>0 ) then
  phi(:,1) = zs(1:ifull) + bet(1)*t(1:ifull,1)
  do k = 2,kl
    phi(:,k) = phi(:,k-1) + bet(k)*t(1:ifull,k) + betm(k)*t(1:ifull,k-1)
  end do
  phi(:,:) = phi(:,:) + phi_nh(:,:)
end if
if ( abs(iaero)>=2 .and. nud_aero>0 ) then
  do ntr = 1,naero
    do k = kbotdav,ktopdav
      xtg(1:ifull,k,ntr) = max(xtg(1:ifull,k,ntr)+xtgbb(:,k,ntr)*vertwgt(k), 0.)
    end do
  end do
end if

return
end subroutine getspecdata

!---------------------------------------------------------------------------------
! Slow 2D spectral downscaling - MPI version
! This option is an exact treatment of the filter
subroutine slowspecmpi(cin,pslbb,ubb,vbb,tbb,qbb,xtgbb,lblock,klt,kln,klx)

use aerosol_arrays    ! Aerosol arrays
use cc_mpi            ! CC MPI routines
use newmpar_m         ! Grid parameters
use parm_m            ! Model configuration
use vecsuv_m          ! Map to cartesian coordinates
      
integer, intent(in) :: klt, kln, klx
integer k, n
real, intent(in) :: cin
real, dimension(ifull), intent(inout) :: pslbb
real, dimension(ifull,kl), intent(inout) :: ubb, vbb
real, dimension(ifull,kl), intent(inout) :: tbb, qbb
real, dimension(ifull,kl,naero), intent(inout) :: xtgbb
real, dimension(ifull,kln:klx) :: wbb
real, dimension(ifull_g,klt) :: tt ! large common array
real, dimension(ifull) :: da, db
real cq
logical, intent(in) :: lblock

cq = sqrt(4.5)*cin

if ( nud_p>0 .and. lblock ) then
  ! Create global copy of host data on each processor.  This can require a lot of memory
  call ccmpi_gatherall(pslbb(:), tt(:,1))
  ! Apply 2D filter
  call slowspecmpi_work(cq,tt(:,1),pslbb,1)
end if

if ( nud_uv==3 ) then
  call ccmpi_gatherall(ubb(:,kln:klx),tt)
  call slowspecmpi_work(cq,tt,ubb(:,kln:klx),klt)
else if ( nud_uv>0 ) then
  ! vectors are processed as Cartesian coordinates (X,Y,Z),
  ! avoiding complications along panel boundaries
  do k = kln,klx
    da(:) = ubb(:,k)
    db(:) = vbb(:,k)
    ubb(:,k) = ax(1:ifull)*da(:) + bx(1:ifull)*db(:)
    vbb(:,k) = ay(1:ifull)*da(:) + by(1:ifull)*db(:)
    wbb(:,k) = az(1:ifull)*da(:) + bz(1:ifull)*db(:)
  end do
  call ccmpi_gatherall(ubb(:,kln:klx),tt)
  call slowspecmpi_work(cq,tt,ubb(:,kln:klx),klt)
  call ccmpi_gatherall(vbb(:,kln:klx),tt)
  call slowspecmpi_work(cq,tt,vbb(:,kln:klx),klt)
  call ccmpi_gatherall(wbb(:,kln:klx),tt)
  call slowspecmpi_work(cq,tt,wbb(:,kln:klx),klt)
  ! Convert Cartesian vectors back to Conformal Cubic vectors
  do k = kln,klx
    da(:) = ax(1:ifull)*ubb(:,k) + ay(1:ifull)*vbb(:,k) + az(1:ifull)*wbb(:,k)
    db(:) = bx(1:ifull)*ubb(:,k) + by(1:ifull)*vbb(:,k) + bz(1:ifull)*wbb(:,k)
    ubb(:,k) = da(:)
    vbb(:,k) = db(:)
  end do
end if

if ( nud_t>0 ) then
  call ccmpi_gatherall(tbb(:,kln:klx),tt)
  call slowspecmpi_work(cq,tt,tbb(:,kln:klx),klt)
end if

if ( nud_q>0 ) then
  call ccmpi_gatherall(qbb(:,kln:klx),tt)
  call slowspecmpi_work(cq,tt,qbb(:,kln:klx),klt)
end if

if ( abs(iaero)>=2 .and. nud_aero>0 ) then
  do n = 1,naero
    call ccmpi_gatherall(xtgbb(:,kln:klx,n),tt)
    call slowspecmpi_work(cq,tt,xtgbb(:,kln:klx,n),klt)
  end do
end if      

return
end subroutine slowspecmpi

subroutine slowspecmpi_work(cq,tt,tbb,klt)

use cc_mpi            ! CC MPI routines
use map_m             ! Grid map arrays
use newmpar_m         ! Grid parameters
use xyzinfo_m         ! Grid coordinate arrays
      
integer, intent(in) :: klt
integer iq, iqg, k, n, j, i, kltp1
real, dimension(ifull,klt), intent(out) :: tbb
real, dimension(ifull_g,klt), intent(in) :: tt
real, dimension(ifull_g) :: sm        ! large working array
real, dimension(ifull_g) :: xa, ya, za ! large working array
real, dimension(klt+1) :: local_sum
real, intent(in) :: cq
real, dimension(ifull_g,klt+1) :: tt_l
#ifdef GPU
real, dimension(ifull,klt+1) :: tbb_l
#endif

! evaluate the 2D convolution
call START_LOG(nestcalc_begin)

kltp1 = klt + 1

! discrete normalisation factor
sm = 1./em_g**2
xa = real(x_g)
ya = real(y_g)
za = real(z_g)

do k = 1,klt
  tt_l(:,k) = tt(:,k)*sm
end do
tt_l(:,klt+1) = sm(:)

#ifdef GPU

! GPU version
!$acc parallel loop collapse(2) copyin(cq,klt,kltp1,xa,ya,za,tt_l) &
!$acc   copyout(tbb_l) private(iq,iqg,n,j,i)
do k = 1,klt+1
  do iq = 1,ifull
    n = 1 + (iq-1)/(ipan*jpan)  ! In range 1 .. npan
    j = 1 + ( iq - (n-1)*(ipan*jpan) - 1) / ipan
    i = iq - (j-1)*ipan - (n-1)*ipan*jpan
    iqg = i+ioff + (j+joff-1)*il_g + (n-noff)*il_g**2
    ! apply low band pass filter
    tbb_l(iq,k) = drpdr_fast(iqg,cq,xa,ya,za,tt_l(:,k))
  end do  
end do
!$acc end parallel loop
do k = 1,klt
  tbb(:,k) = tbb_l(:,k)/tbb_l(:,klt+1)
end do

#else

! CPU version
!$omp parallel do schedule(static) private(iqg,iq,local_sum,i,j,n)
do iq = 1,ifull
  n = 1 + (iq-1)/(ipan*jpan)  ! In range 1 .. npan
  j = 1 + ( iq - (n-1)*(ipan*jpan) - 1) / ipan
  i = iq - (j-1)*ipan - (n-1)*ipan*jpan
  iqg = i+ioff + (j+joff-1)*il_g + (n-noff)*il_g**2
  ! apply low band pass filter
  local_sum(1:kltp1) = drpdr_fast(iqg,cq,xa,ya,za,tt_l)
  tbb(iq,1:klt) = local_sum(1:klt)/local_sum(kltp1)
end do
!$omp end parallel do

#endif

call END_LOG(nestcalc_end)

return
end subroutine slowspecmpi_work
!---------------------------------------------------------------------------------

!---------------------------------------------------------------------------------
! Four pass spectral downscaling
subroutine specfastmpi(cin,psls,uu,vv,tt,qgg,xtgg,lblock,klt,kln,klx)
      
use aerosol_arrays     ! Aerosol arrays
use cc_mpi             ! CC MPI routines
use newmpar_m          ! Grid parameters
use parm_m             ! Model configuration
      
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
  ! general version
  call spechost(cin,psls,uu,vv,tt,qgg,xtgg,lblock,klt,kln,klx)
end if

return
end subroutine specfastmpi
!---------------------------------------------------------------------------------


!---------------------------------------------------------------------------------
! This is the main routine for the scale-selective filter
! (see spechost_n for a reduced memory version)
subroutine spechost(cin,pslbb,ubb,vbb,tbb,qbb,xtgbb,lblock,klt,kln,klx)

use aerosol_arrays     ! Aerosol arrays
use cc_mpi             ! CC MPI routines
use newmpar_m          ! Grid parameters
use parm_m             ! Model configuration
use vecsuv_m           ! Map to cartesian coordinates
      
integer, intent(in) :: klt, kln, klx
integer i, k, ppass, ibase
real, intent(in) :: cin
real, dimension(ifull), intent(inout) :: pslbb
real, dimension(ifull,kl), intent(inout) :: ubb, vbb
real, dimension(ifull,kl), intent(inout) :: tbb, qbb
real, dimension(ifull,kl,naero), intent(inout) :: xtgbb
real, dimension(ifull,kln:klx) :: wbb
real, dimension(ipan*jpan,klt) :: qt
real, dimension(ifull) :: da, db
logical, intent(in) :: lblock

if ( nud_p>0 .and. lblock ) then
  call START_LOG(nestwin_begin)
  call ccmpi_gathermap_send2(pslbb)        ! gather data onto global sparse array (1)
  call ccmpi_gathermap_recv2(1)
  call END_LOG(nestwin_end)
end if

if ( nud_uv==3 ) then
  call START_LOG(nestwin_begin)
  call ccmpi_gathermap_send3(ubb(:,kln:klx))        ! gather data onto global sparse array (1)
  call END_LOG(nestwin_end)
else if ( nud_uv>0 ) then
  do k = kln,klx
    da(:) = ubb(:,k)
    db(:) = vbb(:,k)
    ubb(:,k) = ax(1:ifull)*da(:) + bx(1:ifull)*db(:)
    vbb(:,k) = ay(1:ifull)*da(:) + by(1:ifull)*db(:)
    wbb(:,k) = az(1:ifull)*da(:) + bz(1:ifull)*db(:)
  end do
  call START_LOG(nestwin_begin)
  call ccmpi_gathermap_send3(ubb(:,kln:klx))   ! gather data onto global sparse array (1)
  call END_LOG(nestwin_end)
end if

if ( nud_p>0 .and. lblock ) then
  do ppass = pprocn,pprocx
    call copyglobalpack(1,0,1)            ! copy sparse array data (1) to (0)
    call fastspecmpi_work(cin,qt,1,ppass) ! filter sparse array (0)
    ibase=ipan*jpan*(ppass-pprocn)
    pslbb(1+ibase:ipan*jpan+ibase) = qt(1:ipan*jpan,1)
  end do
end if

if ( nud_uv==3 ) then
  call START_LOG(nestwin_begin)
  call ccmpi_gathermap_recv3(klt,klt)        ! gather data onto global sparse array (1)
  call END_LOG(nestwin_end)
  do ppass = pprocn,pprocx
    call copyglobalpack(klt,0,klt)               ! copy sparse array data (1) to (0)
    call fastspecmpi_work(cin,qt,klt,ppass)      ! filter sparse array (0)
    ibase = ipan*jpan*(ppass-pprocn)
    ubb(1+ibase:ipan*jpan+ibase,kln:kln+klt-1) = qt(1:ipan*jpan,1:klt)
  end do
else if ( nud_uv>0 ) then
  call START_LOG(nestwin_begin)
  call ccmpi_gathermap_recv3(klt,klt)          ! gather data onto global sparse array (1)
  call ccmpi_gathermap_send3(vbb(:,kln:klx))   ! gather data onto global sparse array (1)
  call END_LOG(nestwin_end)
  do ppass = pprocn,pprocx
    call copyglobalpack(klt,0,klt)          ! copy sparse array data (1) to (0)
    call fastspecmpi_work(cin,qt,klt,ppass) ! filter sparse array (0)
    ibase = ipan*jpan*(ppass-pprocn)
    ubb(1+ibase:ipan*jpan+ibase,kln:kln+klt-1) = qt(1:ipan*jpan,1:klt)
  end do
  call START_LOG(nestwin_begin)
  call ccmpi_gathermap_recv3(klt,klt)          ! gather data onto global sparse array (1)
  call ccmpi_gathermap_send3(wbb(:,kln:klx))   ! gather data onto global sparse array (1)
  call END_LOG(nestwin_end)
  do ppass = pprocn,pprocx
    call copyglobalpack(klt,0,klt)          ! copy sparse array data (1) to (0)
    call fastspecmpi_work(cin,qt,klt,ppass) ! filter sparse array (0)
    ibase = ipan*jpan*(ppass-pprocn)
    vbb(1+ibase:ipan*jpan+ibase,kln:kln+klt-1) = qt(1:ipan*jpan,1:klt)
  end do
  call START_LOG(nestwin_begin)
  call ccmpi_gathermap_recv3(klt,klt)   ! gather data onto global sparse array (1)
  call END_LOG(nestwin_end)
end if

if ( nud_t>0 ) then
  call START_LOG(nestwin_begin)
  call ccmpi_gathermap_send3(tbb(:,kln:klx))   ! gather data onto global sparse array (1)
  call END_LOG(nestwin_end)
end if

if ( nud_uv>0 .and. nud_uv/=3 ) then
  do ppass = pprocn,pprocx
    call copyglobalpack(klt,0,klt)          ! copy sparse array data (1) to (0)
    call fastspecmpi_work(cin,qt,klt,ppass) ! filter sparse array (0)
    ibase = ipan*jpan*(ppass-pprocn)
    wbb(1+ibase:ipan*jpan+ibase,kln:kln+klt-1) = qt(1:ipan*jpan,1:klt)
  end do
  do k = kln,klx
    da(:) = ax(1:ifull)*ubb(:,k) + ay(1:ifull)*vbb(:,k) + az(1:ifull)*wbb(:,k)
    db(:) = bx(1:ifull)*ubb(:,k) + by(1:ifull)*vbb(:,k) + bz(1:ifull)*wbb(:,k)
    ubb(:,k) = da(:)
    vbb(:,k) = db(:)
  end do
endif

if ( nud_t>0 ) then
  call START_LOG(nestwin_begin)
  call ccmpi_gathermap_recv3(klt,klt)   ! gather data onto global sparse array (1)
  call END_LOG(nestwin_end)
end if  

if ( nud_q>0 ) then
  call START_LOG(nestwin_begin)
  call ccmpi_gathermap_send3(qbb(:,kln:klx))   ! gather data onto global sparse array (1)
  call END_LOG(nestwin_end)
end if

if ( nud_t>0 ) then
  do ppass = pprocn,pprocx
    call copyglobalpack(klt,0,klt)          ! copy sparse array data (1) to (0)
    call fastspecmpi_work(cin,qt,klt,ppass) ! filter sparse array (0)
    ibase = ipan*jpan*(ppass-pprocn)
    tbb(1+ibase:ipan*jpan+ibase,kln:kln+klt-1) = qt(1:ipan*jpan,1:klt)
  end do
end if

if ( nud_q>0 ) then
  call START_LOG(nestwin_begin)
  call ccmpi_gathermap_recv3(klt,klt)   ! gather data onto global sparse array (1)
  call END_LOG(nestwin_end)
end if  

if ( abs(iaero)>=2 .and. nud_aero>0 ) then
  call START_LOG(nestwin_begin)  
  call ccmpi_gathermap_send3(xtgbb(:,kln:klx,1)) ! gather data onto global sparse array (1)
  call END_LOG(nestwin_end)
end if  

if ( nud_q>0 ) then
  do ppass = pprocn,pprocx
    call copyglobalpack(klt,0,klt)          ! copy sparse array data (1) to (0)
    call fastspecmpi_work(cin,qt,klt,ppass) ! filter sparse array (0)
    ibase = ipan*jpan*(ppass-pprocn)
    qbb(1+ibase:ipan*jpan+ibase,kln:kln+klt-1) = qt(1:ipan*jpan,1:klt)
  end do
end if

if ( abs(iaero)>=2 .and. nud_aero>0 ) then
  do i = 1,naero
    call START_LOG(nestwin_begin)  
    call ccmpi_gathermap_recv3(klt,klt) ! gather data onto global sparse array (1)
    call END_LOG(nestwin_end)
    if ( i<naero ) then
      call START_LOG(nestwin_begin)  
      call ccmpi_gathermap_send3(xtgbb(:,kln:klx,i+1)) ! gather data onto global sparse array (1)
      call END_LOG(nestwin_end)
    end if    
    do ppass = pprocn,pprocx
      call copyglobalpack(klt,0,klt)            ! copy sparse array data (1) to (0)
      call fastspecmpi_work(cin,qt,klt,ppass)   ! filter sparse array (0)
      ibase = ipan*jpan*(ppass-pprocn)
      xtgbb(1+ibase:ipan*jpan+ibase,kln:kln+klt-1,i) = qt(1:ipan*jpan,1:klt)
    end do
  end do
end if      

return
end subroutine spechost
!---------------------------------------------------------------------------------

! This version of spechost is for one panel per processor (reduced memory)
subroutine spechost_n(cin,pslbb,ubb,vbb,tbb,qbb,xtgbb,lblock,klt,kln,klx)

use aerosol_arrays     ! Aerosol arrays
use cc_mpi             ! CC MPI routines
use newmpar_m          ! Grid parameters
use parm_m             ! Model configuration
use vecsuv_m           ! Map to cartesian coordinates
      
integer, intent(in) :: klt,kln,klx
integer k,n
real, intent(in) :: cin
real, dimension(ifull), intent(inout) :: pslbb
real, dimension(ifull,kl), intent(inout) :: ubb, vbb
real, dimension(ifull,kl), intent(inout) :: tbb, qbb
real, dimension(ifull,kl,naero), intent(inout) :: xtgbb
real, dimension(ifull,kln:klx) :: wbb
real, dimension(ifull,klt) :: qt
real, dimension(ifull) :: da, db
logical, intent(in) :: lblock

if ( nud_p>0 .and. lblock ) then
  call START_LOG(nestwin_begin)  
  call ccmpi_gathermap_send2(pslbb) ! gather data onto global sparse array (0)
  call ccmpi_gathermap_recv2(0)     ! gather data onto global sparse array (0)
  call END_LOG(nestwin_end)
end if

if ( nud_uv==3 ) then
  call START_LOG(nestwin_begin)  
  call ccmpi_gathermap_send3(ubb(:,kln:klx))   ! gather data onto global sparse array (0)
  call END_LOG(nestwin_end)
else if ( nud_uv>0 ) then
  do k = kln,klx
    da(:) = ubb(:,k)
    db(:) = vbb(:,k)
    ubb(:,k) = ax(1:ifull)*da(:) + bx(1:ifull)*db(:)
    vbb(:,k) = ay(1:ifull)*da(:) + by(1:ifull)*db(:)
    wbb(:,k) = az(1:ifull)*da(:) + bz(1:ifull)*db(:)
  end do
  call START_LOG(nestwin_begin)
  call ccmpi_gathermap_send3(ubb(:,kln:klx))   ! gather data onto global sparse array (0)
  call END_LOG(nestwin_end)
end if

if ( nud_p>0 .and. lblock ) then
  call fastspecmpi_work(cin,qt(:,1),1,pprocn) ! filter sparse array (0)
  pslbb(:) = qt(:,1)
end if

if ( nud_uv==3 ) then
  call START_LOG(nestwin_begin)  
  call ccmpi_gathermap_recv3(klt,0)   ! gather data onto global sparse array (0)
  call END_LOG(nestwin_end)
  call fastspecmpi_work(cin,qt,klt,pprocn) ! filter sparse array (0)
  ubb(:,kln:klx) = qt(:,:)
else if ( nud_uv>0 ) then
  call START_LOG(nestwin_begin)
  call ccmpi_gathermap_recv3(klt,0)   ! gather data onto global sparse array (0)
  call ccmpi_gathermap_send3(vbb(:,kln:klx))   ! gather data onto global sparse array (0)
  call END_LOG(nestwin_end)
  call fastspecmpi_work(cin,qt,klt,pprocn) ! filter sparse array (0)
  ubb(:,kln:klx) = qt(:,:)
  call START_LOG(nestwin_begin)
  call ccmpi_gathermap_recv3(klt,0)   ! gather data onto global sparse array (0)
  call ccmpi_gathermap_send3(wbb(:,kln:klx))   ! gather data onto global sparse array (0)
  call END_LOG(nestwin_end)
  call fastspecmpi_work(cin,qt,klt,pprocn) ! filter sparse array (0)
  vbb(:,kln:klx) = qt(:,:)
  call START_LOG(nestwin_begin)
  call ccmpi_gathermap_recv3(klt,0)   ! gather data onto global sparse array (0)
  call END_LOG(nestwin_end)
end if

if ( nud_t>0 ) then
  call START_LOG(nestwin_begin)  
  call ccmpi_gathermap_send3(tbb(:,kln:klx))   ! gather data onto global sparse array (0)
  call END_LOG(nestwin_end)
end if

if ( nud_uv>0 .and. nud_uv/=3 ) then
  call fastspecmpi_work(cin,qt,klt,pprocn) ! filter sparse array (0)
  wbb(:,kln:klx) = qt(:,:)
  do k = kln,klx
    da(:) = ax(1:ifull)*ubb(:,k) + ay(1:ifull)*vbb(:,k) + az(1:ifull)*wbb(:,k)
    db(:) = bx(1:ifull)*ubb(:,k) + by(1:ifull)*vbb(:,k) + bz(1:ifull)*wbb(:,k)
    ubb(:,k) = da(:)
    vbb(:,k) = db(:)
  end do
end if

if ( nud_t>0 ) then
  call START_LOG(nestwin_begin)  
  call ccmpi_gathermap_recv3(klt,0)   ! gather data onto global sparse array (0)
  call END_LOG(nestwin_end)
end if

if ( nud_q>0 ) then
  call START_LOG(nestwin_begin)  
  call ccmpi_gathermap_send3(qbb(:,kln:klx))   ! gather data onto global sparse array (0)
  call END_LOG(nestwin_end)
end if

if ( nud_t>0 ) then
  call fastspecmpi_work(cin,qt,klt,pprocn) ! filter sparse array (0)
  tbb(:,kln:klx) = qt(:,:)
end if

if ( nud_q>0 ) then
  call START_LOG(nestwin_begin)  
  call ccmpi_gathermap_recv3(klt,0)   ! gather data onto global sparse array (0)
  call END_LOG(nestwin_end)
end if

if ( abs(iaero)>=2 .and. nud_aero>0 ) then
  call START_LOG(nestwin_begin)  
  call ccmpi_gathermap_send3(xtgbb(:,kln:klx,1)) ! gather data onto global sparse array (0)
  call END_LOG(nestwin_end) 
end if

if ( nud_q>0 ) then
  call fastspecmpi_work(cin,qt,klt,pprocn) ! filter sparse array (0)
  qbb(:,kln:klx) = qt(:,:)
end if

if ( abs(iaero)>=2 .and. nud_aero>0 ) then
  do n = 1,naero
    call START_LOG(nestwin_begin)  
    call ccmpi_gathermap_recv3(klt,0) ! gather data onto global sparse array (0)
    call END_LOG(nestwin_end)
    if ( n<naero ) then
      call START_LOG(nestwin_begin)  
      call ccmpi_gathermap_send3(xtgbb(:,kln:klx,n+1)) ! gather data onto global sparse array (0)
      call END_LOG(nestwin_end)
    end if  
    call fastspecmpi_work(cin,qt,klt,pprocn)   ! filter sparse array (0)
    xtgbb(:,kln:klx,n) = qt(:,:)
  end do
end if      

return
end subroutine spechost_n

! determine if panel is 'left' or 'right'.  This
! basically reorders the convolution so that the
! final convolution leaves the local processors
! data in the correct orientation.

subroutine fastspecmpi_work(cin,tt,klt,ppass)

use cc_mpi             ! CC MPI routines
use newmpar_m          ! Grid parameters
      
integer, intent(in) :: klt,ppass
real, intent(in) :: cin
real, dimension(ipan*jpan,klt), intent(inout) :: tt
real cq

cq = sqrt(4.5)*cin ! filter length scale

! computations for the local processor group
select case(ppass)
  case(1, 2, 3)
    call speclocal_left(cq,ppass,tt,klt)
  case(0, 4, 5)
    call speclocal_right(cq,ppass,tt,klt)
end select
     
return
end subroutine fastspecmpi_work
      
!---------------------------------------------------------------------------------
! This code runs between the local processors
      
subroutine speclocal_left(cq,ppass,qt,klt)

use cc_acc             ! CC OpenACC routines
use cc_mpi             ! CC MPI routines
use map_m              ! Grid map arrays
use newmpar_m          ! Grid parameters
use parm_m             ! Model configuration
use xyzinfo_m          ! Grid coordinate arrays

integer, intent(in) :: ppass, klt
integer :: j, k, n, ipass
integer :: jpoff, ibase
integer :: me, ns, ne, os, oe
integer :: til, a, b, c, sn, sy, jj, nn
integer :: ibeg, iend, kltp1
integer, dimension(0:3) :: astr, bstr, cstr
integer, dimension(0:3) :: maps
real, intent(in) :: cq
real, dimension(ipan*jpan,klt), intent(out) :: qt
real, dimension(il_g*ipan*(klt+1)*3) :: dd    ! subset of sparse array
real, dimension(ipan*jpan*(klt+1),0:2) :: ff
real, dimension(il_g) :: at, asum             ! subset of sparse array
real, dimension(klt+1) :: local_sum
real, dimension(4*il_g,max(ipan,jpan)) :: xa, ya, za ! subset of shared array
real, dimension(4*il_g,klt+1,max(ipan,jpan)) :: at_l         ! subset of sparse array
#ifdef GPU
integer :: async_counter
real, dimension(jpan,klt+1,ipan) :: qt_l
#endif
 
! matched for panels 1,2 and 3
      
maps = (/ il_g, il_g, 4*il_g, 3*il_g /)
til = il_g**2
kltp1 = klt + 1
astr = 0
bstr = 0
cstr = 0
qt = 0.

ns = joff + 1
ne = joff + jpan
os = ioff + 1
oe = ioff + ipan

call START_LOG(nestcalc_begin)

do ipass = 0,2
  me = maps(ipass)
  call getiqa(astr,bstr,cstr,me,ipass,ppass,il_g)

#ifndef GPU
  !$omp parallel do schedule(static) private(j,jj,sn,sy,a,b,c,ibeg,iend,asum,k,at)
#endif
  do j = 1,jpan
    ! pack data from sparse arrays
    jj = j + ns - 1
    do sn = 1,me,il_g
      sy = (sn-1)/il_g
      a = astr(sy)
      b = bstr(sy)
      c = cstr(sy)
      ibeg = a*sn + b*jj + c
      iend = a*(sn+il_g-1) + b*jj + c
      xa(sn:sn+il_g-1,j) = real(x_g(ibeg:iend:a))
      ya(sn:sn+il_g-1,j) = real(y_g(ibeg:iend:a))
      za(sn:sn+il_g-1,j) = real(z_g(ibeg:iend:a))
      asum = 1./em_g(ibeg:iend:a)**2
      do k = 1,klt
        call getglobalpack_v(at,ibeg,iend,k)
        at_l(sn:sn+il_g-1,k,j) = at*asum
      end do
      at_l(sn:sn+il_g-1,klt+1,j) = asum      
    end do
  end do
#ifndef GPU
  !$omp end parallel do
#endif
    
  ! start convolution

#ifdef GPU

  async_counter = mod(ipass,async_length)
  !$acc parallel loop collapse(3) copyin(xa(1:me,1:jpan),ya(1:me,1:jpan),za(1:me,1:jpan),at_l(1:me,:,1:jpan)) &
  !$acc   copyout(ff(:,ipass)) private(j,n,nn,ibase,k) async(async_counter)
  do j = 1,jpan
    do k = 1,klt+1
      do n = 1,ipan
        nn = n + os - 1
        ibase = n + (j-1)*ipan
        ff(ibase+(k-1)*ipan*jpan,ipass) = drpdr_fast(nn,cq,xa(1:me,j),ya(1:me,j),za(1:me,j),at_l(1:me,k,j))   
      end do  
    end do 
  end do 
  !$acc end parallel loop

#else

  !$omp parallel do collapse(2) schedule(static) private(j,n,nn,ibase,k,local_sum)
  do j = 1,jpan
    do n = 1,ipan
      nn = n + os - 1
      ibase = n + (j-1)*ipan
      local_sum(1:klt+1) = drpdr_fast(nn,cq,xa(1:me,j),ya(1:me,j),za(1:me,j),at_l(:,:,j))   
      do k = 1,klt+1
        ff(ibase+(k-1)*ipan*jpan,ipass) = local_sum(k)
      end do  
    end do 
  end do 
  !$omp end parallel do

#endif

end do ! ipass

#ifdef GPU
!$acc wait
#endif

call END_LOG(nestcalc_end)

call START_LOG(nestwin_begin)
call ccmpi_gathermap_wait
call END_LOG(nestwin_end)

call START_LOG(nestcomm_begin)

! gather data for final pass
call ccmpi_allgatherx(dd(1:il_g*ipan*(klt+1)*3),ff(1:ipan*jpan*(klt+1),0:2),comm_cols)

do ipass = 0,2
  me = maps(ipass)
  call getiqa(astr,bstr,cstr,me,ipass,ppass,il_g)

  ! unpacking grid
  a = astr(0)
  b = bstr(0)
  c = cstr(0)

  ! unpack to sparse arrays
  do n = 1,ipan
    nn = n + os - 1
    do k = 1,klt
      do jpoff = 0,il_g-1,jpan
        sy = jpoff/jpan
        ibase = n + ipan*jpan*(k-1) + ipan*jpan*(klt+1)*ipass + ipan*jpan*(klt+1)*3*sy
        at(1+jpoff:jpan+jpoff) = dd(ibase:ibase+ipan*(jpan-1):ipan)
      end do
      ibeg = a*nn + b*1 + c
      iend = a*nn + b*il_g + c
      call setglobalpack_v(at(1:il_g),ibeg,iend,k)
    end do  
    do jpoff = 0,il_g-1,jpan
      sy = jpoff/jpan
      ibase = n + ipan*jpan*klt + ipan*jpan*(klt+1)*ipass + ipan*jpan*(klt+1)*3*sy
      at(1+jpoff:jpan+jpoff) = dd(ibase:ibase+ipan*(jpan-1):ipan)
    end do
    ibeg = a*nn + b*1 + c
    iend = a*nn + b*il_g + c
    call setglobalpack_v(at(1:il_g),ibeg,iend,0)
  end do  

end do ! ipass

call END_LOG(nestcomm_end)

ns = ioff + 1
ne = ioff + ipan
os = joff + 1
oe = joff + jpan

ipass = 3
me = maps(ipass)
call getiqa(astr,bstr,cstr,me,ipass,ppass,il_g)

call START_LOG(nestcalc_begin)

#ifndef GPU
!$omp parallel do schedule(static) private(j,jj,sn,sy,a,b,c,ibeg,iend,asum,k,at)
#endif
do j = 1,ipan    
  ! pack from sparse arrays
  jj = j + ns - 1
  do sn = 1,me,il_g
    sy = (sn-1)/il_g
    a = astr(sy)
    b = bstr(sy)
    c = cstr(sy)
    ibeg = a*sn + b*jj + c
    iend = a*(sn+il_g-1) + b*jj + c
    xa(sn:sn+il_g-1,j) = real(x_g(ibeg:iend:a))
    ya(sn:sn+il_g-1,j) = real(y_g(ibeg:iend:a))
    za(sn:sn+il_g-1,j) = real(z_g(ibeg:iend:a))
    call getglobalpack_v(asum,ibeg,iend,0) 
    do k = 1,klt
      call getglobalpack_v(at,ibeg,iend,k) 
      at_l(sn:sn+il_g-1,k,j) = at
    end do
    at_l(sn:sn+il_g-1,klt+1,j) = asum
  end do
end do
#ifndef GPU
!$omp end parallel do
#endif
  
! start convolution
#ifdef GPU

! GPU version
!$acc parallel loop collapse(3) copyin(xa(1:me,1:ipan),ya(1:me,1:ipan),za(1:me,1:ipan),at_l(1:me,:,1:ipan)) &
!$acc   copyout(qt_l) private(j,n,nn)
do j = 1,ipan
  do k = 1,klt+1  
    do n = 1,jpan
      nn = n + os - 1
      qt_l(n,k,j) = drpdr_fast(nn,cq,xa(1:me,j),ya(1:me,j),za(1:me,j),at_l(1:me,k,j))
    end do  
  end do
end do
!$acc end parallel loop
do j = 1,ipan
  do k = 1,klt
    do n = 1,jpan
      qt(j+ipan*(n-1),k) = qt_l(n,k,j)/qt_l(n,klt+1,j)
    end do
  end do
end do

#else

! CPU version
!$omp parallel do collapse(2) schedule(static) private(j,n,nn,local_sum)
do j = 1,ipan
  do n = 1,jpan
    nn = n + os - 1
    local_sum(:) = drpdr_fast(nn,cq,xa(1:me,j),ya(1:me,j),za(1:me,j),at_l(:,:,j))
    qt(j+ipan*(n-1),1:klt) = local_sum(1:klt)/local_sum(kltp1) ! = dot_product(ra(1:me)*at(1:me,k))/dot_product(ra(1:me)*asum(1:me))
  end do
end do
!$omp end parallel do

#endif

call END_LOG(nestcalc_end)

return  
end subroutine speclocal_left

subroutine speclocal_right(cq,ppass,qt,klt)

use cc_acc             ! CC OpenACC routines
use cc_mpi             ! CC MPI routines
use map_m              ! Grid map arrays
use newmpar_m          ! Grid parameters
use parm_m             ! Model configuration
use xyzinfo_m          ! Grid coordinate arrays

integer, intent(in) :: ppass, klt
integer j, k, n, ipass
integer jpoff, ibase
integer me, ns, ne, os, oe
integer til, a, b, c, sn, sy, jj, nn
integer ibeg, iend, kltp1
integer, dimension(0:3) :: astr, bstr, cstr
integer, dimension(0:3) :: maps
real, intent(in) :: cq
real, dimension(ipan*jpan,klt), intent(out) :: qt
real, dimension(il_g) :: at, asum
real, dimension(il_g*jpan*(klt+1)*3) :: dd
real, dimension(ipan*jpan*(klt+1),0:2) :: ff
real, dimension(klt+1) :: local_sum
real, dimension(4*il_g,max(ipan,jpan)) :: xa, ya, za
real, dimension(4*il_g,klt+1,max(ipan,jpan)) :: at_l
#ifdef GPU
integer async_counter
real, dimension(ipan,klt+1,jpan) :: qt_l
#endif
      
! matched for panels 0, 4 and 5
      
maps = (/ il_g, il_g, 4*il_g, 3*il_g /)
til = il_g**2
kltp1 = klt + 1
astr = 0
bstr = 0
cstr = 0
qt = 0.

ns = ioff + 1
ne = ioff + ipan
os = joff + 1
oe = joff + jpan
  
call START_LOG(nestcalc_begin)

do ipass = 0,2
  me = maps(ipass)
  call getiqa(astr,bstr,cstr,me,ipass,ppass,il_g)

#ifndef GPU
  !$omp parallel do schedule(static) private(j,jj,sn,sy,a,b,c,ibeg,iend,asum,k,at)
#endif
  do j = 1,ipan      
    ! pack data from sparse arrays
    jj = j + ns - 1
    do sn = 1,me,il_g
      sy = (sn-1)/il_g
      a = astr(sy)
      b = bstr(sy)
      c = cstr(sy)
      ibeg = a*sn + b*jj + c
      iend = a*(sn+il_g-1) + b*jj + c
      xa(sn:sn+il_g-1,j) = real(x_g(ibeg:iend:a))
      ya(sn:sn+il_g-1,j) = real(y_g(ibeg:iend:a))
      za(sn:sn+il_g-1,j) = real(z_g(ibeg:iend:a))
      asum = 1./em_g(ibeg:iend:a)**2
      do k = 1,klt
        call getglobalpack_v(at,ibeg,iend,k) 
        at_l(sn:sn+il_g-1,k,j) = at*asum
      end do
      at_l(sn:sn+il_g-1,klt+1,j) = asum
    end do
  end do
#ifndef GPU
  !$omp end parallel do
#endif
  
  ! start convolution

#ifdef GPU

  ! GPU version
  async_counter = mod(ipass,async_length)
  !$acc parallel loop collapse(3) copyin(xa(1:me,1:ipan),ya(1:me,1:ipan),za(1:me,1:ipan),at_l(1:me,:,1:ipan)) &
  !$acc   copyout(ff(:,ipass)) private(j,n,nn,k) async(async_counter)
  do j = 1,ipan
    do k = 1,klt+1  
      do n = 1,jpan
        nn = n + os - 1
        ff(n+(j-1)*jpan+(k-1)*ipan*jpan,ipass) = drpdr_fast(nn,cq,xa(1:me,j),ya(1:me,j),za(1:me,j),at_l(1:me,k,j))
      end do  
    end do
  end do
  !$acc end parallel loop

#else

  ! CPU version
  !$omp parallel do collapse(2) schedule(static) private(j,n,nn,k,local_sum)
  do j = 1,ipan
    do n = 1,jpan
      nn = n + os - 1
      local_sum(1:klt+1) = drpdr_fast(nn,cq,xa(1:me,j),ya(1:me,j),za(1:me,j),at_l(:,:,j))
      do k = 1,klt+1
        ff(n+(j-1)*jpan+(k-1)*ipan*jpan,ipass) = local_sum(k)
      end do  
    end do
  end do
  !$omp end parallel do

#endif

end do ! ipass

#ifdef GPU
!$acc wait
#endif

call END_LOG(nestcalc_end)

call START_LOG(nestwin_begin)
call ccmpi_gathermap_wait
call END_LOG(nestwin_end)

call START_LOG(nestcomm_begin)

! gather data for final pass
call ccmpi_allgatherx(dd(1:il_g*jpan*(klt+1)*3),ff(1:ipan*jpan*(klt+1),0:2),comm_rows)

do ipass = 0,2
  me = maps(ipass)
  call getiqa(astr,bstr,cstr,me,ipass,ppass,il_g)

  ! unpacking grid
  a = astr(0)
  b = bstr(0)
  c = cstr(0)

  ! unpack data to sparse arrays
  do n = 1,jpan
    nn = n + os - 1
    do k = 1,klt
      do jpoff = 0,il_g-1,ipan
        sy = jpoff/ipan
        ibase = n + ipan*jpan*(k-1) + ipan*jpan*(klt+1)*ipass + ipan*jpan*(klt+1)*3*sy
        at(1+jpoff:ipan+jpoff) = dd(ibase:ibase+jpan*(ipan-1):jpan)
      end do
      ibeg = a*nn + b*1 + c
      iend = a*nn + b*il_g + c
      call setglobalpack_v(at(1:il_g),ibeg,iend,k)
    end do
    do jpoff = 0,il_g-1,ipan
      sy = jpoff/ipan
      ibase = n + ipan*jpan*klt + ipan*jpan*(klt+1)*ipass + ipan*jpan*(klt+1)*3*sy
      at(1+jpoff:ipan+jpoff) = dd(ibase:ibase+jpan*(ipan-1):jpan)
    end do
    ibeg = a*nn + b*1 + c
    iend = a*nn + b*il_g + c
    call setglobalpack_v(at(1:il_g),ibeg,iend,0)
  end do  

end do ! ipass

call END_LOG(nestcomm_end)    
    
ns = joff + 1
ne = joff + jpan
os = ioff + 1
oe = ioff + ipan

ipass = 3
me = maps(ipass)
call getiqa(astr,bstr,cstr,me,ipass,ppass,il_g)

call START_LOG(nestcalc_begin)

#ifndef GPU
!$omp parallel do schedule(static) private(j,jj,sn,sy,a,b,c,ibeg,iend,asum,k,at)
#endif
do j = 1,jpan    
  ! pack data from sparse arrays
  jj = j + ns - 1
  do sn = 1,me,il_g
    sy = (sn-1)/il_g
    a = astr(sy)
    b = bstr(sy)
    c = cstr(sy)
    ibeg = a*sn + b*jj + c
    iend = a*(sn+il_g-1) + b*jj + c
    xa(sn:sn+il_g-1,j) = real(x_g(ibeg:iend:a))
    ya(sn:sn+il_g-1,j) = real(y_g(ibeg:iend:a))
    za(sn:sn+il_g-1,j) = real(z_g(ibeg:iend:a))
    ! v version is faster for getglobalpack  
    call getglobalpack_v(asum,ibeg,iend,0)
    do k = 1,klt
      call getglobalpack_v(at,ibeg,iend,k)
      at_l(sn:sn+il_g-1,k,j) = at
    end do
    at_l(sn:sn+il_g-1,klt+1,j) = asum    
  end do
end do
#ifndef GPU
!$omp end parallel do
#endif

! start convolution

#ifdef GPU

! GPU version
!$acc parallel loop collapse(3) copyin(xa(1:me,1:jpan),ya(1:me,1:jpan),za(1:me,1:jpan),at_l(1:me,:,1:jpan)) &
!$acc   copyout(qt_l) private(j,n,nn)
do j = 1,jpan
  do k = 1,klt+1  
    do n = 1,ipan
      nn = n + os - 1
      qt_l(n,k,j) = drpdr_fast(nn,cq,xa(1:me,j),ya(1:me,j),za(1:me,j),at_l(1:me,k,j))
    end do  
  end do
end do
!$acc end parallel loop
do j = 1,jpan
  do k = 1,klt
    do n = 1,ipan
      qt(n+ipan*(j-1),k) = qt_l(n,k,j)/qt_l(n,klt+1,j)
    end do
  end do
end do

#else

! CPU version
!$omp parallel do collapse(2) schedule(static) private(j,n,nn,local_sum)
do j = 1,jpan
  do n = 1,ipan
    nn = n + os - 1
    local_sum = drpdr_fast(nn,cq,xa(1:me,j),ya(1:me,j),za(1:me,j),at_l(:,:,j))
    qt(n+ipan*(j-1),1:klt) = local_sum(1:klt)/local_sum(klt+1)
  end do
end do
!$omp end parallel do

#endif

call END_LOG(nestcalc_end)

return  
end subroutine speclocal_right

!---------------------------------------------------------------------------------
! Map from 1D convolution to global index
subroutine getiqa(a,b,c,ne,ipass,ppass,il_g)
      
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
                wlev,wrtemp,minsal,mloexpmelt
use mlodynamicsarrays_m                             ! Ocean dynamics data
use newmpar_m                                       ! Grid parameters
use parm_m                                          ! Model configuration
use soil_m                                          ! Soil and surface data
use vecsuv_m                                        ! Map to cartesian coordinates

integer, intent(in) :: wl
integer k,ka,kb,kc,kln,klx,klt,kbb
real, dimension(ifull), intent(in) :: sfh
real, dimension(ifull,wlev), intent(in) :: sstb,sssb
real, dimension(ifull,wlev,2), intent(in) :: suvb
real, dimension(ifull,1) :: diffh_l
real, dimension(ifull,kblock) :: diff_l,diffs_l
real, dimension(ifull,kblock) :: diffu_l,diffv_l
real, dimension(ifull) :: dz, old
logical lblock
logical, dimension(ifull,wlev) :: wtr
real nudgewgt
 
ka = 0
kb = 0
kc = min(kbotmlo, ktopmlo+wl-1)

if ( mloalpha==0 ) then
  nudgewgt = 1.
else if ( nud_period == -1 ) then
  nudgewgt = (mtimeb-mtimea)/(120.*real(mloalpha))  ! mloalpha=12 implies 24 hours
else
  nudgewgt = real(nud_period)/(120.*real(mloalpha)) ! mloalpha=12 implies 24 hours
end if

wtr(:,:) = .false.
do k = 1,wlev
  dz(:) = 0.  
  call mloexpdep(1,dz(:),k,0)
  where ( dz(:)>1.e-4 )
    wtr(:,k) = .true.  
  end where  
end do


if ( nud_sfh/=0 ) then
  old = sfh
  call mloexport("eta",old,0,0)
  where ( wtr(:,1) )
    diffh_l(:,1) = sfh - old
  elsewhere
    diffh_l(:,1) = 0.
  end where
end if
      
do kbb = ktopmlo,kc,kblock
                 
  kln = kbb
  klx = min(kbb+kblock-1, kc)
  klt = klx - kln + 1
  lblock = (kbb==ktopmlo)
      
  if ( nud_sst/=0 ) then
    do k = kln,klx
      kb = k - kln + 1
      old = sstb(:,k)
      call mloexport("temp",old,k,0)
      where ( wtr(:,k) )
        diff_l(:,kb) = sstb(:,k) - old
      elsewhere
        diff_l(:,kb) = 0.
      end where
    end do
  end if
  
  if ( nud_sss/=0 ) then
    do k = kln,klx
      kb = k - kln + 1
      old = sssb(:,k)
      call mloexport("sal",old,k,0)
      where ( wtr(:,k) .and. old>minsal .and. sssb(:,k)>minsal )
        diffs_l(:,kb) = sssb(:,k) - old
      elsewhere
        diffs_l(:,kb) = 0.
      end where
    end do
  end if

  if ( nud_ouv/=0 ) then
    do k = kln,klx
      kb = k - kln + 1
      old = suvb(:,k,1)
      call mloexport("u",old,k,0)
      where ( wtr(:,k) )
        diffu_l(:,kb) = suvb(:,k,1) - old
      elsewhere
        diffu_l(:,kb) = 0.
      end where
      old = suvb(:,k,2)
      call mloexport("v",old,k,0)
      where ( wtr(:,k) )
        diffv_l(:,kb) = suvb(:,k,2) - old
      elsewhere
        diffv_l(:,kb) = 0.
      end where
    end do
  end if


  ! Note that nmlo=1 or nmlo=-1 will use the 2D filter

  if ( (nud_uv/=9.and.abs(nmlo)/=1) .or. namip/=0 ) then
    if (myid==0) then
      if (klt==1) then
        write(6,*) "MLO 1D filter - Single level         ",kbb,min(kbb+kblock-1,kc)
      else
        write(6,*) "MLO 1D filter - Multiple level       ",kbb,min(kbb+kblock-1,kc)
      end if
    end if
    call mlofilterfast(diff_l(:,1:klt),diffs_l(:,1:klt),diffu_l(:,1:klt),diffv_l(:,1:klt),diffh_l(:,1),lblock,klt)
  else
    if (myid==0) then
      if (klt==1) then
        write(6,*) "MLO 2D filter - Single level         ",kbb,min(kbb+kblock-1,kc)
      else
        write(6,*) "MLO 2D filter - Multiple level       ",kbb,min(kbb+kblock-1,kc)
      end if
    end if
    call mlofilter(diff_l(:,1:klt),diffs_l(:,1:klt),diffu_l(:,1:klt),diffv_l(:,1:klt),diffh_l(:,1),lblock,klt)
  end if


  if ( nud_sst/=0 ) then
    !call mloexpmelt(timelt)
    do k = kln,klx
      ka = min(wl, k)
      kb = k - kln + 1
      old = sstb(:,ka)
      call mloexport("temp",old,k,0)
      old = old + min( max( diff_l(:,kb)*nudgewgt, 271.16-wrtemp-old ), 373.16-wrtemp-old )
      call mloimport("temp",old,k,0)
    end do
    if ( klx==kc ) then
      do k = kc+1,kbotmlo
        old = sstb(:,ka)
        call mloexport("temp",old,k,0)
        old = old + min( max( diff_l(:,kb)*nudgewgt, 271.16-wrtemp-old ), 373.16-wrtemp-old )
        call mloimport("temp",old,k,0)
      end do
    end if
  end if

  if ( nud_sss/=0 ) then
    do k = kln,klx
      ka = min(wl, k)
      kb = k-kln+1
      old = sssb(:,ka)
      call mloexport("sal",old,k,0)
      where ( old>minsal )
        diffs_l(:,kb) = max( diffs_l(:,kb), (minsal-old)/nudgewgt )
        old = old + diffs_l(:,kb)*nudgewgt
      end where  
      call mloimport("sal",old,k,0)
    end do
    if ( klx==kc ) then
      do k = kc+1,kbotmlo
        old = sssb(:,ka)
        call mloexport("sal",old,k,0)
        where ( old>minsal )
          diffs_l(:,kb) = max( diffs_l(:,kb), (minsal-old)/nudgewgt )  
          old = old + diffs_l(:,kb)*nudgewgt ! kb saved from above loop
        end where  
        call mloimport("sal",old,k,0)
      end do
    end if
  end if

  if ( nud_ouv/=0 ) then
    do k = kln,klx
      ka = min(wl, k)
      kb = k - kln + 1
      old = suvb(:,ka,1)
      call mloexport("u",old,k,0)
      old = old + diffu_l(:,kb)*nudgewgt
      call mloimport("u",old,k,0)
      if ( allocated(oldu1) ) then
        oldu1(:,k) = oldu1(:,k) + diffu_l(:,kb)*nudgewgt
        oldu2(:,k) = oldu2(:,k) + diffu_l(:,kb)*nudgewgt
      end if
    end do
    if ( klx==kc ) then
      do k = kc+1,kbotmlo
        old = suvb(:,ka,1)
        call mloexport("u",old,k,0)
        old = old + diffu_l(:,kb)*nudgewgt ! kb saved from above loop
        call mloimport("u",old,k,0)
        if ( allocated(oldu1) ) then
          oldu1(:,k) = oldu1(:,k) + diffu_l(:,kb)*nudgewgt
          oldu2(:,k) = oldu2(:,k) + diffu_l(:,kb)*nudgewgt
        end if
      end do
    end if
    do k = kln,klx
      ka = min(wl, k)
      kb = k - kln + 1
      old = suvb(:,ka,2)
      call mloexport("v",old,k,0)
      old = old + diffv_l(:,kb)*nudgewgt
      call mloimport("v",old,k,0)
      if ( allocated(oldv1) ) then
        oldv1(:,k) = oldv1(:,k) + diffv_l(:,kb)*nudgewgt
        oldv2(:,k) = oldv2(:,k) + diffv_l(:,kb)*nudgewgt
      end if
    end do
    if ( klx==kc ) then
      do k = kc+1,kbotmlo
        old = suvb(:,ka,2)
        call mloexport("v",old,k,0)
        old = old + diffv_l(:,kb)*nudgewgt
        call mloimport("v",old,k,0)
        if ( allocated(oldv1) ) then
          oldv1(:,k) = oldv1(:,k) + diffv_l(:,kb)*nudgewgt
          oldv2(:,k) = oldv2(:,k) + diffv_l(:,kb)*nudgewgt
        end if
      end do
    end if
  end if
      
end do
     
if ( nud_sfh/=0 ) then
  old = sfh
  call mloexport("eta",old,0,0)
  old = old + diffh_l(:,1)*nudgewgt
  call mloimport("eta",old,0,0)
end if

return
end subroutine mlofilterhub
      
!---------------------------------------------------------------------------------
! 2D Filter for MLO 
subroutine mlofilter(diff_l,diffs_l,diffu_l,diffv_l,diffh_l,lblock,kd)

use cc_mpi                  ! CC MPI routines
use newmpar_m               ! Grid parameters
use parm_m                  ! Model configuration
use vecsuv_m                ! Map to cartesian coordinates

integer, intent(in) :: kd
integer k
real, dimension(ifull,1), intent(inout) :: diffh_l
real, dimension(ifull,kd), intent(inout) :: diff_l, diffs_l
real, dimension(ifull,kd), intent(inout) :: diffu_l, diffv_l
real, dimension(ifull,kd) :: diffw_l
real, dimension(ifull_g,kd) :: diff_g ! large common array
real, dimension(ifull) :: xa_l, xb_l
logical, intent(in) :: lblock

if (nud_sst/=0) then
  call ccmpi_gatherall(diff_l(:,1:kd),diff_g(:,1:kd))
  call mlofilterhost(diff_g,diff_l,kd)
end if

if (nud_sss/=0) then
  call ccmpi_gatherall(diffs_l(:,1:kd),diff_g(:,1:kd))
  call mlofilterhost(diff_g,diffs_l,kd)
end if

if (nud_ouv/=0) then
  do k = 1,kd
    xa_l = diffu_l(:,k)
    xb_l = diffv_l(:,k)
    diffu_l(:,k) = ax(1:ifull)*xa_l + bx(1:ifull)*xb_l
    diffv_l(:,k) = ay(1:ifull)*xa_l + by(1:ifull)*xb_l
    diffw_l(:,k) = az(1:ifull)*xa_l + bz(1:ifull)*xb_l
  end do        
  call ccmpi_gatherall(diffu_l(:,1:kd),diff_g(:,1:kd))
  call mlofilterhost(diff_g,diffu_l,kd)
  call ccmpi_gatherall(diffv_l(:,1:kd),diff_g(:,1:kd))
  call mlofilterhost(diff_g,diffv_l,kd)
  call ccmpi_gatherall(diffw_l(:,1:kd),diff_g(:,1:kd))
  call mlofilterhost(diff_g,diffw_l,kd)
  do k = 1,kd
    xa_l = ax(1:ifull)*diffu_l(:,k) + ay(1:ifull)*diffv_l(:,k) + az(1:ifull)*diffw_l(:,k)
    xb_l = bx(1:ifull)*diffu_l(:,k) + by(1:ifull)*diffv_l(:,k) + bz(1:ifull)*diffw_l(:,k)
    diffu_l(:,k) = xa_l
    diffv_l(:,k) = xb_l
  end do
end if

if (nud_sfh/=0.and.lblock) then
  call ccmpi_gatherall(diffh_l(:,1),diff_g(:,1))
  call mlofilterhost(diff_g,diffh_l,1)
end if
      
return
end subroutine mlofilter

subroutine mlofilterhost(diff_g,dd,kd)

use cc_mpi             ! CC MPI routines
use const_phys         ! Physical constants
use map_m              ! Grid map arrays
use newmpar_m          ! Grid parameters
use parm_m             ! Model configuration
use parmgeom_m         ! Coordinate data
use vecsuv_m           ! Map to cartesian coordinates
use xyzinfo_m          ! Grid coordinate arrays

integer, intent(in) :: kd
integer iqq, iqqg, k, n, j, i, kdp1
real, dimension(ifull_g,kd), intent(in) :: diff_g ! large common array
real, dimension(ifull,kd), intent(out) :: dd
real, dimension(ifull_g) :: sm
real, dimension(ifull_g) :: xa, ya, za
real, dimension(kd+1) :: local_sum
real cq
real, dimension(ifull_g,kd+1) :: diff_g_l ! large common array
#ifdef GPU
real, dimension(ifull,kd+1) :: dd_l
#endif

! eventually will be replaced with mbd once full ocean coupling is complete
cq = sqrt(4.5)*.1*real(mbd_mlo)/(pi*schmidt)

call START_LOG(nestcalc_begin)

kdp1 = kd + 1
dd(:,:) = 0.
sm(:) = 1./em_g(:)**2
xa(:) = real(x_g(:))
ya(:) = real(y_g(:))
za(:) = real(z_g(:))
do k = 1,kd
  diff_g_l(:,k) = diff_g(:,k)*sm(:)
end do  
diff_g_l(:,kd+1) = sm(:)

#ifdef GPU

! GPU version
!$acc parallel loop collapse(2) copyin(cq,xa,ya,za,diff_g_l) &
!$acc   copyout(dd_l) private(iqq,iqqg,n,j,i)
do k = 1,kd+1
  do iqq = 1,ifull
    n = 1 + (iqq-1)/(ipan*jpan)  ! In range 1 .. npan
    j = 1 + ( iqq - (n-1)*(ipan*jpan) - 1) / ipan
    i = iqq - (j-1)*ipan - (n-1)*ipan*jpan
    iqqg = i+ioff + (j+joff-1)*il_g + (n-noff)*il_g**2
    dd_l(iqq,k) = drpdr_fast(iqqg,cq,xa,ya,za,diff_g_l(:,k))
  end do  
end do
!$acc end parallel loop
do k = 1,kd
  do iqq = 1,ifull
    dd(iqq,k) = dd_l(iqq,k)/max(dd_l(iqq,kd+1),1.e-8)
  end do
end do

#else

! CPU version
!$omp parallel do schedule(static) private(iqqg,iqq,local_sum,i,j,n)
do iqq = 1,ifull
  n = 1 + (iqq-1)/(ipan*jpan)  ! In range 1 .. npan
  j = 1 + ( iqq - (n-1)*(ipan*jpan) - 1) / ipan
  i = iqq - (j-1)*ipan - (n-1)*ipan*jpan
  iqqg = i+ioff + (j+joff-1)*il_g + (n-noff)*il_g**2
  local_sum(1:kd+1) = drpdr_fast(iqqg,cq,xa,ya,za,diff_g_l)
  dd(iqq,1:kd) = local_sum(1:kd)/max(local_sum(kd+1),1.e-8)
end do
!$omp end parallel do

#endif

call END_LOG(nestcalc_end)

return
end subroutine mlofilterhost

! 1D filer for ocean
subroutine mlofilterfast(diff_l,diffs_l,diffu_l,diffv_l,diffh_l,lblock,kd)

use cc_mpi                  ! CC MPI routines
use const_phys              ! Physical constants
use newmpar_m               ! Grid parameters
use parm_m                  ! Model configuration
use parmgeom_m              ! Coordinate data

integer, intent(in) :: kd
real, dimension(ifull,1), intent(inout) :: diffh_l
real, dimension(ifull,kd), intent(inout) :: diff_l, diffs_l
real, dimension(ifull,kd), intent(inout) :: diffu_l, diffv_l
real cq
logical, intent(in) :: lblock
      
! eventually will be replaced with mbd once full ocean coupling is complete
cq = sqrt(4.5)*.1*real(mbd_mlo)/(pi*schmidt)
      
if ( pprocn==pprocx ) then
  call mlospechost_n(cq,diff_l,diffs_l,diffu_l,diffv_l,diffh_l,lblock,kd)
else
  call mlospechost(cq,diff_l,diffs_l,diffu_l,diffv_l,diffh_l,lblock,kd)
end if

return
end subroutine mlofilterfast

subroutine mlospechost(cq,diff_l,diffs_l,diffu_l,diffv_l,diffh_l,lblock,kd)

use cc_mpi             ! CC MPI routines
use newmpar_m          ! Grid parameters
use parm_m             ! Model configuration
use vecsuv_m           ! Map to cartesian coordinates
      
integer, intent(in) :: kd
integer k,ppass,ibase
real, intent(in) :: cq
real, dimension(ifull,1), intent(inout) :: diffh_l
real, dimension(ifull,kd), intent(inout) :: diff_l,diffs_l
real, dimension(ifull,kd), intent(inout) :: diffu_l,diffv_l
real, dimension(ifull,kd) :: diffw_l
real, dimension(ipan*jpan,kd) :: qp
real, dimension(ifull) :: xa_l,xb_l
logical, intent(in) :: lblock

if (nud_sfh/=0.and.lblock) then
  call START_LOG(nestwin_begin)
  call ccmpi_gathermap_send2(diffh_l(:,1))    ! gather data onto global sparse array (1)
  call ccmpi_gathermap_recv2(1)               ! gather data onto global sparse array (1)
  call END_LOG(nestwin_end)
end if  

if (nud_sst/=0) then
  call START_LOG(nestwin_begin)
  call ccmpi_gathermap_send3(diff_l(:,:))     ! gather data onto global sparse array (1)
  call END_LOG(nestwin_end)
end if

if (nud_sfh/=0.and.lblock) then
  do ppass=pprocn,pprocx
    call copyglobalpack(1,0,1)                ! copy sparse array (1) to (0)
    call mlofastspec_work(cq,qp,1,ppass)      ! filter sparse array (0)
    ibase = ipan*jpan*(ppass-pprocn)
    diffh_l(1+ibase:ipan*jpan+ibase,1) = qp(1:ipan*jpan,1)
  end do
end if

if (nud_sst/=0) then
  call START_LOG(nestwin_begin)
  call ccmpi_gathermap_recv3(kd,kd)           ! gather data onto global sparse array (1)
  call END_LOG(nestwin_end)
end if  

if (nud_sss/=0) then
  call START_LOG(nestwin_begin)  
  call ccmpi_gathermap_send3(diffs_l(:,:))    ! gather data onto global sparse array (1)
  call END_LOG(nestwin_end)
end if

if (nud_sst/=0) then
  do ppass=pprocn,pprocx
    call copyglobalpack(kd,0,kd)              ! copy sparse array (1) to (0)
    call mlofastspec_work(cq,qp,kd,ppass)     ! filter sparse array (0)
    ibase = ipan*jpan*(ppass-pprocn)
    diff_l(1+ibase:ipan*jpan+ibase,1:kd) = qp(1:ipan*jpan,1:kd)
  end do
end if

if (nud_sss/=0) then
  call START_LOG(nestwin_begin)  
  call ccmpi_gathermap_recv3(kd,kd)           ! gather data onto global sparse array (1)
  call END_LOG(nestwin_end)
end if

if (nud_ouv/=0) then
  do k = 1,kd
    xa_l = diffu_l(:,k)
    xb_l = diffv_l(:,k)
    diffu_l(:,k) = ax(1:ifull)*xa_l + bx(1:ifull)*xb_l
    diffv_l(:,k) = ay(1:ifull)*xa_l + by(1:ifull)*xb_l
    diffw_l(:,k) = az(1:ifull)*xa_l + bz(1:ifull)*xb_l
  end do
  call START_LOG(nestwin_begin)
  call ccmpi_gathermap_send3(diffu_l(:,:))    ! gather data onto global sparse array (1)
  call END_LOG(nestwin_end)
end if

if (nud_sss/=0) then
  do ppass = pprocn,pprocx
    call copyglobalpack(kd,0,kd)              ! copy sparse array (1) to (0)
    call mlofastspec_work(cq,qp,kd,ppass)     ! filter sparse array (0)
    ibase = ipan*jpan*(ppass-pprocn)
    diffs_l(1+ibase:ipan*jpan+ibase,1:kd) = qp(1:ipan*jpan,1:kd)
  end do
end if

if (nud_ouv/=0) then
  call START_LOG(nestwin_begin)
  call ccmpi_gathermap_recv3(kd,kd)           ! gather data onto global sparse array (1)
  call ccmpi_gathermap_send3(diffv_l(:,:))    ! gather data onto global sparse array (1)
  call END_LOG(nestwin_end)
  do ppass = pprocn,pprocx
    call copyglobalpack(kd,0,kd)              ! copy sparse array (1) to (0)
    call mlofastspec_work(cq,qp,kd,ppass)     ! filter sparse array (0)
    ibase = ipan*jpan*(ppass-pprocn)
    diffu_l(1+ibase:ipan*jpan+ibase,1:kd) = qp(1:ipan*jpan,1:kd)
  end do
  call START_LOG(nestwin_begin)
  call ccmpi_gathermap_recv3(kd,kd)           ! gather data onto global sparse array (1)
  call ccmpi_gathermap_send3(diffw_l(:,:))    ! gather data onto global sparse array (1)
  call END_LOG(nestwin_end)
  do ppass = pprocn,pprocx
    call copyglobalpack(kd,0,kd)              ! copy sparse array (1) to (0)
    call mlofastspec_work(cq,qp,kd,ppass)     ! filter sparse array (0)
    ibase = ipan*jpan*(ppass-pprocn)
    diffv_l(1+ibase:ipan*jpan+ibase,1:kd) = qp(1:ipan*jpan,1:kd)
  end do
  call START_LOG(nestwin_begin)
  call ccmpi_gathermap_recv3(kd,kd)            ! gather data onto global sparse array (1)
  call END_LOG(nestwin_end)
  do ppass = pprocn,pprocx
    call copyglobalpack(kd,0,kd)               ! copy sparse array (1) to (0)
    call mlofastspec_work(cq,qp,kd,ppass)      ! filter sparse array (0)
    ibase = ipan*jpan*(ppass-pprocn)
    diffw_l(1+ibase:ipan*jpan+ibase,1:kd) = qp(1:ipan*jpan,1:kd)
  end do
  do k = 1,kd
    xa_l = ax(1:ifull)*diffu_l(:,k) + ay(1:ifull)*diffv_l(:,k) + az(1:ifull)*diffw_l(:,k)
    xb_l = bx(1:ifull)*diffu_l(:,k) + by(1:ifull)*diffv_l(:,k) + bz(1:ifull)*diffw_l(:,k)
    diffu_l(:,k) = xa_l
    diffv_l(:,k) = xb_l
  end do
end if

return
end subroutine mlospechost

!---------------------------------------------------------------------------------
! memory reduced version of mlospechost
subroutine mlospechost_n(cq,diff_l,diffs_l,diffu_l,diffv_l,diffh_l,lblock,kd)

use cc_mpi             ! CC MPI routines
use newmpar_m          ! Grid parameters
use parm_m             ! Model configuration
use vecsuv_m           ! Map to cartesian coordinates
      
integer, intent(in) :: kd
integer k
real, intent(in) :: cq
real, dimension(ifull,1), intent(inout) :: diffh_l
real, dimension(ifull,kd), intent(inout) :: diff_l,diffs_l
real, dimension(ifull,kd), intent(inout) :: diffu_l,diffv_l
real, dimension(ifull,kd) :: diffw_l
real, dimension(ifull) :: xa_l,xb_l
logical, intent(in) :: lblock

if (nud_sfh/=0.and.lblock) then
  call START_LOG(nestwin_begin)  
  call ccmpi_gathermap_send2(diffh_l(:,1)) ! gather data onto global sparse array (0)
  call ccmpi_gathermap_recv2(0) ! gather data onto global sparse array (0)
  call END_LOG(nestwin_end)
end if

if (nud_sst/=0) then
  call START_LOG(nestwin_begin)  
  call ccmpi_gathermap_send3(diff_l(:,:)) ! gather data onto global sparse array (0)
  call END_LOG(nestwin_end)
end if

if (nud_sfh/=0.and.lblock) then
  call mlofastspec_work(cq,diffh_l(:,1),1,pprocn) ! filter sparse array (0)
end if

if (nud_sst/=0) then
  call START_LOG(nestwin_begin)  
  call ccmpi_gathermap_recv3(kd,0) ! gather data onto global sparse array (0)
  call END_LOG(nestwin_end)
end if

if (nud_sss/=0) then
  call START_LOG(nestwin_begin)  
  call ccmpi_gathermap_send3(diffs_l(:,:)) ! gather data onto global sparse array (0)
  call END_LOG(nestwin_end)
end if

if (nud_sst/=0) then
  call mlofastspec_work(cq,diff_l,kd,pprocn)      ! filter sparse array (0)
end if

if (nud_sss/=0) then
  call START_LOG(nestwin_begin)  
  call ccmpi_gathermap_recv3(kd,0)             ! gather data onto global sparse array (0)
  call END_LOG(nestwin_end)
end if

if (nud_ouv/=0) then
  do k = 1,kd
    xa_l = diffu_l(:,k)
    xb_l = diffv_l(:,k)
    diffu_l(:,k) = ax(1:ifull)*xa_l + bx(1:ifull)*xb_l
    diffv_l(:,k) = ay(1:ifull)*xa_l + by(1:ifull)*xb_l
    diffw_l(:,k) = az(1:ifull)*xa_l + bz(1:ifull)*xb_l
  end do
  call START_LOG(nestwin_begin)
  call ccmpi_gathermap_send3(diffu_l(:,:)) ! gather data onto global sparse array (0)
  call END_LOG(nestwin_end)
end if

if (nud_sss/=0) then
  call mlofastspec_work(cq,diffs_l,kd,pprocn)      ! filter sparse array (0)
end if

if (nud_ouv/=0) then
  call START_LOG(nestwin_begin)
  call ccmpi_gathermap_recv3(kd,0) ! gather data onto global sparse array (0)
  call ccmpi_gathermap_send3(diffv_l(:,:)) ! gather data onto global sparse array (0)
  call END_LOG(nestwin_end)
  call mlofastspec_work(cq,diffu_l,kd,pprocn)      ! filter sparse array (0)
  call START_LOG(nestwin_begin)
  call ccmpi_gathermap_recv3(kd,0)! gather data onto global sparse array (0)
  call ccmpi_gathermap_send3(diffw_l(:,:)) ! gather data onto global sparse array (0)
  call END_LOG(nestwin_end)
  call mlofastspec_work(cq,diffv_l,kd,pprocn)      ! filter sparse array (0)
  call START_LOG(nestwin_begin)
  call ccmpi_gathermap_recv3(kd,0) ! gather data onto global sparse array (0)
  call END_LOG(nestwin_end)
  call mlofastspec_work(cq,diffw_l,kd,pprocn)      ! filter sparse array (0)
  do k = 1,kd
    xa_l = ax(1:ifull)*diffu_l(:,k) + ay(1:ifull)*diffv_l(:,k) + az(1:ifull)*diffw_l(:,k)
    xb_l = bx(1:ifull)*diffu_l(:,k) + by(1:ifull)*diffv_l(:,k) + bz(1:ifull)*diffw_l(:,k)
    diffu_l(:,k) = xa_l
    diffv_l(:,k) = xb_l
  end do
endif

return
end subroutine mlospechost_n

subroutine mlofastspec_work(cq,diff_g,kd,ppass)

use cc_mpi
use newmpar_m

integer, intent(in) :: kd,ppass
real, intent(in) :: cq
real, dimension(ipan*jpan,kd), intent(inout) :: diff_g

! computations for the local processor group
select case(ppass)
  case(1,2,3)
    call mlospeclocal_left(cq,ppass,diff_g,kd)
  case(0,4,5)
    call mlospeclocal_right(cq,ppass,diff_g,kd)
end select

return
end subroutine mlofastspec_work
      
!---------------------------------------------------------------------------------
! This version is for asymmetric decomposition
subroutine mlospeclocal_left(cq,ppass,qp,kd)

use cc_acc             ! CC OpenACC routines
use cc_mpi             ! CC MPI routines
use map_m              ! Grid map arrays
use newmpar_m          ! Grid parameters
use parm_m             ! Model configuration
use xyzinfo_m          ! Grid coordinate arrays
     
integer, intent(in) :: ppass, kd
integer j, n, ipass, ns, ne, os, oe
integer jpoff, ibase
integer me, k, til, sn, sy, a, b, c, jj, nn
integer ibeg, iend, kdp1
integer, dimension(0:3) :: astr, bstr, cstr
integer, dimension(0:3) :: maps
real, intent(in) :: cq
real, dimension(ipan*jpan,kd), intent(out) :: qp
real, dimension(il_g) :: ap, asum      
real, dimension(il_g*ipan*(kd+1)*3) :: zz
real, dimension(ipan*jpan*(kd+1),0:2) :: yy
real, dimension(kd+1) :: local_sum
real, dimension(4*il_g,max(ipan,jpan)) :: xa, ya, za
real, dimension(4*il_g,kd+1,max(ipan,jpan)) :: ap_l
#ifdef GPU
integer async_counter
real, dimension(jpan,kd+1,ipan) :: qp_l
#endif
      
maps = (/ il_g, il_g, 4*il_g, 3*il_g /)
til = il_g**2
kdp1 = kd + 1
astr = 0
bstr = 0
cstr = 0

ns = joff + 1
ne = joff + jpan
os = ioff + 1
oe = ioff + ipan
 
call START_LOG(nestcalc_begin)

do ipass = 0,2
  me = maps(ipass)
  call getiqa(astr,bstr,cstr,me,ipass,ppass,il_g)

#ifndef GPU
  !$omp parallel do schedule(static) private(j,jj,sn,sy,a,b,c,ibeg,iend,asum,k,ap)
#endif
  do j = 1,jpan      
    ! pack data from sparse arrays
    jj = j + ns - 1
    do sn = 1,me,il_g
      sy = (sn-1)/il_g
      a = astr(sy)
      b = bstr(sy)
      c = cstr(sy)
      ibeg = a*sn + b*jj + c
      iend = a*(sn+il_g-1) + b*jj + c
      xa(sn:sn+il_g-1,j) = real(x_g(ibeg:iend:a))
      ya(sn:sn+il_g-1,j) = real(y_g(ibeg:iend:a))
      za(sn:sn+il_g-1,j) = real(z_g(ibeg:iend:a))
      asum = 1./em_g(ibeg:iend:a)**2
      do k = 1,kd
        ! v version is faster for getglobalpack  
        call getglobalpack_v(ap,ibeg,iend,k) 
        ap_l(sn:sn+il_g-1,k,j) = ap*asum
      end do
      ap_l(sn:sn+il_g-1,kd+1,j) = asum      
    end do
  end do
#ifndef GPU
  !$omp end parallel do
#endif
    
  ! start convolution

#ifdef GPU

  async_counter = mod(ipass,async_length)
  !$acc parallel loop collapse(3) copyin(xa(1:me,1:jpan),ya(1:me,1:jpan),za(1:me,1:jpan),ap_l(1:me,:,1:jpan)) &
  !$acc   copyout(yy(:,ipass)) private(j,n,nn,k) async(async_counter)
  do j = 1,jpan
    do k = 1,kd+1  
      do n = 1,ipan
        nn = n + os - 1
        yy(n+(j-1)*ipan+(k-1)*ipan*jpan,ipass) = drpdr_fast(nn,cq,xa(1:me,j),ya(1:me,j),za(1:me,j),ap_l(1:me,k,j))
      end do  
    end do
  end do
  !$acc end parallel loop

#else

  !$omp parallel do collapse(2) schedule(static) private(j,n,nn,k,local_sum)
  do j = 1,jpan
    do n = 1,ipan
      nn = n + os - 1
      local_sum(1:kd+1) = drpdr_fast(nn,cq,xa(1:me,j),ya(1:me,j),za(1:me,j),ap_l(:,:,j))
      do k = 1,kd+1
        yy(n+(j-1)*ipan+(k-1)*ipan*jpan,ipass) = local_sum(k)
      end do  
    end do
  end do
  !$omp end parallel do

#endif

end do ! ipass

#ifdef GPU
!$acc wait
#endif

call END_LOG(nestcalc_end)

call START_LOG(nestwin_begin)
call ccmpi_gathermap_wait
call END_LOG(nestwin_end)

call START_LOG(nestcomm_begin)

! gather data on host processors
call ccmpi_allgatherx(zz(1:il_g*ipan*(kd+1)*3),yy(1:ipan*jpan*(kd+1),0:2),comm_cols)

do ipass = 0,2
  me = maps(ipass)
  call getiqa(astr,bstr,cstr,me,ipass,ppass,il_g)
  
  ! unpack grid
  a = astr(0)
  b = bstr(0)
  c = cstr(0)

  ! unpack data to sparse arrays
  do n = 1,ipan
    nn = n + os - 1
    do k = 1,kd
      do jpoff = 0,il_g-1,jpan
        sy = jpoff/jpan
        ibase = n + ipan*jpan*(k-1) + ipan*jpan*(kd+1)*ipass + ipan*jpan*(kd+1)*3*sy
        ap(1+jpoff:jpan+jpoff) = zz(ibase:ibase+ipan*(jpan-1):ipan)
      end do
      ibeg = a*nn + b*1 + c
      iend = a*nn + b*il_g + c
      call setglobalpack_v(ap(1:il_g),ibeg,iend,k)
    end do  
    do jpoff = 0,il_g-1,jpan
      sy = jpoff/jpan
      ibase = n + ipan*jpan*kd + ipan*jpan*(kd+1)*ipass + ipan*jpan*(kd+1)*3*sy
      ap(1+jpoff:jpan+jpoff) = zz(ibase:ibase+ipan*(jpan-1):ipan)
    end do
    ibeg = a*nn + b*1 + c
    iend = a*nn + b*il_g + c
    call setglobalpack_v(ap(1:il_g),ibeg,iend,0)
  end do  
          
end do ! ipass

call END_LOG(nestcomm_end)
    
ns = ioff + 1
ne = ioff + ipan
os = joff + 1
oe = joff + jpan

ipass = 3
me = maps(ipass)
call getiqa(astr,bstr,cstr,me,ipass,ppass,il_g)

call START_LOG(nestcalc_begin)

#ifndef GPU
!$omp parallel do schedule(static) private(j,jj,sn,sy,a,b,c,ibeg,iend,asum,k,ap)
#endif
do j = 1,ipan    
  ! pack data from sparse arrays
  jj = j + ns - 1
  do sn = 1,me,il_g
    sy = (sn-1)/il_g
    a = astr(sy)
    b = bstr(sy)
    c = cstr(sy)
    ibeg = a*sn + b*jj + c
    iend = a*(sn+il_g-1) + b*jj + c
    xa(sn:sn+il_g-1,j) = real(x_g(ibeg:iend:a))
    ya(sn:sn+il_g-1,j) = real(y_g(ibeg:iend:a))
    za(sn:sn+il_g-1,j) = real(z_g(ibeg:iend:a))
    ! v version is faster for getglobalpack  
    call getglobalpack_v(asum,ibeg,iend,0) 
    do k = 1,kd
      call getglobalpack_v(ap,ibeg,iend,k)
      ap_l(sn:sn+il_g-1,k,j) = ap
    end do
    ap_l(sn:sn+il_g-1,kd+1,j) = asum    
  end do
end do
#ifndef GPU
!$omp end parallel do
#endif

! start convolution
#ifdef GPU

!$acc parallel loop collapse(3) copyin(xa(1:me,1:ipan),ya(1:me,1:ipan),za(1:me,1:ipan),ap_l(1:me,:,1:ipan)) &
!$acc   copyout(qp_l) private(j,n,nn)
do j = 1,ipan
  do k = 1,kd+1  
    do n = 1,jpan
      nn = n + os - 1
      qp_l(n,k,j) = drpdr_fast(nn,cq,xa(1:me,j),ya(1:me,j),za(1:me,j),ap_l(1:me,k,j))
    end do  
  end do
end do
!$acc end parallel loop
do j = 1,ipan
  do k = 1,kd
    do n = 1,jpan
      qp(j+ipan*(n-1),k) = qp_l(n,k,j)/max(qp_l(n,kd+1,j),1.e-8)
    end do
  end do
end do

#else

!$omp parallel do collapse(2) schedule(static) private(j,n,nn,local_sum)
do j = 1,ipan
  do n = 1,jpan
    nn = n + os - 1
    local_sum = drpdr_fast(nn,cq,xa(1:me,j),ya(1:me,j),za(1:me,j),ap_l(:,:,j))
    qp(j+ipan*(n-1),1:kd) = local_sum(1:kd)/max(local_sum(kd+1),1.e-8)
  end do
end do
!$omp end parallel do

#endif

call END_LOG(nestcalc_end)
      
return  
end subroutine mlospeclocal_left

subroutine mlospeclocal_right(cq,ppass,qp,kd)

use cc_acc             ! CC OpenACC routines
use cc_mpi             ! CC MPI routines
use map_m              ! Grid map arrays
use newmpar_m          ! Grid parameters
use parm_m             ! Model configuration
use xyzinfo_m          ! Grid coordinate arrays
     
integer, intent(in) :: ppass, kd
integer j, n, ipass, ns, ne, os, oe
integer jpoff, ibase
integer me, k, til, sn, sy, a, b, c, jj, nn
integer ibeg, iend, kdp1
integer, dimension(0:3) :: astr, bstr, cstr
integer, dimension(0:3) :: maps
real, intent(in) :: cq
real, dimension(ipan*jpan,kd), intent(out) :: qp
real, dimension(il_g) :: ap, asum      
real, dimension(il_g*jpan*(kd+1)*3) :: zz
real, dimension(ipan*jpan*(kd+1),0:2) :: yy
real, dimension(kd+1) :: local_sum
real, dimension(4*il_g,max(ipan,jpan)) :: xa, ya, za
real, dimension(4*il_g,kd+1,max(ipan,jpan)) :: ap_l
#ifdef GPU
integer async_counter
real, dimension(ipan,kd+1,jpan) :: qp_l
#endif
      
maps = (/ il_g, il_g, 4*il_g, 3*il_g /)
til = il_g**2
kdp1 = kd + 1
astr = 0
bstr = 0
cstr = 0

ns = ioff + 1
ne = ioff + ipan
os = joff + 1
oe = joff + jpan
     
call START_LOG(nestcalc_begin)

do ipass = 0,2
  me = maps(ipass)
  call getiqa(astr,bstr,cstr,me,ipass,ppass,il_g)

#ifndef GPU
  !$omp parallel do schedule(static) private(j,jj,sn,sy,a,b,c,ibeg,iend,asum,k,ap)
#endif
  do j = 1,ipan      
    ! pack data from sparse arrays
    jj = j + ns - 1
    do sn = 1,me,il_g
      sy = (sn-1)/il_g
      a = astr(sy)
      b = bstr(sy)
      c = cstr(sy)
      ibeg = a*sn + b*jj + c
      iend = a*(sn+il_g-1) + b*jj + c
      xa(sn:sn+il_g-1,j) = real(x_g(ibeg:iend:a))
      ya(sn:sn+il_g-1,j) = real(y_g(ibeg:iend:a))
      za(sn:sn+il_g-1,j) = real(z_g(ibeg:iend:a))
      asum = 1./em_g(ibeg:iend:a)**2
      do k = 1,kd
        call getglobalpack_v(ap,ibeg,iend,k)           
        ap_l(sn:sn+il_g-1,k,j) = ap*asum
      end do
      ap_l(sn:sn+il_g-1,kd+1,j) = asum      
    end do
  end do
#ifndef GPU
  !$omp end parallel do
#endif
    
  ! start convolution

#ifdef GPU

! GPU version
  async_counter = mod(ipass,async_length)
  !$acc parallel loop collapse(3) copyin(xa(1:me,1:ipan),ya(1:me,1:ipan),za(1:me,1:ipan),ap_l(1:me,:,1:ipan)) &
  !$acc   copyout(yy(:,ipass)) private(j,n,nn,k) async(async_counter)
  do j = 1,ipan
    do k = 1,kd+1  
      do n = 1,jpan
        nn = n + os - 1
        yy(n+(j-1)*jpan+(k-1)*ipan*jpan,ipass) = drpdr_fast(nn,cq,xa(1:me,j),ya(1:me,j),za(1:me,j),ap_l(1:me,k,j))
      end do  
    end do    
  end do
  !$acc end parallel loop

#else

  ! CPU version
  !$omp parallel do collapse(2) schedule(static) private(j,n,nn,k,local_sum)
  do j = 1,ipan
    do n = 1,jpan
      nn = n + os - 1
      local_sum(1:kd+1) = &
        drpdr_fast(nn,cq,xa(1:me,j),ya(1:me,j),za(1:me,j),ap_l(:,:,j))
      do k = 1,kd+1
        yy(n+(j-1)*jpan+(k-1)*ipan*jpan,ipass) = local_sum(k)
      end do  
    end do    
  end do
  !$omp end parallel do

#endif

end do ! ipass

#ifdef GPU
!$acc wait
#endif

call END_LOG(nestcalc_end)

call START_LOG(nestwin_begin)
call ccmpi_gathermap_wait
call END_LOG(nestwin_end)

call START_LOG(nestcomm_begin)

! gather data on host processors
call ccmpi_allgatherx(zz(1:il_g*jpan*(kd+1)*3),yy(1:ipan*jpan*(kd+1),0:2),comm_rows)

do ipass = 0,2
  me = maps(ipass)
  call getiqa(astr,bstr,cstr,me,ipass,ppass,il_g)

  ! unpack grid
  a = astr(0)
  b = bstr(0)
  c = cstr(0)
  
  ! unpack data to sparse arrays
  do n = 1,jpan
    nn = n + os - 1
    do k = 1,kd
      do jpoff = 0,il_g-1,ipan
        sy = jpoff/ipan
        ibase = n + ipan*jpan*(k-1) + ipan*jpan*(kd+1)*ipass + ipan*jpan*(kd+1)*3*sy
        ap(1+jpoff:ipan+jpoff) = zz(ibase:ibase+jpan*(ipan-1):jpan)
      end do
      ibeg = a*nn + b*1 + c
      iend = a*nn + b*il_g + c
      call setglobalpack_v(ap(1:il_g),ibeg,iend,k)
    end do  
    do jpoff = 0,il_g-1,ipan
      sy = jpoff/ipan
      ibase = n + ipan*jpan*kd + ipan*jpan*(kd+1)*ipass + ipan*jpan*(kd+1)*3*sy
      ap(1+jpoff:ipan+jpoff) = zz(ibase:ibase+jpan*(ipan-1):jpan)
    end do
    ibeg = a*nn + b*1 + c
    iend = a*nn + b*il_g + c
    call setglobalpack_v(ap(1:il_g),ibeg,iend,0)
  end do  
          
end do ! ipass

call END_LOG(nestcomm_end)
    
ns = joff + 1
ne = joff + jpan
os = ioff + 1
oe = ioff + ipan

ipass = 3
me = maps(ipass)
call getiqa(astr,bstr,cstr,me,ipass,ppass,il_g)

call START_LOG(nestcalc_begin)

#ifndef GPU
!$omp parallel do schedule(static) private(j,jj,sn,sy,a,b,c,ibeg,iend,asum,k,ap)
#endif
do j = 1,jpan    
  ! pack data from sparse arrays
  jj = j + ns - 1
  do sn = 1,me,il_g
    sy = (sn-1)/il_g
    a = astr(sy)
    b = bstr(sy)
    c = cstr(sy)
    ibeg = a*sn + b*jj + c
    iend = a*(sn+il_g-1) + b*jj + c
    xa(sn:sn+il_g-1,j) = real(x_g(ibeg:iend:a))
    ya(sn:sn+il_g-1,j) = real(y_g(ibeg:iend:a))
    za(sn:sn+il_g-1,j) = real(z_g(ibeg:iend:a))
    ! v version is faster for getglobalpack  
    call getglobalpack_v(asum,ibeg,iend,0)
    do k = 1,kd
      call getglobalpack_v(ap,ibeg,iend,k)
      ap_l(sn:sn+il_g-1,k,j) = ap
    end do
    ap_l(sn:sn+il_g-1,kd+1,j) = asum    
  end do
end do
#ifndef GPU
!$omp end parallel do
#endif
  
! start convolution

#ifdef GPU

! GPU version
!$acc parallel loop collapse(3) copyin(xa(1:me,1:jpan),ya(1:me,1:jpan),za(1:me,1:jpan),ap_l(1:me,:,1:jpan)) &
!$acc   copyout(qp_l) private(j,n,nn)
do j = 1,jpan
  do k = 1,kd+1  
    do n = 1,ipan
      nn = n + os - 1
      qp_l(n,k,j) = drpdr_fast(nn,cq,xa(1:me,j),ya(1:me,j),za(1:me,j),ap_l(1:me,k,j))
    end do  
  end do  
end do
!$acc end parallel loop
do k = 1,kd
  do j = 1,jpan
    do n = 1,ipan
      qp(n+ipan*(j-1),k) = qp_l(n,k,j)/max(qp_l(n,kd+1,j),1.e-8)
    end do
  end do  
end do        

#else

! CPU version
!$omp parallel do collapse(2) schedule(static) private(j,n,nn,local_sum)
do j = 1,jpan
  do n = 1,ipan
    nn = n + os - 1
    local_sum = drpdr_fast(nn,cq,xa(1:me,j),ya(1:me,j),za(1:me,j),ap_l(:,:,j))
    qp(n+ipan*(j-1),1:kd) = local_sum(1:kd)/max(local_sum(kd+1), 1.e-8)  
  end do  
end do
!$omp end parallel do

#endif

call END_LOG(nestcalc_end)
      
return  
end subroutine mlospeclocal_right


! Relaxtion method for mlo
subroutine mlonudge(new,sssb,suvb,sfh,wl)

use mlo, only : mloimport,mloexport, &
                minsal,wlev,wrtemp,  &
                mloexpmelt               ! Ocean physics and prognostic arrays
use newmpar_m                            ! Grid parameters
use parm_m                               ! Model configuration
      
integer, intent(in) :: wl
integer k, ka
real, dimension(ifull), intent(in) :: sfh
real, dimension(ifull,wlev), intent(in) :: new,sssb
real, dimension(ifull,wlev,2), intent(in) :: suvb
real, dimension(ifull) :: old
real wgt
      
wgt=dt/real(nud_hrs*3600)

if (nud_sst/=0) then
  !call mloexpmelt(timelt)  
  do k=ktopmlo,kbotmlo
    ka=min(k,wl)
    old=new(:,ka)
    call mloexport("temp",old,k,0)
    !old=old*(1.-wgt)+max(new(:,ka),275.16-wrtemp)*wgt
    old=old*(1.-wgt)+max(new(:,ka),271.16-wrtemp)*wgt
    call mloimport("temp",old,k,0)
  end do
end if
      
if (nud_sss/=0) then
  do k=ktopmlo,kbotmlo
    ka=min(k,wl)
    old=sssb(:,ka)
    call mloexport("sal",old,k,0)
    where ( old>minsal .and. sssb(:,ka)>minsal )
      old=old*(1.-wgt)+sssb(:,ka)*wgt
      old=max(old,0.)
    end where  
    call mloimport("sal",old,k,0)
  end do
end if
      
if (nud_ouv/=0) then
  do k=ktopmlo,kbotmlo
    ka=min(k,wl)
    old=suvb(:,ka,1)
    call mloexport("u",old,k,0)
    old=old*(1.-wgt)+suvb(:,ka,1)*wgt
    call mloimport("u",old,k,0)
    ka=min(k,wl)
    old=suvb(:,ka,2)
    call mloexport("v",old,k,0)
    old=old*(1.-wgt)+suvb(:,ka,2)*wgt
    call mloimport("v",old,k,0) 
  end do
end if

if (nud_sfh/=0) then
  old=sfh
  call mloexport("eta",old,0,0)
  old=old*(1.-wgt)+sfh*wgt
  call mloimport("eta",old,0,0)
end if
      
return
end subroutine mlonudge

!--------------------------------------------------------------
! Initialise RMA windows used by 1D convolutions to retrieve
! data from other processor ranks
subroutine specinit
      
use cc_mpi
use newmpar_m
use parm_m

integer ncount,ipass,ppass,me
integer n,j,jj,sn,sy,kx,ky
integer iqg,ng,jg,ig
integer a,b,c,ns
integer iproc
integer, dimension(0:3) :: maps
integer, dimension(0:3) :: astr,bstr,cstr
logical, dimension(0:nproc-1) :: lproc, lproc_t
logical, save :: first = .true.
      
! can be called from nestinb or amipsst
if ( .not.first ) return
first = .false.

if ( myid==0 ) then
  write(6,*) "Create map for spectral nudging data"
end if

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

! Construct a map of processes that this process requires
if ( myid==0 .and. nmaxpr==1 ) then
  write(6,*) "-> Create map of processes required by this process"
end if
ncount = count(lproc)
allocate(specmap_req(ncount))
ncount = 0
do iproc = 0,nproc-1
  if ( lproc(iproc) ) then
    ncount = ncount + 1
    specmap_req(ncount) = iproc
  end if
end do

! Create create a list of processes receiving data (node aware)
! (use lproc_t for node communication as lproc will define memory
!  allocated to sparse arrays below)
if ( myid==0 .and. nmaxpr==1 ) then
  write(6,*) "-> Create map of processes received by this process"
end if
call ccmpi_allreduce(lproc,lproc_t,"or",comm_node) ! collect required data for node on node_myid==0
! redistribute node messages across cores on the node
allocate( specmap_indx(0:nproc-1) )
specmap_indx = 0
specmap_indxlen = 0
ncount = 0
do iproc = 0,nproc-1
  if ( lproc_t(iproc) ) then
    specmap_indxlen = specmap_indxlen + 1
    specmap_indx(iproc) = specmap_indxlen
    ncount = ncount + 1
    if ( mod(ncount,node_nproc)/=node_myid ) then
      lproc_t(iproc) = .false.
    end if
  end if
end do
! construct list of required messages for this process on behalf of the node
ncount = count(lproc_t)
allocate(specmap_recv(ncount))
ncount = 0
do iproc = 0,nproc-1
  if ( lproc_t(iproc) ) then
    ncount = ncount + 1
    specmap_recv(ncount) = iproc
  end if
end do

! Construct a map of processes that is required by this rank
if ( myid==0 .and. nmaxpr==1 ) then
  write(6,*) "-> Create map of processes sent by this process"
end if
call ccmpi_alltoall(lproc_t,comm_world) ! global transpose
ncount = count(lproc_t(0:nproc-1))
allocate( specmap_send(ncount) )
ncount = 0
do iproc = 0,nproc-1
  if ( lproc_t(iproc) ) then
    ncount = ncount + 1
    specmap_send(ncount) = iproc
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

if ( myid==0 .and. nmaxpr==1 ) then
  write(6,*) "-> Allocating buffers on each process"
end if
ncount = count(lproc)
allocate(specmap_ext(ncount)) ! for allocating internal buffers on each process
ncount = 0
do iproc = 0,nproc-1
  if ( lproc(iproc) ) then
    ncount = ncount + 1
    specmap_ext(ncount) = iproc
  end if
end do
   
ky = max( kl, ol )
if ( npta==1 ) then
  ! face version (nproc>=6)
  kx = min(ky,kblock)
else
  ! general version
  kx = 2*min(ky,kblock) ! extra memory for copy
end if
call allocateglobalpack(kx,ky)
deallocate(specmap_ext)
      
return
end subroutine specinit

subroutine setdavvertwgt

use daviesnudge                  ! Far-field nudging
use newmpar_m                    ! Grid parameters
use parm_m                       ! Model configuration
use sigs_m                       ! Atmosphere sigma levels

integer klow, khigh, k
real slow, shigh

if ( mbd==0 ) return

slow = sig(kbotdav) - sigramplow
shigh = sig(ktopdav) + sigramphigh

do klow = 1,kl-1
  if ( slow>=sig(klow) ) exit
end do
do khigh = kl,2,-1
  if ( shigh<=sig(khigh) ) exit
end do

do k = 1,kbotdav-1
  vertwgt(k) = 0.
end do
do k = kbotdav,klow-1
  vertwgt(k) = (sig(kbotdav)-sig(k))/(sig(kbotdav)-slow)
end do
do k = klow,khigh
  vertwgt(k) = 1.
end do
do k = khigh+1,ktopdav
  vertwgt(k) = (sig(k)-sig(ktopdav))/(shigh-sig(ktopdav))
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

pure function iabsdate(kdate_r,kdate) result(ans)

use parm_m

integer, intent(in) :: kdate_r, kdate
integer iyear_r,iyear0,imonth_r,iday_r
integer months,nl
integer newdate_r, diffyear
integer, dimension(-1:13) :: mdays
integer ans

iyear_r = kdate_r/10000
iyear0  = kdate/10000                ! year of kdate
newdate_r = kdate_r - 10000*iyear_r
imonth_r = newdate_r/100
newdate_r = newdate_r - 100*imonth_r
iday_r = newdate_r

! calculate number of months since start of kdate year
diffyear = iyear_r - iyear0
months = diffyear*12 + imonth_r - 1  

if ( leap==0 ) then ! 365 day calendar
  mdays = (/-31,0,31,59,90,120,151,181,212,243,273,304,334,365,396/)
else if ( leap==1 ) then ! 365/366 day calendar  
  mdays = (/-31,0,31,59,90,120,151,181,212,243,273,304,334,365,396/)
  nl = 0
  if ( mod(iyear0,4)==0 ) nl = 1
  if ( mod(iyear0,100)==0 ) nl = 0
  if ( mod(iyear0,400)==0 ) nl = 1
  mdays(2:13) = mdays(2:13) + nl
else if ( leap==2 ) then ! 360 day calendar
  mdays = (/-30,0,30,60,90,120,150,180,210,240,270,300,330,360,390/)  
end if

! Accumulate days month by month, up to last completed month
ans = mdays(months)

! Add days from this current month
ans = ans + iday_r

end function iabsdate


#ifdef GPU

! GPU version
pure function drpdr_fast(nn,cq,xa,ya,za,at) result(out_sum)
!$acc routine seq

integer, intent(in) :: nn
integer i, ilen
real, intent(in) :: cq
real, dimension(:), intent(in) :: xa, ya, za
real, dimension(:), intent(in) :: at
real ra
real out_sum, rb
real at_t, e, t1, t2
complex local_sum

ilen = size(xa)
local_sum = (0.,0.)

do i = 1,ilen
  ra = xa(nn)*xa(i) + ya(nn)*ya(i) + za(nn)*za(i)
  ra = acos(max(min(ra, 1.), -1.))
  rb = exp(-min((cq*ra)**2,40.))
  at_t = rb*at(i)
  t1 = at_t + real(local_sum)
  e  = t1 - at_t
  t2 = ((real(local_sum) - e) + (at_t - (t1 - e))) + aimag(local_sum)
  local_sum = cmplx( t1 + t2, t2 - ((t1 + t2) - t1) )
end do  

out_sum = real(local_sum)

return
end function drpdr_fast

#else

! CPU version
pure function drpdr_fast(nn,cq,xa,ya,za,at) result(out_sum)

integer, intent(in) :: nn
integer i, k, ilen, kx
real, intent(in) :: cq
real, dimension(:), intent(in) :: xa, ya, za
real, dimension(:,:), intent(in) :: at
real, dimension(size(at,2)) :: at_k
real, dimension(size(at,2)) :: out_sum
real, dimension(size(xa)) :: rb
real at_t, e, t1, t2
real ra
complex, dimension(size(at,2)) :: local_sum

ilen = size(xa)
kx = size(at,2)

local_sum(1:kx) = (0.,0.)
do i = 1,ilen
  ra = xa(nn)*xa(i) + ya(nn)*ya(i) + za(nn)*za(i)
  ra = acos(max(min(ra, 1.), -1.))
  rb(i) = exp(-min((cq*ra)**2,50.))  
end do

do i = 1,ilen
  at_k(:) = at(i,:)  
  do k = 1,kx
    at_t = rb(i)*at_k(k)
    t1 = at_t + real(local_sum(k))
    e  = t1 - at_t
    t2 = ((real(local_sum(k)) - e) + (at_t - (t1 - e))) + aimag(local_sum(k))
    local_sum(k) = cmplx( t1 + t2, t2 - ((t1 + t2) - t1) )
  end do  
end do  

out_sum(1:kx) = real(local_sum(1:kx))

return
end function drpdr_fast
    
#endif

subroutine nestin_exit

use onthefly_m

call onthefly_exit

end subroutine nestin_exit

end module nesting
