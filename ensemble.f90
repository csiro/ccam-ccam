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
    
! CCAM ensemble initialisation and breeding

! This module is to assist with the creation and breeding of
! ensemble members for forecasting applications
    
! ensemble_mode = 0       off
! ensemble_mode = 1       control
! ensemble_mode = 2       breeding

! ensemble_period = 720   (12 hour update by default)
    
! ensemble_rsfactor = 0.1 (rescaling of initial perturbation)
    
! ensembleoutfile         (filename prefix when creating ensemble members)
    

module ensemble

private
public update_ensemble
public ensembleoutfile

integer, save :: mtimea = -1  ! previous mesonest time (mins)
integer, save :: mtimeb = -1  ! timer for updating analysis data
real, save :: refsum = -1.

real, dimension(:), allocatable, save :: tssa
real, dimension(:), allocatable, save :: pslb, tssb, fraciceb
real, dimension(:), allocatable, save :: sicedepb
real, dimension(:,:), allocatable, save :: tb, ub, vb, qb

interface rmse
  module procedure rmse_2d
  module procedure rmse_3d
end interface rmse

interface calcmean
  module procedure calcmean_2d
  module procedure calcmean_3d
end interface calcmean

contains
    
subroutine update_ensemble

use aerosolldr                   ! LDR prognostic aerosols
use arrays_m                     ! Atmosphere dyamics prognostic arrays
use cc_mpi                       ! CC MPI routines
use dates_m                      ! Date data
use mlo                          ! Ocean physics and prognostic arrays
use newmpar_m                    ! Grid parameters
use onthefly_m                   ! Input interpolation routines
use parm_m                       ! Model configuration
use pbl_m                        ! Boundary layer arrays
use soil_m                       ! Soil and surface data
use soilsnow_m                   ! Soil, snow and surface data

implicit none

integer, dimension(ifull) :: dumm
integer, save :: num = 0
integer kdate_r, ktime_r, k
integer kdhour, kdmin, kddate, khour_r, khour, kmin_r, kmin
integer mtimeb_old
real timerm, cona, conb, wgt
real, dimension(ifull) :: dd, old, new, delta
real, dimension(ifull) :: zsb, timelt
real, dimension(:), allocatable :: psl_pos, psl_neg
real, dimension(:), allocatable, save :: u_rms, v_rms
real, dimension(:), allocatable, save :: t_rms
real, dimension(ifull,3) :: duma
real, dimension(ifull,6) :: dumd
real, dimension(ifull,kl) :: ee
real, dimension(:,:), allocatable :: u_pos, v_pos, t_pos
real, dimension(:,:), allocatable :: u_neg, v_neg, t_neg
real, dimension(ifull,ms,3) :: dumg
real, dimension(ifull,kl,5) :: dumv
real, dimension(ifull,3,3) :: dums
real, dimension(ifull,wlev,8) :: dumo
real, dimension(ifull,kl,naero) :: dumr
real, save :: psl_rms

if ( mtimer>mtimeb ) then

  ! allocate arrays on first call     
  if ( .not.allocated(tb) ) then
    allocate( tssa(ifull) )
    allocate( tb(ifull,kl), ub(ifull,kl), vb(ifull,kl), qb(ifull,kl) )
    allocate( pslb(ifull), tssb(ifull), fraciceb(ifull) )
    allocate( sicedepb(ifull) )
    allocate( u_rms(kl), v_rms(kl), t_rms(kl) )
    psl_rms = 1.
    u_rms(:) = 1.
    v_rms(:) = 1.
    t_rms(:) = 1.
    if ( abs(io_in)==1 ) then
      call onthefly(3,kdate_r,ktime_r,                            &
                    pslb,zsb,tssb,sicedepb,fraciceb,tb,ub,vb,qb,  &
                    dumg(:,:,1),dumg(:,:,2),dumg(:,:,3),          &
                    duma(:,1),dumv(:,:,1),dumv(:,:,2),            &
                    dumv(:,:,3),dumv(:,:,4),dumv(:,:,5),          &
                    dums(:,:,1),dums(:,:,2),dums(:,:,3),          &
                    duma(:,2),duma(:,3),dumm,dumo,dumd,dumr)
      call retopo(pslb,zsb,zs,tb,qb)
    else
      write(6,*) 'ERROR: Scale-selective filter requires abs(io_in)=1'
      call ccmpi_abort(-1)
    endif   ! (abs(io_in)==1)
    if ( ensemble_mode==2 ) then
      ! calculate RMSE of perturbation
      psl_rms = -1.
      u_rms(:) = -1.
      v_rms(:) = -1.
      t_rms(:) = -1.
      dd(:) = psl(1:ifull) - pslb
      call rmse(dd,psl_rms)
      psl_rms = ensemble_rsfactor*psl_rms
      ee(:,:) = u(1:ifull,:) - ub
      call rmse(ee,u_rms)
      u_rms = ensemble_rsfactor*u_rms
      ee(:,:) = v(1:ifull,:) - vb
      call rmse(ee,v_rms)
      v_rms = ensemble_rsfactor*v_rms
      ee(:,:) = t(1:ifull,:) - tb
      call rmse(ee,t_rms)
      t_rms = ensemble_rsfactor*t_rms
    end if
    mtimeb = 0
  end if
  
  mtimea = mtimeb
  tssa(:) = tssb(:)

  ! following (till end of subr) reads in next bunch of data in readiness
  ! read tb etc  - for globpea, straight into tb etc
  mtimeb_old = mtimeb
  do while ( mtimeb-mtimeb_old  < ensemble_period )
    if ( abs(io_in)==1 ) then
      call onthefly(3,kdate_r,ktime_r,                            &
                    pslb,zsb,tssb,sicedepb,fraciceb,tb,ub,vb,qb,  &
                    dumg(:,:,1),dumg(:,:,2),dumg(:,:,3),          &
                    duma(:,1),dumv(:,:,1),dumv(:,:,2),            &
                    dumv(:,:,3),dumv(:,:,4),dumv(:,:,5),          &
                    dums(:,:,1),dums(:,:,2),dums(:,:,3),          &
                    duma(:,2),duma(:,3),dumm,dumo,dumd,dumr)
    else
      write(6,*) 'ERROR: Scale-selective filter requires abs(io_in)=1'
      call ccmpi_abort(-1)
    endif   ! (abs(io_in)==1)
    tssb(:) = abs(tssb(:))  ! moved here Mar '03
    qb(1:ifull,:) = max( qb(1:ifull,:), 0. )

    ! calculate time for next filter call
    khour_r = ktime_r/100
    khour   = ktime/100
    kdhour = khour_r - khour
    kmin_r = ktime_r - 100*khour_r
    kmin   = ktime   - 100*khour
    kdmin  = kmin_r - kmin
    kddate = iabsdate(kdate_r,kdate) - iabsdate(kdate,kdate)
    mtimeb = 1440*kddate + 60*kdhour + kdmin
  end do
  
  if ( myid==0 ) then
    write(6,*) "ensemble_period, mtimeb, mtimeb_old ",ensemble_period,mtimeb,mtimeb_old    
  end if

  ! adjust input data for change in orography
  call retopo(pslb,zsb,zs,tb,qb)
  
end if ! (mtimer>mtimeb)

! Apply filter to model data using previously loaded host data
if ( mtimer>=mtimeb .and. mod(nint(ktau*dt),60)==0 ) then

  ! create or breed ensemble members
  select case(ensemble_mode)
    case(1) ! control
      if ( myid==0 ) then
        write(6,*) "Create ensemble initial conditions"  
      end if    
      allocate( psl_pos(ifull), psl_neg(ifull) )
      allocate( u_pos(ifull,kl), u_neg(ifull,kl) )
      allocate( v_pos(ifull,kl), v_neg(ifull,kl) )
      allocate( t_pos(ifull,kl), t_neg(ifull,kl) )
      dd(:) = psl(1:ifull) - pslb
      psl_pos(:) = max( min( pslb + dd, 0.3 ), -1.4 )
      psl_neg(:) = max( min( pslb - dd, 0.3 ), -1.4 )
      ee(:,:) = u(1:ifull,:) - ub
      u_pos(:,:) = max( min( ub + ee, 350. ), -350. )
      u_neg(:,:) = max( min( ub - ee, 350. ), -350. )
      ee(:,:) = v(1:ifull,:) - vb
      v_pos(:,:) = max( min( vb + ee, 350. ), -350. )
      v_neg(:,:) = max( min( vb - ee, 350. ), -350. )
      ee(:,:) = t(1:ifull,:) - tb
      t_pos(:,:) = max( min( tb + ee, 400. ), 100. )
      t_neg(:,:) = max( min( tb - ee, 400. ), 100. )
      num = num + 1 ! odd = pos
      call saveoutput(num,psl_pos,u_pos,v_pos,t_pos,qg)
      num = num + 1 ! even = neg
      call saveoutput(num,psl_neg,u_neg,v_neg,t_neg,qg)
      deallocate( psl_pos, psl_neg )
      deallocate( u_pos, u_neg )
      deallocate( v_pos, v_neg )
      deallocate( t_pos, t_neg )
      psl(1:ifull) = pslb(:)
      u(1:ifull,:) = ub(:,:)
      v(1:ifull,:) = vb(:,:)
      t(1:ifull,:) = tb(:,:)
      qg(1:ifull,:) = qb(:,:)
    case(2) ! breed
      if ( myid==0 ) then
        write(6,*) "Update ensemble breeding"  
      end if    
      dd(:) = psl(1:ifull) - pslb
      call rmse(dd,psl_rms)
      psl(1:ifull) = pslb + dd
      ee(:,:) = u(1:ifull,:) - ub
      call rmse(ee,u_rms)
      u(1:ifull,:) = ub + ee
      ee(:,:) = v(1:ifull,:) - vb
      call rmse(ee,v_rms)
      v(1:ifull,:) = vb + ee
      ee(:,:) = t(1:ifull,:) - tb
      call rmse(ee,t_rms)
      t(1:ifull,:) = tb + ee
  end select  
  
end if     ! (mtimer==mtimec).and.(mod(nint(ktau*dt),60)==0)


! use time interpolated values for tss
timerm = real(ktau)*dt/60.
cona = (real(mtimeb)-timerm)/real(mtimeb-mtimea)
conb = (timerm-real(mtimea))/real(mtimeb-mtimea)

! specify sea-ice if not AMIP or if Mixed-Layer-Ocean
if ( namip==0 ) then  ! namip SSTs/sea-ice take precedence
  if ( nmlo==0 ) then
    ! check whether present ice points should change to/from sice points
    where ( fraciceb(:)>0. .and. fracice(:)<1.e-20 .and. .not.land(1:ifull) )
      ! N.B. if already a sice point, keep present tice (in tggsn)
      tggsn(:,1) = min( 271.2, tssb(:) )
    end where
    ! update tss 
    where ( .not.land(1:ifull) )
      tss(:) = cona*tssa(:) + conb*tssb(:)
      tgg(:,1) = tss(:)
      sicedep(:) = sicedepb(:)
      fracice(:) = fraciceb(:)
    end where  ! (.not.land(iq))
  else
    ! nudge Mixed-Layer-Ocean
    if ( nud_hrs<=0 ) then
      write(6,*) "ERROR: ensemble_mode/=0 has been selected with in-line ocean model (nmlo/=0)"  
      write(6,*) "nud_hrs>0 must be specified for relaxiation of SSTs"
      call ccmpi_abort(-1)
    end if
    old = cona*tssa + conb*tssb
    call mloexpmelt(timelt)
    timelt = min(timelt,271.2)
    new = (cona*tssa + conb*tssb)*(1.-fraciceb(:)) + timelt*fraciceb(:) - wrtemp
    wgt = dt/real(nud_hrs*3600)
    call mloexport(0,old,1,0)
    delta = new - old
    do k = 1,wlev
      call mloexport(0,old,k,0)
      old = wgt*delta + old
      call mloimport(0,old,k,0)
    end do  
  end if ! (nmlo==0) ..else..
end if   ! (namip==0)


return
end subroutine update_ensemble

subroutine saveoutput(num,psl_in,u_in,v_in,t_in,q_in)

use filnames_m
use newmpar_m
use outcdf

implicit none

integer, intent(in) :: num
real, dimension(:), intent(in) :: psl_in
real, dimension(:,:), intent(in) :: u_in, v_in, t_in, q_in
character(len=3) :: ensemblenum
character(len=1024) :: ensemblefilename

write(ensemblenum,'(I3.3)') num
ensemblefilename=trim(ensembleoutfile)//'.'//ensemblenum
call outfile(21,ensemblefilename,psl_in,u_in,v_in,t_in,q_in)

return
end subroutine saveoutput

subroutine rmse_2d(datin,ref_rms)

implicit none

real, intent(inout) :: ref_rms
real, dimension(:), intent(inout) :: datin
real, dimension(size(datin)) :: datsq
real :: datsum, datrms

! calculate RMSE
datsq = datin*datin
call calcmean(datsq,datsum)
datrms = sqrt(datsum)

! update reference RMSE if required
if ( ref_rms<0. ) then
  ref_rms = datrms
end if

! adjust
datin = datin*ref_rms/datrms

return
end subroutine rmse_2d

subroutine rmse_3d(datin,ref_rms)

implicit none

integer :: k
real, dimension(:,:), intent(inout) :: datin
real, dimension(size(datin,2)), intent(inout) :: ref_rms
real, dimension(size(datin,1),size(datin,2)) :: datsq
real, dimension(size(datin,2)) :: datrms, datsum

! calculate RMSE
do k = 1,size(datin,2)
  datsq(:,k) = datin(:,k)*datin(:,k)
end do  
call calcmean(datsq,datsum)
do k = 1,size(datin,2)
  datrms(k) = sqrt(datsum(k))
end do  

! update reference RMSE if required
if ( ref_rms(1)<0. ) then
  ref_rms(:) = datrms(:)
end if

! adjust
do k = 1,size(datin,2)
  datin(:,k) = datin(:,k)*ref_rms(k)/datrms(k)
end do  

return
end subroutine rmse_3d

subroutine calcmean_2d(datin,sumout)

use cc_mpi
use xyzinfo_m 

implicit none

integer :: i
real, dimension(:), intent(in) :: datin
real, intent(out) :: sumout
real :: e, t1, t2
complex :: local_sum, global_sum

if ( refsum<0. ) then
  local_sum = (0.,0.)
  do i = 1,size(wts)
    t1 = wts(i) + real(local_sum)
    e  = t1 - wts(i)
    t2 = ((real(local_sum) - e) + (wts(i) - (t1 - e))) + aimag(local_sum)
    local_sum = cmplx( t1 + t2, t2 - ((t1 + t2) - t1) )
  end do  
  call ccmpi_allreduce(local_sum,global_sum,"sumdr",comm_world)
  refsum = real(global_sum)
end if

local_sum = (0.,0.)
do i = 1,size(datin)
  t1 = datin(i)*wts(i)/refsum + real(local_sum)
  e  = t1 - datin(i)*wts(i)/refsum
  t2 = ((real(local_sum) - e) + (datin(i)*wts(i)/refsum - (t1 - e))) + aimag(local_sum)
  local_sum = cmplx( t1 + t2, t2 - ((t1 + t2) - t1) )
end do  

! MPI
call ccmpi_allreduce(local_sum,global_sum,"sumdr",comm_world)

sumout = real(global_sum)

return
end subroutine calcmean_2d

subroutine calcmean_3d(datain,sumout)

use cc_mpi
use xyzinfo_m 

implicit none

integer :: i
real, dimension(:,:), intent(in) :: datain
real, dimension(size(datain,2)), intent(out) :: sumout
real, dimension(size(datain,2)) :: e, t1, t2, datloc
complex, dimension(size(datain,2)) :: local_sum, global_sum

if ( refsum<0. ) then
  local_sum(1) = (0.,0.)
  do i = 1,size(wts)
    t1(1) = wts(i) + real(local_sum(1))
    e(1)  = t1(1) - wts(i)
    t2(1) = ((real(local_sum(1)) - e(1)) + (wts(i) - (t1(1) - e(1)))) + aimag(local_sum(1))
    local_sum(1) = cmplx( t1(1) + t2(1), t2(1) - ((t1(1) + t2(1)) - t1(1)) )
  end do  
  call ccmpi_allreduce(local_sum(1),global_sum(1),"sumdr",comm_world)
  refsum = real(global_sum(1))
end if

local_sum(:) = (0.,0.)
do i = 1,size(datain,1)
  datloc(:) = datain(i,:)*wts(i)/refsum  
  t1(:) = datloc(:) + real(local_sum(:))
  e(:)  = t1(:) - datloc(:)
  t2(:) = ((real(local_sum(:)) - e(:)) + (datloc(:) - (t1(:) - e(:)))) + aimag(local_sum(:))
  local_sum(:) = cmplx( t1(:) + t2(:), t2(:) - ((t1(:) + t2(:)) - t1(:)) )
end do  

! MPI
call ccmpi_allreduce(local_sum,global_sum,"sumdr",comm_world)

sumout(:) = real(global_sum(:))

return
end subroutine calcmean_3d

pure function iabsdate(kdate_r,kdate) result(ans)

use parm_m

implicit none

integer, intent(in) :: kdate_r, kdate
integer iyear,iyear0,month,iday
integer months,nl
integer newdate_r, diffyear
integer, dimension(0:13) :: mdays
integer ans

mdays = (/0,31,59,90,120,151,181,212,243,273,304,334,365,396/)

iyear  = kdate_r/10000
iyear0 = kdate/10000                ! year of kdate
newdate_r = kdate_r - 10000*iyear
month = newdate_r/100
newdate_r = newdate_r - 100*month
iday = newdate_r

! calculate number of months since start of kdate year
diffyear = iyear - iyear0
months = diffyear*12 + month - 1  

if ( leap==1 ) then
  nl = 0
  if ( mod(iyear0,4)==0 ) nl = 1
  if ( mod(iyear0,100)==0 ) nl = 0
  if ( mod(iyear0,400)==0 ) nl = 1
  mdays(2:13) = mdays(2:13) + nl
end if

! Accumulate days month by month, up to last completed month
ans = mdays(months)

! Add days from this current month
ans = ans + iday

end function iabsdate

end module ensemble
