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
    
! This is the river routing which links to mlo.f90 and mlodynamics.f90
! ocean/lake model

! This version currently assumes that the orography is sufficently resolved
! to determine the river flow.  Plan to read in effective gradients between
! grid boxes from high resolution river flow datasets.    
    
module river

implicit none

private
public rvrinit, rvrrouter

integer, dimension(:,:), allocatable, save :: xp
real, dimension(:), allocatable, save ::  ee
real, dimension(:,:), allocatable, save :: idp

integer, parameter :: basinmd = 0         ! basin mode (0=soil, 2=pile-up, 3=leak)

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine initialises river arrays
!
subroutine rvrinit

use arrays_m
use cc_mpi
use indices_m
use nsibd_m
use riverarrays_m
use soil_m
use xyzinfo_m

implicit none

include 'const_phys.h'
include 'newmpar.h'

integer n, iq
real, dimension(ifull) :: r
real(kind=8), dimension(ifull+iextra,3) :: xyzbc

! setup indices and grid spacing
allocate( xp(ifull,8), idp(ifull,8) ) 
xp(1:ifull,1) = in
xp(1:ifull,2) = ie
xp(1:ifull,3) = is
xp(1:ifull,4) = iw
xp(1:ifull,5) = ine
xp(1:ifull,6) = ise
xp(1:ifull,7) = isw
xp(1:ifull,8) = inw

xyzbc(1:ifull,1) = x(1:ifull)
xyzbc(1:ifull,2) = y(1:ifull)
xyzbc(1:ifull,3) = z(1:ifull)
call boundsr8(xyzbc,corner=.true.)

! JLM suggests using x, y and z for calculating these distances
do n = 1,8
  r(:) = real(sum(xyzbc(1:ifull,:)*xyzbc(xp(:,n),:),2))
  r(:) = acos(max( min( r(:), 1. ), -1. ))*rearth
  idp(:,n) = 1./r(:)
end do

! prep land-sea mask
allocate( ee(ifull+iextra) )
ee(1:ifull+iextra) = 0.
where ( .not.land(1:ifull) )
  ee(1:ifull) = 1.
end where
call bounds(ee,corner=.true.)

! define outflow
outflowmask(1:ifull) = .false.
do iq = 1,ifull
  if ( isoilm_in(iq) == -1 ) then ! ee=1 implies water, isoilm_in=-1 implies inland water
    if ( zs(in(iq))-zs(iq)<-0.1 .and. ee(in(iq))<=0.5 ) outflowmask(iq) = .true.
    if ( zs(ie(iq))-zs(iq)<-0.1 .and. ee(ie(iq))<=0.5 ) outflowmask(iq) = .true.
    if ( zs(is(iq))-zs(iq)<-0.1 .and. ee(is(iq))<=0.5 ) outflowmask(iq) = .true.
    if ( zs(iw(iq))-zs(iq)<-0.1 .and. ee(iw(iq))<=0.5 ) outflowmask(iq) = .true.
    if ( zs(ine(iq))-zs(iq)<-0.1 .and. ee(ine(iq))<=0.5 ) outflowmask(iq) = .true.
    if ( zs(ise(iq))-zs(iq)<-0.1 .and. ee(ise(iq))<=0.5 ) outflowmask(iq) = .true.
    if ( zs(isw(iq))-zs(iq)<-0.1 .and. ee(isw(iq))<=0.5 ) outflowmask(iq) = .true.
    if ( zs(inw(iq))-zs(iq)<-0.1 .and. ee(inw(iq))<=0.5 ) outflowmask(iq) = .true.
  end if
end do

return
end subroutine rvrinit

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine calculates the river routing.
!
subroutine rvrrouter

use arrays_m
use cable_ccam 
use cc_mpi
use indices_m
use map_m
use nsibd_m
use riverarrays_m
use soil_m
use soilsnow_m

implicit none

include 'newmpar.h'
include 'const_phys.h'
include 'parm.h'
include 'soilv.h'

integer i, k
integer nit
real, dimension(ifull+iextra) :: netflx, newwat, zsadj
real, dimension(ifull) :: vel, soilsink
real, dimension(ifull) :: tmpry, tmprysave, deltmpry, ll
real, dimension(ifull) :: rate
real, dimension(ifull,8) :: slope, flow
real, dimension(ifull,8) :: fta, ftb, ftx, fty

! To speed up the code, we use a (semi-)implicit solution rather than an iterative approach
! This avoids additional MPI calls.

! Note that unlike Miller et al (1994), this scheme currently has multiple outflows from a given
! grid box (i.e., following Gordon in Mk3.5).  So if the slope of the grid box is divergent, then
! on average water leaves the grid box, whereas if the slope of the grid box is convergent, then
! on average water enters the grid box.

tmpry(:) = 0. ! for cray compiler

!--------------------------------------------------------------------
! update processor boundaries
call bounds(watbdy,corner=.true.)
newwat(1:ifull+iextra) = watbdy(1:ifull+iextra)

!--------------------------------------------------------------------
! predictor-corrector for water level and salinity
do nit = 1,2
  
  if ( nit == 2 ) call bounds(newwat,corner=.true.)

  ! calculate slopes
  do i = 1,8
    zsadj(1:ifull+iextra) = 0.001*0.5*(newwat(1:ifull+iextra)+watbdy(1:ifull+iextra))
    where ( ee(1:ifull)*ee(xp(1:ifull,i))>0.5 .or. edgetest(i) )
      slope(1:ifull,i) = 0. ! no orographic slope within ocean bounds
    elsewhere
      slope(1:ifull,i) = (zs(1:ifull)/grav+zsadj(1:ifull)-zs(xp(1:ifull,i))/grav-zsadj(xp(1:ifull,i)))*idp(1:ifull,i)
    end where
  end do

  newwat(1:ifull) = watbdy(1:ifull)

  ! Basic expression

  ! m = mass/area
  ! flow = m * vel / dx
  ! m(t+1)-m(t) = dt*sum(inflow)-dt*sum(outflow)

  ! outflow
  ! compute net outgoing flux for a grid box so that total water is conserved
  do i = 1,8
    vel(1:ifull) = min( 0.35*sqrt(max(slope(1:ifull,i),0.)/0.00005), 5. ) ! from Miller et al (1994)
    fta(1:ifull,i) = -dt*vel*idp(1:ifull,i) ! outgoing flux
  end do
  netflx(1:ifull) = sum( abs(fta(1:ifull,1:8)), 2 ) ! MJT notes - this will never trigger for sensible values of dt
  call bounds(netflx,corner=.true.)
  
  ! water outflow
  do i = 1,8
    where ( netflx(1:ifull)>1.E-10 )
      ftx(1:ifull,i) = -fta(1:ifull,i)/netflx(1:ifull) ! max fraction of total outgoing flux
      flow(1:ifull,i) = watbdy(1:ifull)*min( fta(1:ifull,i), ftx(1:ifull,i) ) ! (kg/m^2)
    elsewhere
      flow(1:ifull,i) = 0.
    end where
  end do
  newwat(1:ifull) = newwat(1:ifull) + sum(flow(1:ifull,1:8),2)

  ! inflow
  ! water inflow
  do i = 1,8
    vel(1:ifull) = min( 0.35*sqrt(max(-slope(1:ifull,i),0.)/0.00005), 5. ) ! from Miller et al (1994)
    ftb(1:ifull,i) = dt*vel*idp(1:ifull,i) ! incomming flux
    where ( netflx(xp(1:ifull,i))>1.E-10 )
      fty(1:ifull,i) = ftb(1:ifull,i)/netflx(xp(1:ifull,i)) ! max fraction of flux from outgoing cell
      flow(1:ifull,i) = watbdy(xp(1:ifull,i))*min( ftb(1:ifull,i), fty(1:ifull,i) ) ! (kg/m^2)
      flow(1:ifull,i) = flow(1:ifull,i)*(em(1:ifull)/em(xp(1:ifull,i)))**2 ! change in gridbox area
    elsewhere
      flow(1:ifull,i) = 0.
    end where
  end do
  newwat(1:ifull) = newwat(1:ifull) + sum(flow(1:ifull,1:8),2)

end do

watbdy(1:ifull) = max( newwat(1:ifull), 0. )

!--------------------------------------------------------------------
! Water losses over land basins

  
! Method for land basins
select case(basinmd)
  case(0)
    ! add water to soil moisture 
    ! estimate rate that water leaves river into soil
    rate(1:ifull) = min( watbdy(1:ifull)/100., 1. ) ! MJT suggestion
    if ( nsib==6 .or. nsib==7 ) then
      ! CABLE
      tmpry(1:ifull) = 0.
      where ( all(slope(1:ifull,1:8)<1.e-4,dim=2) .and. land(1:ifull) )
        tmpry(1:ifull) = watbdy(1:ifull)
      end where
      tmprysave(1:ifull) = tmpry(1:ifull)
      call cableinflow(tmpry,rate)
      soilsink(1:ifull) = (tmpry(1:ifull)-tmprysave(1:ifull))*(1.-sigmu(1:ifull))
      newwat(1:ifull) = newwat(1:ifull) + soilsink(1:ifull)
    else
      ! Standard land surface model
      deltmpry(1:ifull) = 0.
      do k = 1,ms
        where ( all(slope(1:ifull,1:8)<1.e-4,dim=2) .and. land(1:ifull) )
          ll(1:ifull) = max( sfc(isoilm(1:ifull))-wb(1:ifull,k), 0. )*1000.*zse(k)
          ll(1:ifull) = ll(1:ifull)*rate(1:ifull)
          ll(1:ifull) = min( tmpry(1:ifull)+deltmpry(1:ifull), ll(1:ifull) )
          wb(1:ifull,k) = wb(1:ifull,k) + ll(1:ifull)/(1000.*zse(k))
          deltmpry(1:ifull) = deltmpry(1:ifull) - ll(1:ifull)
        end where
      end do
      soilsink(1:ifull) = deltmpry(1:ifull)*(1.-sigmu(1:ifull))
      newwat(1:ifull) = newwat(1:ifull) + soilsink(1:ifull)
    end if
  case(2)
    ! pile-up water
  case(3)
    ! leak
    ! estimate rate that water leaves river into soil
    rate(1:ifull) = min( dt/(192.*3600.), 1. ) ! MJT suggestion
    if ( nsib==6 .or. nsib==7 ) then
      ! CABLE
      tmpry(1:ifull) = watbdy(1:ifull)
      call cableinflow(tmpry,rate)
      soilsink(1:ifull) = (tmpry(1:ifull)-watbdy(1:ifull))*(1.-sigmu(1:ifull))
      newwat(1:ifull) = newwat(1:ifull) + soilsink(1:ifull)
    else
      ! Standard land surface model
      deltmpry = 0.
      do k = 1,ms
        where ( land(1:ifull) )
          ll(1:ifull) = max( sfc(isoilm(1:ifull))-wb(1:ifull,k), 0. )*1000.*zse(k)
          ll(1:ifull) = ll(1:ifull)*rate(1:ifull)
          ll(1:ifull) = min( watbdy(1:ifull)+deltmpry(1:ifull), ll(1:ifull) )
          wb(1:ifull,k) = wb(1:ifull,k) + ll(1:ifull)/(1000.*zse(k))
          deltmpry(1:ifull) = deltmpry(1:ifull) - ll(1:ifull)
        end where
      end do
      soilsink(1:ifull) = deltmpry(1:ifull)*(1.-sigmu(1:ifull))
      newwat(1:ifull) = newwat(1:ifull) + soilsink(1:ifull)
    end if
  case default
    write(6,*) "ERROR: Unsupported basinmd ",basinmd
    call ccmpi_abort(-1)
end select

watbdy(1:ifull) = max( newwat(1:ifull), 0. )

! MLO (or other ocean model) will remove watbdy from ocean points when it updates its
! river inflows.

return
end subroutine rvrrouter

function edgetest(i) result(ans)

use cc_mpi

implicit none

include 'newmpar.h'

integer, intent(in) :: i
integer iq, n
logical, dimension(ifull) :: ans

ans(1:ifull) = .false.

select case(i)
  case(5)
    if ( edge_n .and. edge_e ) then
      do n = 1,npan
        iq = n*ipan*jpan
        ans(iq) = .true.
      end do
    end if
  case(6)
    if ( edge_s .and. edge_e ) then
      do n = 1,npan
        iq = ipan + (n-1)*ipan*jpan
        ans(iq) = .true.
      end do
    end if
  case(7)
    if ( edge_s .and. edge_w ) then
      do n = 1,npan
        iq = 1 + (n-1)*ipan*jpan
        ans(iq) = .true.
      end do
    end if
  case(8)
    if ( edge_n .and. edge_w ) then
      do n = 1,npan
        iq = 1 - ipan + n*ipan*jpan
        ans(iq) = .true.
      end do
    end if
end select
  
return
end function edgetest

end module river
