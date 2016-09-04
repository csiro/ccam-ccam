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

!------------------------------------------------------------------------------
    
! This is the river routing which links to mlo.f90 and mlodynamics.f90
! ocean/lake model

! This version currently assumes that the orography is sufficently resolved
! to determine the river flow.  Plan to read in effective gradients between
! grid boxes from high resolution river flow datasets.    
    
module river

implicit none

private
public basinmd, rivermd, rivercoeff
public rvrinit, rvrrouter

integer, dimension(:,:), allocatable, save :: xp
logical, dimension(:,:), allocatable, save :: river_inflow

integer, save :: basinmd = 0         ! basin mode (0=soil, 2=pile-up, 3=leak)
integer, save :: rivermd = 0         ! river mode (0=Miller, 1=Manning)
real, save :: rivercoeff = 0.02      ! river roughness coeff (Miller=0.02, A&B=0.035)

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine initialises river arrays
!
subroutine rvrinit

use arrays_m
use cc_mpi
use indices_m
use newmpar_m
use nsibd_m
use parm_m
use riverarrays_m
use soil_m
use xyzinfo_m

implicit none

include 'const_phys.h'

integer n, iq, iq_g, xp_g, nbasins, nbasins_g, nbasins_old
integer iqout
real(kind=8), dimension(ifull+iextra,3) :: xyzbc
real, dimension(ifull+iextra) ::  ee, river_outloc
real minzs, testzs, r, slope, vel
logical, dimension(ifull+iextra) :: basinmask

! setup indices
allocate( xp(ifull,8) ) 
xp(1:ifull,1) = in
xp(1:ifull,2) = ie
xp(1:ifull,3) = is
xp(1:ifull,4) = iw
xp(1:ifull,5) = ine
xp(1:ifull,6) = ise
xp(1:ifull,7) = isw
xp(1:ifull,8) = inw


! prep land-sea mask
ee(1:ifull+iextra) = 0.
where ( .not.land(1:ifull) )
  ee(1:ifull) = 1.
end where
call bounds(ee,corner=.true.)


allocate( river_inflow(ifull,8) )

!*******************************
! MOVE SECTION BELOW TO OCNBATH?

! calculate outflow
river_outloc(:) = -1
do iq = 1,ifull
  if ( isoilm_in(iq)/=0 ) then
    minzs = zs(iq)
    do n = 1,8
      testzs = zs(xp(iq,n))
      if ( testzs<minzs ) then
        iq_g = iq2iqg(iq)
        select case(n)
          case(1)
            xp_g = in_g(iq_g)
          case(2)
            xp_g = ie_g(iq_g)  
          case(3)
            xp_g = is_g(iq_g)  
          case(4)
            xp_g = iw_g(iq_g)  
          case(5)
            xp_g = ine_g(iq_g)  
          case(6)
            xp_g = ise_g(iq_g)  
          case(7)
            xp_g = isw_g(iq_g)  
          case(8)
            xp_g = inw_g(iq_g)  
        end select            
        river_outloc(iq) = real(xp_g)
        minzs = testzs
      end if
    end do
  end if
end do
call bounds(river_outloc)

! calculate inflows
river_inflow(:,:) = .false.
do iq = 1,ifull
  iq_g = iq2iqg(iq)
  do n = 1,8
    if ( nint(river_outloc(xp(iq,n)))==iq_g ) then
      river_inflow(iq,n) = .true.
    end if
  end do
end do

! determine locations and number of basins
basinmask(:) = ee(:)<=0.5 .and. river_outloc(:)==-1
nbasins = count( basinmask(:) )
call ccmpi_allreduce(nbasins,nbasins_g,'sum',comm_world)
if ( myid==0 ) then
  write(6,*) "Remaining river basins ",nbasins_g  
end if

! loop to remove basins
nbasins_old = nbasins_g + 1
do while ( nbasins_g<nbasins_old ) 

  ! fix basins
  do iq = 1,ifull
    if ( basinmask(iq) ) then
      minzs = 9.e20
      iq_g = iq2iqg(iq)
      do n = 1,8
        if ( .not.river_inflow(iq,n) .and. .not.basinmask(xp(iq,n)) ) then
          testzs = zs(xp(iq,n))
          if ( testzs<minzs ) then
            select case(n)
              case(1)
                xp_g = in_g(iq_g)
              case(2)
                xp_g = ie_g(iq_g)  
              case(3)
                xp_g = is_g(iq_g)  
              case(4)
                xp_g = iw_g(iq_g)  
              case(5)
                xp_g = ine_g(iq_g)  
              case(6)
                xp_g = ise_g(iq_g)  
              case(7)
                xp_g = isw_g(iq_g)  
              case(8)
                xp_g = inw_g(iq_g)  
            end select            
            river_outloc(iq) = real(xp_g)
            minzs = testzs
          end if
        end if
      end do
    end if
  end do
  call bounds(river_outloc)

  ! re-calculate inflows
  river_inflow(:,:) = .false.
  do iq = 1,ifull
    iq_g = iq2iqg(iq)
    do n = 1,8
      if ( nint(river_outloc(xp(iq,n)))==iq_g ) then
        river_inflow(iq,n) = .true.
      end if
    end do
  end do

  ! calculate number of remaining basins
  nbasins_old = nbasins_g
  basinmask(:) = ee(:)<=0.5 .and. river_outloc(:)==-1
  nbasins = count( basinmask(:) )
  call ccmpi_allreduce(nbasins,nbasins_g,'sum',comm_world)
  if ( myid==0 ) then
    write(6,*) "Remaining river basins ",nbasins_g  
  end if
  
end do

! MOVE SECTION ABOVE TO OCNBATH?
!*******************************


!--------------------------------------------------------------------
! re-calculate inflows
river_inflow(:,:) = .false.
do iq = 1,ifull
  iq_g = iq2iqg(iq)
  do n = 1,8
    if ( nint(river_outloc(xp(iq,n)))==iq_g ) then
      river_inflow(iq,n) = .true.
    end if
  end do
end do


!--------------------------------------------------------------------
! calculate river outflow direction
river_outdir(:) = -1
do iq = 1,ifull
  if ( nint(river_outloc(iq))>0 ) then
    iq_g = iq2iqg(iq)
    do n = 1,8
      select case(n)
        case(1)
          xp_g = in_g(iq_g)
        case(2)
          xp_g = ie_g(iq_g)  
        case(3)
          xp_g = is_g(iq_g)  
        case(4)
          xp_g = iw_g(iq_g)  
        case(5)
          xp_g = ine_g(iq_g)  
        case(6)
          xp_g = ise_g(iq_g)  
        case(7)
          xp_g = isw_g(iq_g)  
        case(8)
          xp_g = inw_g(iq_g)  
      end select            
      if ( xp_g==nint(river_outloc(iq)) ) then
        river_outdir(iq) = n
        exit
      end if
    end do
    if ( river_outdir(iq)==-1 ) then
      write(6,*) "ERROR: Internal error in rvrinit."
      write(6,*) "Cannot match river_outdir to river_outloc"
      call ccmpi_abort(-1)
    end if
  end if
end do


!--------------------------------------------------------------------
! define outflow mask for in-land water
! (typically used for lake overflow or in-land sea overflow)
outflowmask(:) = .false.
do iq = 1,ifull
  if ( isoilm_in(iq)==-1 .and. river_outdir(iq)>0 ) then
    if ( ee(xp(iq,river_outdir(iq)))<=0.5 ) then
      outflowmask(iq) = .true.
    end if
  end if
end do


!--------------------------------------------------------------------
! Update coordinates for distance calculations
xyzbc(1:ifull,1) = x(1:ifull)
xyzbc(1:ifull,2) = y(1:ifull)
xyzbc(1:ifull,3) = z(1:ifull)
call boundsr8(xyzbc,corner=.true.)


!--------------------------------------------------------------------
! calculate outflow flux
river_vel(:) = 0.
river_slope(:) = 0.
river_dx(:) = 1.e-9
do iq = 1,ifull
  if ( river_outdir(iq)>0 ) then  
    iqout = xp(iq,river_outdir(iq))
    if ( ee(iq)*ee(iqout)<=0.5 ) then
      ! JLM suggests using x, y and z for calculating these distances
      r = real(sum(xyzbc(iq,:)*xyzbc(iqout,:)))
      r = acos(max( min( r, 1. ), -1. ))*rearth
      r = max(r, 1.e-9)
      slope = max( (zs(iq)-zs(iqout))/grav, 0. )/r
      ! river_vel is output diagnostic.  Always set to Miller
      river_vel(iq) = max( min( rivercoeff*sqrt(slope), 5. ), 0.15 ) ! from Miller et al (1994)
      river_slope(iq) = slope
      river_dx(iq) = r
    end if
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
use newmpar_m
use nsibd_m
use parm_m
use riverarrays_m
use soil_m
use soilsnow_m
use soilv_m

implicit none

include 'const_phys.h'

integer i, k
real, dimension(ifull+iextra) :: outflow
real, dimension(ifull) :: inflow
real, dimension(ifull) :: tmpry, tmprysave, deltmpry, ll
real, dimension(ifull) :: rate, soilsink

tmpry(:) = 0. ! for cray compiler


! Basic expression

! m = mass/area
! flux = m * vel * dt / dx
! vel = sqrt(slope) * K
! r=Length, K=coeff, slope=delzs/r
! m(t+1)-m(t) = sum(inflow)-outflow

!--------------------------------------------------------------------
! calculate outflow
select case(rivermd)
  case(1) ! Miller
    outflow(1:ifull) = (dt/river_dx(1:ifull))*river_vel(1:ifull)*watbdy(1:ifull) ! (kg/m^2)
  case(2) ! Manning ( approximated )
    outflow(1:ifull) = rivercoeff*sqrt(river_slope(1:ifull))*max(watbdy(1:ifull),0.)**(2./3.)
    outflow(1:ifull) = max( 0.15, min( 5., outflow(1:ifull) ) )
    outflow(1:ifull) = (dt/river_dx(1:ifull))*outflow(1:ifull)*watbdy(1:ifull) ! (kg/m^2)
  case default
    write(6,*) "ERROR: Unknown option for rivermd=",rivermd
    call ccmpi_abort(-1)
end select
outflow(1:ifull) = max( 0., min( watbdy(1:ifull), outflow(1:ifull) ) )
call bounds(outflow)

!--------------------------------------------------------------------
! calculate inflow
inflow(:) = 0.
do i = 1,8
  where ( river_inflow(1:ifull,i) )
    ! adjust for change in grid-box area when calculating inflows
    inflow(1:ifull) = inflow(1:ifull) + outflow(xp(1:ifull,i))*(em(1:ifull)/em(xp(1:ifull,i)))**2
  end where
end do

watbdy(1:ifull) = watbdy(1:ifull) - outflow(1:ifull) + inflow(1:ifull)


!--------------------------------------------------------------------
! Water losses over land basins

  
! Method for land basins
select case(basinmd)
  case(0)
    ! add water to soil moisture 
    ! estimate rate that water leaves river into soil
    rate(1:ifull) = 1. ! MJT suggestion
    where ( river_outdir(1:ifull)==-1 .and. land(1:ifull) )
      tmpry(1:ifull) = watbdy(1:ifull)
    elsewhere
      tmpry(1:ifull) = 0.
    end where
    if ( nsib==6 .or. nsib==7 ) then
      ! CABLE
      tmprysave(1:ifull) = tmpry(1:ifull)
      call cableinflow(tmpry,rate)
      deltmpry(1:ifull) = tmpry(1:ifull) - tmprysave(1:ifull)
      soilsink(1:ifull) = deltmpry(1:ifull)*(1.-sigmu(1:ifull))
      watbdy(1:ifull) = watbdy(1:ifull) + soilsink(1:ifull)
    else
      ! Standard land surface model
      deltmpry(1:ifull) = 0.
      do k = 1,ms
        where ( river_outdir(1:ifull)==-1 .and. land(1:ifull) )
          ll(1:ifull) = max( sfc(isoilm(1:ifull))-wb(1:ifull,k), 0. )*1000.*zse(k)
          ll(1:ifull) = ll(1:ifull)*rate(1:ifull)
          ll(1:ifull) = min( tmpry(1:ifull)+deltmpry(1:ifull), ll(1:ifull) )
          wb(1:ifull,k) = wb(1:ifull,k) + ll(1:ifull)/(1000.*zse(k))
          deltmpry(1:ifull) = deltmpry(1:ifull) - ll(1:ifull)
        end where
      end do
      soilsink(1:ifull) = deltmpry(1:ifull)*(1.-sigmu(1:ifull))
      watbdy(1:ifull) = watbdy(1:ifull) + soilsink(1:ifull)
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
      watbdy(1:ifull) = watbdy(1:ifull) + soilsink(1:ifull)
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
      watbdy(1:ifull) = watbdy(1:ifull) + soilsink(1:ifull)
    end if
  case default
    write(6,*) "ERROR: Unsupported basinmd ",basinmd
    call ccmpi_abort(-1)
end select

watbdy(1:ifull) = max( watbdy(1:ifull), 0. )

! MLO (or other ocean model) will remove watbdy from ocean points when it updates its
! river inflows.

return
end subroutine rvrrouter

end module river
