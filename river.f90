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

integer, save :: basinmd = 0         ! basin mode (0=soil)
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

integer n, iq, iq_g, xp_g, xpb_g, xp_i, xp_j, xp_n
integer iqout
real(kind=8), dimension(ifull+iextra,3) :: xyzbc
real, dimension(ifull+iextra) ::  ee
real, dimension(ifull+iextra,3) :: river_outloc
real minzs, testzs, r, slope, vel

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
river_outloc(:,1) = -1.
river_outloc(:,2) = 1.
river_outloc(:,3) = 0.
do iq = 1,ifull
  if ( isoilm_in(iq)/=0 ) then
    minzs = zs(iq)
    do n = 1,8
      testzs = zs(xp(iq,n))
      if ( testzs<minzs .and. ee(iq)*ee(xp(iq,n))<=0.5 ) then
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
        xp_n = xp_g/(il_g*il_g)
        xp_g = xp_g - xp_n*il_g*il_g
        xp_j = (xp_g-1)/il_g + 1
        xp_g = xp_g - (xp_j-1)*il_g
        xp_i = xp_g
        river_outloc(iq,:) = real((/ xp_i, xp_j, xp_n /))
        minzs = testzs
      end if
    end do
  end if
end do
call bounds(river_outloc,corner=.true.)


! MOVE SECTION ABOVE TO OCNBATH?
!*******************************


!--------------------------------------------------------------------
! calculate inflows
river_inflow(:,:) = .false.
do iq = 1,ifull
  iq_g = iq2iqg(iq)
  do n = 1,8
    iqout = xp(iq,n)  
    xp_i = nint(river_outloc(iqout,1))
    xp_j = nint(river_outloc(iqout,2))
    xp_n = nint(river_outloc(iqout,3))
    xp_g = xp_i + (xp_j-1)*il_g + xp_n*il_g*il_g
    if ( xp_g==iq_g ) then
      river_inflow(iq,n) = .true.
    end if
  end do
end do


!--------------------------------------------------------------------
! calculate river outflow direction
river_outdir(:) = -1
do iq = 1,ifull
  xp_i = nint(river_outloc(iq,1))
  xp_j = nint(river_outloc(iq,2))
  xp_n = nint(river_outloc(iq,3))
  xp_g = xp_i + (xp_j-1)*il_g + xp_n*il_g*il_g
  if ( xp_g>0 ) then
    iq_g = iq2iqg(iq)
    do n = 1,8
      select case(n)
        case(1)
          xpb_g = in_g(iq_g)
        case(2)
          xpb_g = ie_g(iq_g)  
        case(3)
          xpb_g = is_g(iq_g)  
        case(4)
          xpb_g = iw_g(iq_g)  
        case(5)
          xpb_g = ine_g(iq_g)  
        case(6)
          xpb_g = ise_g(iq_g)  
        case(7)
          xpb_g = isw_g(iq_g)  
        case(8)
          xpb_g = inw_g(iq_g)  
      end select            
      iqout = xp(iq,n)
      if ( xpb_g==xp_g .and. ee(iq)*ee(iqout)<=0.5 ) then
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
real, dimension(ifull) :: soilsink

tmpry(:) = 0. ! for cray compiler


! Basic expression

! m = mass/area
! flux = m * vel * dt / dx
! vel = sqrt(slope) * K
! r=Length, K=coeff, slope=delzs/r
! m(t+1)-m(t) = sum(inflow)-outflow

!--------------------------------------------------------------------
! calculate outflow
outflow(1:ifull) = 0.
select case(rivermd)
  case(0) ! Miller
    where ( river_outdir(1:ifull)>0 )  
      outflow(1:ifull) = (dt/river_dx(1:ifull))*river_vel(1:ifull)*watbdy(1:ifull) ! (kg/m^2)
    end where
  case(1) ! Manning ( approximated )
    where ( river_outdir(1:ifull)>0 )  
      outflow(1:ifull) = rivercoeff*sqrt(river_slope(1:ifull))*max(watbdy(1:ifull),0.)**(2./3.)
      outflow(1:ifull) = max( 0.15, min( 5., outflow(1:ifull) ) )
      outflow(1:ifull) = (dt/river_dx(1:ifull))*outflow(1:ifull)*watbdy(1:ifull) ! (kg/m^2)
    end where
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
    where ( river_outdir(1:ifull)==-1 .and. land(1:ifull) )
      tmpry(1:ifull) = watbdy(1:ifull)
    elsewhere
      tmpry(1:ifull) = 0.
    end where
    if ( nsib==6 .or. nsib==7 ) then
      ! CABLE
      tmprysave(1:ifull) = tmpry(1:ifull)
      call cableinflow(tmpry)
      soilsink(1:ifull) = (tmpry(1:ifull)-tmprysave(1:ifull))*(1.-sigmu(1:ifull))
      watbdy(1:ifull) = watbdy(1:ifull) + soilsink(1:ifull)
    else
      ! Standard land surface model
      deltmpry(1:ifull) = 0.
      do k = 1,ms
        where ( river_outdir(1:ifull)==-1 .and. land(1:ifull) )
          ll(1:ifull) = max( sfc(isoilm(1:ifull))-wb(1:ifull,k), 0. )*1000.*zse(k)
          ll(1:ifull) = min( tmpry(1:ifull)+deltmpry(1:ifull), ll(1:ifull) )
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
