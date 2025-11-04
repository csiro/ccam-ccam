! Conformal Cubic Atmospheric Model
    
! Copyright 2015-2025 Commonwealth Scientific Industrial Research Organisation (CSIRO)
    
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

! This version uses FAM to determine the river flow direction.  The river
! velocity is based on Miller or a modified Manning approach.  Future work
! will focus on wetlands.
    
module river

implicit none

private
public basinmd, rivermd, rivercoeff, wt_transport
public rvrinit, rvrrouter
public water_table_transport

integer, dimension(:,:), allocatable, save :: xp
logical, dimension(:,:), allocatable, save :: river_inflow

integer, save :: basinmd      = 1    ! basin mode (0=soil, 1=redistribute)
integer, save :: rivermd      = 0    ! river mode (0=Miller, 1=Manning, 2=Wilms)
integer, save :: wt_transport = 0    ! water table transport (0=off, 1=on)
real, save :: rivercoeff = 0.02      ! river roughness coeff (Miller=0.02, A&B=0.035)
real, parameter :: rhow = 1000.      ! density of water (kg/m^3)

! Arrays for rivermd=2 (Wilms)
integer, save :: total_keys = 0
integer, save :: count_inflows = 0
integer, dimension(:,:), allocatable, save :: indices_list
integer, dimension(:,:), allocatable, save :: river_inflow_keys
integer, dimension(:,:,:), allocatable, save :: river_outloc
real, dimension(:,:), allocatable, save :: riverwidth

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine initialises river arrays
!
subroutine rvrinit(river_accin)

use arrays_m
use cc_mpi
use const_phys
use indices_m
use map_m
use newmpar_m
use nsibd_m
use parm_m
use riverarrays_m
use soil_m
use xyzinfo_m

integer, dimension(ifull), intent(in) :: river_accin
integer, dimension(ifull+iextra) :: river_outloc_cc, river_acc
integer, dimension(ifull) :: xp_i, xp_j, xp_n
integer n, iq, iq_g, xp_g, xpb_g, i, j, ridx
integer iqout, maxacc, testacc
real(kind=8) :: dr
real, dimension(ifull+iextra) ::  ee
real, dimension(ifull+iextra,3) :: r_outloc
real minzs, testzs, slope
real grid_area

! for cray compiler
xp_g = 0
xpb_g = 0


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



! update river acc halo
r_outloc(:,:) = 0.
r_outloc(1:ifull,1) = real(river_accin(1:ifull))
call bounds(r_outloc(:,1),corner=.true.)
river_acc(:) = nint(r_outloc(:,1))



! calculate outflow
river_outloc_cc(:) = -1
do iq = 1,ifull
  if ( isoilm_in(iq)/=0 ) then
 
    ! use FAM method
    maxacc = river_acc(iq)
    do n = 1,8
      iqout = xp(iq,n)
      testacc = river_acc(iqout)
      if ( testacc>maxacc .and. ee(iq)*ee(iqout)<=0.5 ) then
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
        river_outloc_cc(iq) = xp_g
        maxacc = testacc
      end if
    end do
    
    ! No data for FAM, try steepest descent
    if ( maxacc==river_acc(iq) ) then
      minzs = zs(iq)
      do n = 1,8
        iqout = xp(iq,n)
        testzs = zs(iqout)
        if ( testzs<minzs .and. ee(iq)*ee(iqout)<=0.5 .and.  &
            (river_acc(iq)==river_acc(iqout).or.river_acc(iqout)==0) ) then
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
          river_outloc_cc(iq) = xp_g
          minzs = testzs
        end if
      end do
    end if
    
  end if
end do

! outflow calculation
river_outdir(:) = -1
do iq = 1,ifull
  iq_g = iq2iqg(iq)
  xp_g = river_outloc_cc(iq)
  if ( xp_g>0 ) then
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
      if ( xpb_g==xp_g ) then
        river_outdir(iq) = n
        exit
      end if
    end do
    if ( river_outdir(iq)==-1 ) then
      write(6,*) "ERROR: Internal error in rvrinit."
      write(6,*) "Cannot match river_outdir to river_outloc_cc"
      call ccmpi_abort(-1)
    end if
  end if
  
end do


! update river_outloc_cc halo
r_outloc(1:ifull+iextra,1) = -1.
r_outloc(1:ifull+iextra,2) = -1.
r_outloc(1:ifull+iextra,3) = -1.
where ( river_outloc_cc(1:ifull)>0 )
  xp_n(:) = river_outloc_cc(1:ifull)/(il_g*il_g)
  river_outloc_cc(1:ifull) = river_outloc_cc(1:ifull) - xp_n(:)*il_g*il_g
  xp_j(:) = (river_outloc_cc(1:ifull)-1)/il_g + 1
  river_outloc_cc(1:ifull) = river_outloc_cc(1:ifull) - (xp_j(:)-1)*il_g
  xp_i(:) = river_outloc_cc(1:ifull)
  r_outloc(1:ifull,1) = real(xp_i(:))
  r_outloc(1:ifull,2) = real(xp_j(:))
  r_outloc(1:ifull,3) = real(xp_n(:))
end where
call bounds(r_outloc,corner=.true.)
where ( nint(r_outloc(1:ifull+iextra,1))>0 )
  river_outloc_cc(1:ifull+iextra) = nint(r_outloc(1:ifull+iextra,1)) + (nint(r_outloc(1:ifull+iextra,2))-1)*il_g + &
                                    nint(r_outloc(1:ifull+iextra,3))*il_g*il_g
elsewhere
  river_outloc_cc(1:ifull+iextra) = -1
end where


!--------------------------------------------------------------------
! calculate inflow directions
allocate( river_inflow(ifull,8) )
river_inflow(:,:) = .false.
do iq = 1,ifull
  iq_g = iq2iqg(iq)
  do n = 1,8
    iqout = xp(iq,n)  
    xp_g = river_outloc_cc(iqout)
    if ( xp_g==iq_g ) then
      river_inflow(iq,n) = .true.
    end if
  end do
end do


!--------------------------------------------------------------------
! define outflow mask for in-land water
! (typically used for lake overflow or in-land sea overflow)
outflowmask(:) = .false.
do iq = 1,ifull
  ridx = max( river_outdir(iq), 1 )
  if ( isoilm_in(iq)==-1 .and. river_outdir(iq)>0 ) then
    if ( ee(xp(iq,ridx))<=0.5 ) then
      outflowmask(iq) = .true.
    end if
  end if
end do


!--------------------------------------------------------------------
! calculate outflow velocity, slope and dx
! JLM suggests using x, y and z for calculating distances,
! instead of using the map factor.
river_vel(:) = 0.
river_dx(:) = 1.e-9
do iq = 1,ifull
  if ( river_outdir(iq)>0 ) then  
    iqout = xp(iq,river_outdir(iq))
    dr = x(iq)*x(iqout) + y(iq)*y(iqout) + z(iq)*z(iqout)
    dr = acos(max( min( dr, 1._8 ), -1._8 ))*real(rearth,8)
    dr = max(dr, 1.e-9_8)
    slope = max( (zs(iq)-zs(iqout))/grav, 0. )/real(dr)
    ! river_vel is output diagnostic that is always calculated according to Miller
    river_vel(iq) = max( min( rivercoeff*sqrt(slope), 5. ), 0.15 ) ! from Miller et al (1994)
    river_dx(iq) = real(dr)
  end if
end do


! initialise arrays for Wilms
if ( rivermd==2 ) then
  total_keys = count( land(1:ifull) )
  allocate( indices_list(total_keys,2) )
  total_keys = 0
  do iq = 1,ifull
    if ( land(iq) ) then
      total_keys = total_keys + 1  
      indices_list(total_keys,1) = mod(iq-1,il)+1
      indices_list(total_keys,2) = (iq-1)/il + 1
      ! iq = indices_list(:,1) + (indices_list(:,2)-1)*il
    end if
  end do  
  count_inflows = 0
  do iq = 1,ifull
    if ( any( river_inflow(iq,:) ) ) then
      count_inflows = count_inflows + 1  
    end if
  end do  
  allocate( river_inflow_keys(2,count_inflows) )
  count_inflows = 0
  do j = 1,jl
    do i = 1,il
      iq = i + (j-1)*il  
      if ( any( river_inflow(iq,:) ) ) then
        count_inflows = count_inflows + 1
        river_inflow_keys(1,count_inflows) = i
        river_inflow_keys(2,count_inflows) = j
      end if    
    end do    
  end do
  allocate( riverwidth(il,jl) )
  !the width through which water may move is equal to sqrt of area of cell (yes, i know...)
  !should try linking this to the waterheight.
  do j = 1,jl
    do i = 1,il
      iq = i + (j-1)*il
      grid_area = (ds/em(iq))**2
      riverwidth(i,j) = sqrt(grid_area)
    end do
  end do  
  allocate( river_outloc(2,il,jl) )
  do j = 1,jl
    do i = 1,il
      iq = i + (j-1)*il
      select case(river_outdir(iq))
        case(1) !n
          river_outloc(1,i,j) = i
          river_outloc(2,i,j) = j + 1
        case(2) !e
          river_outloc(1,i,j) = i + 1
          river_outloc(2,i,j) = j
        case(3) !s
          river_outloc(1,i,j) = i
          river_outloc(2,i,j) = j - 1
        case(4) !w
          river_outloc(1,i,j) = i - 1
          river_outloc(2,i,j) = j
        case(5) !ne
          river_outloc(1,i,j) = i + 1
          river_outloc(2,i,j) = j + 1
        case(6) !se
          river_outloc(1,i,j) = i + 1
          river_outloc(2,i,j) = j - 1
        case(7) !sw
          river_outloc(1,i,j) = i - 1
          river_outloc(2,i,j) = j - 1
        case(8) !nw
          river_outloc(1,i,j) = i - 1
          river_outloc(2,i,j) = j + 1
        case default
          river_outloc(1,i,j) = 0
          river_outloc(2,i,j) = 0
      end select    
    end do
  end do  
end if

return
end subroutine rvrinit

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine calculates the river routing.
!
subroutine rvrrouter

use arrays_m
use cc_mpi
use const_phys
use indices_m
use map_m
use newmpar_m
use nsibd_m
use parm_m
use riverarrays_m
use sflux_m
use soil_m
use soilsnow_m
use soilv_m

integer i, k, iq, iqout
real alph_p, delpos, delneg
real, dimension(ifull+iextra) :: outflow
real, dimension(ifull) :: inflow, vel, river_slope
real, dimension(ifull) :: tmpry, tmprysave, deltmpry, ll
real, dimension(ifull) :: soilsink
real, dimension(ifull) :: watbdy_mask, newwatbdy_mask, diffwatbdy
logical, dimension(ifull) :: basin_mask

tmpry(:) = 0. ! for cray compiler


if ( rivermd==2 ) then

  ! call Josefine Wilms method
  call wilms(il,jl)

else
  
  ! Basic expression

  ! m = mass/area
  ! vel = sqrt(slope) * K
  ! dx=Length, K=coeff, slope=delzs/dx
  ! outflow = dt/dx*vel*m
  ! m(t+1)-m(t) = sum(inflow)-outflow

  ! Approximating Manning's formula
  !
  ! h*W = (m/1000)*dx
  ! W = width of the river and h = height of the river cross-section
  ! vel = K*sqrt(slope)*(A/P)^(2/3)
  ! A = area and P = perimeter of river cross-section
  ! vel = K*sqrt(slope)*(h*W/(2*W+2*h))^(2/3)
  ! Assume W>>h
  ! vel = K*sqrt(slope)*(h/2)^(2/3)
  ! Assume dx=W (approx)
  ! vel = K*sqrt(slope)*(m/2000)^(2/3)

  !--------------------------------------------------------------------
  ! calculate slope
  river_slope(:) = 0.
  do iq = 1,ifull
    if ( river_outdir(iq)>0 ) then  
      iqout = xp(iq,river_outdir(iq))
      river_slope(iq) = max( (zs(iq)-zs(iqout))/grav, 0. )/river_dx(iq)
    end if
  end do

  !--------------------------------------------------------------------
  ! calculate outflow
  outflow(1:ifull) = 0.
  select case(rivermd)
    case(0) ! Miller
      where ( river_outdir(1:ifull)>0 )  
        vel(1:ifull) = rivercoeff*sqrt(river_slope(1:ifull))
        vel(1:ifull) = max( 0.15, min( 5., vel(1:ifull) ) )
        outflow(1:ifull) = (dt/river_dx(1:ifull))*vel(1:ifull)*watbdy(1:ifull) ! (kg/m^2)
      end where
    case(1) ! Manning ( approximated )
      where ( river_outdir(1:ifull)>0 )  
        vel(1:ifull) = rivercoeff*sqrt(river_slope(1:ifull))
        vel(1:ifull) = max( 0.15, min( 5., vel(1:ifull) ) )*max(watbdy(1:ifull)/(2.*1000.),0.)**(2./3.)
        vel(1:ifull) = max( 0.15, vel(1:ifull) )
        outflow(1:ifull) = (dt/river_dx(1:ifull))*vel(1:ifull)*watbdy(1:ifull) ! (kg/m^2)
      end where
    case default
      write(6,*) "ERROR: Unknown option for rivermd=",rivermd
      call ccmpi_abort(-1)
  end select
  outflow(1:ifull) = max( 0., min( watbdy(1:ifull), outflow(1:ifull) ) )
  call bounds(outflow,corner=.true.)


  river_discharge(1:ifull) = outflow(1:ifull)*ds**2/(em(1:ifull)**2*dt*rhow)


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
  
end if


!--------------------------------------------------------------------
! Water losses over land basins

  
! Method for land basins
select case(basinmd)
  case(0)
    ! add water to soil moisture 
    where ( river_outdir(1:ifull)==-1 .and. land(1:ifull) ) ! basin
      tmpry(1:ifull) = watbdy(1:ifull)
    elsewhere ( watbdy(1:ifull)>1000. .and. land(1:ifull) ) ! water exceeds a threshold
      tmpry(1:ifull) = max(watbdy(1:ifull)-1000.,0.)
    elsewhere
      tmpry(1:ifull) = 0.
    end where
    if ( nsib==6 .or. nsib==7 ) then
      ! CABLE
      tmprysave(1:ifull) = tmpry(1:ifull)
      call cableinflow(tmpry)
      soilsink(1:ifull) = (tmpry(1:ifull)-tmprysave(1:ifull))*(1.-sigmu(1:ifull))
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
    end if
    watbdy(1:ifull) = watbdy(1:ifull) + soilsink(1:ifull) ! soilsink is -ve
  case(1)
    ! add water to ocean runoff
    basin_mask(1:ifull) = river_outdir(1:ifull)==-1 .and. land(1:ifull)
    where ( basin_mask(1:ifull) .or. .not.land(1:ifull) )
      watbdy_mask(1:ifull) = watbdy(1:ifull)        ! only includes ocean runoff and basins
    elsewhere
      watbdy_mask(1:ifull) = 0.
    end where
    where ( basin_mask(1:ifull) )
      newwatbdy_mask(1:ifull) = 0.                         ! remove water from basins  
    elsewhere ( watbdy_mask(1:ifull)>1.e-20 )
      newwatbdy_mask(1:ifull) = watbdy_mask(1:ifull) + 0.1 ! small pertubation to allow for increase
    elsewhere
      newwatbdy_mask(1:ifull) = watbdy_mask(1:ifull)  
    end where
    diffwatbdy(1:ifull) = newwatbdy_mask(1:ifull) - watbdy_mask(1:ifull)
    call ccglobal_posneg(diffwatbdy,delpos,delneg)
    if ( delpos>1.e-30 .and. delneg<-1.e-30 ) then
      alph_p = -delneg/delpos
      alph_p = min( sqrt(alph_p), alph_p )
      newwatbdy_mask(1:ifull) = watbdy_mask(1:ifull) + max(0.,diffwatbdy(1:ifull))*alph_p  &
                              + min(0.,diffwatbdy(1:ifull))/max(1.,alph_p)
      where ( basin_mask(1:ifull) .or. .not.land(1:ifull) )
        watbdy(1:ifull) = newwatbdy_mask(1:ifull)
      end where
    end if
  case default
    write(6,*) "ERROR: Unsupported basinmd ",basinmd
    call ccmpi_abort(-1)
end select

watbdy(1:ifull) = max( watbdy(1:ifull), 0. ) ! for rounding errors

return
end subroutine rvrrouter

subroutine water_table_transport

use arrays_m                               ! Atmosphere dyamics prognostic arrays
use cable_ccam                             ! CABLE interface
use cc_mpi                                 ! CC MPI routines
use indices_m                              ! Grid index arrays
use map_m                                  ! Grid map arrays
use newmpar_m                              ! Grid parameters
use nsibd_m                                ! Land-surface arrays
use parm_m                                 ! Model configuration
use soil_m                                 ! Soil and surface data
use soilv_m                                ! Soil parameters
use riverarrays_m                          ! River rarrays
use xyzinfo_m                              ! Grid coordinate arrays

integer iq, n
real, dimension(ifull+iextra,2) :: dumw
real, dimension(ifull+iextra) :: wth_ave, gwwb_min
real, dimension(ifull) :: flux, flux_m, flux_c, wconst
logical, save :: first_call = .true.

if ( wt_transport==1 ) then

  if ( cable_gw_model/=1 ) then
    write(6,*) "ERROR: wt_transport==1 requires cable_gw_model==1"
    call ccmpi_abort(-1)
  end if

  ! calculate average wt height and minimum gw amount
  call calc_wt_ave( wth_ave, gwwb_min, gwdz )

  ! initialise
  if ( first_call ) then
    first_call = .false.
    ! estimate (saturated) hydraulic conductivity (m/s)
    k0(:) = 0.
    do iq = 1,ifull
      if ( land(iq) ) then
        k0(iq) = hyds(isoilm(iq))
      end if
    end do
    ! broadcast
    dumw(1:ifull,1) = k0(1:ifull)
    dumw(1:ifull,2) = gwdz(1:ifull)
    call bounds(dumw(:,1:2))
    k0(ifull+1:ifull+iextra) = dumw(ifull+1:ifull+iextra,1)
    gwdz(ifull+1:ifull+iextra) = dumw(ifull+1:ifull+iextra,2)
  end if

  ! update halo
  dumw(1:ifull,1) = wth_ave(1:ifull)
  dumw(1:ifull,2) = gwwb_min(1:ifull)
  call bounds(dumw(:,1:2))
  wth_ave(ifull+1:ifull+iextra) = dumw(ifull+1:ifull+iextra,1)
  gwwb_min(ifull+1:ifull+iextra) = dumw(ifull+1:ifull+iextra,2)

  ! calculate gradients
  ! gwwb_min also accounts for land-sea mask
  flux(:) = 0.
  flux_m(:) = 0.
  flux_c(:) = 0.
  wconst(:) = sqrt(0.5*tan(3.14159/8.))    ! octagon for eight directions
  call add_flux(flux,flux_m,flux_c,gwwb_min,wth_ave,zs,em,wconst,x,y,z,in)
  call add_flux(flux,flux_m,flux_c,gwwb_min,wth_ave,zs,em,wconst,x,y,z,ie)
  call add_flux(flux,flux_m,flux_c,gwwb_min,wth_ave,zs,em,wconst,x,y,z,is)
  call add_flux(flux,flux_m,flux_c,gwwb_min,wth_ave,zs,em,wconst,x,y,z,iw)
  wconst(:) = sqrt(0.5*tan(3.14159/8.))
  if ( edge_n .and. edge_e ) then
    do n = 1,npan
      iq = ipan + (jpan-1)*il + (n-1)*il**2
      wconst(iq) = 0.
    end do
  end if
  call add_flux(flux,flux_m,flux_c,gwwb_min,wth_ave,zs,em,wconst,x,y,z,ine)
  wconst(:) = sqrt(0.5*tan(3.14159/8.))
  if ( edge_s .and. edge_e ) then
    do n = 1,npan
      iq = ipan + (n-1)*il**2
      wconst(iq) = 0.
    end do
  end if
  call add_flux(flux,flux_m,flux_c,gwwb_min,wth_ave,zs,em,wconst,x,y,z,ise)
  wconst(:) = sqrt(0.5*tan(3.14159/8.))
  if ( edge_s .and. edge_w ) then
    do n = 1,npan
      iq = 1 + (n-1)*il**2
      wconst(iq) = 0.
    end do
  end if
  call add_flux(flux,flux_m,flux_c,gwwb_min,wth_ave,zs,em,wconst,x,y,z,isw)
  wconst(:) = sqrt(0.5*tan(3.14159/8.))
  if ( edge_n .and. edge_w ) then
    do n = 1,npan
      iq = 1 + (jpan-1)*il + (n-1)*il**2
      wconst(iq) = 0.
    end do
  end if
  call add_flux(flux,flux_m,flux_c,gwwb_min,wth_ave,zs,em,wconst,x,y,z,inw)

  ! losing streams (exchange between rivers and GW)

  ! distribute GWwb flux
  call calc_wt_flux( flux, flux_m, flux_c, dt )

end if

return
end subroutine water_table_transport

subroutine add_flux(flux,flux_m,flux_c,gwwb_min,wth,zs,em,wconst,x,y,z,dir)

use const_phys                             ! Physical constants
use newmpar_m                              ! Grid parameters
use parm_m                                 ! Model configuration
use riverarrays_m                          ! River rarrays

integer iq
integer, dimension(ifull), intent(in) :: dir
real wth_del, wth_ave, w, t, dx, f, k0_ave, slope
real wth_max, flux_add, dr
real, dimension(ifull), intent(inout) :: flux, flux_m, flux_c
real, dimension(ifull+iextra), intent(in) :: wth, gwwb_min, zs
real, dimension(ifull+iextra), intent(in) :: em
real, dimension(ifull), intent(in) :: wconst
real(kind=8), dimension(ifull+iextra), intent(in) :: x, y, z
real(kind=8) dotprod

! Calculation based on Fan et al Water Table Observations doi: 10.1029/2006JD008111

do iq = 1,ifull
  dotprod = x(iq)*x(dir(iq)) + y(iq)*y(dir(iq)) + z(iq)*z(dir(iq))
  dr = real( acos( max( min( dotprod, 1._8 ), -1._8 ) ) )*rearth ! distance between grid cells (m)
  dx = 2.*ds/(em(iq)+em(dir(iq)))             ! width of grid cell (m)
  w = dx*wconst(iq)                           ! width of cell boundary (m)
  k0_ave = 0.5*(k0(iq)+k0(dir(iq)))           ! m/s
  slope = abs(zs(iq)-zs(dir(iq)))/(grav*dr)   ! m/m
  wth_del = wth(iq) - wth(dir(iq))            ! m
  wth_ave = 0.5*(wth(iq)+wth(dir(iq)))        ! m
  wth_max = 0.5*(zs(iq)+zs(dir(iq)))/grav     ! m
  f = 120./(1.+150.*slope)
  f = max( f, 5. )                            ! m
  t = k0_ave*f*exp(min(wth_ave-wth_max,0.)/f) ! m2/s
  flux_add = w*t*wth_del/dr                   ! m3/s
  ! limit flux based on avaliable GW
  flux_add = max( min( flux_add, gwwb_min(iq)*gwdz(iq)/dt ), -gwwb_min(dir(iq))*gwdz(dir(iq))/dt )
  if ( gwdz(iq)>0. .and. gwdz(dir(iq))>0. ) then
    flux(iq) = flux(iq) + flux_add/gwdz(iq)                ! m2/s
    flux_m(iq) = flux_m(iq) + w*t/dr/gwdz(iq)              ! m/2
    flux_c(iq) = flux_c(iq) - w*t*wth(dir(iq))/dr/gwdz(iq) ! m2/s
  end if
end do

return
end subroutine add_flux

subroutine wilms(INDEX1,INDEX2)

use arrays_m
use cc_mpi
use const_phys
use indices_m
use map_m
use newmpar_m, only : il
use nsibd_m
use parm_m
use riverarrays_m, only : watbdy_cc => watbdy, river_dx_cc => river_dx, river_discharge_cc => river_discharge
use sflux_m
use soil_m
use soilsnow_m
use soilv_m

implicit none

integer, intent(in) :: INDEX1, INDEX2
integer idx1, idx2, count_keys
integer key1, key2, val1, val2
integer iq, nb1, nb2, count_nb
real, dimension(INDEX1,INDEX2) :: watbdy, discharge_minimum, inflow
real, dimension(INDEX1,INDEX2) :: river_discharge, discharge_factor, ocean_discharge
real, dimension(INDEX1,INDEX2) :: grid_area
real, dimension(0:INDEX1+1,0:INDEX2+1) :: z_dem, watbdy_height, outflow
real river_dx, slope
real, parameter :: n = 0.025    ! manning constant
real, parameter :: rho = 1000.  ! density of water

! convert CCAM grid to INDEX1, INDEX2
watbdy(1:INDEX1,1:INDEX2)    = reshape( watbdy_cc(1:INDEX1*INDEX2),  (/ INDEX1, INDEX2 /) )
grid_area(1:INDEX1,1:INDEX2) = reshape( (ds/em(1:INDEX1*INDEX2))**2, (/ INDEX1, INDEX2 /) )
watbdy(1:INDEX1,1:INDEX2)    = watbdy(1:INDEX1,1:INDEX2)*grid_area(1:INDEX1,1:INDEX2)
! copy halo data into z_dem (zs halo already updated during CCAM initialisation in indata.f90)
call ccreshape(z_dem,zs)

discharge_minimum = 0.0
inflow = 0.0
outflow = 0.0
river_discharge = 0.0
discharge_factor = 0.0
ocean_discharge = 0.0

!add runoff (allready as a volume) to waterbody for this day
do idx1=1,INDEX1
    do idx2=1,INDEX2
        if (z_dem(idx1,idx2) < 0) z_dem(idx1,idx2) = 0.0
        ! runoff has already been included in sflux.f90
        ! watbdy(idx1,idx2) =  watbdy(idx1,idx2) + runoffavg(day,idx1,idx2)
        !calculate the new waterbody height by dividing the volume of water through the area of that cell
        watbdy_height(idx1,idx2) = watbdy(idx1,idx2)/grid_area(idx1,idx2) 
    end do!do idx2
end do!do idx1


! need to update halo between CCAM processes since it uses val1 and val2
call updatehalo(watbdy_height)


!calculate volume of water that leaves a cell:
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
do count_keys=1,total_keys

    key1 = indices_list(count_keys,1) 
    key2 = indices_list(count_keys,2) 
                
    discharge_minimum(key1,key2) = 0.15*riverwidth(key1,key2)*watbdy_height(key1,key2) 


    discharge_factor(key1, key2) = (watbdy_height(key1,key2) * riverwidth(key1,key2) * &
                                  &(watbdy_height(key1,key2) * riverwidth(key1,key2) / &
                                  & (riverwidth(key1,key2) + 2. * watbdy_height(key1,key2))) ** ( &
                                  2.0 / 3.0)) / n
                                                    
    val1 = river_outloc(1,key1,key2)
    val2 = river_outloc(2,key1,key2)
    
    ! replace the following with the calculation for river_dx
    !temp_lon1 = x(key1)
    !temp_lat1 = y(key2)
    !temp_lon2 = x(val1)
    !temp_lat2 = y(val2)
    !call haversine(temp_lon1,temp_lat1,temp_lon2,temp_lat2,river_dx)
    river_dx = river_dx_cc(key1+(key2-1)*il)
    
   
    !calculate the slope
    slope = (z_dem(key1,key2) + watbdy_height(key1,key2)/rho - z_dem(val1,val2) - &
            & watbdy_height(val1,val2) / rho) / river_dx  
                
                
    !calculate the volume of water that leaves the cell
    if (slope>0.) then
        river_discharge(key1,key2) = discharge_factor(key1,key2) * sqrt(slope)
    else 
        river_discharge(key1,key2) = discharge_minimum(key1,key2)
    endif
                
    outflow(key1,key2) = river_discharge(key1,key2)*dt !outflow in m^3 for the entire timestep
                
    !if the amount that should flow out is larger than volume of water cell has to give, set outflow
    !to the amount the cell has to give and recalculate the discharge
    if (outflow(key1,key2) > watbdy(key1,key2)) then 
        outflow(key1,key2) = watbdy(key1,key2)
        river_discharge(key1,key2) = outflow(key1,key2)/dt
    endif
                
                
    !now we know how much water leaves the cell
end do!count_keys
            

! need to update halo between CCAM processes since it uses val1 and val2
call updatehalo(outflow)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!proceed to check what the inflows are
do count_keys=1,count_inflows !now we calculate the inflows
               
    key1 = river_inflow_keys(1,count_keys)
    key2 = river_inflow_keys(2,count_keys)
    !in_arr = (/ key1, key2 /)
    !call neighbours (in_arr, neighbour_list)
    
    iq = key1 + (key2-1)*il    
    
    do count_nb = 1,8
        
        !nb1 = neighbour_list(1,count_nb)
        !nb2 = neighbour_list(2,count_nb)
        !            
        !rv_out1 = river_outloc(1, nb1, nb2)
        !rv_out2 = river_outloc(2, nb1, nb2)
        !   
        !if (rv_out1==key1 .and. rv_out2==key2) then
        !   
        !    inflow(key1, key2) = inflow(key1, key2) + outflow(nb1, nb2)
        !                
        !endif

        if ( river_inflow(iq,count_nb) ) then
            select case(count_nb)
                case(1) !n
                    nb1 = key1
                    nb2 = key2 + 1
                case(2) !e    
                    nb1 = key1 + 1
                    nb2 = key2
                case(3) !s
                    nb1 = key1
                    nb2 = key2 - 1
                case(4) !w
                    nb1 = key1 - 1
                    nb2 = key2
                case(5) !ne
                    nb1 = key1 + 1
                    nb2 = key2 + 1
                case(6) !se
                    nb1 = key1 + 1
                    nb2 = key2 - 1
                case(7) !sw
                    nb1 = key1 - 1
                    nb2 = key2 - 1
                case(8) !nw
                    nb1 = key1 - 1
                    nb2 = key2 + 1
            end select      
                
            inflow(key1, key2) = inflow(key1, key2) + outflow(nb1, nb2)

        end if     
            
    end do !(count_nb = 1,8)

end do !(do count_keys=1,count_inflows)
            
!update waterbody   
do idx1=1,INDEX1
    do idx2=1,INDEX2
        watbdy(idx1,idx2) =  watbdy(idx1,idx2) + inflow(idx1,idx2)-outflow(idx1,idx2)
    end do!do idx2
end do!do idx1
   
! The following lines are handled by the basinmd options in the subroutine rvrrouter (above)
!!now handle instances where a grid point doesn't have an output but still receives water...this may include some coastal points
!!too.
!do count_keys=1,count_inflows 
!    key1 = river_inflow_keys(1,count_keys)
!    key2 = river_inflow_keys(2,count_keys)
!    if (river_outloc(1, key1, key2) == 0) then !if no output loc, dump into the ocean
!        watbdy_ocn(key1, key2)=watbdy(key1, key2) 
!        watbdy(key1, key2) = 0.0    ! now set watbdy of this cell to 0
!                    
!    endif
!end do !(do count_keys=1,count_inflows )

! The following lines are handled by sflux.f90
!!now handle instances where the grid point is located on the coast
!!each of the coast_keys are located on the coast of Africa.  Determined with python's basemap.
!do coast_keys = 1, 2450
!    key1 = coast_list(coast_keys,1)
!    key2 = coast_list(coast_keys,2)
!    !if we did not give water to this coastal point during the previous iteration we'll do the following.  
!    !However, if we did give water from a point with no outflow the job has already been done since each grid point donates only
!    !to one neighbour.
!    if (watbdy_ocn(key1, key2)==0.0) then 
!        watbdy_ocn(key1, key2) = watbdy(key1, key2)
!        ocean_discharge(key1, key2) = watbdy_ocn(key1, key1)/dt !determine how fast it flows into the ocean.
!        watbdy(key1, key2) = 0.0 ! I'm giving all of the water away to the ocean so need to set the watbdy to zero for the point
!                                 ! that gave its water. (Probably not such a great idea?)
!                    
!    endif !(if (watbdy_ocn(key1, key2)==0.0))
!end do !(do coast_keys = 1, 2450)


! convert arrays back for CCAM
watbdy(1:INDEX1,1:INDEX2)           = watbdy(1:INDEX1,1:INDEX2)/grid_area(1:INDEX1,1:INDEX2)
watbdy_cc(1:INDEX1*INDEX2)          = reshape( watbdy(1:INDEX1,1:INDEX2),          (/ INDEX1*INDEX2 /) )
river_discharge_cc(1:INDEX1*INDEX2) = reshape( river_discharge(1:INDEX1,1:INDEX2), (/ INDEX1*INDEX2 /) )

return
end subroutine wilms

subroutine ccreshape(outdata,indata)

use indices_m
use newmpar_m

implicit none

integer i, j, iq
real, dimension(ifull+iextra), intent(in) :: indata
real, dimension(0:il+1,0:jl+1), intent(out) :: outdata

outdata(1:il,1:jl) = reshape( indata(1:il*jl), (/ il, jl /) )
do i = 1,il
  iq = i
  outdata(i,0) = indata(is(iq))
  iq = i + (jl-1)*il
  outdata(i,jl+1) = indata(in(iq))
end do
do j = 1,jl
  iq = 1 + (j-1)*il
  outdata(0,j) = indata(iw(iq))
  iq = il + (j-1)*il
  outdata(il+1,j) = indata(ie(iq))
end do
iq = 1
outdata(0,0) = indata(isw(iq))
iq = il
outdata(il+1,0) = indata(ise(iq))
iq = 1 + (jl-1)*il
outdata(0,jl+1) = indata(inw(iq))
iq = il + (jl-1)*il
outdata(il+1,jl+1) = indata(ine(iq))

return
end subroutine ccreshape

subroutine updatehalo(val)

use cc_mpi
use newmpar_m

implicit none

real, dimension(0:il+1,0:jl+1), intent(inout) :: val
real, dimension(ifull+iextra) :: val_cc

val_cc(1:il*jl) = reshape( val(1:il,1:jl), (/ il*jl /) )
call bounds(val_cc,corner=.true.)

call ccreshape(val,val_cc)

return
end subroutine updatehalo

end module river
