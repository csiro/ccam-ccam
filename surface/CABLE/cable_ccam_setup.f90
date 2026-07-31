! Conformal Cubic Atmospheric Model
    
! Copyright 2026 Commonwealth Scientific Industrial Research Organisation (CSIRO)
    
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

! Subroutines to configue CABLE in CCAM
    
module cable_ccam_setup

use cable_ccam_common

implicit none

private
public loadcbmparm, cbmparm

!! Expected CABLE input data version
!real, save :: cable_version = 6608. ! expected version id for input data
!                                    ! 6608 includes update to ice albedo
!                                    ! 3939 had fixes for soil albedo

contains

subroutine loadcbmparm(fveg,fphen,casafile, &
                       ivs,svs,vlin,casapoint,greenup,fall,phendoy1)

use cc_mpi
use newmpar_m
use parm_m

integer n
integer, dimension(ifull,maxtile), intent(out) :: ivs
real, dimension(ifull,maxtile), intent(out) :: svs,vlin
real, dimension(ifull,maxtile), intent(out) :: casapoint
real, dimension(ifull) :: svs_sum
real cableformat
integer, dimension(271,mxvt), intent(out) :: greenup, fall, phendoy1
character(len=*), intent(in) :: fveg,fphen,casafile

! read CABLE biome and LAI data
if ( myid==0 ) then
  write(6,*) "Reading tiled surface data for CABLE"
  call vegta(ivs,svs,vlin,fveg,cableformat)
else
  call vegtb(ivs,svs,vlin,cableformat)
end if
svs_sum = sum(svs,dim=2)
do n = 1,maxtile
  svs(:,n) = svs(:,n)/svs_sum(:)
end do

if ( abs(cableformat-1.)<1.e-20 ) then
  if ( myid==0 ) write(6,*) "Procesing CSIRO PFTs"    
else
  if ( myid==0 ) write(6,*) "Processing IGBP and converting to CSIRO PFTs"
  call convertigbp(ivs,svs,vlin)
end if

if ( ccycle==0 ) then
  casapoint(:,:) = 0.
  greenup(:,:)  = -50
  fall(:,:)     = 367
  phendoy1(:,:) = 2  
else
  call casa_readpoint(casafile,casapoint)         ! read point sources
  call casa_readphen(fphen,greenup,fall,phendoy1) ! read MODIS leaf phenology
end if

if ( any(vlin>10.) ) then
  write(6,*) "ERROR: LAI is out of range"
  write(6,*) "vlin ",maxval(vlin)
  call ccmpi_abort(-1)
end if

return
end subroutine loadcbmparm

! legacy code for IGBP vegetation classes    
subroutine convertigbp(ivs,svs,vlin)

use cc_mpi
use const_phys
use latlong_m
use newmpar_m
use soil_m

integer, dimension(ifull,maxtile), intent(inout) :: ivs
integer, dimension(1) :: pos
integer iq, n, ipos, iv
real, dimension(ifull,maxtile), intent(inout) :: svs,vlin
real, dimension(18) :: newlai
real, dimension(18) :: newgrid
real fc3, fc4, ftu, fg3, fg4, clat, nsum
real xp

do iq = 1,ifull
  if ( land(iq) ) then
    newgrid(:)  = 0.
    newlai(:) = 0.      
    clat = rlatt(iq)*180./pi
    ! grass
    if (abs(clat)>50.5) then
      fg3=0.
      fg4=0.
    else if (abs(clat)>49.5) then
      xp=abs(clat)-49.5
      fg3=(1.-xp)*0.9
      fg4=(1.-xp)*0.1
    else if (abs(clat)>40.5) then
      fg3=0.9
      fg4=0.1
    else if (abs(clat)>39.5) then
      xp=abs(clat)-39.5
      fg3=(1.-xp)*0.8+xp*0.9
      fg4=(1.-xp)*0.2+xp*0.1
    else if (abs(clat)>30.5) then
      fg3=0.8
      fg4=0.2
    else if (abs(clat)>29.5) then
      xp=abs(clat)-29.5
      fg3=(1.-xp)*0.5+xp*0.8
      fg4=(1.-xp)*0.5+xp*0.2
    else if (abs(clat)>25.5) then
      fg3=0.5
      fg4=0.5
    else if (abs(clat)>24.5) then
      xp=abs(clat)-24.5
      fg3=(1.-xp)*0.05+xp*0.5
      fg4=(1.-xp)*0.95+xp*0.5
    else
      fg3=0.05
      fg4=0.95
    end if
    ftu=1.-fg3-fg4
    ! crops
    if (abs(clat)>40.5) then
      fc3=1.
    else if (abs(clat)>39.5) then
      xp=abs(clat)-39.5
      fc3=(1.-xp)*0.9+xp
    else if (abs(clat)>30.5) then
      fc3=0.9
    else if (abs(clat)>29.5) then
      xp=abs(clat)-29.5
      fc3=(1.-xp)*0.7+xp*0.9
    else
      fc3=0.7
    end if
    fc4=1.-fc3
    do n = 1,maxtile
      select case (ivs(iq,n))
        case (1,2,3,4,11)
          newgrid(ivs(iq,n))=newgrid(ivs(iq,n))+svs(iq,n)
          newlai(ivs(iq,n))=newlai(ivs(iq,n))+svs(iq,n)*vlin(iq,n)
        case (5)
          if (abs(clat)>25.5) then
            newgrid(1)=newgrid(1)+svs(iq,n)*0.5
            newlai(1)=newlai(1)+svs(iq,n)*0.5*vlin(iq,n)
            newgrid(4)=newgrid(4)+svs(iq,n)*0.5
            newlai(4)=newlai(4)+svs(iq,n)*0.5*vlin(iq,n)
          else if (abs(clat)>24.5) then
            xp=abs(clat)-24.5
            newgrid(1)=newgrid(1)+svs(iq,n)*0.5*xp
            newlai(1)=newlai(1)+svs(iq,n)*0.5*vlin(iq,n)*xp
            newgrid(4)=newgrid(4)+svs(iq,n)*(1.-0.5*xp)
            newlai(4)=newlai(4)+svs(iq,n)*vlin(iq,n)*(1.-0.5*xp)
          else
            newgrid(4)=newgrid(4)+svs(iq,n)
            newlai(4)=newlai(4)+svs(iq,n)*vlin(iq,n)
          end if
        case (6)
          newgrid(5)=newgrid(5)+svs(iq,n)*0.8
          newlai(5)=newlai(5)+svs(iq,n)*0.8*vlin(iq,n)
          newgrid(6)=newgrid(6)+svs(iq,n)*0.2*fg3
          newlai(6)=newlai(6)+svs(iq,n)*0.2*fg3*vlin(iq,n)
          newgrid(7)=newgrid(7)+svs(iq,n)*0.2*fg4
          newlai(7)=newlai(7)+svs(iq,n)*0.2*fg4*vlin(iq,n)
          newgrid(8)=newgrid(8)+svs(iq,n)*0.2*ftu
          newlai(8)=newlai(8)+svs(iq,n)*0.2*ftu*vlin(iq,n)
        case (7)
          newgrid(5)=newgrid(5)+svs(iq,n)*0.2
          newlai(5)=newlai(5)+svs(iq,n)*0.2*vlin(iq,n)
          newgrid(6)=newgrid(6)+svs(iq,n)*0.8*fg3
          newlai(6)=newlai(6)+svs(iq,n)*0.8*fg3*vlin(iq,n)
          newgrid(7)=newgrid(7)+svs(iq,n)*0.8*fg4
          newlai(7)=newlai(7)+svs(iq,n)*0.8*fg4*vlin(iq,n)
          newgrid(8)=newgrid(8)+svs(iq,n)*0.8*ftu
          newlai(8)=newlai(8)+svs(iq,n)*0.8*ftu*vlin(iq,n)
        case (8)
          if (abs(clat)>40.5) then
            newgrid(1)=newgrid(1)+svs(iq,n)*0.4
            newlai(1)=newlai(1)+svs(iq,n)*0.4*vlin(iq,n)
          else if (abs(clat)>39.5) then
            xp=abs(clat)-39.5
            newgrid(1)=newgrid(1)+svs(iq,n)*0.4*xp
            newlai(1)=newlai(1)+svs(iq,n)*vlin(iq,n)*0.4*xp
            newgrid(18)=newgrid(18)+svs(iq,n)*0.4*(1.-xp)
            newlai(18)=newlai(18)+svs(iq,n)*vlin(iq,n)*0.4*(1.-xp)
          else
            newgrid(18)=newgrid(18)+svs(iq,n)*0.4
            newlai(18)=newlai(18)+svs(iq,n)*0.4*vlin(iq,n)
          end if
          newgrid(6)=newgrid(6)+svs(iq,n)*0.6*fg3
          newlai(6)=newlai(6)+svs(iq,n)*0.6*fg3*vlin(iq,n)
          newgrid(7)=newgrid(7)+svs(iq,n)*0.6*fg4
          newlai(7)=newlai(7)+svs(iq,n)*0.6*fg4*vlin(iq,n)
          newgrid(8)=newgrid(8)+svs(iq,n)*0.6*ftu
          newlai(8)=newlai(8)+svs(iq,n)*0.6*ftu*vlin(iq,n)
        case (9)
          if (abs(clat)>40.5) then
            newgrid(1)=newgrid(1)+svs(iq,n)*0.1
            newlai(1)=newlai(1)+svs(iq,n)*0.1*vlin(iq,n)
          else if (abs(clat)>39.5) then
            xp=abs(clat)-39.5
            newgrid(1)=newgrid(1)+svs(iq,n)*0.1*xp
            newlai(1)=newlai(1)+svs(iq,n)*vlin(iq,n)*0.1*xp
            newgrid(18)=newgrid(18)+svs(iq,n)*0.1*(1.-xp)
            newlai(18)=newlai(18)+svs(iq,n)*vlin(iq,n)*0.1*(1.-xp)
          else
            newgrid(18)=newgrid(18)+svs(iq,n)*0.1
            newlai(18)=newlai(18)+svs(iq,n)*0.1*vlin(iq,n)
          end if
          newgrid(6)=newgrid(6)+svs(iq,n)*0.9*fg3
          newlai(6)=newlai(6)+svs(iq,n)*0.9*fg3*vlin(iq,n)
          newgrid(7)=newgrid(7)+svs(iq,n)*0.9*fg4
          newlai(7)=newlai(7)+svs(iq,n)*0.9*fg4*vlin(iq,n)
          newgrid(8)=newgrid(8)+svs(iq,n)*0.9*ftu
          newlai(8)=newlai(8)+svs(iq,n)*0.9*ftu*vlin(iq,n)
        case (10)
          newgrid(6)=newgrid(6)+svs(iq,n)*fg3
          newlai(6)=newlai(6)+svs(iq,n)*fg3*vlin(iq,n)
          newgrid(7)=newgrid(7)+svs(iq,n)*fg4
          newlai(7)=newlai(7)+svs(iq,n)*fg4*vlin(iq,n)
          newgrid(8)=newgrid(8)+svs(iq,n)*ftu
          newlai(8)=newlai(8)+svs(iq,n)*ftu*vlin(iq,n)
        case (12,14)
          newgrid(9)=newgrid(9)+svs(iq,n)*fc3
          newlai(9)=newlai(9)+svs(iq,n)*fc3*vlin(iq,n)
          newgrid(10)=newgrid(10)+svs(iq,n)*fc4
          newlai(10)=newlai(10)+svs(iq,n)*fc4*vlin(iq,n)
        case (13)
          newgrid(15)=newgrid(15)+svs(iq,n)
          newlai(15)=newlai(15)+svs(iq,n)*vlin(iq,n)
        case (15)
          newgrid(17)=newgrid(17)+svs(iq,n)
          newlai(17)=newlai(17)+svs(iq,n)*vlin(iq,n)
        case (16)
          newgrid(14)=newgrid(14)+svs(iq,n)
          newlai(14)=newlai(14)+svs(iq,n)*vlin(iq,n)
        case (17)
          newgrid(16)=newgrid(16)+svs(iq,n)
          newlai(16)=newlai(16)+svs(iq,n)*vlin(iq,n)
        case DEFAULT
          write(6,*) "ERROR: Land-type/lsmask mismatch at myid,iq,ivs,land=",myid,iq,ivs(iq,n),land(iq)
          call ccmpi_abort(-1)
      end select
    end do
    where ( newgrid(:)>0. )
      newlai(:) = newlai(:)/newgrid(:)
    end where
    ipos = count(newgrid(:)>0.)
    do while ( ipos>maxtile )
      pos = minloc(newgrid(:), newgrid(:)>0.)
      newgrid(pos(1)) = 0.
      nsum = sum(newgrid(:))
      newgrid(:) = newgrid(:)/nsum
      ipos = count(newgrid(:)>0.)
    end do    
    !do while ( any(newgrid(:)<minfrac.and.newgrid(:)>0.) )
    !  pos = minloc(newgrid(:), newgrid(:)>0.)
    !  newgrid(pos(1)) = 0.
    !  nsum = sum(newgrid(:))
    !  newgrid(:) = newgrid(:)/nsum
    !end do

    n = 0
    ivs(iq,:) = 0
    svs(iq,:) = 0.
    vlin(iq,:) = 0.
    do iv = 1,18
      if ( newgrid(iv)>0. ) then
        n = n + 1
        ivs(iq,n)  = iv
        svs(iq,n)  = newgrid(iv)
        vlin(iq,n) = newlai(iv)
      end if
    end do

  end if
end do

return
end subroutine convertigbp

subroutine cbmparm(ivs,svs,vlin,casapoint,greenup,fall,phendoy1,fcasapft)

use carbpools_m
use cc_mpi
use const_phys
use infile
use latlong_m
use newmpar_m
use nsibd_m
use parm_m
use pbl_m
use riverarrays_m
use sigs_m
use soil_m
use soilsnow_m
use soilv_m
use vegpar_m
  
integer, dimension(ifull,maxtile), intent(in) :: ivs
integer, dimension(271,mxvt), intent(in) :: greenup, fall, phendoy1
integer(kind=4), dimension(:), allocatable, save :: Iwood
integer(kind=4), dimension(:,:), allocatable, save :: disturbance_interval
integer i,iq,n,k,ipos,ilat,ivp,is,ie,iis,iie
integer landcount
integer(kind=4) mp_POP
real ivmax, landsum
real, dimension(ifull,maxtile), intent(in) :: svs,vlin
real, dimension(ifull,5), intent(in) :: casapoint
real, dimension(ifull,2) :: albsoilsn
real, dimension(ifull) :: dummy_pack, albsoil
real, dimension(0:maxtile) :: stat_count, global_stat_count
real(kind=8), dimension(:), allocatable :: dummy_unpack
integer :: tile, popcount
character(len=*), intent(in) :: fcasapft

if ( myid==0 ) write(6,*) "Initialising CABLE"

if ( cbm_ms/=ms ) then
  write(6,*) "ERROR: CABLE and CCAM soil levels do not match"
  call ccmpi_abort(-1)
end if

!POP_NPATCH = NPATCH
!POP_NLAYER = NLAYER
!POP_NCOHORT = NCOHORT_MAX
!POP_HEIGHT_BINS = HEIGHT_BINS
!POP_NDISTURB = NDISTURB
!POP_AGEMAX = AGEMAX

mvtype = mxvt
mstype = mxst
icycle = ccycle

if ( cbm_ms /= ms ) then
  write(6,*) "ERROR: Mismatch in number of soil layers between CABLE and CCAM"
  write(6,*) "cbm_ms = ",cbm_ms
  write(6,*) "ms = ",ms
  call ccmpi_abort(-1)
end if

if ( myid==0 .and. nmaxpr==1 ) write(6,*) "-> Define CABLE and CASA CNP arrays"

! default values (i.e., no land)  
ivegt = 0
albsoilsn = 0.08  
albsoil   = 0.08
albvisdir = 0.08
albvisdif = 0.08
albnirdir = 0.08
albnirdif = 0.08
zolnd     = 0.

! calculate CABLE vector length
allocate( tdata(ntiles) )
mp_max = 0
do tile = 1,ntiles
  tdata(tile)%mp = 0
  tdata(tile)%np = 0
  is = 1 + (tile-1)*imax
  ie = tile*imax
  do iq = is,ie
    if ( land(iq) ) then
      landcount = count(svs(iq,:)>0.)
      tdata(tile)%mp = tdata(tile)%mp + landcount
      if ( landcount==0 ) then
        write(6,*) "ERROR: No CABLE tiles assigned to land point: myid,iq,tile",myid,iq,tile
        call ccmpi_abort(-1)
      end if
      landsum = sum(svs(iq,:))
      if ( landsum<0.99 .or. landsum>1.01 ) then
        write(6,*) "ERROR: CABLE tiles do not sum to 1."
        write(6,*) "myid,iq,tile,landsum ",myid,iq,tile,landsum
        call ccmpi_abort(-1)
      end if
    end if
  end do
  
  if ( tile==1 ) then
    mp_global = tdata(1)%mp  
    tdata(1)%toffset = 0
  else  
    mp_global = mp_global + tdata(tile)%mp
    tdata(tile)%toffset = tdata(tile-1)%toffset + tdata(tile-1)%mp
  end if
  
  mp_max = max( mp_max, tdata(tile)%mp )
  
end do

ktau_gl = 900
kend_gl = 999

! if CABLE is present on this processor, then start allocating arrays
! Write messages here in case myid==0 has no land-points (mp_global==0)
if ( myid==0 ) then
  if ( nmaxpr==1 ) write(6,*) "-> Allocating CABLE and CASA CNP arrays"
  if ( soil_struc==1 ) write(6,*) "-> Using SLI soil model"
  if ( ccycle==0 ) then
    write(6,*) "-> Using CABLE without carbon cycle"
  else if ( ccycle==1 ) then
    write(6,*) "-> Using CASA C"
  else if ( ccycle==2 ) then
    write(6,*) "-> Using CASA CN"
  else
    write(6,*) "-> Using CASA CNP"
  end if
end if

do tile = 1,ntiles
  allocate( tdata(tile)%tind(maxtile,2) )
  tdata(tile)%tind(:,1) = 1
  tdata(tile)%tind(:,2) = 0
  allocate( tdata(tile)%pind(maxtile,2) )
  tdata(tile)%pind(:,1) = 1
  tdata(tile)%pind(:,2) = 0
  tdata(tile)%maxnb = 0
end do

allocate( air(ntiles), bgc(ntiles), canopy(ntiles), met(ntiles) )
allocate( bal(ntiles), rad(ntiles), rough(ntiles), soil(ntiles) )
allocate( ssnow(ntiles), sum_flux(ntiles), veg(ntiles), climate(ntiles) )
allocate( casabiome(ntiles), casapool(ntiles), casaflux(ntiles), casamet(ntiles) )
allocate( casabal(ntiles), phen(ntiles) )
allocate( POP(ntiles) )

if ( mp_global>0 ) then

  do tile = 1,ntiles
    if ( tdata(tile)%mp>0 ) then
      allocate( tdata(tile)%sv(tdata(tile)%mp) )
      allocate( tdata(tile)%vl2(tdata(tile)%mp) )
      allocate( tdata(tile)%cveg(tdata(tile)%mp) )
      allocate( tdata(tile)%old_sv(tdata(tile)%mp) )
      allocate( tdata(tile)%old_cveg(tdata(tile)%mp) )
      call alloc_cbm_var(air(tile), tdata(tile)%mp)
      call alloc_cbm_var(bgc(tile), tdata(tile)%mp)
      call alloc_cbm_var(canopy(tile), tdata(tile)%mp)
      call alloc_cbm_var(met(tile), tdata(tile)%mp)
      call alloc_cbm_var(bal(tile), tdata(tile)%mp)
      call alloc_cbm_var(rad(tile), tdata(tile)%mp)
      call alloc_cbm_var(rough(tile), tdata(tile)%mp)
      call alloc_cbm_var(soil(tile), tdata(tile)%mp)
      call alloc_cbm_var(ssnow(tile), tdata(tile)%mp)
      call alloc_cbm_var(sum_flux(tile), tdata(tile)%mp)
      call alloc_cbm_var(veg(tile), tdata(tile)%mp)
    end if
    climate(tile)%nyear_average = 20
    climate(tile)%nday_average = 31
  end do ! tile

  ! Cable configuration
  select case( cable_potev )
    case(0)
      cable_user%ssnow_POTEV = "P-M" ! Penman Monteith
    case default
      ! default for offline and ACCESS
      cable_user%ssnow_POTEV = "HDM" ! default Humidity Deficit
  end select    
  cable_user%MetType = "defaul"    ! only 6 characters for "default"
  cable_user%diag_soil_resp = "ON"
  cable_user%leaf_respiration = "ON"
  cable_user%run_diag_level = "NONE"
  cable_user%consistency_check = .false.
  cable_user%logworker = .false.
  cable_user%l_new_roughness_soil = cable_roughness==1
  cable_user%l_rev_corr = .true.
  cable_user%gw_model = cable_gw_model==1
  cable_user%soil_thermal_fix = cable_thermfix==1
  cable_user%call_climate = .false.
  cable_user%phenology_switch = "MODIS"
  !cable_user%finite_gm = .false.
  cable_user%call_pop = cable_pop==1
  cable_user%or_evap = .false.
  select case ( soil_struc )
    case(1)
      cable_user%soil_struc = "sli"  
    case default
      cable_user%soil_struc = "default"
  end select
  select case ( fwsoil_switch )
    case(3)
      ! default for ACCESS  
      cable_user%fwsoil_switch = "Haverd2013"
    case(2)
      cable_user%fwsoil_switch = "Lai and Ktaul 2000"  
    case(1)
      cable_user%fwsoil_switch = "non-linear extrapola" ! only 20 characters for "non-linear extrapolation"
    case default
      ! default for offline
      cable_user%fwsoil_switch = "standard"      
  end select
  select case ( gs_switch )  
    case(1)
      ! default for ACCESS  
      cable_user%gs_switch = "medlyn"  
    case default
      ! default for offline
      cable_user%gs_switch = "leuning"
  end select
  cable_user%litter = cable_litter==1
  select case ( progvcmax )
    case(2)
      cable_user%vcmax = "Walker2014"        
    case default
      cable_user%vcmax = "standard"    
  end select
  select case ( smrf_switch )
    case(1)
      cable_user%SMRF_NAME = "CASA-CNP"
    case(2)
      cable_user%SMRF_NAME = "SOILN"
    case(3)
      cable_user%SMRF_NAME = "TRIFFID"
    case(4)
      cable_user%SMRF_NAME = "Trudinger2" ! 2016
    case(5)
      cable_user%SMRF_NAME = "DAMM"
    case default
      cable_user%SMRF_NAME = "Trudinger2" ! 2016
  end select
  select case ( strf_switch )
    case(1)
      cable_user%STRF_NAME = "CASA-CNP"
    case(2)
      cable_user%STRF_NAME = "K1995"
    case(3)
      cable_user%STRF_NAME = "PnET-CN"
    case(4)
      cable_user%STRF_NAME = "LT1994"
    case(5)
      cable_user%STRF_NAME = "DAMM"
    case default
      cable_user%STRF_NAME = "LT1994"
  end select
  cable_user%l_new_runoff_speed = cable_runoff==1
  cable_user%l_new_reduce_soilevp = cable_soilevap==1
  !kwidth_gl = nint(dt) ! MJT notes - what happens when the timestep is less than a second?
  !if ( kwidth_gl==0) then
  !  write(6,*) "ERROR: Time-step dt<1 is too small for CABLE"
  !  call ccmpi_abort(-1)
  !end if

  do tile = 1,ntiles
    is = 1 + (tile-1)*imax
    ie = tile*imax
      
    ! soil layer thickness
    if ( tdata(tile)%mp>0 ) then
      soil(tile)%zse        = real(zse,8) ! soil layer thickness
      soil(tile)%zshh(1)    = 0.5_8*soil(tile)%zse(1)
      soil(tile)%zshh(2:cbm_ms) = 0.5_8*(soil(tile)%zse(1:cbm_ms-1) + soil(tile)%zse(2:cbm_ms))
      soil(tile)%zshh(cbm_ms+1) = 0.5_8*soil(tile)%zse(cbm_ms)
    end if  

    ! allocate map between gridbox tile and vegetation tile
    allocate(tdata(tile)%tmap(imax,maxtile))
    tdata(tile)%tmap(:,:) = .false.
  
    ! define indices for vegetation tiles
    ipos = 0
    do n = 1,maxtile
      tdata(tile)%tind(n,1) = ipos + 1
      do iq = is,ie
        if ( land(iq) ) then
          if ( svs(iq,n)>0. ) then
            ipos = ipos + 1
            tdata(tile)%tmap(iq-is+1,n) = .true.
          end if
        end if
      end do
      tdata(tile)%tind(n,2) = ipos
    end do
  
    if ( ipos/=tdata(tile)%mp ) then
      write(6,*) "ERROR: Internal memory allocation error for CABLE set-up"
      call ccmpi_abort(-1)
    end if
    
    ! Maximum number of vegetation tiles in a gridbox tile
    tdata(tile)%maxnb = 0
    do n = 1,maxtile
      if ( tdata(tile)%tind(n,2)>=tdata(tile)%tind(n,1) ) then
        tdata(tile)%maxnb = n
      end if
    end do
    
    do n = 1,maxtile
      call cable_pack(svs(:,n),tdata(tile)%sv,tile,nb=n)
      call cable_pack(ivs(:,n),tdata(tile)%cveg,tile,nb=n)
      dummy_pack(is:ie) = max( vlin(is:ie,n), 0.01 )
      call cable_pack(dummy_pack(is:ie),tdata(tile)%vl2,tile,nb=n)
    end do
  
    call cable_pack(isoilm,soil(tile)%isoilm,tile)
    call cable_pack(slope_ave,soil(tile)%slope,tile)
    call cable_pack(slope_std,soil(tile)%slope_std,tile)

    ! slope is used by the gw model
    if ( tdata(tile)%mp>0 ) then
      soil(tile)%slope=min(max(soil(tile)%slope,1.e-8_8),0.95_8)
      soil(tile)%slope_std=min(max(soil(tile)%slope_std,1.e-8_8),0.95_8)
    end if

  end do ! tile

  ! Create a simple map of dominant land-use
  do iq = 1,ifull
    if ( land(iq) ) then
      ivmax = -1.
      do n = 1,maxtile
        if ( svs(iq,n)>ivmax ) then
          ivmax = svs(iq,n)
          ivegt(iq) = ivs(iq,n) ! diagnostic (CSIRO pft)
        end if
      end do
    end if
  end do
  
  ! Load canopy parameters
  call cable_biophysic_parm

  do tile = 1,ntiles
    ! Fixes for LAI
    if ( tdata(tile)%mp>0 ) then  
      where ( veg(tile)%iveg>=14 .and. veg(tile)%iveg<=17 )
        tdata(tile)%vl2(:) = 1.e-8
      end where
    end if
    
  end do ! tile

  ! Load soil parameters
  call cable_soil_parm
  
  do tile = 1,ntiles
    is = 1 + (tile-1)*imax
    ie = tile*imax
    if ( tdata(tile)%mp>0 ) then
      ! store bare soil albedo and define snow free albedo
      call cable_pack(albvisnir(is:ie,1),soil(tile)%albsoil(:,1),tile)
      call cable_pack(albvisnir(is:ie,2),soil(tile)%albsoil(:,2),tile)
      soil(tile)%albsoil(:,3) = 0.05_8
    end if
  end do
  
  ! Average albedo (for radiation)
  where ( land(1:ifull) )
    albsoil(1:ifull) = 0.5*sum(albvisnir(1:ifull,:),2)
  end where
  where ( albsoil(1:ifull)<=0.14 .and. land(1:ifull) )
    !sfact=0.5 for alb <= 0.14
    albsoilsn(1:ifull,1) = (1.00/1.50)*albsoil(1:ifull)
    albsoilsn(1:ifull,2) = (2.00/1.50)*albsoil(1:ifull)
  elsewhere ( albsoil(1:ifull)<=0.2 .and. land(1:ifull) )
    !sfact=0.62 for 0.14 < alb <= 0.20
    albsoilsn(1:ifull,1) = (1.24/1.62)*albsoil(1:ifull)
    albsoilsn(1:ifull,2) = (2.00/1.62)*albsoil(1:ifull)
  elsewhere ( land(1:ifull) )
    !sfact=0.68 for 0.2 < alb
    albsoilsn(1:ifull,1) = (1.36/1.68)*albsoil(1:ifull)
    albsoilsn(1:ifull,2) = (2.00/1.68)*albsoil(1:ifull)
  end where

  ! MJT patch
  do tile = 1,ntiles
    is = 1 + (tile-1)*imax
    ie = tile*imax
    if ( tdata(tile)%mp>0 ) then
      call cable_pack(albvisnir(is:ie,1),soil(tile)%albsoil(:,1),tile)
      call cable_pack(albvisnir(is:ie,2),soil(tile)%albsoil(:,2),tile)
      call cable_pack(albsoilsn(is:ie,1),ssnow(tile)%albsoilsn(:,1),tile) ! overwritten by CABLE
      call cable_pack(albsoilsn(is:ie,2),ssnow(tile)%albsoilsn(:,2),tile) ! overwritten by CABLE
      call cable_pack(albsoil(is:ie),rad(tile)%albedo_T(:),tile)
      dummy_pack(is:ie) = rlatt(is:ie)*180./pi
      call cable_pack(dummy_pack(is:ie),rad(tile)%latitude(:),tile)
      dummy_pack(is:ie) = rlongg(is:ie)*180./pi
      call cable_pack(dummy_pack(is:ie),rad(tile)%longitude(:),tile)

      !veg(tile)%vcmax_shade = veg(tile)%vcmax
      !veg(tile)%ejmax_shade = veg(tile)%ejmax
      !veg(tile)%vcmax_sun   = veg(tile)%vcmax
      !veg(tile)%ejmax_sun   = veg(tile)%ejmax
 
      ssnow(tile)%albsoilsn(:,3)=0.05_8    
      ssnow(tile)%t_snwlr=0.05_8
      ssnow(tile)%pudsmx=0._8
      ssnow(tile)%zdelta = 0._8
      ssnow(tile)%le = 0._8
  
      canopy(tile)%oldcansto=0._8  
      canopy(tile)%ghflux=0._8
      canopy(tile)%sghflux=0._8
      canopy(tile)%ga=0._8
      canopy(tile)%dgdtg=0._8
      canopy(tile)%fhs_cor=0._8
      canopy(tile)%fes_cor=0._8
      canopy(tile)%fns_cor=0._8
      canopy(tile)%ga_cor=0._8
      canopy(tile)%us=0.01_8
      ssnow(tile)%wb_lake=0._8 ! not used when mlo.f90 is active
      ssnow(tile)%fland=1._8
      ssnow(tile)%ifland=soil(tile)%isoilm
    
      ! Initialise sum flux variables
      sum_flux(tile)%sumpn=0._8
      sum_flux(tile)%sumrp=0._8
      sum_flux(tile)%sumrs=0._8
      sum_flux(tile)%sumrd=0._8
      sum_flux(tile)%sumrpw=0._8
      sum_flux(tile)%sumrpr=0._8
      sum_flux(tile)%dsumpn=0._8
      sum_flux(tile)%dsumrp=0._8
      sum_flux(tile)%dsumrs=0._8
      sum_flux(tile)%dsumrd=0._8
  
      bal(tile)%evap_tot=0._8
      bal(tile)%precip_tot=0._8
      bal(tile)%ebal_tot=0._8
      bal(tile)%rnoff_tot=0._8
    end if ! mp>0 
  end do   ! tile

  
  ! Initialise CABLE carbon pools
  if ( ccycle==0 ) then
    if ( cable_pop==1 ) then
      write(6,*) "ERROR: cable_pop=1 requires ccycle>0"
      call ccmpi_abort(-1)
    end if    
  else if ( ccycle>=1 .and. ccycle<=3 ) then
    ! CASA CNP
    do tile = 1,ntiles
      is = 1 + (tile-1)*imax
      ie = tile*imax
      if ( tdata(tile)%mp>0 ) then
        call alloc_casavariable(casabiome(tile),casapool(tile),casaflux(tile),casamet(tile), &
                                casabal(tile),tdata(tile)%mp)
        call alloc_phenvariable(phen(tile),tdata(tile)%mp)
        
        casamet(tile)%lat = rad(tile)%latitude
        call cable_pack(casapoint(:,1),casamet(tile)%isorder,tile)
        dummy_pack(is:ie) = casapoint(is:ie,2)/365.*1.E-3
        call cable_pack(dummy_pack(is:ie),casaflux(tile)%Nmindep,tile)
        dummy_pack(is:ie) = casapoint(is:ie,3)/365.
        call cable_pack(dummy_pack(is:ie),casaflux(tile)%Nminfix,tile)
        dummy_pack(is:ie) = casapoint(is:ie,4)/365.
        call cable_pack(dummy_pack(is:ie),casaflux(tile)%Pdep,tile)
        dummy_pack(is:ie) = casapoint(is:ie,5)/365.
        call cable_pack(dummy_pack(is:ie),casaflux(tile)%Pwea,tile)

        where ( veg(tile)%iveg==9 .or. veg(tile)%iveg==10 ) ! crops
          ! P fertilizer =13 Mt P globally in 1994
          casaflux(tile)%Pdep = casaflux(tile)%Pdep + 0.7_8/365._8
          ! N fertilizer =86 Mt N globally in 1994
          casaflux(tile)%Nminfix = casaflux(tile)%Nminfix + 4.3_8/365._8
        end where

        if ( any( casamet(tile)%isorder<1 .or. casamet(tile)%isorder>12 ) ) then
          write(6,*) "ERROR: Invalid isorder in CASA-CNP"
          call ccmpi_abort(-1)
        end if
        
      end if ! mp>0
    end do   ! tile

    ! Load CASA parameters
    call casa_readbiome(veg,casabiome,casapool,casaflux, &
                        casamet,phen,fcasapft)
    
    do tile = 1,ntiles
      is = 1 + (tile-1)*imax
      ie = tile*imax
      if ( tdata(tile)%mp>0 ) then
        do i = 1,tdata(tile)%mp
          ilat = nint((rad(tile)%latitude(i)+55.25)*2.) + 1
          ilat = min( 271, max( 1, ilat ) )
          ivp = veg(tile)%iveg(i)
          phen(tile)%phen(i)       = 1._8
          phen(tile)%aphen(i)      = 0._8
          phen(tile)%phase(i)      = phendoy1(ilat,ivp)
          phen(tile)%doyphase(i,1) = greenup(ilat,ivp)                ! DOY for greenup
          phen(tile)%doyphase(i,2) = phen(tile)%doyphase(i,1) + 14    ! DOY for steady LAI
          phen(tile)%doyphase(i,3) = fall(ilat,ivp)                   ! DOY for leaf senescence
          phen(tile)%doyphase(i,4) = phen(tile)%doyphase(i,3) + 14    ! DOY for minimal LAI season
          if ( phen(tile)%doyphase(i,2) > 365 ) then
            phen(tile)%doyphase(i,2) = phen(tile)%doyphase(i,2) - 365
          end if  
          if ( phen(tile)%doyphase(i,4) > 365 ) then
            phen(tile)%doyphase(i,4) = phen(tile)%doyphase(i,4) - 365
          end if
        end do ! i = 1,mp

        casamet(tile)%tairk         = 0._8
        casamet(tile)%tsoil         = 0._8
        casamet(tile)%moist         = 0._8
    
        casaflux(tile)%cgpp         = 0._8
        casaflux(tile)%Crsoil       = 0._8
        casaflux(tile)%crgplant     = 0._8
        casaflux(tile)%crmplant     = 0._8
        casaflux(tile)%clabloss     = 0._8
        casaflux(tile)%frac_sapwood = 1._8
        casaflux(tile)%sapwood_area = 0._8
        casaflux(tile)%stemnpp      = 0._8
        casaflux(tile)%Cnpp         = 0._8
        !casaflux(tile)%fHarvest     = 0._8
        !casaflux(tile)%NHarvest     = 0._8
        !casaflux(tile)%CHarvest     = 0._8
        !casaflux(tile)%fcrop        = 0._8

        canopy(tile)%fnee = 0._8
        canopy(tile)%fpn = 0._8
        canopy(tile)%frday = 0._8
        canopy(tile)%frp = 0._8
        canopy(tile)%frpw = 0._8
        canopy(tile)%frpr = 0._8
        canopy(tile)%frs = 0._8
   
        casabal(tile)%LAImax = 0._8
        casabal(tile)%Cleafmean = 0._8
        casabal(tile)%Crootmean = 0._8
      end if
 
      ! Update gridbox average carbon pools
      cplant(is:ie,:)=0.
      clitter(is:ie,:)=0.
      csoil(is:ie,:)=0.
      niplant(is:ie,:)=0.
      nilitter(is:ie,:)=0.
      nisoil(is:ie,:)=0.
      pplant(is:ie,:)=0.
      plitter(is:ie,:)=0.
      psoil(is:ie,:)=0.
      do k = 1,mplant
        call cable_update(cplant(is:ie,k),casapool(tile)%cplant(:,k),tile)
        call cable_update(niplant(is:ie,k),casapool(tile)%nplant(:,k),tile)
        call cable_update(pplant(is:ie,k),casapool(tile)%pplant(:,k),tile)
      end do
      do k = 1,mlitter
        call cable_update(clitter(is:ie,k),casapool(tile)%clitter(:,k),tile)
        call cable_update(nilitter(is:ie,k),casapool(tile)%nlitter(:,k),tile)
        call cable_update(plitter(is:ie,k),casapool(tile)%plitter(:,k),tile)
      end do
      do k = 1,msoil
        call cable_update(csoil(is:ie,k),casapool(tile)%csoil(:,k),tile)
        call cable_update(nisoil(is:ie,k),casapool(tile)%nsoil(:,k),tile)
        call cable_update(psoil(is:ie,k),casapool(tile)%psoil(:,k),tile)
      end do
    
    end do ! tile
    
    ! POP
    if ( cable_pop==1 ) then
        
      mp_pop_max = 0  
      do tile = 1,ntiles
        is = 1 + (tile-1)*imax
        ie = tile*imax
        mp_POP = count(casamet(tile)%iveg2==forest.or.casamet(tile)%iveg2==shrub)
        allocate( Iwood(mp_POP), disturbance_interval(mp_POP,2) )
        allocate( tdata(tile)%pmap(imax,maxtile) )
        tdata(tile)%pmap = .false.
        ipos = 0
        do n = 1,tdata(tile)%maxnb
          iis = tdata(tile)%tind(n,1)
          iie = tdata(tile)%tind(n,2)
          tdata(tile)%pind(n,1) = ipos + 1
          do i = iis,iie
            if ( casamet(tile)%iveg2(i)==forest .or. casamet(tile)%iveg2(i)==shrub ) then
              ipos = ipos + 1
              Iwood(ipos) = i
              tdata(tile)%pmap(i,n) = .true.
              disturbance_interval(ipos,:) = veg(tile)%disturbance_interval(i,:)
            end if
          end do
          tdata(tile)%pind(n,2) = ipos
        end do
        tdata(tile)%np = ipos
        mp_pop_max = max( mp_pop_max, ipos )
        call POP_init(POP(tile), disturbance_interval(:,:), mp_POP, Iwood(1:mp_POP)) 
        deallocate( Iwood, disturbance_interval )
      end do ! tile 

    end if ! cable_pop==1
    
  else  
    write(6,*) "ERROR: Unknown option ccycle ",ccycle
    call ccmpi_abort(-1)
  end if ! ccycle>0

  
  ! Calculate LAI and veg fraction diagnostics
  ! (needs to occur after CASA-CNP in case prognostic LAI is required)
  vlai(:) = 0.
  sigmf(:) = 0.
  allocate( dummy_unpack(mp_max) )
  do tile = 1,ntiles
    if ( tdata(tile)%mp>0 ) then
      call setlai(tdata(tile)%sv,tdata(tile)%vl2,casamet(tile),veg(tile),tdata(tile)%mp)
      call cable_update(vlai(:),veg(tile)%vlai,tile)
      dummy_unpack(1:tdata(tile)%mp) = 1._8 - exp(-0.4_8*veg(tile)%vlai(:))
      call cable_update(sigmf(:),dummy_unpack(1:tdata(tile)%mp),tile)
    end if  
  end do ! tile   
  deallocate( dummy_unpack )

  ! MJT suggestion to get an approx inital albedo (before cable is called)
  where ( land(1:ifull) )
    albvisnir(:,1) = albsoilsn(:,1)*(1.-sigmf) + 0.03*sigmf
    albvisnir(:,2) = albsoilsn(:,2)*(1.-sigmf) + 0.40*sigmf
  end where
  albvisdir(:) = albvisnir(:,1) ! To be updated by CABLE
  albvisdif(:) = albvisnir(:,1) ! To be updated by CABLE
  albnirdir(:) = albvisnir(:,2) ! To be updated by CABLE
  albnirdif(:) = albvisnir(:,2) ! To be updated by CABLE
  
else

  do tile = 1,ntiles
    allocate( tdata(tile)%sv(0) )
    allocate( tdata(tile)%vl2(0) )
    allocate( tdata(tile)%cveg(0) )
    allocate( tdata(tile)%old_sv(0) )
    allocate( tdata(tile)%old_cveg(0) )
    call alloc_cbm_var(air(tile), 0)
    call alloc_cbm_var(bgc(tile), 0)
    call alloc_cbm_var(canopy(tile), 0)
    call alloc_cbm_var(met(tile), 0)
    call alloc_cbm_var(bal(tile), 0)
    call alloc_cbm_var(rad(tile), 0)
    call alloc_cbm_var(rough(tile), 0)
    call alloc_cbm_var(soil(tile), 0)
    call alloc_cbm_var(ssnow(tile), 0)
    call alloc_cbm_var(sum_flux(tile), 0)
    call alloc_cbm_var(veg(tile), 0)
    climate(tile)%nyear_average = 20
    climate(tile)%nday_average = 31
  end do
  call cable_biophysic_parm
  call cable_soil_parm
  if ( ccycle>=1 .and. ccycle<=3 ) then
    do tile = 1,ntiles  
      call alloc_casavariable(casabiome(tile),casapool(tile),casaflux(tile),casamet(tile), &
                              casabal(tile),0)
      call alloc_phenvariable(phen(tile),0)
    end do  
    call casa_readbiome(veg,casabiome,casapool,casaflux,casamet,phen,fcasapft)
    if ( cable_pop==1 ) then
      mp_POP = 0
      mp_pop_max = 0
      do tile = 1,ntiles
        allocate( Iwood(0), disturbance_interval(0,2) )
        call POP_init(POP(tile), disturbance_interval, 0, Iwood) 
        deallocate( Iwood, disturbance_interval )
      end do
    end if
  end if
  
end if
  
! statistics
stat_count(:) = 0.
global_stat_count(:) = 0.
do iq = 1,ifull
  if ( land(iq) ) then  
    landcount = count( svs(iq,:)>0. )
    stat_count(landcount) = stat_count(landcount) + 1.
    stat_count(0) = stat_count(0) + 1.
  end if  
end do  
call ccmpi_reduce(stat_count,global_stat_count,"sum",0,comm_world)
if ( myid==0 ) then
  write(6,*) "CABLE statistics:"
  do n = 1,maxtile
    write(6,'(A,I1.1,A,F5.1,A)') "   Percentage of gridpoints with ",n," tile(s) is ", &
        100.*global_stat_count(n)/global_stat_count(0),"%"
  end do  
end if
  
if ( myid==0 .and. nmaxpr==1 ) write(6,*) "Finished defining CABLE and CASA CNP arrays"

return
end subroutine cbmparm

subroutine casa_readbiome(veg,casabiome,casapool,casaflux,casamet,phen,fcasapft)

use cc_mpi
use newmpar_m

type(veg_parameter_type), dimension(ntiles), intent(in) :: veg
type(casa_biome), dimension(ntiles), intent(inout) :: casabiome
type(casa_pool), dimension(ntiles), intent(inout) :: casapool
type(casa_flux), dimension(ntiles), intent(inout) :: casaflux
type(casa_met), dimension(ntiles), intent(inout) :: casamet
type(phen_variable), dimension(ntiles), intent(inout) :: phen
character(len=*), intent(in) :: fcasapft

real(kind=8), dimension(mxvt,mplant) :: ratiocnplant
real(kind=8), dimension(mxvt,msoil) :: ratiocnsoil,ratiocnsoilmax,ratiocnsoilmin
real(kind=8), dimension(mso,msoil) :: rationpsoil

real(kind=8), dimension(mxvt) :: leafage,woodage,frootage,metage
real(kind=8), dimension(mxvt) :: strage,cwdage,micage,slowage,passage,slax
real(kind=8), dimension(mxvt) :: xfherbivore,xxkleafcoldmax,xxkleafdrymax
real(kind=8), dimension(mxvt) :: xratioNPleafmin,xratioNPleafmax,xratioNPwoodmin,xratioNPwoodmax
real(kind=8), dimension(mxvt) :: xratioNPfrootmin,xratioNPfrootmax,xfNminloss,xfNminleach,xnfixrate
real(kind=8), dimension(mxvt) :: xnsoilmin,xplab,xpsorb,xpocc
real(kind=8), dimension(mxvt) :: clabileage
real(kind=8), dimension(mxvt) :: xxnpmax,xq10soil,xxkoptlitter,xxkoptsoil,xprodptase
real(kind=8), dimension(mxvt) :: xcostnpup,xmaxfinelitter,xmaxcwd,xnintercept,xnslope
real(kind=8), dimension(mxvt) :: xla_to_sa,xvcmax_scalar,xdisturbance_interval
real(kind=8), dimension(mxvt) :: xDAMM_EnzPool,xDAMM_KMO2,xDAMM_KMcp,xDAMM_Ea,xDAMM_alpha
real(kind=8), dimension(mso) :: xxkplab,xxkpsorb,xxkpocc
real(kind=8), dimension(mso) :: xkmlabp,xpsorbmax,xfPleach

integer i, iso, nv, ierr
integer nv0,nv1,nv2,nv3,nv4,nv5,nv6,nv7,nv8,nv9,nv10,nv11,nv12,nv13
integer :: fflag=0
integer tile

integer, dimension(mxvt) :: ivt2
real(kind=8), dimension(mxvt) :: kroot
real(kind=8), dimension(mxvt) :: rootdepth
real(kind=8), dimension(mxvt) :: kuptake
real(kind=8), dimension(mxvt) :: krootlen
real(kind=8), dimension(mxvt) :: kminn
real(kind=8), dimension(mxvt) :: kuplabp
real(kind=8), dimension(mxvt,mplant) :: fracnpptop
real(kind=8), dimension(mxvt,mplant) :: rmplant
real(kind=8), dimension(mxvt,mplant) :: ftransnptol
real(kind=8), dimension(mxvt,mplant) :: fracligninplant
real(kind=8), dimension(mxvt) :: glaimax
real(kind=8), dimension(mxvt) :: glaimin
real(kind=8), dimension(mxvt) :: xkleafcoldexp
real(kind=8), dimension(mxvt) :: xkleafdryexp
real(kind=8), dimension(mxvt,mplant) :: rationcplantmin
real(kind=8), dimension(mxvt,mplant) :: rationcplantmax
real(kind=8), dimension(mxvt,mplant) :: ftranspptol
real(kind=8), dimension(mxvt) :: tkshed

if ( trim(fcasapft) /= '' ) fflag = 1
call ccmpi_bcast(fflag,0,comm_world)

if ( fflag == 1  ) then
  if ( myid == 0 ) then
    write(6,*) "-> Using user defined CASA PFT parameter tables"
    open(86,file=fcasapft,status='old',action='read',iostat=ierr)
    write(6,*) "ERROR: CASA PFT not supported for user file ",trim(fcasapft)
    call ccmpi_abort(-1)

    write(6,*) "-> Reading ",trim(fcasapft)

    do i=1,3
      read(86,*)
    enddo

    do nv=1,mxvt
      read(86,*,iostat=ierr) nv0,ivt2(nv)
      if ( ierr/=0 .or. nv0/=nv  ) then
        write(6,*) "ERROR: Cannot read casapftfile file ",trim(fcasapft)
        write(6,*) "Formatting error in line",3+0*(2+mxvt)+nv,"reading nv0"
        if ( nv0/=nv ) then
          write(6,*) "pft index is not sequential"
        end if
        call ccmpi_abort(-1)
       end if
    enddo

    read(86,*)
    read(86,*)
    do nv=1,mxvt
      read(86,*,iostat=ierr) nv1,kroot(nv),rootdepth(nv),      &
                  kuptake(nv),krootlen(nv),         &
                  kminn(nv), kuplabp(nv),           &
                  xfherbivore(nv),leafage(nv),woodage(nv),frootage(nv), &
                  metage(nv),strage(nv),cwdage(nv),  &
                  micage(nv),slowage(nv),passage(nv),clabileage(nv),slax(nv)
      if ( ierr/=0 .or. nv1/=nv ) then
        write(6,*) "ERROR: Cannot read casapftfile file ",trim(fcasapft)
        write(6,*) "Formatting error in line",3+1*(2+mxvt)+nv,"reading nv1"
        if ( nv1/=nv ) then
          write(6,*) "pft index is not sequential"
        end if
        call ccmpi_abort(-1)
       end if
    enddo

    read(86,*)
    read(86,*)
    do nv=1,mxvt
      read(86,*,iostat=ierr) nv2, &
                  fracnpptop(nv,leaf),fracnpptop(nv,wood), &
                  fracnpptop(nv,xroot),rmplant(nv,leaf),   &
                  rmplant(nv,wood),rmplant(nv,xroot)
      if ( ierr/=0 .or. nv2/=nv ) then
        write(6,*) "ERROR: Cannot read casapftfile file ",trim(fcasapft)
        write(6,*) "Formatting error in line",3+2*(2+mxvt)+nv,"reading nv2"
        if ( nv2/=nv ) then
          write(6,*) "pft index is not sequential"
        end if
        call ccmpi_abort(-1)
       end if
    enddo

    read(86,*)
    read(86,*)
    do nv=1,mxvt
      read(86,*,iostat=ierr) nv3, ratiocnplant(nv,leaf),ratiocnplant(nv,wood),   &
           ratiocnplant(nv,xroot),                                         &
           ftransnptol(nv,leaf), ftransnptol(nv,wood), &
           ftransnptol(nv,xroot),                                &
           fracligninplant(nv,leaf),                             &
           fracligninplant(nv,wood),                             &
           fracligninplant(nv,xroot),                            &
           ratiocnsoil(nv,mic),ratiocnsoil(nv,slow),ratiocnsoil(nv,pass),  &
           ratiocnsoilmin(nv,mic),ratiocnsoilmin(nv,slow),ratiocnsoilmin(nv,pass),  &
           ratiocnsoilmax(nv,mic),ratiocnsoilmax(nv,slow),ratiocnsoilmax(nv,pass),  &
           glaimax(nv),glaimin(nv)
      if ( ierr/=0 .or. nv3/=nv ) then
        write(6,*) "ERROR: Cannot read casapftfile file ",trim(fcasapft)
        write(6,*) "Formatting error in line",3+3*(2+mxvt)+nv,"reading nv3"
        if ( nv3/=nv ) then
          write(6,*) "pft index is not sequential"
        end if
        call ccmpi_abort(-1)
       end if
    enddo

    read(86,*)
    read(86,*)
    do nv=1,mxvt
      read(86,*,iostat=ierr) nv4, cleaf(nv),cwood(nv),cfroot(nv),cmet(nv),   &
                  cstr(nv),ccwd(nv), cmic(nv), cslow(nv),cpass(nv)
      if ( ierr/=0 .or. nv4/=nv ) then
        write(6,*) "ERROR: Cannot read casapftfile file ",trim(fcasapft)
        write(6,*) "Formatting error in line",3+4*(2+mxvt)+nv,"reading nv4"
        if ( nv4/=nv ) then
          write(6,*) "pft index is not sequential"
        end if
        call ccmpi_abort(-1)
       end if
    enddo
    
    if ( cable_pop==1 ) then
      if ( cwood(1)>0.02_8 ) then
        write(6,*) "ERROR: casapftfile is not configured for POP"
        call ccmpi_abort(-1)
      end if
    else
      if ( cwood(1)<0.02_8 ) then
        write(6,*) "ERROR: casapftfile is configured for POP"
        call ccmpi_abort(-1)          
      end if
    end if

    read(86,*)
    read(86,*)
    do nv=1,mxvt
      read(86,*,iostat=ierr) nv5, &
           tkshed(nv),xxkleafcoldmax(nv),xkleafcoldexp(nv), &
           xxkleafdrymax(nv),xkleafdryexp(nv)
      if ( ierr/=0 .or. nv5/=nv ) then
        write(6,*) "ERROR: Cannot read casapftfile file ",trim(fcasapft)
        write(6,*) "Formatting error in line",3+5*(2+mxvt)+nv,"reading nv5"
        if ( nv5/=nv ) then
          write(6,*) "pft index is not sequential"
        end if
        call ccmpi_abort(-1)
       end if
    enddo

    read(86,*)
    read(86,*)
    do nv=1,mxvt
      read(86,*,iostat=ierr) nv6, &
        rationcplantmin(nv,leaf),rationcplantmax(nv,leaf), &
        rationcplantmin(nv,wood),rationcplantmax(nv,wood), &
        rationcplantmin(nv,xroot),rationcplantmax(nv,xroot), &
        xfnminloss(nv), xfnminleach(nv),xnfixrate(nv)
      if ( ierr/=0 .or. nv6/=nv ) then
        write(6,*) "ERROR: Cannot read casapftfile file ",trim(fcasapft)
        write(6,*) "Formatting error in line",3+6*(2+mxvt)+nv,"reading nv6"
        if ( nv6/=nv ) then
          write(6,*) "pft index is not sequential"
        end if
        call ccmpi_abort(-1)
       end if
    enddo

    read(86,*)
    read(86,*)
    do nv=1,mxvt
      read(86,*,iostat=ierr) nv7,nleaf(nv),nwood(nv),nfroot(nv), &
                  nmet(nv),nstr(nv), ncwd(nv), &
                  nmic(nv),nslow(nv),npass(nv),xnsoilmin(nv)
      if ( ierr/=0 .or. nv7/=nv ) then
        write(6,*) "ERROR: Cannot read casapftfile file ",trim(fcasapft)
        write(6,*) "Formatting error in line",3+7*(2+mxvt)+nv,"reading nv7"
        if ( nv7/=nv ) then
          write(6,*) "pft index is not sequential"
        end if
        call ccmpi_abort(-1)
       end if
    enddo

    read(86,*)
    read(86,*)
    do nv=1,mxvt
      read(86,*,iostat=ierr) nv8,xrationpleafmin(nv),xrationpleafmax(nv),    &
           xrationpwoodmin(nv),xrationpwoodmax(nv),              &
           xrationpfrootmin(nv),xrationpfrootmax(nv),            &
           ftranspptol(nv,leaf), ftranspptol(nv,wood), &
           ftranspptol(nv,xroot)
      if ( ierr/=0 .or. nv8/=nv ) then
        write(6,*) "ERROR: Cannot read casapftfile file ",trim(fcasapft)
        write(6,*) "Formatting error in line",3+8*(2+mxvt)+nv,"reading nv8"
        if ( nv8/=nv ) then
          write(6,*) "pft index is not sequential"
        end if
        call ccmpi_abort(-1)
       end if
    enddo

    read(86,*)
    read(86,*)
    do iso=1,mso
      read(86,*,iostat=ierr) nv9,xkmlabp(iso),xpsorbmax(iso),xfpleach(iso), &
                  rationpsoil(iso,mic),rationpsoil(iso,slow),rationpsoil(iso,pass), &
                  xxkplab(iso),xxkpsorb(iso),xxkpocc(iso)
      if ( ierr/=0 .or. nv9/=iso ) then
        write(6,*) "ERROR: Cannot read casapftfile file ",trim(fcasapft)
        write(6,*) "Formatting error in line",3+9*(2+mxvt)+nv,"reading nv9"
        if ( nv9/=iso ) then
          write(6,*) "pft index is not sequential"
        end if
        call ccmpi_abort(-1)
       end if
    enddo
    read(86,*)
    read(86,*)
    do nv=1,mxvt
      read(86,*,iostat=ierr) nv10, &
           xpleaf(nv),xpwood(nv),xpfroot(nv),xpmet(nv),xpstr(nv),xpcwd(nv), &
           xpmic(nv),xpslow(nv),xppass(nv),xplab(nv),xpsorb(nv),xpocc(nv)
      if ( ierr/=0 .or. nv10/=nv ) then
        write(6,*) "ERROR: Cannot read casapftfile file ",trim(fcasapft)
        write(6,*) "Formatting error in line",3+10*(2+mxvt)+nv,"reading nv10"
        if ( nv10/=nv ) then
          write(6,*) "pft index is not sequential"
        end if
        call ccmpi_abort(-1)
       end if
    enddo

    read(86,*)
    read(86,*)
    do nv=1,mxvt
      read(86,*,iostat=ierr) nv11, &
           xxnpmax(nv),xq10soil(nv),xxkoptlitter(nv),xxkoptsoil(nv),xprodptase(nv), &
           xcostnpup(nv),xmaxfinelitter(nv),xmaxcwd(nv),xnintercept(nv),xnslope(nv)
      if ( ierr/=0 .or. nv11/=nv ) then
        write(6,*) "ERROR: Cannot read casapftfile file ",trim(fcasapft)
        write(6,*) "Formatting error in line",3+11*(2+mxvt)+nv,"reading nv11"
        if ( nv11/=nv ) then
          write(6,*) "pft index is not sequential"
        end if
        call ccmpi_abort(-1)
       end if
    enddo

    read(86,*)
    read(86,*)
    do nv=1,mxvt
      read(86,*,iostat=ierr) nv12, &
           xla_to_sa(nv),xdisturbance_interval(nv),xvcmax_scalar(nv)
      if ( ierr/=0 .or. nv12/=nv ) then
        write(6,*) "ERROR: Cannot read casapftfile file ",trim(fcasapft)
        write(6,*) "Formatting error in line",3+12*(2+mxvt)+nv,"reading nv12"
        if ( nv12/=nv ) then
          write(6,*) "pft index is not sequential"
        end if
        call ccmpi_abort(-1)
       end if
    enddo

    read(86,*)
    read(86,*)
    do nv=1,mxvt
      read(86,*,iostat=ierr) nv13, &
           xdamm_enzpool(nv), xdamm_kmo2(nv),xdamm_kmcp(nv), xdamm_ea(nv), xdamm_alpha(nv)
      if ( ierr/=0 .or. nv13/=nv ) then
        write(6,*) "ERROR: Cannot read casapftfile file ",trim(fcasapft)
        write(6,*) "Formatting error in line",3+13*(2+mxvt)+nv,"reading nv13"
        if ( nv13/=nv ) then
          write(6,*) "pft index is not sequential"
        end if
        call ccmpi_abort(-1)
       end if
    enddo

    close(86)
  end if

  call ccmpi_bcast(ivt2,0,comm_world)
  call ccmpi_bcastr8(kroot,0,comm_world)
  call ccmpi_bcastr8(rootdepth,0,comm_world)
  call ccmpi_bcastr8(kuptake,0,comm_world)
  call ccmpi_bcastr8(krootlen,0,comm_world)
  call ccmpi_bcastr8(kminn,0,comm_world)
  call ccmpi_bcastr8(kuplabp,0,comm_world)
  call ccmpi_bcastr8(xfherbivore,0,comm_world)
  call ccmpi_bcastr8(leafage,0,comm_world)
  call ccmpi_bcastr8(woodage,0,comm_world)
  call ccmpi_bcastr8(frootage,0,comm_world)
  call ccmpi_bcastr8(metage,0,comm_world)
  call ccmpi_bcastr8(strage,0,comm_world)
  call ccmpi_bcastr8(cwdage,0,comm_world)
  call ccmpi_bcastr8(micage,0,comm_world)
  call ccmpi_bcastr8(slowage,0,comm_world)
  call ccmpi_bcastr8(passage,0,comm_world)
  call ccmpi_bcastr8(clabileage,0,comm_world)
  call ccmpi_bcastr8(slax,0,comm_world)

  call ccmpi_bcastr8(fracnpptop,0,comm_world)
  call ccmpi_bcastr8(rmplant,0,comm_world)

  call ccmpi_bcastr8(ratiocnplant,0,comm_world)
  call ccmpi_bcastr8(ftransnptol,0,comm_world)
  call ccmpi_bcastr8(fracligninplant,0,comm_world)
  call ccmpi_bcastr8(ratiocnsoil,0,comm_world)
  call ccmpi_bcastr8(ratiocnsoilmin,0,comm_world)
  call ccmpi_bcastr8(ratiocnsoilmax,0,comm_world)
  call ccmpi_bcastr8(glaimax,0,comm_world)
  call ccmpi_bcastr8(glaimin,0,comm_world)

  call ccmpi_bcastr8(cleaf,0,comm_world)
  call ccmpi_bcastr8(cwood,0,comm_world)
  call ccmpi_bcastr8(cfroot,0,comm_world)
  call ccmpi_bcastr8(cmet,0,comm_world)
  call ccmpi_bcastr8(cstr,0,comm_world)
  call ccmpi_bcastr8(ccwd,0,comm_world)
  call ccmpi_bcastr8(cmic,0,comm_world)
  call ccmpi_bcastr8(cslow,0,comm_world)
  call ccmpi_bcastr8(cpass,0,comm_world)

  call ccmpi_bcastr8(tkshed,0,comm_world)
  call ccmpi_bcastr8(xxkleafcoldmax,0,comm_world)
  call ccmpi_bcastr8(xkleafcoldexp,0,comm_world)
  call ccmpi_bcastr8(xxkleafdrymax,0,comm_world)
  call ccmpi_bcastr8(xkleafdryexp,0,comm_world)

  call ccmpi_bcastr8(rationcplantmin,0,comm_world)
  call ccmpi_bcastr8(rationcplantmax,0,comm_world)
  call ccmpi_bcastr8(xfnminloss,0,comm_world)
  call ccmpi_bcastr8(xfnminleach,0,comm_world)
  call ccmpi_bcastr8(xnfixrate,0,comm_world)

  call ccmpi_bcastr8(nleaf,0,comm_world)
  call ccmpi_bcastr8(nwood,0,comm_world)
  call ccmpi_bcastr8(nfroot,0,comm_world)
  call ccmpi_bcastr8(nmet,0,comm_world)
  call ccmpi_bcastr8(nstr,0,comm_world)
  call ccmpi_bcastr8(ncwd,0,comm_world)
  call ccmpi_bcastr8(nmic,0,comm_world)
  call ccmpi_bcastr8(nslow,0,comm_world)
  call ccmpi_bcastr8(npass,0,comm_world)
  call ccmpi_bcastr8(xnsoilmin,0,comm_world)

  call ccmpi_bcastr8(xrationpleafmin,0,comm_world)
  call ccmpi_bcastr8(xrationpleafmax,0,comm_world)
  call ccmpi_bcastr8(xrationpwoodmin,0,comm_world)
  call ccmpi_bcastr8(xrationpwoodmax,0,comm_world)
  call ccmpi_bcastr8(xrationpfrootmin,0,comm_world)
  call ccmpi_bcastr8(xrationpfrootmax,0,comm_world)
  call ccmpi_bcastr8(ftranspptol,0,comm_world)

  call ccmpi_bcastr8(xkmlabp,0,comm_world)
  call ccmpi_bcastr8(xpsorbmax,0,comm_world)
  call ccmpi_bcastr8(xfpleach,0,comm_world)
  call ccmpi_bcastr8(rationpsoil,0,comm_world)
  call ccmpi_bcastr8(xxkplab,0,comm_world)
  call ccmpi_bcastr8(xxkpsorb,0,comm_world)
  call ccmpi_bcastr8(xxkpocc,0,comm_world)

  call ccmpi_bcastr8(xpleaf,0,comm_world)
  call ccmpi_bcastr8(xpwood,0,comm_world)
  call ccmpi_bcastr8(xpfroot,0,comm_world)
  call ccmpi_bcastr8(xpmet,0,comm_world)
  call ccmpi_bcastr8(xpstr,0,comm_world)
  call ccmpi_bcastr8(xpcwd,0,comm_world)
  call ccmpi_bcastr8(xpmic,0,comm_world)
  call ccmpi_bcastr8(xpslow,0,comm_world)
  call ccmpi_bcastr8(xppass,0,comm_world)
  call ccmpi_bcastr8(xplab,0,comm_world)
  call ccmpi_bcastr8(xpsorb,0,comm_world)
  call ccmpi_bcastr8(xpocc,0,comm_world)

  call ccmpi_bcastr8(xxnpmax,0,comm_world)
  call ccmpi_bcastr8(xq10soil,0,comm_world)
  call ccmpi_bcastr8(xxkoptlitter,0,comm_world)
  call ccmpi_bcastr8(xxkoptsoil,0,comm_world)
  call ccmpi_bcastr8(xprodptase,0,comm_world)
  call ccmpi_bcastr8(xcostnpup,0,comm_world)
  call ccmpi_bcastr8(xmaxfinelitter,0,comm_world)
  call ccmpi_bcastr8(xmaxcwd,0,comm_world)
  call ccmpi_bcastr8(xnintercept,0,comm_world)
  call ccmpi_bcastr8(xnslope,0,comm_world)

  call ccmpi_bcastr8(xla_to_sa,0,comm_world)
  call ccmpi_bcastr8(xdisturbance_interval,0,comm_world)
  call ccmpi_bcastr8(xvcmax_scalar,0,comm_world)

  call ccmpi_bcastr8(xdamm_enzpool,0,comm_world)
  call ccmpi_bcastr8(xdamm_kmo2,0,comm_world)
  call ccmpi_bcastr8(xdamm_kmcp,0,comm_world)
  call ccmpi_bcastr8(xdamm_ea,0,comm_world)
  call ccmpi_bcastr8(xdamm_alpha,0,comm_world)
  
  do tile = 1,ntiles
    if ( tdata(tile)%mp>0 ) then
      casabiome(tile)%ivt2 = ivt2
      casabiome(tile)%kroot = kroot
      casabiome(tile)%rootdepth = rootdepth
      casabiome(tile)%kuptake = kuptake
      casabiome(tile)%krootlen = krootlen
      casabiome(tile)%kminn = kminn
      casabiome(tile)%kuplabp = kuplabp
      casabiome(tile)%fracnpptop = fracnpptop
      casabiome(tile)%rmplant = rmplant
      casabiome(tile)%ftransnptol = ftransnptol
      casabiome(tile)%fracligninplant = fracligninplant
      casabiome(tile)%glaimax = glaimax
      casabiome(tile)%glaimin = glaimin
      phen(tile)%tkshed = tkshed
      casabiome(tile)%xkleafcoldexp = xkleafcoldexp
      casabiome(tile)%xkleafdryexp = xkleafdryexp
      casabiome(tile)%rationcplantmin = rationcplantmin
      casabiome(tile)%rationcplantmax = rationcplantmax
      casabiome(tile)%ftranspptol = ftranspptol
    end if
  end do  
else
  if ( myid == 0 ) then
    write(6,*) "Using default CASA PFT parameter tables"
  end if
  if ( mp_global>0 ) then

    leafage =(/ 2.0_8, 1.5_8, 1.0_8, 1.0_8, 1.0_8, 0.8_8, 0.8_8, 1.0_8,      0.8_8,      0.8_8, 1.0_8, &
                1.0_8, 1.0_8, 1.0_8, 1.0_8, 1.0_8, 1.0_8 /)
    woodage =(/ 70._8, 60._8, 80._8, 40._8, 40._8, 1.0_8, 1.0_8, 1.0_8,      1.0_8,      1.0_8, 1.0_8, &
                1.0_8, 1.0_8, 5.0_8, 1.0_8, 1.0_8, 1.0_8 /)
    frootage=(/ 18._8, 10._8, 10._8, 10._8, 5.0_8, 3.0_8, 3.0_8, 3.0_8, 0.884227_8, 0.884227_8, 1.0_8, &
                1.0_8, 1.0_8, 4.0_8, 1.0_8, 1.0_8, 1.0_8 /)
    metage=0.04_8
    strage=0.23_8
    cwdage=0.824_8
    micage=0.137_8
    slowage=5._8
    passage=222.22_8
    clabileage=0.2_8
    slax = (/ 0.007178742_8, 0.015319264_8, 0.023103346_8, 0.026464564_8, 0.009919764_8, 0.029590494_8, &
              0.022417511_8, 0.026704118_8, 0.029590494_8, 0.022417511_8,        0.02_8,        0.02_8, &
              0.02_8,        0.024470894_8,        0.02_8,        0.02_8,        0.02_8 /)

    xfherbivore   =(/ 0.068_8, 0.406_8, 0.068_8, 0.134_8, 0.022_8, 0.109_8, 0.109_8, 0.109_8, 0.140_8, &
                      0.140_8, 0.000_8, 0.000_8, 0.000_8, 0.010_8, 0.000_8, 0.000_8, 0.000_8 /)
    xxkleafcoldmax=(/   0.2_8,  0.1_8,  0.1_8,  0.6_8,       1._8,  0.2_8,    0.2_8,   0.2_8,   0.3_8, &
                        0.3_8,  0.1_8,  0.1_8,  0.1_8,      0.1_8,  0.1_8,    0.1_8,   0.1_8 /)
    xxkleafdrymax =(/   0.1_8,  0.1_8,  0.1_8,   1._8,      0.1_8,  0.1_8,    0.1_8,   0.1_8,   0.1_8, &
                        0.1_8,  0.1_8,  0.1_8,  0.1_8,      0.1_8,  0.1_8,    0.1_8,   0.1_8 /)
    xratioNPleafmin =(/ 10.92308_8, 15.95339_8, 9.254839_8, 12.73848_8, 12.07217_8, 13.51473_8,    14.05_8, &
                        12.57800_8, 15.12262_8,      10._8,      13._8,      10._8,      10._8,  16.2336_8, &
                             10._8,      10._8,      10._8 /)
    xratioNPleafmax =(/ 12.07288_8,  17.6327_8, 10.22903_8, 14.07938_8, 13.34292_8, 14.93733_8, 15.52895_8, &
                          13.902_8, 16.71447_8,      10._8,      13._8,      10._8,      10._8,  17.9424_8, &
                             10._8,      10._8,      10._8 /)
    xratioNPwoodmin =(/ 20.30167_8, 15.89425_8, 17.48344_8, 19.08018_8, 22.46035_8,      15._8,      15._8, &
                           15.96_8,    20.52_8,      15._8,      15._8,      15._8,      15._8,  17.5275_8, &
                             15._8,      15._8,      15._8 /)
    xratioNPwoodmax =(/ 22.43869_8, 17.56733_8,  19.3238_8, 21.08862_8,  24.8246_8,      15._8,      15._8, &
                           17.64_8,    20.52_8,      15._8,      15._8,      15._8,      15._8,  19.3725_8, &
                             15._8,      15._8,      15._8 /)
    xratioNPfrootmin=(/ 20.29341_8, 15.87155_8, 17.39767_8,  19.0601_8, 22.49363_8, 15.63498_8, 16.08255_8, &
                        14.49241_8, 22.69109_8,      15._8,      15._8,      15._8,      15._8, 22.13268_8, &
                             15._8,      15._8,      15._8 /)
    xratioNPfrootmax=(/ 22.42955_8, 17.54224_8,   19.229_8, 21.06643_8, 24.86138_8, 17.28077_8, 17.77545_8, &
                        16.01793_8, 25.07962_8,      15._8,      15._8,      15._8,      15._8, 24.46244_8, &
                             15._8,      15._8,      15._8 /)
    xfNminloss=0.05_8
    xfNminleach=0.05_8
    xnfixrate=(/ 0.08_8, 2.6_8, 0.21_8, 1.64_8, 0.37_8, 0.95_8, 0.95_8, 0.95_8, 4._8, 4._8, 0._8, &
                   0._8,  0._8, 0.35_8,   0._8,   0._8,   0._8 /)
    xnsoilmin=1000._8
  
    ratiocnplant(:,leaf)=(/  49.8_8,  23.1_8,  59.3_8,  31.4_8,  37.6_8, 34.8_8,  44._8,  49.2_8,  21.6_8, &
                              25._8,   30._8,   30._8,   30._8,   50._8,  40._8,  40._8,   40._8 /)
    ratiocnplant(:,wood)=(/ 238.1_8, 134.9_8, 243.8_8, 156.2_8, 142.1_8, 150._8, 150._8, 147.3_8,  150._8, &
                             125._8,  150._8,  150._8,  150._8,  150._8, 150._8, 135._8,  150._8 /)
    ratiocnplant(:,xroot)=(/ 73.7_8,  61.2_8,   75._8,  63.2_8,  67.1_8, 64.5_8, 62.7_8,   69._8,  60.7_8, &
                              71._8,   71._8,   71._8,   71._8,   71._8,  71._8,  71._8,   71._8 /)
    ratiocnsoil(:,mic)=8._8
    ratiocnsoil(:,slow)=(/ 16.1_8, 12.8_8, 24.8_8,  30._8, 19.3_8, 13.1_8, 13.1_8, 13.1_8, 13.2_8, 13.2_8, &
                           13.1_8, 13.1_8, 13.1_8, 26.8_8,  20._8,  20._8,  20._8 /)
    ratiocnsoil(:,pass)=(/ 16.1_8, 12.8_8, 24.8_8,  30._8, 19.3_8, 13.1_8, 13.1_8, 13.1_8, 13.2_8, 13.2_8, &
                           13.1_8, 13.1_8, 13.1_8, 26.8_8,  20._8,  20._8,  20._8 /)
    ratiocnsoilmin(:,mic)=3._8
    ratiocnsoilmin(:,slow)=12._8
    ratiocnsoilmin(:,pass)=7._8
    ratiocnsoilmax(:,mic)=15._8
    ratiocnsoilmax(:,slow)=30._8
    ratiocnsoilmax(:,pass)=15._8
   
    xplab  =(/   26.737_8,   19.947_8,   29.107_8,   30.509_8,   23.206_8,   25.538_8,   25.538_8,   25.538_8, &
                 27.729_8,   27.729_8,       0._8,       0._8,       0._8,   21.038_8,       0._8,       0._8, &
                  0.103_8 /)
    xpsorb =(/   126.73_8,   92.263_8,  134.639_8,  132.012_8,   173.47_8,  186.207_8,  186.207_8,  186.207_8, &
                155.518_8,  155.518_8,       0._8,       0._8,       0._8,   255.79_8,       0._8,       0._8, &
                  1.176_8 /)
    xpocc  =(/  138.571_8,  120.374_8,   138.22_8,  148.083_8,  114.496_8,  145.163_8,  145.163_8,  145.163_8, &
                158.884_8,  158.884_8,       0._8,       0._8,       0._8,  108.897_8,       0._8,       0._8, &
                  0.688_8 /)

    xkmlabp  =(/ 74.5408_8,  68.1584_8,   77.952_8, 64.41918_8, 64.41918_8, 70.5856_8,  64.5888_8, 54.1692_8, &
                  9.7704_8,    28.29_8,   63.963_8,   32.402_8 /)
    xpsorbmax=(/ 745.408_8, 788.0815_8, 1110.816_8,  744.847_8,  744.847_8, 816.146_8, 746.8081_8, 722.256_8, &
                 293.112_8,   311.19_8, 373.1175_8, 615.6381_8 /)
    xfPleach =0.0005_8
    ratioNPsoil(:,mic)=4._8
    ratioNPsoil(:,slow)=(/ 5._8, 5._8, 5._8, 15._8, 5._8, 5._8, 5._8, 5._8, 7._8, 7._8, 7._8, 7._8 /)
    ratioNPsoil(:,pass)=(/ 5._8, 5._8, 5._8, 15._8, 5._8, 5._8, 5._8, 5._8, 7._8, 7._8, 7._8, 7._8 /)
  
    xxnpmax = (/ 1.510856726_8, 1.27916225_8, 1.591076159_8, 1.186066584_8, 1.358075681_8,  1.45621905_8, &
                  1.45621905_8, 1.45621905_8, 1.210382326_8, 1.210382326_8,  1.45621905_8, 1.365993164_8, &
                 1.210382326_8,         1._8, 1.399652677_8,          1._8,          1._8 /)
    xq10soil = 1.72_8
    xxkoptlitter = 0.4_8
    xxkoptsoil = (/ 0.33_8, 0.6_8, 0.15_8, 0.6_8, 0.16_8, 0.4_8, 0.3_8, 0.2_8, 0.2_8, 0.25_8,  1._8, &
                    0.65_8, 0.5_8,   2._8, 0.5_8,   1._8,  1._8 /)
    xprodptase = (/  0.5_8, 0.2_8,  0.5_8, 0.5_8,  0.5_8, 0.5_8, 0.5_8, 0.5_8, 0.5_8,  0.5_8, 0.5_8, &
                      4._8, 0.5_8,  0.5_8, 0.5_8,  0.5_8, 0.5_8 /)
    xcostnpup = (/   40._8, 25._8,  40._8, 40._8,  40._8, 40._8, 40._8, 40._8, 40._8,  40._8, 40._8, &
                     40._8, 40._8,  40._8, 40._8,  40._8, 40._8 /)
    xmaxfinelitter = (/ 1524._8, 384._8, 1527._8, 887._8, 157._8, 361._8, 225._8, 913._8, 660._8, 100._8, &
                         100._8, 100._8,  100._8,  83._8, 100._8, 100._8, 100._8 /)
    xmaxcwd = (/ 1795._8, 613._8, 1918._8, 1164._8, 107._8, 420._8, 228._8, 573._8, 811._8, 100._8, &
                  100._8, 100._8,  100._8,   23._8, 100._8, 100._8, 100._8 /)
    xnintercept = (/ 6.32_8, 4.19_8,  6.32_8,  5.73_8, 14.71_8,  6.42_8,    2._8, 14.71_8, 4.71_8, 14.71_8, &
                    14.71_8,   7._8, 14.71_8, 14.71_8, 14.71_8, 14.71_8, 14.71_8 /)
    xnslope = (/ 18.15_8, 26.19_8, 18.15_8, 29.81_8, 23.15_8, 40.96_8, 8._8, 23.15_8, 59.23_8, 23.15_8, &
                 23.15_8,   10._8, 23.15_8, 23.15_8, 23.15_8, 23.15_8, 23.15_8 /)
    xla_to_sa = (/ 5000._8, 5000._8, 5000._8, 5000._8, 5000._8, 5000._8, 5000._8, 5000._8, 5000._8, 5000._8, &
                   5000._8, 5000._8, 5000._8, 5000._8, 5000._8, 5000._8, 5000._8 /)
    xvcmax_scalar = (/ 0.92_8, 1.10_8, 0.92_8, 0.92_8, 1.25_8, 1.25_8, 1.25_8, 1.25_8, 1.0_8, 1.0_8, &
                        1.0_8,  1.0_8,  1.0_8,  1.0_8,  1.0_8,  1.0_8,  1.0_8 /)
    xdisturbance_interval = (/ 100._8, 100._8, 100._8, 100._8, 100._8, 100._8, 100._8, 100._8, 100._8, &
                               100._8, 100._8, 100._8, 100._8, 100._8, 100._8, 100._8, 100._8 /)
    xDAMM_EnzPool = (/ 10.0_8, 10.0_8, 10.0_8, 10.0_8, 10.0_8, 10.0_8, 10.0_8, 10.0_8, 10.0_8, 10.0_8, &
                       10.0_8, 10.0_8, 10.0_8, 10.0_8, 10.0_8, 10.0_8, 10.0_8 /)
    xDAMM_KMO2 = (/ 0.01_8, 0.01_8, 0.01_8, 0.01_8, 0.01_8, 0.01_8, 0.01_8, 0.01_8, 0.01_8, 0.01_8, &
                    0.01_8, 0.01_8, 0.01_8, 0.01_8, 0.01_8, 0.01_8, 0.01_8 /)
    xDAMM_KMcp = (/  0.1_8,  0.1_8,  0.1_8,  0.1_8,  0.1_8,  0.1_8,  0.1_8,  0.1_8,  0.1_8,  0.1_8, &
                     0.1_8,  0.1_8,  0.1_8,  0.1_8,  0.1_8,  0.1_8,  0.1_8 /)
    xDAMM_Ea = (/ 62._8, 62._8, 62._8, 62._8, 62._8, 62._8, 62._8, 62._8, 62._8, 62._8, 62._8, 62._8, &
                  62._8, 62._8, 62._8, 62._8, 62._8 /)
    xDAMM_alpha = (/ 10.6_8, 10.4_8, 10.6_8, 10.6_8, 10.6_8, 10.6_8, 10.6_8, 10.6_8, 10.6_8, 10.6_8, &
                     10.6_8, 10.6_8, 10.6_8, 10.6_8, 10.6_8, 10.6_8, 10.6_8 /)
  
    xxkplab = 0.001369863_8
    xxkpsorb = (/ 1.8356191E-05_8, 2.0547975E-05_8, 1.3698650E-05_8, 1.4794542E-05_8, 2.1369894E-05_8, &
                  2.3835651E-05_8, 1.9452083E-05_8, 2.1095921E-05_8, 2.7123327E-05_8, 2.1095921E-05_8, &
                  2.7123327E-05_8, 2.1095921E-05_8 /)
    xxkpocc = 2.73973E-05_8

    do tile = 1,ntiles
      if ( tdata(tile)%mp>0 ) then
        casabiome(tile)%ivt2     =(/        3,        3,        3,        3,        2,        1,   1,   2, &
                                            1,        1,        0,        0,        0,        1,   0,   0, &
                                            0 /)
        casabiome(tile)%kroot    =(/      5.5_8,      3.9_8,      5.5_8,      3.9_8,      2.0_8,      5.5_8, 5.5_8, 5.5_8, &
                                          5.5_8,      5.5_8,      5.5_8,      5.5_8,      5.5_8,      2.0_8, 2.0_8, 5.5_8, &
                                          5.5_8 /)
        casabiome(tile)%rootdepth=(/      1.5_8,      1.5_8,      1.5_8,      1.5_8,      0.5_8,      0.5_8, 0.5_8, 0.5_8, &
                                          0.5_8,      0.5_8,      0.5_8,      0.5_8,      0.5_8,      0.5_8, 0.5_8, 1.5_8, &
                                          0.5_8 /)
        casabiome(tile)%kuptake  =(/      2.0_8,      1.9_8,      2.0_8,      2.0_8,      1.8_8,      2.0_8, 2.0_8, 2.0_8, &
                                          1.6_8,      1.6_8,      1.6_8,      1.8_8,      1.8_8,      1.8_8, 1.8_8, 1.8_8, &
                                          1.8_8 /)
        casabiome(tile)%krootlen =(/ 14.87805_8, 14.38596_8, 14.02597_8, 18.94737_8, 32.30769_8,      84._8, 84._8, 84._8, &
                                        120.5_8,    120.5_8,       0._8,       0._8,       0._8, 30.76923_8,  0._8,  0._8, &
                                           0._8 /)
        casabiome(tile)%kminN=2.0_8
        casabiome(tile)%kuplabP=0.5_8
        casabiome(tile)%fracnpptoP(:,leaf) =(/ 0.25_8, 0.20_8, 0.40_8, 0.35_8, 0.35_8, 0.35_8, 0.35_8, 0.50_8, 0.50_8, &
                                               0.50_8, 0.50_8, 0.50_8, 0.50_8, 0.25_8, 0.50_8, 0.60_8, 0.50_8 /)
        casabiome(tile)%fracnpptoP(:,wood) =(/ 0.40_8, 0.35_8, 0.30_8, 0.25_8, 0.25_8, 0.00_8, 0.00_8, 0.10_8, 0.00_8, &
                                               0.00_8, 0.00_8, 0.00_8, 0.00_8, 0.25_8, 0.00_8, 0.40_8, 0.00_8 /)
        casabiome(tile)%fracnpptoP(:,xroot)=(/ 0.35_8, 0.45_8, 0.30_8, 0.40_8, 0.40_8, 0.65_8, 0.65_8, 0.40_8, 0.50_8, &
                                               0.50_8, 0.50_8, 0.50_8, 0.50_8, 0.50_8, 0.50_8, 0.00_8, 0.50_8 /)
        casabiome(tile)%rmplant(:,leaf)    =0.1_8
        casabiome(tile)%rmplant(:,wood)    =(/ 2.0_8, 1.0_8, 1.5_8, 0.8_8, 0.5_8, 0.5_8, 0.4_8, 1.8_8, 2.0_8, 1.0_8, &
                                               1.0_8, 1.0_8, 1.0_8, 2.0_8, 1.0_8, 1.0_8, 1.0_8 /)
        casabiome(tile)%rmplant(:,xroot)   =(/ 10._8, 2.0_8, 7.5_8, 2.5_8, 4.5_8, 4.5_8, 4.0_8, 15._8, 25._8, 10._8, &
                                               10._8, 10._8, 10._8, 10._8, 10._8, 10._8, 10._8 /)
        casabiome(tile)%ftransNPtoL(:,leaf) =0.5_8
        casabiome(tile)%ftransNPtoL(:,wood) =0.95_8
        casabiome(tile)%ftransNPtoL(:,xroot)=0.9_8
        casabiome(tile)%fracligninplant(:,leaf) =(/ 0.25_8, 0.20_8, 0.20_8, 0.20_8, 0.20_8, 0.10_8, 0.10_8, 0.10_8, &
                                                    0.10_8, 0.10_8, 0.15_8, 0.15_8, 0.15_8, 0.15_8, 0.15_8, 0.25_8, &
                                                    0.10_8 /)
        casabiome(tile)%fracligninplant(:,wood) =0.4_8
        casabiome(tile)%fracligninplant(:,xroot)=(/ 0.25_8, 0.20_8, 0.20_8, 0.20_8, 0.20_8, 0.10_8, 0.10_8, 0.10_8, &
                                                    0.10_8, 0.10_8, 0.15_8, 0.15_8, 0.15_8, 0.15_8, 0.15_8, 0.25_8, &
                                                    0.10_8 /)
        casabiome(tile)%glaimax=(/ 10._8, 10._8, 10._8, 10._8, 10._8, 3._8, 3._8, 3._8, 6._8, 6._8,  5._8, 5._8, &
                                    5._8, 1._8,  6._8,   1._8,  0._8 /)
        casabiome(tile)%glaimin=(/ 1._8,  1._8, .5_8,  .5_8, .1_8, .1_8, .1_8, .1_8, .1_8, .1_8, .05_8, .05_8, &
                                  .05_8, .05_8, 0._8, .05_8, 0._8 /)
        phen(tile)%TKshed=(/ 268._8,   260._8, 263.15_8, 268.15_8, 277.15_8, 275.15_8, 275.15_8, 275.15_8, 278.15_8, &
                           278.15_8, 277.15_8, 277.15_8, 277.15_8, 277.15_8, 277.15_8, 277.15_8, 283.15_8 /)
        casabiome(tile)%xkleafcoldexp=3._8
        casabiome(tile)%xkleafdryexp=3._8
        casabiome(tile)%ratioNCplantmin(:,leaf) =(/     0.02_8,     0.04_8, 0.016667_8, 0.028571_8,    0.025_8,  0.02631_8, &
                                                        0.02_8,     0.02_8,     0.04_8,     0.04_8, 0.033333_8,    0.025_8, &
                                                       0.025_8, 0.018182_8,    0.025_8,    0.025_8,    0.025_8 /)
        casabiome(tile)%ratioNCplantmax(:,leaf) =(/    0.024_8,    0.048_8,     0.02_8, 0.034286_8,     0.03_8, 0.031572_8, &
                                                       0.024_8,    0.024_8,    0.048_8,    0.048_8,     0.04_8,     0.03_8, &
                                                        0.03_8, 0.022222_8,     0.03_8,     0.03_8,     0.03_8 /)
        casabiome(tile)%ratioNCplantmin(:,wood) =(/    0.004_8, 0.006667_8,    0.004_8, 0.005714_8, 0.006667_8, 0.006667_8, &
                                                    0.006667_8, 0.006667_8,    0.008_8,    0.008_8, 0.006667_8, 0.006667_8, &
                                                    0.006667_8, 0.006667_8, 0.006667_8, 0.007307_8, 0.006667_8 /)
        casabiome(tile)%ratioNCplantmax(:,wood) =(/   0.0048_8,    0.008_8,   0.0048_8, 0.006857_8,    0.008_8,    0.008_8, &
                                                       0.008_8,    0.008_8,   0.0096_8,   0.0096_8,    0.008_8,    0.008_8, &
                                                       0.008_8,    0.008_8,    0.008_8, 0.008889_8,    0.008_8 /)
        casabiome(tile)%ratioNCplantmin(:,xroot)=(/ 0.012821_8, 0.014706_8, 0.012821_8, 0.014085_8, 0.014085_8, 0.014085_8, &
                                                    0.014085_8, 0.014085_8, 0.014085_8, 0.014085_8, 0.014085_8, 0.014085_8, &
                                                    0.014085_8, 0.014085_8, 0.014085_8, 0.014085_8, 0.014085_8 /)
        casabiome(tile)%ratioNCplantmax(:,xroot)=(/ 0.015385_8, 0.017647_8, 0.015385_8, 0.016901_8, 0.016901_8, 0.016901_8, &
                                                    0.016901_8, 0.016901_8, 0.016901_8, 0.016901_8, 0.016901_8, 0.016901_8, &
                                                    0.016901_8, 0.016901_8, 0.016901_8, 0.016901_8, 0.016901_8 /)
        casabiome(tile)%ftransPPtoL(:,leaf)=0.5_8
        casabiome(tile)%ftransPPtoL(:,wood)=0.95_8
        casabiome(tile)%ftransPPtoL(:,xroot)=0.9_8
       
      end if 
    end do   

  end if

end if

do tile = 1,ntiles
  if ( tdata(tile)%mp>0 ) then
    casabiome(tile)%ratioPcplantmin(:,leaf)  = 1._8/(xratioNPleafmax*ratioCNplant(:,leaf))
    casabiome(tile)%ratioPcplantmax(:,leaf)  = 1._8/(xratioNPleafmin*ratioCNplant(:,leaf))
    casabiome(tile)%ratioPcplantmin(:,wood)  = 1._8/(xratioNPwoodmax*ratioCNplant(:,wood))
    casabiome(tile)%ratioPcplantmax(:,wood)  = 1._8/(xratioNPwoodmin*ratioCNplant(:,wood))
    casabiome(tile)%ratioPcplantmin(:,xroot) = 1._8/(xratioNPfrootmax*ratioCNplant(:,xroot))
    casabiome(tile)%ratioPcplantmax(:,xroot) = 1._8/(xratioNPfrootmin*ratioCNplant(:,xroot))

    casabiome(tile)%ratioNPplantmin(:,leaf)  = xratioNPleafmin
    casabiome(tile)%ratioNPplantmax(:,leaf)  = xratioNPleafmax
    casabiome(tile)%ratioNPplantmin(:,wood)  = xratioNPwoodmin
    casabiome(tile)%ratioNPplantmax(:,wood)  = xratioNPwoodmax
    casabiome(tile)%ratioNPplantmin(:,xroot) = xratioNPfrootmin
    casabiome(tile)%ratioNPplantmax(:,xroot) = xratioNPfrootmax    

    casabiome(tile)%sla                = slax
    casabiome(tile)%fraclabile(:,leaf) = deltcasa*0.6_8    !1/day
    casabiome(tile)%fraclabile(:,xroot)= deltcasa*0.4_8    !1/day
    casabiome(tile)%fraclabile(:,wood) = 0._8
    casabiome(tile)%plantrate(:,leaf)  = deltcasa/(leafage*(1._8-xfherbivore))
    casabiome(tile)%plantrate(:,xroot) = deltcasa/frootage
    casabiome(tile)%plantrate(:,wood)  = deltcasa/woodage
    casabiome(tile)%litterrate(:,metb) = deltcasa/metage
    casabiome(tile)%litterrate(:,str)  = deltcasa/strage
    casabiome(tile)%litterrate(:,cwd)  = deltcasa/cwdage
    casabiome(tile)%soilrate(:,mic)    = deltcasa/micage
    casabiome(tile)%soilrate(:,slow)   = deltcasa/slowage
    casabiome(tile)%soilrate(:,pass)   = deltcasa/passage
    casabiome(tile)%xkleafcoldmax      = deltcasa*xxkleafcoldmax
    casabiome(tile)%xkleafdrymax       = deltcasa*xxkleafdrymax
    casabiome(tile)%rmplant            = deltcasa*casabiome(tile)%rmplant
    casabiome(tile)%kclabrate          = deltcasa/clabileage

    casabiome(tile)%xnpmax(:)          = xxnpmax(:)
    casabiome(tile)%q10soil(:)         = xq10soil(:)
    casabiome(tile)%xkoptlitter(:)     = xxkoptlitter(:)
    casabiome(tile)%xkoptsoil(:)       = xxkoptsoil(:)
    casabiome(tile)%prodptase(:)       = xprodptase(:)/365._8   ! convert from yearly to daily
    casabiome(tile)%costnpup(:)        = xcostnpup(:)
    casabiome(tile)%maxfinelitter(:)   = xmaxfinelitter(:)
    casabiome(tile)%maxcwd(:)          = xmaxcwd(:)
    casabiome(tile)%nintercept(:)      = xnintercept(:)
    casabiome(tile)%nslope(:)          = xnslope(:)    

    !casabiome(tile)%la_to_sa(:)             = xla_to_sa(:)
    !casabiome(tile)%vcmax_scalar(:)         = xvcmax_scalar(:)
    !casabiome(tile)%disturbance_interval(:) = xdisturbance_interval(:)
    !casabiome(tile)%DAMM_EnzPool(:)         = xDAMM_EnzPool(:)p
    !casabiome(tile)%DAMM_KMO2(:)            = xDAMM_KMO2(:)
    !casabiome(tile)%DAMM_KMcp(:)            = xDAMM_KMcp(:)
    !casabiome(tile)%DAMM_Ea(:)              = xDAMM_Ea(:)
    !casabiome(tile)%DAMM_alpha(:)           = xDAMM_alpha(:)

    casabiome(tile)%xkplab = xxkplab
    casabiome(tile)%xkpsorb = xxkpsorb
    casabiome(tile)%xkpocc = xxkpocc

    casamet(tile)%iveg2 = casabiome(tile)%ivt2(veg(tile)%iveg)
    where (casamet(tile)%iveg2==forest.or.casamet(tile)%iveg2==shrub)
      casamet(tile)%lnonwood = 0
      casapool(tile)%cplant(:,wood)  = cwood(veg(tile)%iveg) 
      casapool(tile)%clitter(:,cwd)  = ccwd(veg(tile)%iveg)
      casapool(tile)%nplant(:,wood)  = nwood(veg(tile)%iveg) 
      casapool(tile)%nlitter(:,cwd)  = ncwd(veg(tile)%iveg)
      casapool(tile)%pplant(:,wood)  = xpwood(veg(tile)%iveg)
      casapool(tile)%plitter(:,cwd)  = xpcwd(veg(tile)%iveg)
    elsewhere
      casamet(tile)%lnonwood = 1
      casapool(tile)%cplant(:,wood)  = 0._8
      casapool(tile)%clitter(:,cwd)  = 0._8
      casapool(tile)%nplant(:,wood)  = 0._8
      casapool(tile)%nlitter(:,cwd)  = 0._8
      casapool(tile)%pplant(:,wood)  = 0._8
      casapool(tile)%plitter(:,cwd)  = 0._8
    end where
    if ( cable_pop==1 ) then
      where (casamet(tile)%iveg2==forest.or.casamet(tile)%iveg2==shrub)
        casapool(tile)%cplant(:,wood)  = 0.01_8
        casapool(tile)%nplant(:,wood)  = casabiome(tile)%ratioNCplantmin(veg(tile)%iveg,wood)*casapool(tile)%cplant(:,wood)
        casapool(tile)%pplant(:,wood)  = casabiome(tile)%ratioPCplantmin(veg(tile)%iveg,wood)*casapool(tile)%cplant(:,wood)
      end where
    end if
    casapool(tile)%cplant(:,leaf)     = cleaf(veg(tile)%iveg)
    casapool(tile)%cplant(:,xroot)    = cfroot(veg(tile)%iveg)
    casapool(tile)%clabile            = 0._8
    casapool(tile)%clitter(:,metb)    = cmet(veg(tile)%iveg)
    casapool(tile)%clitter(:,str)     = cstr(veg(tile)%iveg)
    casapool(tile)%csoil(:,mic)       = cmic(veg(tile)%iveg)
    casapool(tile)%csoil(:,slow)      = cslow(veg(tile)%iveg)
    casapool(tile)%csoil(:,pass)      = cpass(veg(tile)%iveg)
    if ( ccycle==1 ) then
      casapool(tile)%ratioNCplant     = 1._8/ratioCNplant(veg(tile)%iveg,:)  
    end if
    casapool(tile)%dclabiledt         = 0._8

    ! initializing glai in case not reading pool file (eg. during spin)
    casamet(tile)%glai = max(casabiome(tile)%glaimin(veg(tile)%iveg), &
        casabiome(tile)%sla(veg(tile)%iveg)*casapool(tile)%cplant(:,leaf))
    casaflux(tile)%fNminloss   = xfNminloss(veg(tile)%iveg)
    casaflux(tile)%fNminleach  = 10._8*xfNminleach(veg(tile)%iveg)*deltcasa
    casapool(tile)%nplant(:,leaf) = nleaf(veg(tile)%iveg)
    casapool(tile)%nplant(:,xroot)= nfroot(veg(tile)%iveg)
    casapool(tile)%nlitter(:,metb)= nmet(veg(tile)%iveg)
    casapool(tile)%nlitter(:,str) = cstr(veg(tile)%iveg)*ratioNCstrfix
    casapool(tile)%nsoil(:,mic)   = nmic(veg(tile)%iveg)
    casapool(tile)%nsoil(:,slow)  = nslow(veg(tile)%iveg)
    casapool(tile)%nsoil(:,pass)  = npass(veg(tile)%iveg) 
    casapool(tile)%nsoilmin       = xnsoilmin(veg(tile)%iveg) 
    casapool(tile)%pplant(:,leaf) = xpleaf(veg(tile)%iveg)
    casapool(tile)%pplant(:,xroot)= xpfroot(veg(tile)%iveg) 
    casapool(tile)%plitter(:,metb)= xpmet(veg(tile)%iveg)
    casapool(tile)%plitter(:,str) = casapool(tile)%nlitter(:,str)/ratioNPstrfix
    casapool(tile)%psoil(:,mic)   = xpmic(veg(tile)%iveg)
    casapool(tile)%psoil(:,slow)  = xpslow(veg(tile)%iveg)
    casapool(tile)%psoil(:,pass)  = xppass(veg(tile)%iveg)
    casapool(tile)%psoillab       = xplab(veg(tile)%iveg)
    casapool(tile)%psoilsorb      = xpsorb(veg(tile)%iveg)
    casapool(tile)%psoilocc       = xpocc(veg(tile)%iveg)
    casaflux(tile)%kmlabp         = xkmlabp(casamet(tile)%isorder)
    casaflux(tile)%psorbmax       = xpsorbmax(casamet(tile)%isorder)
    casaflux(tile)%fpleach        = xfPleach(casamet(tile)%isorder)/365._8

    casapool(tile)%ratioNCplant   = 1._8/ratioCNplant(veg(tile)%iveg,:)
    casapool(tile)%ratioNPplant   = casabiome(tile)%ratioNPplantmin(veg(tile)%iveg,:)
    casapool(tile)%ratioNClitter  = casapool(tile)%nlitter/(casapool(tile)%clitter+1.0e-10_8)
    casapool(tile)%ratioNPlitter  = casapool(tile)%nlitter/(casapool(tile)%plitter+1.0e-10_8)
    casapool(tile)%ratioNCsoil    = 1._8/ratioCNsoil(veg(tile)%iveg,:)
    casapool(tile)%ratioNPsoil    = ratioNPsoil(casamet(tile)%isorder,:)
    casapool(tile)%ratioNCsoilmin = 1._8/ratioCNsoilmax(veg(tile)%iveg,:)
    casapool(tile)%ratioNCsoilmax = 1._8/ratioCNsoilmin(veg(tile)%iveg,:)
    casapool(tile)%ratioNCsoilnew = casapool(tile)%ratioNCsoilmax

    casapool(tile)%ratioPCplant   = casabiome(tile)%ratioPcplantmax(veg(tile)%iveg,:)
    casapool(tile)%ratioPClitter  = casapool(tile)%plitter/(casapool(tile)%clitter(:,:)+1.0e-10_8)
    casapool(tile)%ratioPCsoil    = 1._8/(ratioCNsoil(veg(tile)%iveg,:)*ratioNPsoil(casamet(tile)%isorder,:))

    if ( ccycle<2 ) then
      casapool(tile)%Nplant         = casapool(tile)%Cplant*casapool(tile)%ratioNCplant
      casapool(tile)%Nsoil          = casapool(tile)%ratioNCsoil*casapool(tile)%Csoil
    end if
    if ( ccycle<3 ) then
      casapool(tile)%Psoil          = casapool(tile)%Nsoil/casapool(tile)%ratioNPsoil
      casapool(tile)%psoilsorb      = casaflux(tile)%psorbmax*casapool(tile)%psoillab &
                              /(casaflux(tile)%kmlabp+casapool(tile)%psoillab)
    end if
  
  end if
end do  

end subroutine casa_readbiome

subroutine cable_biophysic_parm

use cc_mpi 
use darcdf_m
use infile
use newmpar_m, only : ntiles
use nsibd_m

integer k, chid_len
integer tile
integer, dimension(2) :: nstart, ncount
integer, dimension(:), allocatable, save :: csiropft
real totdepth
real, dimension(:), allocatable, save :: hc, xfang, leaf_w, leaf_l, canst1
real, dimension(:), allocatable, save :: shelrb, extkn, vcmax, rpcoef
real, dimension(:), allocatable, save :: rootbeta, c4frac, vbeta
real, dimension(:), allocatable, save :: a1gs, d0gs, alpha, convex, cfrd
real, dimension(:), allocatable, save :: gswmin, conkc0, conko0, ekc, eko, g0, g1
real, dimension(:), allocatable, save :: zr, clitt
real, dimension(:,:), allocatable, save :: refl, taul, froot2

if ( lncveg_numpft<1 ) then
    
  ! default biophysical parameter tables
  if ( myid==0 ) then
    write(6,*) "-> Using default CABLE biophysical parameter tables"
  end if
  lncveg_numpft = 18
  allocate( csiropft(lncveg_numpft), hc(lncveg_numpft), xfang(lncveg_numpft), leaf_w(lncveg_numpft), leaf_l(lncveg_numpft) )
  allocate( canst1(lncveg_numpft), shelrb(lncveg_numpft), extkn(lncveg_numpft), refl(lncveg_numpft,2), taul(lncveg_numpft,2) )
  allocate( vcmax(lncveg_numpft), rpcoef(lncveg_numpft), rootbeta(lncveg_numpft), c4frac(lncveg_numpft) )
  allocate( froot2(lncveg_numpft,cbm_ms), vbeta(lncveg_numpft) )
  allocate( a1gs(lncveg_numpft), d0gs(lncveg_numpft), alpha(lncveg_numpft), convex(lncveg_numpft), cfrd(lncveg_numpft) )
  allocate( gswmin(lncveg_numpft), conkc0(lncveg_numpft), conko0(lncveg_numpft), ekc(lncveg_numpft), eko(lncveg_numpft) )
  allocate( g0(lncveg_numpft), g1(lncveg_numpft), zr(lncveg_numpft), clitt(lncveg_numpft) )
  allocate( ivegt_name(lncveg_numpft) )
  
  csiropft=(/ 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 2 /)
  hc    =(/   17.,  35.,  15.5,  20.,   0.6, 0.567, 0.567, 0.567, 0.55, 0.55, 0.567,  0.2, 6.017,  0.2,  0.2,  0.2,  0.2, 17. /)
  xfang =(/  0.01,  0.1,  0.01, 0.25,  0.01,  -0.3,  -0.3,  -0.3, -0.3, -0.3,  -0.3,  0.1,    0.,   0.,   0.,   0.,   0., 0.1 /)
  leaf_w=(/ 0.001, 0.05, 0.001, 0.08, 0.005,  0.01,  0.01,  0.01, 0.01, 0.01,  0.01, 0.03, 0.015, 0.00,   0.,   0.,   0., 0.05 /)
  leaf_l=(/ 0.055, 0.10, 0.040, 0.15, 0.100,  0.30,  0.30,  0.30, 0.30, 0.30,  0.30, 0.30, 0.242, 0.03, 0.03, 0.03, 0.03, 0.10 /)
  canst1=0.1
  shelrb=2.
  extkn=0.001
  refl(:,1)=(/ 0.062,0.076,0.056,0.092,0.100,0.110,0.100,0.117,0.100,0.090,0.108,0.055,0.091,0.238,0.143,0.143,0.159,0.076 /)
  refl(:,2)=(/ 0.302,0.350,0.275,0.380,0.400,0.470,0.400,0.343,0.400,0.360,0.343,0.190,0.310,0.457,0.275,0.275,0.305,0.350 /)
  taul(:,1)=(/ 0.050,0.050,0.045,0.050,0.050,0.070,0.100,0.080,0.100,0.090,0.075,0.023,0.059,0.039,0.023,0.023,0.026,0.050 /)
  taul(:,2)=(/ 0.100,0.250,0.144,0.250,0.240,0.250,0.150,0.124,0.150,0.225,0.146,0.198,0.163,0.189,0.113,0.113,0.113,0.250 /)
  vcmax=(/ 40.E-6,55.E-6,40.E-6,60.E-6,40.E-6,60.E-6,10.E-6,40.E-6,80.E-6,80.E-6,60.E-6,17.E-6,1.E-6,17.E-6,17.E-6,17.E-6, &
           17.E-6,55.E-6 /)
  rpcoef=0.0832
  rootbeta=(/ 0.943,0.962,0.966,0.961,0.964,0.943,0.943,0.943,0.961,0.961,0.943,0.975,0.961,0.961,0.961,0.961,0.961,0.962 /)
  c4frac=(/ 0., 0., 0., 0., 0., 0., 1., 0., 0., 1., 0., 0., 0., 0., 0., 0., 0., 0. /)
  vbeta=(/ 2., 2., 2., 2., 4., 4., 4., 4., 2., 2., 4., 4., 2., 4., 4., 4., 4., 2. /)
  a1gs=(/ 9., 9., 9., 9., 9., 9., 4., 9., 9., 4., 9., 9., 9., 9., 9., 9., 9., 9. /)
  d0gs=1500.
  alpha=(/ 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.05, 0.2, 0.2, 0.05, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2 /)
  convex=(/ 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.8, 0.01, 0.01, 0.8, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01 /)
  cfrd=(/ 0.015, 0.015, 0.015, 0.015, 0.015, 0.015, 0.025, 0.015, 0.015, 0.025, 0.015, 0.015, 0.015, 0.015, 0.015, 0.015, &
          0.015, 0.015 /)
  gswmin=(/ 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.04, 0.01, 0.01, 0.04, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01 /)
  conkc0=302.e-6
  conko0=256.e-3
  ekc=59430.
  eko=36000.
  g0=0.
  g1=(/ 2.346064, 4.114762, 2.346064, 4.447321, 4.694803, 5.248500, 1.616178, 2.222156, 5.789377, 1.616178, 5.248500, 5.248500, &
        0.000000, 5.248500, 5.248500, 5.248500, 5.248500, 2.346064 /)
  zr=(/ 1.8, 3., 2., 2., 2.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.8, 3.1, 3., 1., 1., 1., 1., 3. /)
  clitt=(/ 20., 6., 10., 13., 2., 2., 0.3, 0.3, 0., 0., 2., 2., 0., 0., 0., 0., 0., 6. /) 
  
  !ivegt_name(0) = "water"
  ivegt_name(1) = "evergreen_needleleaf"
  ivegt_name(2) = "evergreen_broadleaf"
  ivegt_name(3) = "deciduous_needleleaf"
  ivegt_name(4) = "deciduous_broadleaf"
  ivegt_name(5) = "shrub"
  ivegt_name(6) = "C3_grassland"
  ivegt_name(7) = "C4_grassland"
  ivegt_name(8) = "tundra"
  ivegt_name(9) = "C3_cropland"
  ivegt_name(10) = "C4_cropland"
  ivegt_name(11) = "wetland"
  ivegt_name(12) = "empty"
  ivegt_name(13) = "empty"
  ivegt_name(14) = "barren"
  ivegt_name(15) = "urban"
  ivegt_name(16) = "lakes"
  ivegt_name(17) = "ice"
  ivegt_name(18) = "evergreen_broadleaf_sava"
  
  do k = 1,lncveg_numpft
    if ( myid==0 ) then
      write(6,*) " -> Define ",trim(ivegt_name(k))
    end if
  end do  
  
else

  ! user defined biophysical parameter tables  
  if ( myid==0 ) then
    write(6,*) "-> Using user defined CABLE biophysical parameter tables"
  end if
  allocate( csiropft(lncveg_numpft), hc(lncveg_numpft), xfang(lncveg_numpft), leaf_w(lncveg_numpft), leaf_l(lncveg_numpft) )
  allocate( canst1(lncveg_numpft), shelrb(lncveg_numpft), extkn(lncveg_numpft), refl(lncveg_numpft,2), taul(lncveg_numpft,2) )
  allocate( vcmax(lncveg_numpft), rpcoef(lncveg_numpft), rootbeta(lncveg_numpft), c4frac(lncveg_numpft) )
  allocate( froot2(lncveg_numpft,cbm_ms), vbeta(lncveg_numpft) )
  allocate( a1gs(lncveg_numpft), d0gs(lncveg_numpft), alpha(lncveg_numpft), convex(lncveg_numpft), cfrd(lncveg_numpft) )
  allocate( gswmin(lncveg_numpft), conkc0(lncveg_numpft), conko0(lncveg_numpft), ekc(lncveg_numpft), eko(lncveg_numpft) )
  allocate( g0(lncveg_numpft), g1(lncveg_numpft), zr(lncveg_numpft), clitt(lncveg_numpft) )
  allocate( ivegt_name(lncveg_numpft) )

  ! set default values to prevent false float invalid when calling Bcast.
  csiropft = 0.
  hc = 0.
  xfang = 0.
  leaf_w = 0.
  leaf_l = 0.
  canst1 = 0.
  shelrb = 0.
  extkn = 0.
  refl = 0.
  taul = 0.
  vcmax = 0.
  rpcoef = 0.
  rootbeta = 0.
  c4frac = 0.
  vbeta = 0.
  a1gs = 0.
  d0gs = 0.
  alpha = 0.
  convex = 0.
  cfrd = 0.
  gswmin = 0.
  conkc0 = 0.
  conko0 = 0.
  ekc = 0.
  eko = 0.
  g0 = 0.
  g1 = 0.
  zr = 0.
  clitt = 0.
  ivegt_name = ""
  
  if ( myid==0 ) then
    chid_len = 0  
    call ccnf_inq_dimlen(ncidveg,'chid',chid_len,failok=.true.)
    do k = 1,lncveg_numpft
      nstart(1) = 1
      nstart(2) = k
      ncount(1) = chid_len
      ncount(2) = 1
      call ccnf_get_vara(ncidveg,'pftname',nstart,ncount,ivegt_name(k))
    end do      
    nstart(1) = 1
    ncount(1) = lncveg_numpft
    call ccnf_get_vara(ncidveg,'csiropft',nstart,ncount,csiropft)
    call ccnf_get_vara(ncidveg,'hc',nstart,ncount,hc)
    call ccnf_get_vara(ncidveg,'xfang',nstart,ncount,xfang)
    call ccnf_get_vara(ncidveg,'leaf_w',nstart,ncount,leaf_w)
    call ccnf_get_vara(ncidveg,'leaf_l',nstart,ncount,leaf_l)
    call ccnf_get_vara(ncidveg,'canst1',nstart,ncount,canst1)
    call ccnf_get_vara(ncidveg,'shelrb',nstart,ncount,shelrb)
    call ccnf_get_vara(ncidveg,'extkn',nstart,ncount,extkn)
    call ccnf_get_vara(ncidveg,'rholeaf-vis',nstart,ncount,refl(:,1))
    call ccnf_get_vara(ncidveg,'rholeaf-nir',nstart,ncount,refl(:,2))
    call ccnf_get_vara(ncidveg,'tauleaf-vis',nstart,ncount,taul(:,1))
    call ccnf_get_vara(ncidveg,'tauleaf-nir',nstart,ncount,taul(:,2))
    call ccnf_get_vara(ncidveg,'vcmax',nstart,ncount,vcmax)
    call ccnf_get_vara(ncidveg,'rpcoef',nstart,ncount,rpcoef)
    call ccnf_get_vara(ncidveg,'rootbeta',nstart,ncount,rootbeta)
    call ccnf_get_vara(ncidveg,'c4frac',nstart,ncount,c4frac)
    call ccnf_get_vara(ncidveg,'vbeta',nstart,ncount,vbeta)
    call ccnf_get_vara(ncidveg,'a1gs',nstart,ncount,a1gs)
    call ccnf_get_vara(ncidveg,'d0gs',nstart,ncount,d0gs)
    call ccnf_get_vara(ncidveg,'alpha',nstart,ncount,alpha)
    call ccnf_get_vara(ncidveg,'convex',nstart,ncount,convex)
    call ccnf_get_vara(ncidveg,'cfrd',nstart,ncount,cfrd)
    call ccnf_get_vara(ncidveg,'gswmin',nstart,ncount,gswmin)
    call ccnf_get_vara(ncidveg,'conkc0',nstart,ncount,conkc0)
    call ccnf_get_vara(ncidveg,'conko0',nstart,ncount,conko0)
    call ccnf_get_vara(ncidveg,'ekc',nstart,ncount,ekc)
    call ccnf_get_vara(ncidveg,'eko',nstart,ncount,eko)
    call ccnf_get_vara(ncidveg,'g0',nstart,ncount,g0)
    call ccnf_get_vara(ncidveg,'g1',nstart,ncount,g1)
    call ccnf_get_vara(ncidveg,'zr',nstart,ncount,zr)
    call ccnf_get_vara(ncidveg,'clitt',nstart,ncount,clitt)
  end if
  call ccmpi_bcast(csiropft,0,comm_world)  
  do k = 1,lncveg_numpft
    if ( myid==0 ) then
      write(6,*) "--> Found ",trim(ivegt_name(k))
    end if
    call ccmpi_bcast(ivegt_name(k),0,comm_world)  
  end do
  call ccmpi_bcast(hc,0,comm_world)
  call ccmpi_bcast(xfang,0,comm_world)
  call ccmpi_bcast(leaf_w,0,comm_world)
  call ccmpi_bcast(leaf_l,0,comm_world)
  call ccmpi_bcast(canst1,0,comm_world)
  call ccmpi_bcast(shelrb,0,comm_world)
  call ccmpi_bcast(extkn,0,comm_world)
  call ccmpi_bcast(refl,0,comm_world)
  call ccmpi_bcast(taul,0,comm_world)
  call ccmpi_bcast(vcmax,0,comm_world)
  call ccmpi_bcast(rpcoef,0,comm_world)
  call ccmpi_bcast(rootbeta,0,comm_world)
  call ccmpi_bcast(c4frac,0,comm_world)
  call ccmpi_bcast(vbeta,0,comm_world)
  call ccmpi_bcast(a1gs,0,comm_world)
  call ccmpi_bcast(d0gs,0,comm_world)
  call ccmpi_bcast(alpha,0,comm_world)
  call ccmpi_bcast(convex,0,comm_world)
  call ccmpi_bcast(cfrd,0,comm_world)
  call ccmpi_bcast(gswmin,0,comm_world)
  call ccmpi_bcast(conkc0,0,comm_world)
  call ccmpi_bcast(conko0,0,comm_world)
  call ccmpi_bcast(ekc,0,comm_world)
  call ccmpi_bcast(eko,0,comm_world)
  call ccmpi_bcast(g0,0,comm_world)
  call ccmpi_bcast(g1,0,comm_world)
  call ccmpi_bcast(zr,0,comm_world)
  call ccmpi_bcast(clitt,0,comm_world)
  
end if

do tile = 1,ntiles
  if ( tdata(tile)%mp>0 ) then  

    ! froot is now calculated from soil depth and the new parameter rootbeta 
    ! according to Jackson et al. 1996, Oceologica, 108:389-411
    totdepth = 0.
    do k = 1,cbm_ms
      totdepth = totdepth + real(soil(tile)%zse(k))*100.
      froot2(:,k) = min(1.,1.-rootbeta(:)**totdepth)
    end do
    do k = cbm_ms-1, 2, -1
      froot2(:,k) = froot2(:,k) - froot2(:,k-1)
    end do
    froot2(:,cbm_ms) = 1.-sum(froot2(:,1:cbm_ms-1),2)
  
    ! Eva's method for ACCESS1.3
    !froot2(:,1)=0.05
    !froot2(:,2)=0.20
    !froot2(:,3)=0.20
    !froot2(:,4)=0.20
    !froot2(:,5)=0.20
    !froot2(:,6)=0.15

    if ( maxval(tdata(tile)%cveg)>lncveg_numpft .or. minval(tdata(tile)%cveg)<1 ) then
      write(6,*) "ERROR: Invalid range of vegetation classes for CABLE"
      write(6,*) "cveg min,max           = ",minval(tdata(tile)%cveg),maxval(tdata(tile)%cveg)
      write(6,*) "Expected range min,max = ",1,lncveg_numpft
      call ccmpi_abort(-1)
    end if

    veg(tile)%meth      = 1
    veg(tile)%iveg      = csiropft(tdata(tile)%cveg)
    veg(tile)%hc        = real(hc(tdata(tile)%cveg),8)
    veg(tile)%xfang     = real(xfang(tdata(tile)%cveg),8)  
    veg(tile)%dleaf     = real(sqrt(max(leaf_w(tdata(tile)%cveg)*leaf_l(tdata(tile)%cveg),1.e-20)),8)
    veg(tile)%canst1    = real(canst1(tdata(tile)%cveg),8)
    veg(tile)%shelrb    = real(shelrb(tdata(tile)%cveg),8)
    veg(tile)%extkn     = real(extkn(tdata(tile)%cveg),8)
    veg(tile)%refl(:,1) = real(refl(tdata(tile)%cveg,1),8)
    veg(tile)%refl(:,2) = real(refl(tdata(tile)%cveg,2),8)  
    veg(tile)%taul(:,1) = real(taul(tdata(tile)%cveg,1),8)
    veg(tile)%taul(:,2) = real(taul(tdata(tile)%cveg,2),8)  
    veg(tile)%vcmax     = real(vcmax(tdata(tile)%cveg),8)
    veg(tile)%ejmax     = real(2.*veg(tile)%vcmax,8)
    veg(tile)%rpcoef    = real(rpcoef(tdata(tile)%cveg),8)
    do k = 1,cbm_ms
      veg(tile)%froot(:,k)=real(froot2(tdata(tile)%cveg,k),8)
    end do
    veg(tile)%frac4     = real(c4frac(tdata(tile)%cveg),8)
    veg(tile)%xalbnir   = 1._8 ! not used
    veg(tile)%vbeta     = real(vbeta(tdata(tile)%cveg),8)
    veg(tile)%a1gs      = real(a1gs(tdata(tile)%cveg),8)   
    veg(tile)%d0gs      = real(d0gs(tdata(tile)%cveg),8)
    veg(tile)%g0        = real(g0(tdata(tile)%cveg),8)
    veg(tile)%g1        = real(g1(tdata(tile)%cveg),8)
    veg(tile)%alpha     = real(alpha(tdata(tile)%cveg),8)
    veg(tile)%convex    = real(convex(tdata(tile)%cveg),8) 
    veg(tile)%cfrd      = real(cfrd(tdata(tile)%cveg),8)
    veg(tile)%gswmin    = real(gswmin(tdata(tile)%cveg),8)
    veg(tile)%conkc0    = real(conkc0(tdata(tile)%cveg),8)
    veg(tile)%conko0    = real(conko0(tdata(tile)%cveg),8)
    veg(tile)%ekc       = real(ekc(tdata(tile)%cveg),8)
    veg(tile)%eko       = real(eko(tdata(tile)%cveg),8)
    veg(tile)%zr        = real(zr(tdata(tile)%cveg),8)
    veg(tile)%clitt     = real(clitt(tdata(tile)%cveg),8)
  
    veg(tile)%gamma     = 3.e-2_8
    veg(tile)%F10       = 0.85_8
    !veg(tile)%ZR        = 5._8
    veg(tile)%disturbance_interval = 100
    veg(tile)%disturbance_intensity = 0._8
  
    ! depeciated
    !veg(tile)%tminvj    = real(tminvj(veg(tile)%iveg),8)
    !veg(tile)%tmaxvj    = real(tmaxvj(veg(tile)%iveg),8)
    !veg(tile)%rp20      = real(rp20(veg(tile)%iveg),8)
    !veg(tile)%rs20      = real(rs20(veg(tile)%iveg),8)
    !veg(tile)%vegcf     = real(vegcf(veg(tile)%iveg),8)
  
    ! patch
    if ( gs_switch==1 ) then
      where ( veg(tile)%g0<1.e-8_8 ) 
        veg(tile)%g0 = 0.01_8
      end where
    end if    
    
  end if
end do  

deallocate( csiropft, hc, xfang, leaf_w, leaf_l )
deallocate( canst1, shelrb, extkn, refl, taul )
deallocate( vcmax, rpcoef, rootbeta, c4frac, froot2 )
deallocate( vbeta )
deallocate( a1gs, d0gs, alpha, convex, cfrd )
deallocate( gswmin, conkc0, conko0, ekc, eko, g0, g1 )
deallocate( zr, clitt )

return
end subroutine cable_biophysic_parm
                   
subroutine cable_soil_parm

use cc_mpi
use darcdf_m
use infile
use newmpar_m
use nsibd_m
use parm_m, only : nmaxpr
use soilv_m

integer isoil, k, chid_len
integer tile
integer, dimension(2) :: nstart, ncount

if ( lncveg_numsoil<1 ) then
    
  ! default soil parameter tables
  if ( myid==0 ) then
    write(6,*) "-> Using default CABLE soil parameter tables"
  end if

  lncveg_numsoil = 9
  
  ! redefine rhos
  rhos=(/ 1600., 1600., 1381., 1373., 1476., 1521., 1373., 1537.,  910., 2600., 2600., 2600., 2600. /)

  if ( myid==0 .and. nmaxpr==1 ) then
    do isoil = 1,mxst
      write(6,"('isoil,ssat,sfc,swilt,hsbh ',i2,3f7.3,e11.4)") isoil,ssat(isoil),sfc(isoil),swilt(isoil),hsbh(isoil)
    end do
  end if
  
  allocate( isoilm_name(lncveg_numsoil) )
  isoilm_name(1) = "Coarse_sand/Loamy_sand"
  isoilm_name(2) = "Medium_clay_loam/silty_clay_loam/silt_loam"
  isoilm_name(3) = "Fine_clay"
  isoilm_name(4) = "Coarse-medium_sandy_loam/loam"
  isoilm_name(5) = "Coarse-fine_sandy_clay"
  isoilm_name(6) = "Medium-fine_silty_clay"
  isoilm_name(7) = "Coarse-medium-fine_sandy_clay_loam"
  isoilm_name(8) = "Organc_peat"  
  isoilm_name(9) = "Permanent_ice"
  do k = 1,lncveg_numsoil
    if ( myid==0 ) then
      write(6,*) " -> Define ",trim(isoilm_name(k))
    end if 
  end do  
  
else

  ! user defined soil parameter tables  
  if ( myid==0 ) then
    write(6,*) "-> Using user defined CABLE soil parameter tables"
  end if
  if ( lncveg_numsoil > mxst ) then
    write(6,*) "ERROR: Number of soil types larger than maximum, mxst,numsoil:",mxst,lncveg_numsoil
    call ccmpi_abort(-1)
  end if
  
  allocate( isoilm_name(lncveg_numsoil) )
  isoilm_name = ""

  silt = 0.
  clay = 0.
  sand = 0.
  swilt = 0.
  sfc = 0.
  ssat = 0.
  bch = 0.
  hyds = 0.
  sucs = 0.
  rhos = 0.
  css = 0.
  swilt(0) = 0.
  sfc(0) = 1.
  ssat(0) = 2.

  if ( myid==0 ) then
    chid_len = 0  
    call ccnf_inq_dimlen(ncidveg,'chid',chid_len,failok=.true.)
    do k = 1,lncveg_numsoil  
      nstart(1) = 1
      nstart(2) = k
      ncount(1) = chid_len
      ncount(2) = 1
      call ccnf_get_vara(ncidveg,'soilname',nstart,ncount,isoilm_name(k))
    end do  
    nstart(1) = 1
    ncount(1) = lncveg_numsoil    
    call ccnf_get_vara(ncidveg,'silt',nstart,ncount,silt(1:lncveg_numsoil))
    call ccnf_get_vara(ncidveg,'clay',nstart,ncount,clay(1:lncveg_numsoil))
    call ccnf_get_vara(ncidveg,'sand',nstart,ncount,sand(1:lncveg_numsoil))
    call ccnf_get_vara(ncidveg,'swilt',nstart,ncount,swilt(1:lncveg_numsoil))
    call ccnf_get_vara(ncidveg,'sfc',nstart,ncount,sfc(1:lncveg_numsoil))
    call ccnf_get_vara(ncidveg,'ssat',nstart,ncount,ssat(1:lncveg_numsoil))
    call ccnf_get_vara(ncidveg,'bch',nstart,ncount,bch(1:lncveg_numsoil))
    call ccnf_get_vara(ncidveg,'hyds',nstart,ncount,hyds(1:lncveg_numsoil))
    call ccnf_get_vara(ncidveg,'sucs',nstart,ncount,sucs(1:lncveg_numsoil))
    call ccnf_get_vara(ncidveg,'rhosoil',nstart,ncount,rhos(1:lncveg_numsoil))
    call ccnf_get_vara(ncidveg,'css',nstart,ncount,css(1:lncveg_numsoil))
  end if
  do k = 1,lncveg_numsoil
    if ( myid==0 ) then
      write(6,*) "--> Found ",trim(isoilm_name(k))
    end if
    call ccmpi_bcast(isoilm_name(k),0,comm_world)  
  end do     
  call ccmpi_bcast(silt,0,comm_world)
  call ccmpi_bcast(clay,0,comm_world)
  call ccmpi_bcast(sand,0,comm_world)
  call ccmpi_bcast(swilt,0,comm_world)
  call ccmpi_bcast(sfc,0,comm_world)
  call ccmpi_bcast(ssat,0,comm_world)
  call ccmpi_bcast(bch,0,comm_world)
  call ccmpi_bcast(hyds,0,comm_world)
  call ccmpi_bcast(sucs,0,comm_world)
  call ccmpi_bcast(rhos,0,comm_world)
  call ccmpi_bcast(css,0,comm_world)
  
  !redo from insoil
  do isoil = 1,mxst
    cnsd(isoil)  = sand(isoil)*0.3+clay(isoil)*0.25+silt(isoil)*0.265
    hsbh(isoil)  = hyds(isoil)*abs(sucs(isoil))*bch(isoil) !difsat*etasat
    ibp2(isoil)  = nint(bch(isoil))+2
    i2bp3(isoil) = 2*nint(bch(isoil))+3
    if ( myid==0 .and. nmaxpr==1 ) then
      write(6,"('isoil,ssat,sfc,swilt,hsbh ',i2,3f7.3,e11.4)") isoil,ssat(isoil),sfc(isoil),swilt(isoil),hsbh(isoil)
    end if
  end do
  cnsd(9) = 2.51

end if

do tile = 1,ntiles
  if ( tdata(tile)%mp>0 ) then

    ! Load CABLE soil data
    soil(tile)%bch       = real(bch(soil(tile)%isoilm),8)
    soil(tile)%css       = real(css(soil(tile)%isoilm),8)
    soil(tile)%rhosoil   = real(rhos(soil(tile)%isoilm),8)
    soil(tile)%cnsd      = real(cnsd(soil(tile)%isoilm),8)
    soil(tile)%hyds      = real(hyds(soil(tile)%isoilm),8)
    soil(tile)%sucs      = real(sucs(soil(tile)%isoilm),8)
    soil(tile)%hsbh      = real(hsbh(soil(tile)%isoilm),8)
    soil(tile)%sfc       = real(sfc(soil(tile)%isoilm),8)
    soil(tile)%ssat      = real(ssat(soil(tile)%isoilm),8)
    soil(tile)%swilt     = real(swilt(soil(tile)%isoilm),8)
    soil(tile)%ibp2      = real(ibp2(soil(tile)%isoilm),8)
    soil(tile)%i2bp3     = real(i2bp3(soil(tile)%isoilm),8)
    soil(tile)%pwb_min   = (soil(tile)%swilt/soil(tile)%ssat)**soil(tile)%ibp2
    soil(tile)%clay      = real(clay(soil(tile)%isoilm),8)
    soil(tile)%sand      = real(sand(soil(tile)%isoilm),8)
    soil(tile)%silt      = real(silt(soil(tile)%isoilm),8)
    soil(tile)%zeta      = 0._8
    soil(tile)%fsatmax   = 0._8
    soil(tile)%nhorizons = 1
    soil(tile)%ishorizon = 1
    do k = 1,cbm_ms
      soil(tile)%swilt_vec(:,k)   = soil(tile)%swilt
      soil(tile)%ssat_vec(:,k)    = soil(tile)%ssat
      soil(tile)%sfc_vec(:,k)     = soil(tile)%sfc
      soil(tile)%rhosoil_vec(:,k) = soil(tile)%rhosoil
      soil(tile)%sucs_vec(:,k)    = 1000._8*abs(soil(tile)%sucs)
      soil(tile)%bch_vec(:,k)     = soil(tile)%bch
      soil(tile)%hyds_vec(:,k)    = 1000._8*soil(tile)%hyds
      soil(tile)%watr(:,k)        = 0.05_8
      soil(tile)%cnsd_vec(:,k)    = soil(tile)%cnsd
      soil(tile)%clay_vec(:,k)    = soil(tile)%clay
      soil(tile)%sand_vec(:,k)    = soil(tile)%sand
      soil(tile)%silt_vec(:,k)    = soil(tile)%silt
      soil(tile)%zse_vec(:,k)     = soil(tile)%zse(k)
      soil(tile)%css_vec(:,k)     = soil(tile)%css
    end do
    soil(tile)%GWhyds_vec    = soil(tile)%hyds*1000._8  
    soil(tile)%GWsucs_vec    = abs(soil(tile)%sucs)*1000._8
    soil(tile)%GWbch_vec     = soil(tile)%bch
    soil(tile)%GWrhosoil_vec = soil(tile)%rhosoil
    soil(tile)%GWssat_vec    = soil(tile)%ssat
    soil(tile)%GWwatr        = soil(tile)%watr(:,cbm_ms) !residual water content of the aquifer [mm3/mm3]
    soil(tile)%GWdz          = 20._8           !thickness of the aquifer   [m]
    soil(tile)%GWdz          = max( 1._8, min( 20._8, soil(tile)%GWdz - sum(soil(tile)%zse,dim=1) ) )
    soil(tile)%drain_dens    = 0.008_8         !  drainage density ( mean dist to rivers/streams )
  
    soil(tile)%heat_cap_lower_limit = 0.01 ! recalculated in cable_soilsnow.F90
  
  end if
end do
    
  !if ( cable_gw_model==1 ) then
  !
  !  ! note 17 hard coded vegetation PFTs  
  !  psi_o(1:3)  = -66000._r_2
  !  psi_o(4)    = -35000._r_2
  !  psi_o(5)    = -83000._r_2
  !  psi_o(6:17) = -74000._r_2
  !  psi_c(1:3)  = -2550000._r_2
  !  psi_c(4)    = -2240000._r_2
  !  psi_c(5)    = -4280000._r_2
  !  psi_c(6:17) = -2750000._r_2      
  !    
  !  soil_depth(1) = REAL(soil%zse(1),r_2)
  !  DO k=2,cbm_ms
  !     soil_depth(k) = soil_depth(k-1) + REAL(soil%zse(k),r_2)
  !  END DO    
  !  
  !  do k=1,cbm_ms
  !    soil%hyds_vec(:,k) = 0.0070556_r_2*10.0_r_2**(-0.884_r_2 + 0.0153_r_2*soil%sand_Vec(:,k)*100.0_r_2)* &
  !      EXP(-gw_params%hkrz*(MAX(0._r_2,soil_depth(k)-gw_params%zdepth)))
  !    soil%sucs_vec(:,k) = 10.0_r_2 * 10.0_r_2**(1.88_r_2 -0.0131_r_2*soil%Sand_Vec(:,k)*100.0_r_2)
  !    soil%bch_vec(:,k) = 2.91_r_2 + 0.159_r_2*soil%Clay_Vec(:,k)*100.0_r_2
  !    soil%ssat_vec(:,k) = 0.489_r_2 - 0.00126_r_2*soil%Sand_Vec(:,k)*100.0_r_2
  !    soil%watr(:,k) = 0.02_r_2 + 0.00018_r_2*soil%Clay_Vec(:,k)*100.0_r_2
  !  end do
  !  !aquifer share non-organic with last layer if not found in param file
  !  soil%GWhyds_vec(:) = soil%hyds_vec(:,cbm_ms)
  !  soil%GWsucs_vec(:) = soil%sucs_vec(:,cbm_ms)
  !  soil%GWbch_vec(:)  = soil%bch_vec(:,cbm_ms)
  !  soil%GWssat_vec(:) = soil%ssat_vec(:,cbm_ms)
  !  soil%GWwatr(:)     = soil%watr(:,cbm_ms)
  !  !include organin impact.  fraction of grid cell where percolation through
  !  !organic macropores dominates
  !  !soil%Org_Vec = MAX(0._r_2,soil%Org_Vec)
  !  !soil%Org_Vec = MIN(1._r_2,soil%Org_Vec)
  !  !do k=1,3  !0-23.3 cm, data really is to 30cm
  !  !  soil%hyds_vec(:,k)  = (1.-soil%Org_Vec(:,k))*soil%hyds_vec(:,k) + &
  !  !    soil%Org_Vec(:,k)*gw_params%org%hyds_vec_organic
  !  !  soil%sucs_vec(:,k) = (1.-soil%Org_Vec(:,k))*soil%sucs_vec(:,k) + &
  !  !    soil%Org_Vec(:,k)*gw_params%org%sucs_vec_organic
  !  !  soil%bch_vec(:,k) = (1.-soil%Org_Vec(:,k))*soil%bch_vec(:,k) +&
  !  !    soil%Org_Vec(:,k)*gw_params%org%clappb_organic
  !  !  soil%ssat_vec(:,k) = (1.-soil%Org_Vec(:,k))*soil%ssat_vec(:,k) + &
  !  !    soil%Org_Vec(:,k)*gw_params%org%ssat_vec_organic
  !  !  soil%watr(:,k)   = (1.-soil%Org_Vec(:,klev))*soil%watr(:,klev) + &
  !  !    soil%Org_Vec(:,k)*gw_params%org%watr_organic
  !  !end do
  !
  !  !vegetation dependent field capacity (point plants get stressed) and
  !  !wilting point
  !  do k = 1,cbm_ms
  !    psi_tmp(:,k) = -psi_c(veg%iveg(:))
  ! end do
  ! soil%sfc_vec = (soil%ssat_vec-soil%watr) * (ABS(psi_tmp(:,:)) &
  !    /(ABS(soil%sucs_vec)))**(-1.0/soil%bch_vec)+soil%watr
  !  do k = 1,cbm_ms
  !    psi_tmp(:,k) = -psi_c(veg%iveg(:))
  !  end do
  !  soil%swilt_vec = (soil%ssat_vec-soil%watr) * (ABS(psi_tmp(:,:)) &
  !    /(ABS(soil%sucs_vec)))**(-1.0/soil%bch_vec)+soil%watr
  !
  !  !set the non-vectored values to srf value
  !  soil%sfc(:) = REAL(soil%sfc_vec(:,1))
  !  soil%swilt(:) = REAL(soil%swilt_vec(:,1))
  !
  !  !convert the units back to what default uses and GW only uses the
  !  !vectored versions
  !  soil%hyds = REAL(soil%hyds_vec(:,1))/1000.0
  !  soil%sucs = REAL(soil%sucs_vec(:,1))/1000.0
  !  soil%ssat = REAL(soil%ssat_vec(:,1))
  !  soil%bch  = REAL(soil%bch_vec(:,1))  
  !  
  !end if 
  
  !if ( cable_user%gw_model ) then
  !  do k = 1,cbm_ms
  !    soil%hyds_vec(:,k) = soil%hyds_vec(:,k) * exp(-gw_parms%hkrz*(znode(k)-gw_params%zdepth) )
  !  end do
  !  soil%hyds(:) = soil%hyds_vec(:,1)
  !end if
  
  !if ( cable_user%soil_thermal_fix ) then
  !  do k = 1,cbm_ms
  !    ssat_bounded = min( ssat_hi, max( ssat_lo, soil%ssat_vec(:,k) ) )  
  !    rho_soil_bulk = min( rhob_hi, max( rhob_lo, (2700.*(1.-ssat_bounded)) ) )
  !    where ( soil%isoilm(:) /= 9 )
  !     soil%rhosoil_vec(:,k) = 2700.
  !      soil%cnsd_vec(:,k) = ( 0.135*(1.-ssat_bounded(:)) + 64.7/rho_soil_bulk(:) ) &
  !                         / ( 1.-0.947*(1.-ssat_bounded(:)) )
  !    end where    
  !  end do
  !  where ( soil%isoilm(:) /= 9 )
  !    soil%rhosoil(:) = soil%rhosoil_vec(:,1)
  !    soil%cnsd(:) = soil%cnsd_vec(:,1)
  !  end where  
  !end if

return
end subroutine cable_soil_parm
                   
! *************************************************************************************
! Load CABLE biome and LAI data
! vegta is for myid==0
subroutine vegta(ivs,svs,vlin,fveg,cableformat)
  
use cc_mpi
use darcdf_m
use infile
use newmpar_m
use parmgeom_m

character(len=*), intent(in) :: fveg
integer, dimension(ifull,maxtile), intent(out) :: ivs
integer, dimension(ifull_g,maxtile) :: ivsg  
integer, dimension(3) :: spos,npos
integer n,iq,ilx,jlx,iad 
integer ncidx,iernc,varid,ndims
real, dimension(ifull,maxtile), intent(out) :: svs, vlin
real, dimension(ifull_g,maxtile) :: svsg, vling
real, dimension(ifull_g) :: savannafrac_g
real rlong0x,rlat0x,schmidtx,dsx,ra,rb
real cablever
real, intent(out) :: cableformat
character(len=47) header  
character(len=7) vname
logical tst

cableformat=0.

if ( lncveg==1 ) then
  spos(1:3) = 1
  npos(1) = il_g
  npos(2) = 6*il_g
  npos(3) = 1
  call ccnf_inq_dimlen(ncidveg,'longitude',ilx)
  call ccnf_inq_dimlen(ncidveg,'latitude',jlx)
  call ccnf_get_attg(ncidveg,'lon0',rlong0x)
  call ccnf_get_attg(ncidveg,'lat0',rlat0x)
  call ccnf_get_attg(ncidveg,'schmidt',schmidtx)
  if (ilx/=il_g.or.jlx/=jl_g.or.abs(rlong0x-rlong0)>1.e-20.or.abs(rlat0x-rlat0)>1.e-20.or. &
      abs(schmidtx-schmidt)>1.e-20) then
    write(6,*) 'wrong data file supplied ',trim(fveg)
    call ccmpi_abort(-1)
  end if
  !call ccnf_get_attg(ncidveg,'cableversion',cablever,ierr=iernc)
  !if (iernc/=0) then
  !  write(6,*) "Missing version of CABLE data"
  !  write(6,*) "Regenerate land-use data with up-to-date version of igbpveg"
  !  call ccmpi_abort(-1)
  !end if
  !if (abs(cablever-cable_version)>1.e-20 .and. abs(cablever-6608.)>1.e-20 .and. &
  !    abs(cablever-3939.)>1.e-20 ) then
  !  write(6,*) "Wrong version of CABLE data"
  !  write(6,*) "Expecting 3939. or 6608. or ",cable_version
  !  write(6,*) "Found     ",cablever
  !  write(6,*) "Please upgrade igbpveg to fix this error"
  !  call ccmpi_abort(-1)
  !end if
  call ccnf_get_attg(ncidveg,'cableformat',cableformat,ierr=iernc)
  if ( iernc/=0 ) then
    cableformat=0.
  end if
  do n = 1,maxtile
    vling(:,n) = 0.
    svsg(:,n) = 0.
    ivsg(:,n) = 0
    write(vname,"(A,I1.1)") "lai",n
    call ccnf_inq_varid(ncidveg,vname,varid,tst)
    if ( .not.tst ) then
      call ccnf_inq_varndims(ncidveg,varid,ndims)
      call ccnf_get_vara(ncidveg,varid,spos(1:ndims),npos(1:ndims),vling(:,n)) 
      write(vname,"(A,I1.1)") "vegt",n
      call ccnf_inq_varid(ncidveg,vname,varid)
      call ccnf_inq_varndims(ncidveg,varid,ndims)
      call ccnf_get_vara(ncidveg,varid,spos(1:ndims),npos(1:ndims),svsg(:,n)) 
      ivsg(:,n)=nint(svsg(:,n))
      write(vname,"(A,I1.1)") "vfrac",n
      call ccnf_inq_varid(ncidveg,vname,varid)
      call ccnf_inq_varndims(ncidveg,varid,ndims)
      call ccnf_get_vara(ncidveg,varid,spos(1:ndims),npos(1:ndims),svsg(:,n))
    end if  
  end do
  vname="savanna"
  call ccnf_inq_varid(ncidveg,vname,varid,tst)
  if ( .not.tst ) then
    call ccnf_inq_varndims(ncidveg,varid,ndims)
    call ccnf_get_vara(ncidveg,varid,spos(1:ndims),npos(1:ndims),savannafrac_g)
    do n = 1,maxtile
      where ( savannafrac_g>0.5*svsg(:,n) .and. ivsg(:,n)==2 )
        ivsg(:,n) = 18
      end where
    end do
  end if
  call ccmpi_distribute(ivs,ivsg)
  call ccmpi_distribute(svs,svsg)
  call ccmpi_distribute(vlin,vling)
else
  open(87,file=fveg,status='old')
  read(87,*) ilx,jlx,rlong0x,rlat0x,schmidtx,dsx,header
  if (ilx/=il_g.or.jlx/=jl_g.or.abs(rlong0x-rlong0)>1.e-20.or.abs(rlat0x-rlat0)>1.e-20.or. &
      abs(schmidtx-schmidt)>1.e-20) then
    write(6,*) 'wrong data file supplied ',trim(fveg)
    call ccmpi_abort(-1)
  end if
  ivsg = 0
  svsg = 0.
  vling = 0.
  do iq = 1,ifull_g
    read(87,*) iad,ra,rb,ivsg(iq,1),svsg(iq,1),vling(iq,1),ivsg(iq,2),svsg(iq,2),vling(iq,2),ivsg(iq,3),svsg(iq,3),vling(iq,3), &
               ivsg(iq,4),svsg(iq,4),vling(iq,4),ivsg(iq,5),svsg(iq,5),vling(iq,5)
  end do
  close(87)
  call ccmpi_distribute(ivs,ivsg)
  call ccmpi_distribute(svs,svsg)
  call ccmpi_distribute(vlin,vling)
end if

call ccmpi_bcast(cableformat,0,comm_world)

return
end subroutine vegta
  
! vegtb is for myid != 0
subroutine vegtb(ivs,svs,vlin,cableformat)
  
use cc_mpi
use newmpar_m
  
integer, dimension(ifull,maxtile), intent(out) :: ivs
real, dimension(ifull,maxtile), intent(out) :: svs, vlin
real, intent(out) :: cableformat

cableformat = 0.

call ccmpi_distribute(ivs)
call ccmpi_distribute(svs)
call ccmpi_distribute(vlin)
  
call ccmpi_bcast(cableformat,0,comm_world)

return
end subroutine vegtb

! *************************************************************************************
! Transfer grid information from CABLE internally, read N&P input from
! integral NETCDF file
subroutine casa_readpoint(casafile,casapoint)

use cc_mpi
use infile
use newmpar_m
use parmgeom_m

integer ncstatus, ncid, varid, tilg
integer, dimension(2) :: spos, npos
real tlat, tlon, tschmidt
real, dimension(:,:), allocatable, save :: dumg
real, dimension(ifull,5), intent(out) :: casapoint
character(len=*), intent(in) :: casafile
logical tst

if ( myid==0 ) then
  allocate( dumg(ifull_g,5) )
  dumg = 0.
  if ( casafile==" " ) then
    write(6,*) "ERROR: casafile is not specified"
    call ccmpi_abort(-1)
  end if    
  write(6,*) "Reading ",trim(casafile)
  call ccnf_open(casafile,ncid,ncstatus)
  call ncmsg('CASA_readpoint',ncstatus)
  ! check dimensions and location
  call ccnf_get_attg(ncid,'lat0',tlat)
  call ccnf_get_attg(ncid,'lon0',tlon)
  call ccnf_get_attg(ncid,'schmidt0',tschmidt)
  if ( abs(rlong0-tlon)>1.e-20 .or. abs(rlat0-tlat)>1.e-20 .or. abs(schmidt-tschmidt)>1.e-20 ) then
    write(6,*) "ERROR: Grid mismatch for ",trim(casafile)
    write(6,*) "rlong0,rlat0,schmidt ",rlong0,rlat0,schmidt
    write(6,*) "tlon,tlat,tschmidt   ",tlon,tlat,tschmidt
    call ccmpi_abort(-1)
  end if
  call ccnf_inq_dimlen(ncid,'longitude',tilg)
  if ( tilg /= il_g ) then
    write (6,*) "ERROR: Grid mismatch for ",trim(casafile)
    write (6,*) "il_g,tilg ",il_g,tilg
    call ccmpi_abort(-1)
  end if
  ! load casa fields
  spos(1:2) = 1
  npos(1) = il_g
  npos(2) = il_g*6
  write(6,*) "Loading soil order"
  call ccnf_inq_varid(ncid,'sorder',varid,tst)
  call ccnf_get_vara(ncid,varid,spos,npos,dumg(:,1))
  write(6,*) "Loading N deposition rate"
  call ccnf_inq_varid(ncid,'ndep',varid,tst)
  call ccnf_get_vara(ncid,varid,spos,npos,dumg(:,2))
  write(6,*) "Loading N fixation rate"
  call ccnf_inq_varid(ncid,'nfix',varid,tst)
  call ccnf_get_vara(ncid,varid,spos,npos,dumg(:,3))
  write(6,*) "Loading P dust deposition"
  call ccnf_inq_varid(ncid,'pdust',varid,tst)
  call ccnf_get_vara(ncid,varid,spos,npos,dumg(:,4))
  write(6,*) "Loading P weathering rate"
  call ccnf_inq_varid(ncid,'pweather',varid,tst)
  call ccnf_get_vara(ncid,varid,spos,npos,dumg(:,5))
  call ccnf_close(ncid)
  call ccmpi_distribute(casapoint,dumg)
  deallocate(dumg)
else
  call ccmpi_distribute(casapoint)
end if

end subroutine casa_readpoint

! *************************************************************************************  
! This subroutine reads the MODIS derived leaf phenology data
subroutine casa_readphen(fphen,greenup,fall,phendoy1)

use cc_mpi
use infile
use newmpar_m
use parm_m

integer, parameter :: nphen = 8
integer ilat, ierr
integer ncid
integer, dimension(271,mxvt), intent(out) :: greenup, fall, phendoy1
integer, dimension(nphen) :: greenupx, fallx, xphendoy1
integer, dimension(nphen) :: ivtx
real :: xlat
logical :: ncfile
character(len=*), intent(in) :: fphen

! initilize for evergreen PFTs
greenup(:,:)  = -50
fall(:,:)     = 367
phendoy1(:,:) = 2

! MJT notes - if the CSIRO PFTs are redefined, then this phenology data will be mismatched

if ( myid==0 ) then
  write(6,*) "-> Reading CASA leaf phenology data ",trim(fphen)
  call ccnf_open(fphen,ncid,ierr)
  if ( ierr==0 ) then
    ncfile = .true.
    if ( nmaxpr==1 ) write(6,*) "-> Found netcdf file ",trim(fphen)
  else
    call ccnf_open(fphen//'.nc',ncid,ierr)
    if ( ierr==0 ) then
      ncfile = .true.
      if ( nmaxpr==1 ) write(6,*) "-> Found netcdf file ",trim(fphen)
    else  
      open(87,file=fphen,status='old',iostat=ierr)
      if ( ierr/=0 ) then
        write(6,*) "ERROR: Cannot open phenfile=",trim(fphen)
        call ccmpi_abort(-1)
      end if
      read(87,*)
      read(87,*) ivtx
    end if
  end if
  
  if ( ncfile ) then
    call ccnf_get_vara(ncid,"greenup",(/1,1/),(/271,mxvt/),greenup)
    call ccnf_get_vara(ncid,"fall",(/1,1/),(/271,mxvt/),fall)
    call ccnf_get_vara(ncid,"phendoy1",(/1,1/),(/271,mxvt/),phendoy1)
    call ccnf_close(ncid)
  else    
    do ilat = 271,1,-1
      read(87,*) xlat,greenupx,fallx,xphendoy1 
      greenup(ilat,ivtx(:))  = greenupx(:)
      fall(ilat,ivtx(:))     = fallx(:)
      phendoy1(ilat,ivtx(:)) = xphendoy1(:)
    end do
    close(87)
  end if  
end if
call ccmpi_bcast(greenup, 0,comm_world)
call ccmpi_bcast(fall,    0,comm_world)
call ccmpi_bcast(phendoy1,0,comm_world)

return
end subroutine casa_readphen

end module cable_ccam_setup

