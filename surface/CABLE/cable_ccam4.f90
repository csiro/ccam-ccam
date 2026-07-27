! Conformal Cubic Atmospheric Model
    
! Copyright 2015-2026 Commonwealth Scientific Industrial Research Organisation (CSIRO)
    
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
    
! This model provides routines for packing and unpacking cable land-surface
! and carbon cycle data
    
module cable_ccam4
 
use cable_def_types_mod, only : air_type, balances_type, canopy_type, climate_type, met_type, radiation_type, &
                                roughness_type, soil_parameter_type, soil_snow_type, sum_flux_type,           &
                                veg_parameter_type, mp, r_2
use casavariable, only : casa_balance, casa_biome, casa_flux, casa_met, casa_pool
use newmpar_m, only : ntiles, imax, ifull
use phenvariable, only : phen_variable
use pop_types, only : pop_type, dp

implicit none

private
public tdata
public cable_pack, cable_unpack, pop_pack, pop_unpack
public cable_update
    
type tiledata
  integer, dimension(:,:), allocatable :: tind
  integer, dimension(:,:), allocatable :: pind
  logical, dimension(:,:), allocatable :: tmap
  logical, dimension(:,:), allocatable :: pmap
  integer :: mp
  integer :: np
  integer :: maxnb
  integer :: toffset
  integer, dimension(:), allocatable :: cveg ! CABLE vegetation index for each point
  integer, dimension(:), allocatable :: old_cveg ! for redistribution
  real, dimension(:), allocatable :: sv      ! area fraction for each point
  real, dimension(:), allocatable :: vl2     ! prescribed leaf area index (LAI) for point
  real, dimension(:), allocatable :: old_sv      ! for redustribution
end type tiledata

type(tiledata), dimension(:), allocatable :: tdata

interface cable_pack
  module procedure cable_pack_r4_2_r4, cable_pack_r4_2_r8
  module procedure cable_pack_i4_2_i4, cable_pack_r4_2_i4
  module procedure cable_pack_i4_2_i8, cable_pack_r4_2_i8
  module procedure cable_pack_r4_2_r4_map, cable_pack_r4_2_r8_map
  module procedure cable_pack_i4_2_i4_map, cable_pack_i4_2_i8_map
#ifndef i8r8
  module procedure cable_pack_r8_2_r8
  module procedure cable_pack_r8_2_r8_map
#endif
end interface

interface cable_unpack
  module procedure cable_unpack_r4_2_r4, cable_unpack_r8_2_r4
  module procedure cable_unpack_i4_2_i4
#ifndef i8r8
  module procedure cable_unpack_r4_2_r8
  module procedure cable_unpack_r8_2_r8
#endif
end interface

interface pop_pack
  module procedure pop_pack_r8_2_r8, pop_pack_r8_2_r8_map
  module procedure pop_pack_i4_2_i4, pop_pack_i4_2_i4_map
  module procedure pop_pack_i4_2_i8, pop_pack_i4_2_i8_map
end interface

interface pop_unpack
  module procedure pop_unpack_r8_2_r8, pop_unpack_i4_2_r8
  module procedure pop_unpack_i8_2_r8
end interface

interface cable_update
  module procedure cable_update_r8_2_r4
#ifndef i8r8
  module procedure cable_update_r8_2_r8
#endif
end interface

contains

subroutine cable_pack_r4_2_r4(indata,outdata,tile,nb)
  real, dimension(:), intent(in) :: indata
  real(kind=4), dimension(:), intent(inout) :: outdata
  integer, intent(in) :: tile
  integer, intent(in), optional :: nb
  integer js, je, is, ie, lnb
  
  if ( size(indata)==ifull ) then
    js = (tile-1)*imax + 1
    je = tile*imax
  else if ( size(indata)==imax ) then
    js = 1
    je = imax
  else
    write(6,*) "Invalid input lengh for cable_pack (tile)"
    stop
  end if  
  
  if ( present(nb) ) then
      
    is = tdata(tile)%tind(nb,1)
    ie = tdata(tile)%tind(nb,2)
    if ( is<=ie ) then
      outdata(is:ie) = pack(real(indata(js:je),4),tdata(tile)%tmap(:,nb))
    end if
          
  else

    do lnb = 1,tdata(tile)%maxnb
      is = tdata(tile)%tind(lnb,1)
      ie = tdata(tile)%tind(lnb,2)
      if ( is<=ie ) then
        outdata(is:ie) = pack(real(indata(js:je),4),tdata(tile)%tmap(:,lnb))
      end if
    end do  
  
  end if
 
end subroutine cable_pack_r4_2_r4

subroutine cable_pack_r4_2_r8(indata,outdata,tile,nb)
  real, dimension(:), intent(in) :: indata
  real(kind=8), dimension(:), intent(inout) :: outdata
  integer, intent(in) :: tile
  integer, intent(in), optional :: nb
  integer js, je, is, ie, lnb

  if ( size(indata)==ifull ) then
    js = (tile-1)*imax + 1
    je = tile*imax
  else if ( size(indata)==imax ) then
    js = 1
    je = imax
  else
    write(6,*) "Invalid input lengh for cable_pack (tile)"
    stop
  end if
  
  if ( present(nb) ) then

    is = tdata(tile)%tind(nb,1)
    ie = tdata(tile)%tind(nb,2)
    if ( is<=ie ) then
      outdata(is:ie) = pack(real(indata(js:je),8),tdata(tile)%tmap(:,nb))
    end if      
      
  else
  
    do lnb = 1,tdata(tile)%maxnb
      is = tdata(tile)%tind(lnb,1)
      ie = tdata(tile)%tind(lnb,2)
      if ( is<=ie ) then
        outdata(is:ie) = pack(real(indata(js:je),8),tdata(tile)%tmap(:,lnb))
      end if
    end do

end if

end subroutine cable_pack_r4_2_r8

subroutine cable_pack_r4_2_r4_map(indata,outdata,tile,nmp)
  real, dimension(:), intent(in) :: indata
  real(kind=4), dimension(:), intent(inout) :: outdata
  integer, intent(in) :: tile
  integer, dimension(imax), intent(in) :: nmp
  integer nb, is, ie, iqt, iq, js, je, iqx

  if ( size(indata)==ifull ) then
    js = (tile-1)*imax + 1
    je = tile*imax
  else if ( size(indata)==imax ) then
    js = 1
    je = imax
  else
    write(6,*) "Invalid input lengh for cable_pack (map)"
    stop
  end if
  
  if ( all( nmp(js)==nmp(js:je) ) ) then
    nb = nmp(js)
    is = tdata(tile)%tind(nb,1)
    ie = tdata(tile)%tind(nb,2)
    if ( is<=ie ) then
      outdata(is:ie) = pack(real(indata(js:je),4),tdata(tile)%tmap(:,nb))
    end if  
  else    
    do iq = js,je
      nb = nmp(iq)
      iqx = iq-js+1
      if ( tdata(tile)%tmap(iqx,nb) ) then
        is = tdata(tile)%tind(nb,1)  
        iqt = is + count( tdata(tile)%tmap(1:iqx,nb) ) - 1
        outdata(iqt) = real(indata(iq),4)
      end if
    end do  
  end if  

end subroutine cable_pack_r4_2_r4_map

subroutine cable_pack_r4_2_r8_map(indata,outdata,tile,nmp)
  real, dimension(:), intent(in) :: indata
  real(kind=8), dimension(:), intent(inout) :: outdata
  integer, intent(in) :: tile
  integer, dimension(imax), intent(in) :: nmp
  integer nb, is, ie, iqt, iq, js, je, iqx

  if ( size(indata)==ifull ) then
    js = (tile-1)*imax + 1
    je = tile*imax
  else if ( size(indata)==imax ) then
    js = 1
    je = imax
  else
    write(6,*) "Invalid input lengh for cable_pack (map)"
    stop
  end if

  if ( all( nmp(js)==nmp(js:je) ) ) then
    nb = nmp(js)
    is = tdata(tile)%tind(nb,1)
    ie = tdata(tile)%tind(nb,2)
    if ( is<=ie ) then
      outdata(is:ie) = pack(real(indata(js:je),8),tdata(tile)%tmap(:,nb))
    end if  
  else    
    do iq = js,je
      nb = nmp(iq)
      iqx = iq-js+1
      if ( tdata(tile)%tmap(iqx,nb) ) then
        is = tdata(tile)%tind(nb,1)  
        iqt = is + count( tdata(tile)%tmap(1:iqx,nb) ) - 1
        outdata(iqt) = real(indata(iq),8)
      end if
    end do   
  end if  
  
end subroutine cable_pack_r4_2_r8_map

#ifndef i8r8
subroutine cable_pack_r8_2_r8(indata,outdata,tile,nb)
  real(kind=8), dimension(:), intent(in) :: indata
  real(kind=8), dimension(:), intent(inout) :: outdata
  integer, intent(in) :: tile
  integer, intent(in), optional :: nb
  integer js, je, is, ie, lnb

  if ( size(indata)==ifull ) then
    js = (tile-1)*imax + 1
    je = tile*imax
  else if ( size(indata)==imax ) then
    js = 1
    je = imax
  else
    write(6,*) "Invalid input lengh for cable_pack (tile)"
    stop
  end if
  
  if ( present(nb) ) then
  
    is = tdata(tile)%tind(nb,1)
    ie = tdata(tile)%tind(nb,2)
    if ( is<=ie ) then
      outdata(is:ie) = pack(indata(js:je),tdata(tile)%tmap(:,nb))
    end if
    
  else

     do lnb = 1,tdata(tile)%maxnb
      is = tdata(tile)%tind(lnb,1)
      ie = tdata(tile)%tind(lnb,2)
      if ( is<=ie ) then
        outdata(is:ie) = pack(indata(js:je),tdata(tile)%tmap(:,lnb))
      end if
    end do

  end if
  
end subroutine cable_pack_r8_2_r8

subroutine cable_pack_r8_2_r8_map(indata,outdata,tile,nmp)
  real(kind=8), dimension(:), intent(in) :: indata
  real(kind=8), dimension(:), intent(inout) :: outdata
  integer, intent(in) :: tile
  integer, dimension(imax), intent(in) :: nmp
  integer nb, is, ie, iqt, iq, js, je, iqx

  if ( size(indata)==ifull ) then
    js = (tile-1)*imax + 1
    je = tile*imax
  else if ( size(indata)==imax ) then
    js = 1
    je = imax
  else
    write(6,*) "Invalid input lengh for cable_pack (map)"
    stop
  end if

  if ( all( nmp(js)==nmp(js:je) ) ) then
    nb = nmp(js)
    is = tdata(tile)%tind(nb,1)
    ie = tdata(tile)%tind(nb,2)
    if ( is<=ie ) then
      outdata(is:ie) = pack(indata(js:je),tdata(tile)%tmap(:,nb))
    end if  
  else    
    do iq = js,je
      nb = nmp(iq)
      iqx = iq-js+1
      if ( tdata(tile)%tmap(iqx,nb) ) then
        is = tdata(tile)%tind(nb,1)  
        iqt = is + count( tdata(tile)%tmap(1:iqx,nb) ) - 1
        outdata(iqt) = indata(iq)
      end if
    end do  
  end if  

end subroutine cable_pack_r8_2_r8_map
#endif

subroutine cable_pack_i4_2_i4(indata,outdata,tile,nb)
  integer, dimension(:), intent(in) :: indata
  integer(kind=4), dimension(:), intent(inout) :: outdata
  integer, intent(in) :: tile
  integer, intent(in), optional :: nb
  integer js, je, is, ie, lnb
  
  if ( size(indata)==ifull ) then
    js = (tile-1)*imax + 1
    je = tile*imax
  else if ( size(indata)==imax ) then
    js = 1
    je = imax
  else
    write(6,*) "Invalid input lengh for cable_pack (tile)"
    stop
  end if  
  
  if ( present(nb) ) then
  
    is = tdata(tile)%tind(nb,1)
    ie = tdata(tile)%tind(nb,2)
    if ( is<=ie ) then
      outdata(is:ie) =  pack(int(indata(js:je),4),tdata(tile)%tmap(:,nb))
    end if
  
  else

    do lnb = 1,tdata(tile)%maxnb
      is = tdata(tile)%tind(lnb,1)
      ie = tdata(tile)%tind(lnb,2)
      if ( is<=ie ) then
        outdata(is:ie) =  pack(int(indata(js:je),4),tdata(tile)%tmap(:,lnb))
      end if
    end do
      
  end if

end subroutine cable_pack_i4_2_i4

subroutine cable_pack_i4_2_i4_map(indata,outdata,tile,nmp)
  integer, dimension(:), intent(in) :: indata
  integer(kind=4), dimension(:), intent(inout) :: outdata
  integer, intent(in) :: tile
  integer, dimension(imax), intent(in) :: nmp
  integer nb, is, ie, iqt, iq, js, je, iqx

  if ( size(indata)==ifull ) then
    js = (tile-1)*imax + 1
    je = tile*imax
  else if ( size(indata)==imax ) then
    js = 1
    je = imax
  else
    write(6,*) "Invalid input lengh for cable_pack (map)"
    stop
  end if

  if ( all( nmp(js)==nmp(js:je) ) ) then
    nb = nmp(js)
    is = tdata(tile)%tind(nb,1)
    ie = tdata(tile)%tind(nb,2)
    if ( is<=ie ) then
      outdata(is:ie) = pack(int(indata(js:je),4),tdata(tile)%tmap(:,nb))
    end if
  else
    do iq = js,je
      nb = nmp(iq)
      iqx = iq-js+1
      if ( tdata(tile)%tmap(iqx,nb) ) then
        is = tdata(tile)%tind(nb,1)  
        iqt = is + count( tdata(tile)%tmap(1:iqx,nb) ) - 1
        outdata(iqt) = int(indata(iq),4)
      end if
    end do  
  end if

end subroutine cable_pack_i4_2_i4_map

subroutine cable_pack_i4_2_i8_map(indata,outdata,tile,nmp)
  integer, dimension(:), intent(in) :: indata
  integer(kind=8), dimension(:), intent(inout) :: outdata
  integer, intent(in) :: tile
  integer, dimension(imax), intent(in) :: nmp
  integer nb, is, ie, iqt, iq, js, je, iqx

  if ( size(indata)==ifull ) then
    js = (tile-1)*imax + 1
    je = tile*imax
  else if ( size(indata)==imax ) then
    js = 1
    je = imax
  else
    write(6,*) "Invalid input lengh for cable_pack (map)"
    stop
  end if

  if ( all( nmp(js)==nmp(js:je) ) ) then
    nb = nmp(js)
    is = tdata(tile)%tind(nb,1)
    ie = tdata(tile)%tind(nb,2)
    if ( is<=ie ) then
      outdata(is:ie) = pack(int(indata(js:je),8),tdata(tile)%tmap(:,nb))
    end if  
  else    
    do iq = js,je
      nb = nmp(iq)
      iqx = iq-js+1
      if ( tdata(tile)%tmap(iqx,nb) ) then
        is = tdata(tile)%tind(nb,1)  
        iqt = is + count( tdata(tile)%tmap(1:iqx,nb) ) - 1
        outdata(iqt) = int(indata(iq),8)
      end if
    end do  
  end if

end subroutine cable_pack_i4_2_i8_map

subroutine cable_pack_i4_2_i8(indata,outdata,tile,nb)
  integer, dimension(:), intent(in) :: indata
  integer(kind=8), dimension(:), intent(inout) :: outdata
  integer, intent(in) :: tile
  integer, intent(in), optional :: nb
  integer js, je, is, ie, lnb

  if ( size(indata)==ifull ) then
    js = (tile-1)*imax + 1
    je = tile*imax
  else if ( size(indata)==imax ) then
    js = 1
    je = imax
  else
    write(6,*) "Invalid input lengh for cable_pack (tile)"
    stop
  end if  
  
  if ( present(nb) ) then
  
    is = tdata(tile)%tind(nb,1)
    ie = tdata(tile)%tind(nb,2)
    if ( is<=ie ) then
      outdata(is:ie) =  pack(int(indata(js:je),8),tdata(tile)%tmap(:,nb))
    end if
    
  else

    do lnb = 1,tdata(tile)%maxnb
      is = tdata(tile)%tind(lnb,1)
      ie = tdata(tile)%tind(lnb,2)
      if ( is<=ie ) then
        outdata(is:ie) =  pack(int(indata(js:je),8),tdata(tile)%tmap(:,lnb))
      end if
    end do  
      
  end if

end subroutine cable_pack_i4_2_i8

subroutine cable_pack_r4_2_i4(indata,outdata,tile,nb)
  real, dimension(:), intent(in) :: indata
  integer(kind=4), dimension(:), intent(inout) :: outdata
  integer, intent(in) :: tile
  integer, intent(in), optional :: nb
  integer js, je, is, ie, lnb

  if ( size(indata)==ifull ) then
    js = (tile-1)*imax + 1
    je = tile*imax
  else if ( size(indata)==imax ) then
    js = 1
    je = imax
  else
    write(6,*) "Invalid input lengh for cable_pack (tile)"
    stop
  end if  
  
  if ( present(nb) ) then
  
    is = tdata(tile)%tind(nb,1)
    ie = tdata(tile)%tind(nb,2)
    if ( is<=ie ) then
      outdata(is:ie) =  pack(nint(indata(js:je),4),tdata(tile)%tmap(:,nb))
    end if
  
  else
      
    do lnb = 1,tdata(tile)%maxnb  
      is = tdata(tile)%tind(lnb,1)
      ie = tdata(tile)%tind(lnb,2)
      if ( is<=ie ) then
        outdata(is:ie) =  pack(nint(indata(js:je),4),tdata(tile)%tmap(:,lnb))
      end if
    end do  
      
  end if

end subroutine cable_pack_r4_2_i4

subroutine cable_pack_r4_2_i8(indata,outdata,tile,nb)
  real, dimension(:), intent(in) :: indata
  integer(kind=8), dimension(:), intent(inout) :: outdata
  integer, intent(in) :: tile
  integer, intent(in), optional :: nb
  integer js, je, is, ie, lnb
  
  if ( size(indata)==ifull ) then
    js = (tile-1)*imax + 1
    je = tile*imax
  else if ( size(indata)==imax ) then
    js = 1
    je = imax
  else
    write(6,*) "Invalid input lengh for cable_pack (tile)"
    stop
  end if
  
  if ( present(nb) ) then

    is = tdata(tile)%tind(nb,1)
    ie = tdata(tile)%tind(nb,2)
    if ( is<=ie ) then
      outdata(is:ie) =  pack(nint(indata(js:je),8),tdata(tile)%tmap(:,nb))
    end if  
  
  else

    do lnb = 1,tdata(tile)%maxnb
      is = tdata(tile)%tind(lnb,1)
      ie = tdata(tile)%tind(lnb,2)
      if ( is<=ie ) then
        outdata(is:ie) =  pack(nint(indata(js:je),8),tdata(tile)%tmap(:,lnb))
      end if
    end do  

  end if

end subroutine cable_pack_r4_2_i8

subroutine cable_unpack_r4_2_r4(indata,outdata,tile,nb)
  real(kind=4), dimension(:), intent(in) :: indata
  real, dimension(:), intent(inout) :: outdata
  integer, intent(in) :: tile
  integer, intent(in), optional :: nb
  integer js, je, is, ie, lnb
  
  if ( size(outdata)==ifull ) then
    js = (tile-1)*imax + 1
    je = tile*imax
  else if ( size(outdata)==imax ) then
    js = 1
    je = imax
  else
    write(6,*) "Invalid input lengh for cable_unpack (tile)"
    stop
  end if

  if ( present(nb) ) then
      
    is = tdata(tile)%tind(nb,1)
    ie = tdata(tile)%tind(nb,2)
    if ( is<=ie ) then
      outdata(js:je) = unpack(real(indata(is:ie)),tdata(tile)%tmap(:,nb),outdata(js:je))
    end if  
    
  else

    do lnb = 1,tdata(tile)%maxnb
      is = tdata(tile)%tind(lnb,1)
      ie = tdata(tile)%tind(lnb,2)
      if ( is<=ie ) then
        outdata(js:je) = unpack(real(indata(is:ie)),tdata(tile)%tmap(:,lnb),outdata(js:je))
      end if  
    end do

  end if

end subroutine cable_unpack_r4_2_r4

subroutine cable_unpack_r8_2_r4(indata,outdata,tile,nb)
  real(kind=8), dimension(:), intent(in) :: indata
  real, dimension(:), intent(inout) :: outdata
  integer, intent(in) :: tile
  integer, intent(in), optional :: nb
  integer js, je, is, ie, lnb
  
  if ( size(outdata)==ifull ) then
    js = (tile-1)*imax + 1
    je = tile*imax
  else if ( size(outdata)==imax ) then
    js = 1
    je = imax
  else
    write(6,*) "Invalid input lengh for cable_unpack (tile)"
    stop
  end if
  
  if ( present(nb) ) then

    is = tdata(tile)%tind(nb,1)
    ie = tdata(tile)%tind(nb,2)
    if ( is<=ie ) then
      outdata(js:je) = unpack(real(indata(is:ie)),tdata(tile)%tmap(:,nb),outdata(js:je))
    end if

  else
      
    do lnb = 1,tdata(tile)%maxnb 
      is = tdata(tile)%tind(lnb,1)
      ie = tdata(tile)%tind(lnb,2)
      if ( is<=ie ) then
        outdata(js:je) = unpack(real(indata(is:ie)),tdata(tile)%tmap(:,lnb),outdata(js:je))
      end if
    end do  
      
  end if

end subroutine cable_unpack_r8_2_r4

subroutine cable_unpack_i4_2_i4(indata,outdata,tile,nb)
  integer, dimension(:), intent(in) :: indata
  integer, dimension(:), intent(inout) :: outdata
  integer, intent(in) :: tile
  integer, intent(in), optional :: nb
  integer js, je, is, ie, lnb
  
  if ( size(outdata)==ifull ) then
    js = (tile-1)*imax + 1
    je = tile*imax
  else if ( size(outdata)==imax ) then
    js = 1
    je = imax
  else
    write(6,*) "Invalid input lengh for cable_unpack (tile)"
    stop
  end if
  
  if ( present(nb) ) then

    is = tdata(tile)%tind(nb,1)
    ie = tdata(tile)%tind(nb,2)
    if ( is<=ie ) then
      outdata(js:je) = unpack(indata(is:ie),tdata(tile)%tmap(:,nb),outdata(js:je))
    end if  
  
  else

    do lnb = 1,tdata(tile)%maxnb
      is = tdata(tile)%tind(lnb,1)
      ie = tdata(tile)%tind(lnb,2)
      if ( is<=ie ) then
        outdata(js:je) = unpack(indata(is:ie),tdata(tile)%tmap(:,lnb),outdata(js:je))
      end if  
    end do
      
  end if

end subroutine cable_unpack_i4_2_i4

#ifndef i8r8
subroutine cable_unpack_r4_2_r8(indata,outdata,tile,nb)
  real(kind=4), dimension(:), intent(in) :: indata
  real(kind=8), dimension(:), intent(inout) :: outdata
  integer, intent(in) :: tile
  integer, intent(in), optional :: nb
  integer js, je, is, ie, lnb
  
  if ( size(outdata)==ifull ) then
    js = (tile-1)*imax + 1
    je = tile*imax
  else if ( size(outdata)==imax ) then
    js = 1
    je = imax
  else
    write(6,*) "Invalid input lengh for cable_unpack (tile)"
    stop
  end if

  if ( present(nb) ) then
  
    is = tdata(tile)%tind(nb,1)
    ie = tdata(tile)%tind(nb,2)
    if ( is<=ie ) then
      outdata(js:je) = unpack(real(indata(is:ie),8),tdata(tile)%tmap(:,nb),outdata(js:je))
    end if  

  else

    do lnb = 1,tdata(tile)%maxnb 
      is = tdata(tile)%tind(lnb,1)
      ie = tdata(tile)%tind(lnb,2)
      if ( is<=ie ) then
        outdata(js:je) = unpack(real(indata(is:ie),8),tdata(tile)%tmap(:,lnb),outdata(js:je))
      end if
    end do  
      
  end if
      
end subroutine cable_unpack_r4_2_r8

subroutine cable_unpack_r8_2_r8(indata,outdata,tile,nb)
  real(kind=8), dimension(:), intent(in) :: indata
  real(kind=8), dimension(:), intent(inout) :: outdata
  integer, intent(in) :: tile
  integer, intent(in), optional :: nb
  integer js, je, is, ie, lnb

  if ( size(outdata)==ifull ) then
    js = (tile-1)*imax + 1
    je = tile*imax
  else if ( size(outdata)==imax ) then
    js = 1
    je = imax
  else
    write(6,*) "Invalid input lengh for cable_unpack (tile)"
    stop
  end if
  
  if ( present(nb) ) then
  
    is = tdata(tile)%tind(nb,1)
    ie = tdata(tile)%tind(nb,2)
    if ( is<=ie ) then
      outdata(js:je) = unpack(indata(is:ie),tdata(tile)%tmap(:,nb),outdata(js:je))
    end if  
    
  else

    do lnb = 1,tdata(tile)%maxnb
      is = tdata(tile)%tind(lnb,1)
      ie = tdata(tile)%tind(lnb,2)
      if ( is<=ie ) then
        outdata(js:je) = unpack(indata(is:ie),tdata(tile)%tmap(:,lnb),outdata(js:je))
      end if
    end do  
      
  end if

end subroutine cable_unpack_r8_2_r8
#endif

subroutine pop_pack_r8_2_r8(indata,outdata,tile,nb)
  real(kind=8), dimension(:), intent(in) :: indata
  real(kind=8), dimension(:), intent(inout) :: outdata
  integer, intent(in) :: tile
  integer, intent(in), optional :: nb
  integer is, ie, js, je, lnb

  if ( size(indata)==ifull ) then
    js = (tile-1)*imax + 1
    je = tile*imax
  else if ( size(indata)==imax ) then
    js = 1
    je = imax
  else
    write(6,*) "Invalid input lengh for pop_pack (tile)"
    stop
  end if

  if ( present(nb) ) then

    is = tdata(tile)%pind(nb,1)
    ie = tdata(tile)%pind(nb,2)
    if ( is<=ie ) then
      outdata(is:ie) = pack(indata(js:je),tdata(tile)%pmap(:,nb))
    end if  

  else
      
    do lnb = 1,tdata(tile)%maxnb
      is = tdata(tile)%pind(lnb,1)
      ie = tdata(tile)%pind(lnb,2)
      if ( is<=ie ) then
        outdata(is:ie) = pack(indata(js:je),tdata(tile)%pmap(:,lnb))
      end if  
    end do

  end if
      
end subroutine pop_pack_r8_2_r8

subroutine pop_pack_r8_2_r8_map(indata,outdata,tile,nmp)
  real(kind=8), dimension(:), intent(in) :: indata
  real(kind=8), dimension(:), intent(inout) :: outdata
  integer, intent(in) :: tile
  integer, dimension(ifull), intent(in) :: nmp
  integer nb, is, ie, js, je, iqx, iqt, iq

  if ( size(indata)==ifull ) then
    js = (tile-1)*imax + 1
    je = tile*imax
  else if ( size(indata)==imax ) then
    js = 1
    je = imax
  else
    write(6,*) "Invalid input lengh for pop_pack (map)"
    stop
  end if
  
  if ( all( nmp(js)==nmp(js:je) ) ) then
    nb = nmp(js)  
    is = tdata(tile)%pind(nb,1)
    ie = tdata(tile)%pind(nb,2)
    if ( is<=ie ) then
      outdata(is:ie) =  pack(indata(js:je),tdata(tile)%pmap(:,nb))
    end if  
  else  
    do iq = js,je
      nb = nmp(iq)
      iqx = iq-js+1
      if ( tdata(tile)%pmap(iqx,nb) ) then
        is = tdata(tile)%pind(nb,1)  
        iqt = is + count( tdata(tile)%pmap(1:iqx,nb) ) - 1
        outdata(iqt) = indata(iq)
      end if
    end do 
  end if  

end subroutine pop_pack_r8_2_r8_map

subroutine pop_pack_i4_2_i4(indata,outdata,tile,nb)
  integer, dimension(:), intent(in) :: indata
  integer(kind=4), dimension(:), intent(inout) :: outdata
  integer, intent(in) :: tile
  integer, intent(in), optional :: nb
  integer is, ie, js, je, lnb

  if ( size(indata)==ifull ) then
    js = (tile-1)*imax + 1
    je = tile*imax
  else if ( size(indata)==imax ) then
    js = 1
    je = imax
  else
    write(6,*) "Invalid input lengh for pop_pack (tile)"
    stop
  end if
  
  if ( present(nb) ) then

    is = tdata(tile)%pind(nb,1)
    ie = tdata(tile)%pind(nb,2)
    if ( is<=ie ) then
      outdata(is:ie) =  pack(int(indata(js:je),4),tdata(tile)%pmap(:,nb))
    end if  

  else
      
    do lnb = 1,tdata(tile)%maxnb
      is = tdata(tile)%pind(lnb,1)
      ie = tdata(tile)%pind(lnb,2)
      if ( is<=ie ) then
        outdata(is:ie) =  pack(int(indata(js:je),4),tdata(tile)%pmap(:,lnb))
      end if  
    end do
      
  end if

end subroutine pop_pack_i4_2_i4

subroutine pop_pack_i4_2_i4_map(indata,outdata,tile,nmp)
  integer, dimension(:), intent(in) :: indata
  integer(kind=4), dimension(:), intent(inout) :: outdata
  integer, intent(in) :: tile
  integer, dimension(ifull), intent(in) :: nmp
  integer nb, is, ie, js, je, iqx, iqt, iq

  if ( size(indata)==ifull ) then
    js = (tile-1)*imax + 1
    je = tile*imax
  else if ( size(indata)==imax ) then
    js = 1
    je = imax
  else
    write(6,*) "Invalid input lengh for pop_pack (map)"
    stop
  end if
  
  if ( all( nmp(js)==nmp(js:je) ) ) then
    nb = nmp(js)  
    is = tdata(tile)%pind(nb,1)
    ie = tdata(tile)%pind(nb,2)
    if ( is<=ie ) then
      outdata(is:ie) =  pack(indata(js:je),tdata(tile)%pmap(:,nb))
    end if  
  else  
    do iq = js,je
      nb = nmp(iq)
      iqx = iq-js+1
      if ( tdata(tile)%pmap(iqx,nb) ) then
        is = tdata(tile)%pind(nb,1)  
        iqt = is + count( tdata(tile)%pmap(1:iqx,nb) ) - 1
        outdata(iqt) = indata(iq)
      end if
    end do 
  end if  

end subroutine pop_pack_i4_2_i4_map

subroutine pop_pack_i4_2_i8(indata,outdata,tile,nb)
  integer, dimension(:), intent(in) :: indata
  integer(kind=8), dimension(:), intent(inout) :: outdata
  integer, intent(in) :: tile
  integer, intent(in), optional :: nb
  integer is, ie, js, je, lnb

  if ( size(indata)==ifull ) then
    js = (tile-1)*imax + 1
    je = tile*imax
  else if ( size(indata)==imax ) then
    js = 1
    je = imax
  else
    write(6,*) "Invalid input lengh for pop_pack (tile)"
    stop
  end if
  
  if ( present(nb) ) then

    is = tdata(tile)%pind(nb,1)
    ie = tdata(tile)%pind(nb,2)
    if ( is<=ie ) then
      outdata(is:ie) =  pack(int(indata(js:je),8),tdata(tile)%pmap(:,nb))
    end if  
      
  else

    do lnb = 1,tdata(tile)%maxnb
      is = tdata(tile)%pind(lnb,1)
      ie = tdata(tile)%pind(lnb,2)
      if ( is<=ie ) then
        outdata(is:ie) =  pack(int(indata(js:je),8),tdata(tile)%pmap(:,lnb))
      end if  
    end do    
      
  end if
  

end subroutine pop_pack_i4_2_i8

subroutine pop_pack_i4_2_i8_map(indata,outdata,tile,nmp)
  integer, dimension(:), intent(in) :: indata
  integer(kind=8), dimension(:), intent(inout) :: outdata
  integer, intent(in) :: tile
  integer, dimension(ifull), intent(in) :: nmp
  integer nb, is, ie, js, je, iqx, iqt, iq

  if ( size(indata)==ifull ) then
    js = (tile-1)*imax + 1
    je = tile*imax
  else if ( size(indata)==imax ) then
    js = 1
    je = imax
  else
    write(6,*) "Invalid input lengh for pop_pack (map)"
    stop
  end if
  
  if ( all( nmp(js)==nmp(js:je) ) ) then
    nb = nmp(js)  
    is = tdata(tile)%pind(nb,1)
    ie = tdata(tile)%pind(nb,2)
    if ( is<=ie ) then
      outdata(is:ie) =  pack(indata(js:je),tdata(tile)%pmap(:,nb))
    end if  
  else  
    do iq = js,je
      nb = nmp(iq)
      iqx = iq-js+1
      if ( tdata(tile)%pmap(iqx,nb) ) then
        is = tdata(tile)%pind(nb,1)  
        iqt = is + count( tdata(tile)%pmap(1:iqx,nb) ) - 1
        outdata(iqt) = indata(iq)
      end if
    end do 
  end if  

end subroutine pop_pack_i4_2_i8_map

subroutine pop_unpack_r8_2_r8(indata,outdata,tile,nb)
  real(kind=8), dimension(:), intent(in) :: indata
  real(kind=8), dimension(:), intent(out) :: outdata
  integer, intent(in) :: tile
  integer, intent(in), optional :: nb
  integer is, ie, js, je, lnb

  if ( size(outdata)==ifull ) then
    js = (tile-1)*imax + 1
    je = tile*imax
  else if ( size(outdata)==imax ) then
    js = 1
    je = imax
  else
    write(6,*) "Invalid input lengh for pop_unpack (tile)"
    stop
  end if

  if ( present(nb) ) then

    is = tdata(tile)%pind(nb,1)
    ie = tdata(tile)%pind(nb,2)
    if ( is<=ie ) then
      outdata(js:je) = unpack(indata(is:ie),tdata(tile)%pmap(:,nb),outdata(js:je))
    end if  
      
  else

    do lnb = 1,tdata(tile)%maxnb
      is = tdata(tile)%pind(nb,1)
      ie = tdata(tile)%pind(nb,2)
      if ( is<=ie ) then
        outdata(js:je) = unpack(indata(is:ie),tdata(tile)%pmap(:,nb),outdata(js:je))
      end if  
    end do    
      
  end if
  
end subroutine pop_unpack_r8_2_r8

subroutine pop_unpack_i4_2_r8(indata,outdata,tile,nb)
  integer(kind=4), dimension(:), intent(in) :: indata
  real(kind=8), dimension(ifull), intent(out) :: outdata
  integer, intent(in) :: tile
  integer, intent(in), optional :: nb
  integer is, ie, js, je, lnb

  if ( size(outdata)==ifull ) then
    js = (tile-1)*imax + 1
    je = tile*imax
  else if ( size(outdata)==imax ) then
    js = 1
    je = imax
  else
    write(6,*) "Invalid input lengh for pop_unpack (tile)"
    stop
  end if  
  
  if ( present(nb) ) then

    is = tdata(tile)%pind(nb,1)
    ie = tdata(tile)%pind(nb,2)
    if ( is<=ie ) then
      outdata(js:je) = unpack(real(indata(is:ie),8),tdata(tile)%pmap(:,nb),outdata(js:je))
    end if  
      
  else

    do lnb = 1,tdata(tile)%maxnb
      is = tdata(tile)%pind(lnb,1)
      ie = tdata(tile)%pind(lnb,2)
      if ( is<=ie ) then
        outdata(js:je) = unpack(real(indata(is:ie),8),tdata(tile)%pmap(:,lnb),outdata(js:je))
      end if  
  
    end do    
      
  end if

end subroutine pop_unpack_i4_2_r8

subroutine pop_unpack_i8_2_r8(indata,outdata,tile,nb)
  integer(kind=8), dimension(:), intent(in) :: indata
  real(kind=8), dimension(ifull), intent(out) :: outdata
  integer, intent(in) :: tile
  integer, intent(in), optional :: nb
  integer is, ie, js, je, lnb

  if ( size(outdata)==ifull ) then
    js = (tile-1)*imax + 1
    je = tile*imax
  else if ( size(outdata)==imax ) then
    js = 1
    je = imax
  else
    write(6,*) "Invalid input lengh for pop_unpack (tile)"
    stop
  end if  
  
  if ( present(nb) ) then

    is = tdata(tile)%pind(nb,1)
    ie = tdata(tile)%pind(nb,2)
    if ( is<=ie ) then
      outdata(js:je) = unpack(real(indata(is:ie),8),tdata(tile)%pmap(:,nb),outdata(js:je))
    end if  
      
  else

    do lnb = 1,tdata(tile)%maxnb  
      is = tdata(tile)%pind(lnb,1)
      ie = tdata(tile)%pind(lnb,2)
      if ( is<=ie ) then
        outdata(js:je) = unpack(real(indata(is:ie),8),tdata(tile)%pmap(:,lnb),outdata(js:je))
      end if  
    end do  
      
  end if

end subroutine pop_unpack_i8_2_r8

subroutine cable_update_r8_2_r4(outdata,indata,tile,nb)
  real, dimension(:), intent(inout) :: outdata
  real(kind=8), dimension(:), intent(in) :: indata
  integer, intent(in) :: tile
  integer, intent(in), optional :: nb
  integer js, je, is, ie, lnb
  
  if ( size(outdata)==ifull ) then
    js = (tile-1)*imax + 1
    je = tile*imax
  else if ( size(outdata)==imax ) then
    js = 1
    je = imax
  else
    write(6,*) "Invalid input lengh for cable_update (tile)"
    stop
  end if  
  
  if ( present(nb) ) then

    is = tdata(tile)%tind(nb,1)
    ie = tdata(tile)%tind(nb,2)
    if ( is<=ie ) then
      outdata(js:je) = outdata(js:je) + unpack(tdata(tile)%sv(is:ie)*real(indata(is:ie)),tdata(tile)%tmap(:,nb),0.)
    end if

  else
  
    do lnb = 1,tdata(tile)%maxnb
      is = tdata(tile)%tind(lnb,1)
      ie = tdata(tile)%tind(lnb,2)
      if ( is<=ie ) then
        outdata(js:je) = outdata(js:je) + unpack(tdata(tile)%sv(is:ie)*real(indata(is:ie)),tdata(tile)%tmap(:,lnb),0.)
      end if
    end do  

  end if

end subroutine cable_update_r8_2_r4

#ifndef i8r8
subroutine cable_update_r8_2_r8(outdata,indata,tile,nb)
  real(kind=8), dimension(imax), intent(inout) :: outdata
  real(kind=8), dimension(:), intent(in) :: indata
  integer, intent(in) :: tile
  integer, intent(in), optional :: nb
  integer js, je, is, ie, lnb
  
  if ( size(outdata)==ifull ) then
    js = (tile-1)*imax + 1
    je = tile*imax
  else if ( size(outdata)==imax ) then
    js = 1
    je = imax
  else
    write(6,*) "Invalid input lengh for pop_unpack (tile)"
    stop
  end if  
  
  if ( present(nb) ) then

    is = tdata(tile)%tind(nb,1)
    ie = tdata(tile)%tind(nb,2)
    if ( is<=ie ) then
      outdata(js:je) = outdata(js:je) + unpack(tdata(tile)%sv(is:ie)*real(indata(is:ie)),tdata(tile)%tmap(:,nb),0.)
    end if

  else
 
    do lnb = 1,tdata(tile)%maxnb
      is = tdata(tile)%tind(lnb,1)
      ie = tdata(tile)%tind(lnb,2)
      if ( is<=ie ) then
        outdata(js:je) = outdata(js:je) + unpack(tdata(tile)%sv(is:ie)*real(indata(is:ie)),tdata(tile)%tmap(:,lnb),0.)
      end if
    end do  

  end if

end subroutine cable_update_r8_2_r8
#endif

end module cable_ccam4
