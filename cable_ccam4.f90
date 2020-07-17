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
    
module cable_ccam4
 
use cc_omp, only : ntiles,imax
use cable_def_types_mod, only : air_type, balances_type, canopy_type, climate_type, met_type, radiation_type, &
                                roughness_type, soil_parameter_type, soil_snow_type, sum_flux_type,           &
                                veg_parameter_type, mp
use casavariable, only : casa_balance, casa_biome, casa_flux, casa_met, casa_pool
use phenvariable, only : phen_variable
use pop_types, only : pop_type, dp
use cable_ccam3, only : climate_save_type

implicit none

private
public tdata
public :: cable_pack, cable_unpack, pop_pack, pop_unpack, setp
    
type tiledata
  integer, dimension(:,:), allocatable :: tind
  integer, dimension(:,:), allocatable :: pind
  logical, dimension(:,:), allocatable :: tmap
  logical, dimension(:,:), allocatable :: pmap
  integer :: mp
  integer :: np
  integer :: maxnb
  integer :: toffset
  integer :: poffset
end type tiledata

TYPE(tiledata), dimension(:), allocatable :: tdata

interface cable_pack
  module procedure cable_pack_r4_2_r4, cable_pack_r4_2_r8, &
                   cable_pack_i4_2_i4, cable_pack_r4_2_i4, &
                   cable_pack_i4_2_i8, cable_pack_r4_2_i8

#ifndef i8r8
  module procedure cable_pack_r8_2_r8
#endif
end interface

interface cable_unpack
  module procedure cable_unpack_r4_2_r4, cable_unpack_r8_2_r4
  module procedure cable_unpack_r4_2_r4_tile, cable_unpack_r8_2_r4_tile
#ifndef i8r8
  module procedure cable_unpack_r4_2_r8
  module procedure cable_unpack_r8_2_r8_tile
#endif
end interface

interface pop_pack
  module procedure pop_pack_r8_2_r8_tile, pop_pack_i4_2_i4_tile
end interface

interface pop_unpack
  module procedure pop_unpack_r8_2_r8_tile, pop_unpack_i4_2_r8_tile
  module procedure pop_unpack_i8_2_r8_tile
end interface

interface setp
  module procedure setp_air, setp_bal, setp_canopy, setp_casabal,            &
                   setp_casaflux, setp_casamet, setp_casapool, setp_climate, &
                   setp_met, setp_phen, setp_pop, setp_rad, setp_rough,      &
                   setp_soil, setp_ssnow, setp_sumflux, setp_veg,            &
                   setp_climate_save
end interface

contains

subroutine cable_pack_r4_2_r4(indata,outdata,inb)
  use newmpar_m, only : ifull

  implicit none

  real, dimension(ifull), intent(in) :: indata
  real(kind=4), dimension(:), intent(out) :: outdata
  integer, intent(in), optional :: inb
  integer :: nb, is, ie, js, je, tile

  if ( present(inb) ) then
    nb = inb
    do tile = 1,ntiles
      js = 1 + (tile-1)*imax
      je = tile*imax
      is = tdata(tile)%tind(nb,1)
      ie = tdata(tile)%tind(nb,2)
      if ( is<=ie ) then
        outdata(is:ie) =  pack(indata(js:je),tdata(tile)%tmap(:,nb))
      end if  
    end do
  else
    do tile = 1,ntiles
      js = 1 + (tile-1)*imax
      je = tile*imax
      do nb = 1,tdata(tile)%maxnb
        is = tdata(tile)%tind(nb,1)
        ie = tdata(tile)%tind(nb,2)
        if ( is<=ie ) then
          outdata(is:ie) =  pack(indata(js:je),tdata(tile)%tmap(:,nb))
        end if  
      end do
    end do
  end if

end subroutine cable_pack_r4_2_r4

subroutine cable_pack_r4_2_r8(indata,outdata,inb)
  use newmpar_m, only : ifull

  implicit none

  real, dimension(ifull), intent(in) :: indata
  real(kind=8), dimension(:), intent(out) :: outdata
  integer, intent(in), optional :: inb
  integer :: nb, is, ie, js, je, tile

  if ( present(inb) ) then
    nb = inb
    do tile = 1,ntiles
      js = 1 + (tile-1)*imax
      je = tile*imax
      is = tdata(tile)%tind(nb,1)
      ie = tdata(tile)%tind(nb,2)
      if ( is<=ie ) then
        outdata(is:ie) = pack(real(indata(js:je),8),tdata(tile)%tmap(:,nb))
      end if  
    end do
  else
    do tile = 1,ntiles
      js = 1 + (tile-1)*imax
      je = tile*imax
      do nb = 1,tdata(tile)%maxnb
        is = tdata(tile)%tind(nb,1)
        ie = tdata(tile)%tind(nb,2)
        if ( is<=ie ) then
          outdata(is:ie) = pack(real(indata(js:je),8),tdata(tile)%tmap(:,nb))
        end if  
      end do
    end do
  end if

end subroutine cable_pack_r4_2_r8

#ifndef i8r8
subroutine cable_pack_r8_2_r8(indata,outdata,inb)
  use newmpar_m, only : ifull

  implicit none

  real(kind=8), dimension(ifull), intent(in) :: indata
  real(kind=8), dimension(:), intent(out) :: outdata
  integer, intent(in), optional :: inb
  integer :: nb, is, ie, js, je, tile

  if ( present(inb) ) then
    nb = inb
    do tile = 1,ntiles
      js=1+(tile-1)*imax
      je=tile*imax
      is = tdata(tile)%tind(nb,1)
      ie = tdata(tile)%tind(nb,2)
      if ( is<=ie ) then
        outdata(is:ie) =  pack(indata(js:je),tdata(tile)%tmap(:,nb))
      end if  
    end do
  else
    do tile = 1,ntiles
      js=1+(tile-1)*imax
      je=tile*imax
      do nb = 1,tdata(tile)%maxnb
        is = tdata(tile)%tind(nb,1)
        ie = tdata(tile)%tind(nb,2)
        if ( is<=ie ) then
          outdata(is:ie) =  pack(indata(js:je),tdata(tile)%tmap(:,nb))
        end if  
      end do
    end do
  end if

end subroutine cable_pack_r8_2_r8
#endif

subroutine cable_pack_i4_2_i4(indata,outdata,inb)
  use newmpar_m, only : ifull

  implicit none

  integer, dimension(ifull), intent(in) :: indata
  integer(kind=4), dimension(:), intent(out) :: outdata
  integer, intent(in), optional :: inb
  integer :: nb, is, ie, js, je, tile

  if ( present(inb) ) then
    nb = inb
    do tile = 1,ntiles
      js=1+(tile-1)*imax
      je=tile*imax
      is = tdata(tile)%tind(nb,1)
      ie = tdata(tile)%tind(nb,2)
      if ( is<=ie ) then
        outdata(is:ie) =  pack(indata(js:je),tdata(tile)%tmap(:,nb))
      end if  
    end do
  else
    do tile = 1,ntiles
      js=1+(tile-1)*imax
      je=tile*imax
      do nb = 1,tdata(tile)%maxnb
        is = tdata(tile)%tind(nb,1)
        ie = tdata(tile)%tind(nb,2)
        if ( is<=ie ) then
          outdata(is:ie) =  pack(indata(js:je),tdata(tile)%tmap(:,nb))
        end if  
      end do
    end do
  end if

end subroutine cable_pack_i4_2_i4

subroutine cable_pack_i4_2_i8(indata,outdata,inb)
  use newmpar_m, only : ifull

  implicit none

  integer, dimension(ifull), intent(in) :: indata
  integer(kind=8), dimension(:), intent(out) :: outdata
  integer, intent(in), optional :: inb
  integer :: nb, is, ie, js, je, tile

  if ( present(inb) ) then
    nb = inb
    do tile = 1,ntiles
      js=1+(tile-1)*imax
      je=tile*imax
      is = tdata(tile)%tind(nb,1)
      ie = tdata(tile)%tind(nb,2)
      if ( is<=ie ) then
        outdata(is:ie) =  pack(indata(js:je),tdata(tile)%tmap(:,nb))
      end if  
    end do
  else
    do tile = 1,ntiles
      js=1+(tile-1)*imax
      je=tile*imax
      do nb = 1,tdata(tile)%maxnb
        is = tdata(tile)%tind(nb,1)
        ie = tdata(tile)%tind(nb,2)
        if ( is<=ie ) then
          outdata(is:ie) =  pack(indata(js:je),tdata(tile)%tmap(:,nb))
        end if  
      end do
    end do
  end if

end subroutine cable_pack_i4_2_i8

subroutine cable_pack_r4_2_i4(indata,outdata,inb)
  use newmpar_m, only : ifull

  implicit none

  real, dimension(ifull), intent(in) :: indata
  integer(kind=4), dimension(:), intent(out) :: outdata
  integer, intent(in), optional :: inb
  integer :: nb, is, ie, js, je, tile

  if ( present(inb) ) then
    nb = inb
    do tile = 1,ntiles
      js=1+(tile-1)*imax
      je=tile*imax
      is = tdata(tile)%tind(nb,1)
      ie = tdata(tile)%tind(nb,2)
      if ( is<=ie ) then
        outdata(is:ie) =  nint(pack(indata(js:je),tdata(tile)%tmap(:,nb)))
      end if  
    end do
  else
    do tile = 1,ntiles
      js=1+(tile-1)*imax
      je=tile*imax
      do nb = 1,tdata(tile)%maxnb
        is = tdata(tile)%tind(nb,1)
        ie = tdata(tile)%tind(nb,2)
        if ( is<=ie ) then
          outdata(is:ie) =  nint(pack(indata(js:je),tdata(tile)%tmap(:,nb)))
        end if  
      end do
    end do
  end if

end subroutine cable_pack_r4_2_i4

subroutine cable_pack_r4_2_i8(indata,outdata,inb)
  use newmpar_m, only : ifull

  implicit none

  real, dimension(ifull), intent(in) :: indata
  integer(kind=8), dimension(:), intent(out) :: outdata
  integer, intent(in), optional :: inb
  integer :: nb, is, ie, js, je, tile

  if ( present(inb) ) then
    nb = inb
    do tile = 1,ntiles
      js=1+(tile-1)*imax
      je=tile*imax
      is = tdata(tile)%tind(nb,1)
      ie = tdata(tile)%tind(nb,2)
      if ( is<=ie ) then
        outdata(is:ie) =  nint(pack(indata(js:je),tdata(tile)%tmap(:,nb)))
      end if  
    end do
  else
    do tile = 1,ntiles
      js=1+(tile-1)*imax
      je=tile*imax
      do nb = 1,tdata(tile)%maxnb
        is = tdata(tile)%tind(nb,1)
        ie = tdata(tile)%tind(nb,2)
        if ( is<=ie ) then
          outdata(is:ie) =  nint(pack(indata(js:je),tdata(tile)%tmap(:,nb)))
        end if  
      end do
    end do
  end if

end subroutine cable_pack_r4_2_i8

subroutine cable_unpack_r4_2_r4(indata,outdata)
  use newmpar_m, only : ifull

  implicit none

  real(kind=4), dimension(:), intent(in) :: indata
  real, dimension(ifull), intent(inout) :: outdata
  integer :: is, ie, js, je, tile, nb

  do tile = 1,ntiles
    js = 1 + (tile-1)*imax
    je = tile*imax
    do nb = 1,tdata(tile)%maxnb
      is = tdata(tile)%tind(nb,1)
      ie = tdata(tile)%tind(nb,2)
      if ( is<=ie ) then
        outdata(js:je) = outdata(js:je) + unpack(indata(is:ie),tdata(tile)%tmap(:,nb),0.)
      end if  
    end do
  end do

end subroutine cable_unpack_r4_2_r4

subroutine cable_unpack_r8_2_r4(indata,outdata)
  use newmpar_m, only : ifull

  implicit none

  real(kind=8), dimension(:), intent(in) :: indata
  real, dimension(ifull), intent(inout) :: outdata
  integer :: is, ie, js, je, tile, nb

  do tile =1,ntiles
    js=1+(tile-1)*imax
    je=tile*imax
    do nb = 1,tdata(tile)%maxnb
      is = tdata(tile)%tind(nb,1)
      ie = tdata(tile)%tind(nb,2)
      if ( is<=ie ) then
        outdata(js:je) = outdata(js:je) + real(unpack(indata(is:ie),tdata(tile)%tmap(:,nb),0._8),4)
      end if  
    end do
  end do

end subroutine cable_unpack_r8_2_r4

subroutine cable_unpack_r4_2_r4_tile(indata,outdata,inb)
  use newmpar_m, only : ifull

  implicit none

  real(kind=4), dimension(:), intent(in) :: indata
  real, dimension(ifull), intent(inout) :: outdata
  integer, intent(in) :: inb
  integer :: is, ie, js, je, tile, nb

  nb = inb
  do tile = 1,ntiles
    is = tdata(tile)%tind(nb,1)
    ie = tdata(tile)%tind(nb,2)
    if ( is<=ie ) then
      js=1+(tile-1)*imax
      je=tile*imax
      outdata(js:je) = unpack(indata(is:ie),tdata(tile)%tmap(:,nb),outdata(js:je))
    end if  
  end do

end subroutine cable_unpack_r4_2_r4_tile

subroutine cable_unpack_r8_2_r4_tile(indata,outdata,inb)
  use newmpar_m, only : ifull

  implicit none

  real(kind=8), dimension(:), intent(in) :: indata
  real, dimension(ifull), intent(inout) :: outdata
  integer, intent(in) :: inb
  integer :: is, ie, js, je, tile, nb

  nb = inb
  do tile = 1,ntiles
    is = tdata(tile)%tind(nb,1)
    ie = tdata(tile)%tind(nb,2)
    if ( is<=ie ) then
      js=1+(tile-1)*imax
      je=tile*imax
      outdata(js:je) = real(unpack(indata(is:ie),tdata(tile)%tmap(:,nb),real(outdata(js:je),8)),4)
    end if  
  end do

end subroutine cable_unpack_r8_2_r4_tile

#ifndef i8r8
subroutine cable_unpack_r4_2_r8(indata,outdata)
  use newmpar_m, only : ifull

  implicit none

  real(kind=4), dimension(:), intent(in) :: indata
  real(kind=8), dimension(ifull), intent(inout) :: outdata
  integer :: is, ie, js, je, tile, nb

  do tile =1,ntiles
    js=1+(tile-1)*imax
    je=tile*imax
    do nb = 1,tdata(tile)%maxnb
      is = tdata(tile)%tind(nb,1)
      ie = tdata(tile)%tind(nb,2)
      if ( is<=ie ) then
        outdata(js:je) = outdata(js:je) + unpack(indata(is:ie),tdata(tile)%tmap(:,nb),0._4)
      end if  
    end do
  end do

end subroutine cable_unpack_r4_2_r8

subroutine cable_unpack_r8_2_r8_tile(indata,outdata,inb)
  use newmpar_m, only : ifull

  implicit none

  real(kind=8), dimension(:), intent(in) :: indata
  real(kind=8), dimension(ifull), intent(inout) :: outdata
  integer, intent(in) :: inb
  integer :: is, ie, js, je, tile, nb

  nb = inb
  do tile = 1,ntiles
    is = tdata(tile)%tind(nb,1)
    ie = tdata(tile)%tind(nb,2)
    if ( is<=ie ) then
      js=1+(tile-1)*imax
      je=tile*imax
      outdata(js:je) = unpack(indata(is:ie),tdata(tile)%tmap(:,nb),outdata(js:je))
    end if  
  end do

end subroutine cable_unpack_r8_2_r8_tile
#endif

subroutine pop_pack_r8_2_r8_tile(indata,outdata,inb)
  use newmpar_m, only : ifull
  use TypeDef, only : dp

  implicit none

  real(kind=8), dimension(ifull), intent(in) :: indata
  real(kind=dp), dimension(:), intent(inout) :: outdata
  integer, intent(in) :: inb
  integer :: nb, is, ie, js, je, tile

  nb = inb
  do tile = 1,ntiles
    is = tdata(tile)%pind(nb,1)
    ie = tdata(tile)%pind(nb,2)
    if ( is<=ie ) then
      js=1+(tile-1)*imax
      je=tile*imax
      outdata(is:ie) =  pack(real(indata(js:je),dp),tdata(tile)%pmap(:,nb))
    end if  
  end do

end subroutine pop_pack_r8_2_r8_tile

subroutine pop_pack_i4_2_i4_tile(indata,outdata,inb)
  use newmpar_m, only : ifull

  implicit none

  integer, dimension(ifull), intent(in) :: indata
  integer(kind=4), dimension(:), intent(inout) :: outdata
  integer, intent(in) :: inb
  integer :: nb, is, ie, js, je, tile

  nb = inb
  do tile = 1,ntiles
    is = tdata(tile)%pind(nb,1)
    ie = tdata(tile)%pind(nb,2)
    if ( is<=ie ) then
      js=1+(tile-1)*imax
      je=tile*imax
      outdata(is:ie) =  pack(indata(js:je),tdata(tile)%pmap(:,nb))
    end if  
  end do

end subroutine pop_pack_i4_2_i4_tile

subroutine pop_unpack_r8_2_r8_tile(indata,outdata,inb)
  use newmpar_m, only : ifull
  use TypeDef, only : dp

  implicit none

  real(kind=dp), dimension(:), intent(in) :: indata
  real(kind=8), dimension(ifull), intent(out) :: outdata
  integer, intent(in) :: inb
  integer :: is, ie, js, je, tile, nb

  nb = inb
  do tile =1,ntiles
    is = tdata(tile)%pind(nb,1)
    ie = tdata(tile)%pind(nb,2)
    if ( is<=ie ) then
      js=1+(tile-1)*imax
      je=tile*imax
      outdata(js:je) = unpack(real(indata(is:ie),8),tdata(tile)%pmap(:,nb),0._8)
    end if  
  end do

end subroutine pop_unpack_r8_2_r8_tile

subroutine pop_unpack_i4_2_r8_tile(indata,outdata,inb)
  use newmpar_m, only : ifull

  implicit none

  integer(kind=4), dimension(:), intent(in) :: indata
  real(kind=8), dimension(ifull), intent(out) :: outdata
  integer, intent(in) :: inb
  integer :: is, ie, js, je, tile, nb

  nb = inb
  do tile =1,ntiles
    is = tdata(tile)%pind(nb,1)
    ie = tdata(tile)%pind(nb,2)
    if ( is<=ie ) then
      js=1+(tile-1)*imax
      je=tile*imax
      outdata(js:je) = unpack(real(indata(is:ie),8),tdata(tile)%pmap(:,nb),0._8)
    end if  
  end do

end subroutine pop_unpack_i4_2_r8_tile

subroutine pop_unpack_i8_2_r8_tile(indata,outdata,inb)
  use newmpar_m, only : ifull

  implicit none

  integer(kind=8), dimension(:), intent(in) :: indata
  real(kind=8), dimension(ifull), intent(out) :: outdata
  integer, intent(in) :: inb
  integer :: is, ie, js, je, tile, nb

  nb = inb
  do tile =1,ntiles
    is = tdata(tile)%pind(nb,1)
    ie = tdata(tile)%pind(nb,2)
    if ( is<=ie ) then
      js=1+(tile-1)*imax
      je=tile*imax
      outdata(js:je) = unpack(real(indata(is:ie),8),tdata(tile)%pmap(:,nb),0._8)
    end if  
  end do

end subroutine pop_unpack_i8_2_r8_tile

subroutine setp_air(air,lair,tile)
  implicit none

  type(air_type), intent(in) :: air
  type(air_type), intent(inout) :: lair
  integer, intent(in) :: tile
  integer :: is, ie

  is = tdata(tile)%toffset + 1
  ie = tdata(tile)%toffset + tdata(tile)%mp

  if ( is<=ie ) then
  
    lair%rho => air%rho(is:ie)
    lair%volm => air%volm(is:ie)
    lair%rlam => air%rlam(is:ie)
    lair%qsat => air%qsat(is:ie)
    lair%epsi => air%epsi(is:ie)
    lair%visc => air%visc(is:ie)
    lair%psyc => air%psyc(is:ie)
    lair%dsatdk => air%dsatdk(is:ie)
    lair%cmolar => air%cmolar(is:ie)
    
  end if  
 
end subroutine setp_air

subroutine setp_bal(bal,lbal,tile)
  implicit none

  type(balances_type), intent(in) :: bal
  type(balances_type), intent(inout) :: lbal
  integer, intent(in) :: tile
  integer :: is, ie

  is = tdata(tile)%toffset + 1
  ie = tdata(tile)%toffset + tdata(tile)%mp

  if ( is<=ie ) then
  
    lbal%drybal => bal%drybal(is:ie)
    lbal%ebal => bal%ebal(is:ie)
    lbal%ebal_tot => bal%ebal_tot(is:ie)
    lbal%ebal_cncheck => bal%ebal_cncheck(is:ie)
    lbal%ebal_tot_cncheck => bal%ebal_tot_cncheck(is:ie)
    lbal%ebaltr => bal%ebaltr(is:ie)
    lbal%ebal_tottr => bal%ebal_tottr(is:ie)
    lbal%evap_tot => bal%evap_tot(is:ie)
    lbal%osnowd0 => bal%osnowd0(is:ie)
    lbal%precip_tot => bal%precip_tot(is:ie)
    lbal%rnoff_tot => bal%rnoff_tot(is:ie)
    lbal%wbal => bal%wbal(is:ie)
    lbal%wbal_tot => bal%wbal_tot(is:ie)
    lbal%wbtot0 => bal%wbtot0(is:ie)
    lbal%wetbal => bal%wetbal(is:ie)
    lbal%cansto0 => bal%cansto0(is:ie)
    lbal%owbtot => bal%owbtot(is:ie)
    lbal%evapc_tot => bal%evapc_tot(is:ie)
    lbal%evaps_tot => bal%evaps_tot(is:ie)
    lbal%rnof1_tot => bal%rnof1_tot(is:ie)
    lbal%rnof2_tot => bal%rnof2_tot(is:ie)
    lbal%snowdc_tot => bal%snowdc_tot(is:ie)
    lbal%wbal_tot1 => bal%wbal_tot1(is:ie)
    lbal%delwc_tot => bal%delwc_tot(is:ie)
    lbal%qasrf_tot => bal%qasrf_tot(is:ie)
    lbal%qfsrf_tot => bal%qfsrf_tot(is:ie)
    lbal%qssrf_tot => bal%qssrf_tot(is:ie)
    lbal%Radbal => bal%Radbal(is:ie)
    lbal%EbalSoil => bal%EbalSoil(is:ie)
    lbal%Ebalveg => bal%Ebalveg(is:ie)
    lbal%Radbalsum => bal%Radbalsum(is:ie)
    
  end if  
 
end subroutine setp_bal

subroutine setp_canopy(canopy,lcanopy,tile)
  implicit none

  type(canopy_type), intent(in) :: canopy
  type(canopy_type), intent(inout) :: lcanopy
  integer, intent(in) :: tile
  integer :: is, ie

  is = tdata(tile)%toffset + 1
  ie = tdata(tile)%toffset + tdata(tile)%mp

  if ( is<=ie ) then
  
    lcanopy%cansto => canopy%cansto(is:ie)
    lcanopy%cduv => canopy%cduv(is:ie)
    lcanopy%delwc => canopy%delwc(is:ie)
    lcanopy%dewmm => canopy%dewmm(is:ie)
    lcanopy%fe => canopy%fe(is:ie)
    lcanopy%fh => canopy%fh(is:ie)
    lcanopy%fpn => canopy%fpn(is:ie)
    lcanopy%frp => canopy%frp(is:ie)
    lcanopy%frpw => canopy%frpw(is:ie)
    lcanopy%frpr => canopy%frpr(is:ie)
    lcanopy%frs => canopy%frs(is:ie)
    lcanopy%fnee => canopy%fnee(is:ie)
    lcanopy%frday => canopy%frday(is:ie)
    lcanopy%fnv => canopy%fnv(is:ie)
    lcanopy%fev => canopy%fev(is:ie)
    lcanopy%epot => canopy%epot(is:ie)
    lcanopy%fnpp => canopy%fnpp(is:ie)
    lcanopy%fevw_pot => canopy%fevw_pot(is:ie)
    lcanopy%gswx_T => canopy%gswx_T(is:ie)
    lcanopy%cdtq => canopy%cdtq(is:ie)
    lcanopy%wetfac_cs => canopy%wetfac_cs(is:ie)
    lcanopy%fevw => canopy%fevw(is:ie)
    lcanopy%fhvw => canopy%fhvw(is:ie)
    lcanopy%oldcansto => canopy%oldcansto(is:ie)
    lcanopy%fhv => canopy%fhv(is:ie)
    lcanopy%fns => canopy%fns(is:ie)
    lcanopy%fhs => canopy%fhs(is:ie)
    lcanopy%fhs_cor => canopy%fhs_cor(is:ie)
    lcanopy%ga => canopy%ga(is:ie)
    lcanopy%ghflux => canopy%ghflux(is:ie)
    lcanopy%precis => canopy%precis(is:ie)
    lcanopy%qscrn => canopy%qscrn(is:ie)
    lcanopy%rnet => canopy%rnet(is:ie)
    lcanopy%rniso => canopy%rniso(is:ie)
    lcanopy%segg => canopy%segg(is:ie)
    lcanopy%sghflux => canopy%sghflux(is:ie)
    lcanopy%through => canopy%through(is:ie)
    lcanopy%spill => canopy%spill(is:ie)
    lcanopy%tscrn => canopy%tscrn(is:ie)
    lcanopy%wcint => canopy%wcint(is:ie)
    lcanopy%tv => canopy%tv(is:ie)
    lcanopy%us => canopy%us(is:ie)
    lcanopy%uscrn => canopy%uscrn(is:ie)
    lcanopy%vlaiw => canopy%vlaiw(is:ie)
    lcanopy%rghlai => canopy%rghlai(is:ie)
    lcanopy%fwet => canopy%fwet(is:ie)

    lcanopy%evapfbl => canopy%evapfbl(is:ie,:)
    lcanopy%gswx => canopy%gswx(is:ie,:)
    lcanopy%zetar => canopy%zetar(is:ie,:)
    lcanopy%zetash => canopy%zetash(is:ie,:)

    lcanopy%fess => canopy%fess(is:ie)
    lcanopy%fesp => canopy%fesp(is:ie)
    lcanopy%dgdtg => canopy%dgdtg(is:ie)
    lcanopy%fes => canopy%fes(is:ie)
    lcanopy%fes_cor => canopy%fes_cor(is:ie)
    lcanopy%fevc => canopy%fevc(is:ie)
    lcanopy%ofes => canopy%ofes(is:ie)
  
    !lcanopy%A_sh => canopy%A_sh(is:ie)
    !lcanopy%A_sl => canopy%A_sl(is:ie)
    lcanopy%A_slC => canopy%A_slC(is:ie)
    lcanopy%A_shC => canopy%A_shC(is:ie)
    lcanopy%A_slJ => canopy%A_slJ(is:ie)
    lcanopy%A_shJ => canopy%A_shJ(is:ie)
    !lcanopy%eta_A_cs => canopy%eta_A_cs(is:ie)
    !lcanopy%dAdcs => canopy%dAdcs(is:ie)
    !lcanopy%cs => canopy%cs(is:ie)
    lcanopy%cs_sl => canopy%cs_sl(is:ie)
    lcanopy%cs_sh => canopy%cs_sh(is:ie)
    lcanopy%tlf => canopy%tlf(is:ie)
    lcanopy%dlf => canopy%dlf(is:ie)

    lcanopy%gw => canopy%gw(is:ie,:)
    lcanopy%ancj => canopy%ancj(is:ie,:,:)
    lcanopy%tlfy => canopy%tlfy(is:ie,:)
    lcanopy%ecy => canopy%ecy(is:ie,:)
    lcanopy%ecx => canopy%ecx(is:ie,:)
    lcanopy%ci => canopy%ci(is:ie,:,:)
    lcanopy%fwsoil => canopy%fwsoil(is:ie)
    lcanopy%kthLitt => canopy%kthLitt(is:ie)
    lcanopy%DvLitt => canopy%DvLitt(is:ie)
    
    lcanopy%fns_cor => canopy%fns_cor(is:ie)
    lcanopy%ga_cor => canopy%ga_cor(is:ie)
    
  end if  

end subroutine setp_canopy

subroutine setp_casabal(casabal,lcasabal,tile)
  implicit none

  type(casa_balance), intent(in) :: casabal
  type(casa_balance), intent(inout) :: lcasabal
  integer, intent(in) :: tile
  integer :: is, ie

  is = tdata(tile)%toffset + 1
  ie = tdata(tile)%toffset + tdata(tile)%mp

  if ( is<=ie ) then
  
    lcasabal%FCgppyear => casabal%FCgppyear(is:ie)
    lcasabal%FCnppyear => casabal%FCnppyear(is:ie)
    lcasabal%FCrmleafyear => casabal%FCrmleafyear(is:ie)
    lcasabal%FCrmwoodyear => casabal%FCrmwoodyear(is:ie)
    lcasabal%FCrmrootyear => casabal%FCrmrootyear(is:ie)
    lcasabal%FCrgrowyear => casabal%FCrgrowyear(is:ie)
    lcasabal%FCrpyear => casabal%FCrpyear(is:ie)
    lcasabal%FCrsyear => casabal%FCrsyear(is:ie)
    lcasabal%FCneeyear => casabal%FCneeyear(is:ie)
    lcasabal%dCdtyear => casabal%dCdtyear(is:ie)
    lcasabal%LAImax => casabal%LAImax(is:ie)
    lcasabal%Cleafmean => casabal%Cleafmean(is:ie)
    lcasabal%Crootmean => casabal%Crootmean(is:ie)
    lcasabal%FNdepyear => casabal%FNdepyear(is:ie)
    lcasabal%FNfixyear => casabal%FNfixyear(is:ie)
    lcasabal%FNsnetyear => casabal%FNsnetyear(is:ie)
    lcasabal%FNupyear => casabal%FNupyear(is:ie)
    lcasabal%FNleachyear => casabal%FNleachyear(is:ie)
    lcasabal%FNlossyear => casabal%FNlossyear(is:ie)
    lcasabal%FPweayear => casabal%FPweayear(is:ie)
    lcasabal%FPdustyear => casabal%FPdustyear(is:ie)
    lcasabal%FPsnetyear => casabal%FPsnetyear(is:ie)
    lcasabal%FPupyear => casabal%FPupyear(is:ie)
    lcasabal%FPleachyear => casabal%FPleachyear(is:ie)
    lcasabal%FPlossyear => casabal%FPlossyear(is:ie)

    lcasabal%glaimon => casabal%glaimon(is:ie,:)
    lcasabal%glaimonx => casabal%glaimonx(is:ie,:)
    lcasabal%cplantlast => casabal%cplantlast(is:ie,:)
    lcasabal%nplantlast => casabal%nplantlast(is:ie,:)
    lcasabal%pplantlast => casabal%pplantlast(is:ie,:)
    lcasabal%clitterlast => casabal%clitterlast(is:ie,:)
    lcasabal%nlitterlast => casabal%nlitterlast(is:ie,:)
    lcasabal%plitterlast => casabal%plitterlast(is:ie,:)
    lcasabal%csoillast => casabal%csoillast(is:ie,:)
    lcasabal%nsoillast => casabal%nsoillast(is:ie,:)
    lcasabal%psoillast => casabal%psoillast(is:ie,:)

    lcasabal%nsoilminlast => casabal%nsoilminlast(is:ie)
    lcasabal%psoillablast => casabal%psoillablast(is:ie)
    lcasabal%psoilsorblast => casabal%psoilsorblast(is:ie)
    lcasabal%psoilocclast => casabal%psoilocclast(is:ie)
    lcasabal%cbalance => casabal%cbalance(is:ie)
    lcasabal%nbalance => casabal%nbalance(is:ie)
    lcasabal%pbalance => casabal%pbalance(is:ie)
    lcasabal%sumcbal => casabal%sumcbal(is:ie)
    lcasabal%sumnbal => casabal%sumnbal(is:ie)
    lcasabal%sumpbal => casabal%sumpbal(is:ie)

    lcasabal%clabilelast => casabal%clabilelast(is:ie)
    
  end if  

end subroutine setp_casabal

subroutine setp_casaflux(casaflux,lcasaflux,tile)
  implicit none

  type(casa_flux), intent(in) :: casaflux
  type(casa_flux), intent(inout) :: lcasaflux
  integer, intent(in) :: tile
  integer :: is, ie

  is = tdata(tile)%toffset + 1
  ie = tdata(tile)%toffset + tdata(tile)%mp

  if ( is<=ie ) then
  
    lcasaflux%Cgpp => casaflux%Cgpp(is:ie)
    lcasaflux%Cnpp => casaflux%Cnpp(is:ie)
    lcasaflux%Crp => casaflux%Crp(is:ie)
    lcasaflux%Crgplant => casaflux%Crgplant(is:ie)
    lcasaflux%Nminfix => casaflux%Nminfix(is:ie)
    lcasaflux%Nminuptake => casaflux%Nminuptake(is:ie)
    lcasaflux%Plabuptake => casaflux%Plabuptake(is:ie)
    lcasaflux%Clabloss => casaflux%Clabloss(is:ie)
    lcasaflux%fracClabile => casaflux%fracClabile(is:ie)
    lcasaflux%stemnpp => casaflux%stemnpp(is:ie)
    lcasaflux%frac_sapwood => casaflux%frac_sapwood(is:ie)
    lcasaflux%sapwood_area => casaflux%sapwood_area(is:ie)

    lcasaflux%fracCalloc => casaflux%fracCalloc(is:ie,:)
    lcasaflux%fracNalloc => casaflux%fracNalloc(is:ie,:)
    lcasaflux%fracPalloc => casaflux%fracPalloc(is:ie,:)
    lcasaflux%Crmplant => casaflux%Crmplant(is:ie,:)
    lcasaflux%kplant => casaflux%kplant(is:ie,:)
    lcasaflux%Cplant_turnover => casaflux%Cplant_turnover(is:ie,:)

    lcasaflux%fromPtoL => casaflux%fromPtoL(is:ie,:,:)

    lcasaflux%Cnep => casaflux%Cnep(is:ie)
    lcasaflux%Crsoil => casaflux%Crsoil(is:ie)
    lcasaflux%Nmindep => casaflux%Nmindep(is:ie)
    lcasaflux%Nminloss => casaflux%Nminloss(is:ie)
    lcasaflux%Nminleach => casaflux%Nminleach(is:ie)
    lcasaflux%Nupland => casaflux%Nupland(is:ie)
    lcasaflux%Nlittermin => casaflux%Nlittermin(is:ie)
    lcasaflux%Nsmin => casaflux%Nsmin(is:ie)
    lcasaflux%Nsimm => casaflux%Nsimm(is:ie)
    lcasaflux%Nsnet => casaflux%Nsnet(is:ie)
    lcasaflux%fNminloss => casaflux%fNminloss(is:ie)
    lcasaflux%fNminleach => casaflux%fNminleach(is:ie)
    lcasaflux%Pdep => casaflux%Pdep(is:ie)
    lcasaflux%Pwea => casaflux%Pwea(is:ie)
    lcasaflux%Pleach => casaflux%Pleach(is:ie)
    lcasaflux%Ploss => casaflux%Ploss(is:ie)
    lcasaflux%Pupland => casaflux%Pupland(is:ie)
    lcasaflux%Plittermin => casaflux%Plittermin(is:ie)
    lcasaflux%Psmin => casaflux%Psmin(is:ie)
    lcasaflux%Psimm => casaflux%Psimm(is:ie)
    lcasaflux%Psnet => casaflux%Psnet(is:ie)
    lcasaflux%fPleach => casaflux%fPleach(is:ie)
    lcasaflux%kplab => casaflux%kplab(is:ie)
    lcasaflux%kpsorb => casaflux%kpsorb(is:ie)
    lcasaflux%kpocc => casaflux%kpocc(is:ie)
    lcasaflux%kmlabp => casaflux%kmlabp(is:ie)
    lcasaflux%Psorbmax => casaflux%Psorbmax(is:ie)
    lcasaflux%Cplant_turnover_disturbance => casaflux%Cplant_turnover_disturbance(is:ie)
    lcasaflux%Cplant_turnover_crowding => casaflux%Cplant_turnover_crowding(is:ie)
    lcasaflux%Cplant_turnover_resource_limitation => casaflux%Cplant_turnover_resource_limitation(is:ie)

    lcasaflux%klitter => casaflux%klitter(is:ie,:)
    lcasaflux%ksoil => casaflux%ksoil(is:ie,:)
    lcasaflux%fromLtoS => casaflux%fromLtoS(is:ie,:,:)
    lcasaflux%fromStoS => casaflux%fromStoS(is:ie,:,:)
    lcasaflux%fromLtoCO2 => casaflux%fromLtoCO2(is:ie,:)
    lcasaflux%fromStoCO2 => casaflux%fromStoCO2(is:ie,:)
    lcasaflux%FluxCtolitter => casaflux%FluxCtolitter(is:ie,:)
    lcasaflux%FluxNtolitter => casaflux%FluxNtolitter(is:ie,:)
    lcasaflux%FluxPtolitter => casaflux%FluxPtolitter(is:ie,:)
    lcasaflux%FluxCtosoil => casaflux%FluxCtosoil(is:ie,:)
    lcasaflux%FluxNtosoil => casaflux%FluxNtosoil(is:ie,:)
    lcasaflux%FluxPtosoil => casaflux%FluxPtosoil(is:ie,:)
    lcasaflux%FluxCtoCO2 => casaflux%FluxCtoCO2(is:ie)
    lcasaflux%FluxCtohwp => casaflux%FluxCtohwp(is:ie)
    lcasaflux%FluxNtohwp => casaflux%FluxNtohwp(is:ie)
    lcasaflux%FluxPtohwp => casaflux%FluxPtohwp(is:ie)
    lcasaflux%FluxCtoclear => casaflux%FluxCtoclear(is:ie)
    lcasaflux%FluxNtoclear => casaflux%FluxNtoclear(is:ie)
    lcasaflux%FluxPtoclear => casaflux%FluxPtoclear(is:ie)
    lcasaflux%CtransferLUC => casaflux%CtransferLUC(is:ie)

    !lcasaflux%fHarvest => casaflux%fHarvest(is:ie)
    !lcasaflux%NHarvest => casaflux%NHarvest(is:ie)
    !lcasaflux%CHarvest => casaflux%CHarvest(is:ie)
    !lcasaflux%fcrop => casaflux%fcrop(is:ie)
    
  end if  
  
end subroutine setp_casaflux

subroutine setp_casamet(casamet,lcasamet,tile)
  implicit none

  type(casa_met), intent(in) :: casamet
  type(casa_met), intent(inout) :: lcasamet
  integer, intent(in) :: tile
  integer :: is, ie

  is = tdata(tile)%toffset + 1
  ie = tdata(tile)%toffset + tdata(tile)%mp

  if ( is<=ie ) then
  
    lcasamet%glai => casamet%glai(is:ie)
    lcasamet%Tairk => casamet%Tairk(is:ie)
    lcasamet%precip => casamet%precip(is:ie)
    lcasamet%tsoilavg => casamet%tsoilavg(is:ie)
    lcasamet%moistavg => casamet%moistavg(is:ie)
    lcasamet%btran => casamet%btran(is:ie)

    lcasamet%lnonwood => casamet%lnonwood(is:ie)
    lcasamet%Tsoil => casamet%Tsoil(is:ie,:)
    lcasamet%moist => casamet%moist(is:ie,:)
    lcasamet%iveg2 => casamet%iveg2(is:ie)
    lcasamet%ijgcm => casamet%ijgcm(is:ie)
    lcasamet%isorder => casamet%isorder(is:ie)
    lcasamet%lat => casamet%lat(is:ie)
    lcasamet%lon => casamet%lon(is:ie)
    lcasamet%areacell => casamet%areacell(is:ie)

    lcasamet%Tairkspin => casamet%Tairkspin(is:ie,:)
    lcasamet%cgppspin => casamet%cgppspin(is:ie,:)
    lcasamet%crmplantspin_1 => casamet%crmplantspin_1(is:ie,:)
    lcasamet%crmplantspin_2 => casamet%crmplantspin_2(is:ie,:)
    lcasamet%crmplantspin_3 => casamet%crmplantspin_3(is:ie,:)
    lcasamet%Tsoilspin_1 => casamet%Tsoilspin_1(is:ie,:)
    lcasamet%Tsoilspin_2 => casamet%Tsoilspin_2(is:ie,:)
    lcasamet%Tsoilspin_3 => casamet%Tsoilspin_3(is:ie,:)
    lcasamet%Tsoilspin_4 => casamet%Tsoilspin_4(is:ie,:)
    lcasamet%Tsoilspin_5 => casamet%Tsoilspin_5(is:ie,:)
    lcasamet%Tsoilspin_6 => casamet%Tsoilspin_6(is:ie,:)
    lcasamet%moistspin_1 => casamet%moistspin_1(is:ie,:)
    lcasamet%moistspin_2 => casamet%moistspin_2(is:ie,:)
    lcasamet%moistspin_3 => casamet%moistspin_3(is:ie,:)
    lcasamet%moistspin_4 => casamet%moistspin_4(is:ie,:)
    lcasamet%moistspin_5 => casamet%moistspin_5(is:ie,:)
    lcasamet%moistspin_6 => casamet%moistspin_6(is:ie,:)
    lcasamet%mtempspin => casamet%mtempspin(is:ie,:)
    
  end if  

end subroutine setp_casamet

subroutine setp_casapool(casapool,lcasapool,tile)
  implicit none

  type(casa_pool), intent(in) :: casapool
  type(casa_pool), intent(inout) :: lcasapool
  integer, intent(in) :: tile
  integer :: is, ie

  is = tdata(tile)%toffset + 1
  ie = tdata(tile)%toffset + tdata(tile)%mp

  if ( is<=ie ) then
  
    lcasapool%Clabile => casapool%Clabile(is:ie)
    lcasapool%dClabiledt => casapool%dClabiledt(is:ie)
    lcasapool%Ctot => casapool%Ctot(is:ie)
    lcasapool%Ctot_0 => casapool%Ctot_0(is:ie)

    lcasapool%Cplant => casapool%Cplant(is:ie,:)
    lcasapool%Nplant => casapool%Nplant(is:ie,:)
    lcasapool%Pplant => casapool%Pplant(is:ie,:)
    lcasapool%dCplantdt => casapool%dCplantdt(is:ie,:)
    lcasapool%dNplantdt => casapool%dNplantdt(is:ie,:)
    lcasapool%dPplantdt => casapool%dPplantdt(is:ie,:)
    lcasapool%ratioNCplant => casapool%ratioNCplant(is:ie,:)
    lcasapool%ratioNPplant => casapool%ratioNPplant(is:ie,:)

    lcasapool%Nsoilmin => casapool%Nsoilmin(is:ie)
    lcasapool%Psoillab => casapool%Psoillab(is:ie)
    lcasapool%Psoilsorb => casapool%Psoilsorb(is:ie)
    lcasapool%Psoilocc => casapool%Psoilocc(is:ie)
    lcasapool%dNsoilmindt => casapool%dNsoilmindt(is:ie)
    lcasapool%dPsoillabdt => casapool%dPsoillabdt(is:ie)
    lcasapool%dPsoilsorbdt => casapool%dPsoilsorbdt(is:ie)
    lcasapool%dPsoiloccdt => casapool%dPsoiloccdt(is:ie)

    lcasapool%Clitter => casapool%Clitter(is:ie,:)
    lcasapool%Nlitter => casapool%Nlitter(is:ie,:)
    lcasapool%Plitter => casapool%Plitter(is:ie,:)
    lcasapool%dClitterdt => casapool%dClitterdt(is:ie,:)
    lcasapool%dNlitterdt => casapool%dNlitterdt(is:ie,:)
    lcasapool%dPlitterdt => casapool%dPlitterdt(is:ie,:)
    lcasapool%ratioNClitter => casapool%ratioNClitter(is:ie,:)
    lcasapool%ratioNPlitter => casapool%ratioNPlitter(is:ie,:)

    lcasapool%Csoil => casapool%Csoil(is:ie,:)
    lcasapool%Nsoil => casapool%Nsoil(is:ie,:)
    lcasapool%Psoil => casapool%Psoil(is:ie,:)
    lcasapool%dCsoildt => casapool%dCsoildt(is:ie,:)
    lcasapool%dNsoildt => casapool%dNsoildt(is:ie,:)
    lcasapool%dPsoildt => casapool%dPsoildt(is:ie,:)
    lcasapool%ratioNCsoil => casapool%ratioNCsoil(is:ie,:)
    lcasapool%ratioNCsoilnew => casapool%ratioNCsoilnew(is:ie,:)
    lcasapool%ratioNPsoil => casapool%ratioNPsoil(is:ie,:)
    lcasapool%ratioNCsoilmin => casapool%ratioNCsoilmin(is:ie,:)
    lcasapool%ratioNCsoilmax => casapool%ratioNCsoilmax(is:ie,:)
    lcasapool%ratioPcsoil => casapool%ratioPcsoil(is:ie,:)
    lcasapool%ratioPcplant => casapool%ratioPcplant(is:ie,:)
    lcasapool%ratioPclitter => casapool%ratioPclitter(is:ie,:)
    
  end if  

end subroutine setp_casapool

subroutine setp_climate(climate,lclimate,tile)
  implicit none

  type(climate_type), intent(in) :: climate
  type(climate_type), intent(inout) :: lclimate
  integer, intent(in) :: tile
  integer :: is, ie

  is = tdata(tile)%toffset + 1
  ie = tdata(tile)%toffset + tdata(tile)%mp
  
  lclimate%nyear_average = climate%nyear_average
  lclimate%nday_average = climate%nday_average
  lclimate%nyears = climate%nyears
  lclimate%doy = climate%doy

  if ( is<=ie ) then
  
    lclimate%chilldays => climate%chilldays(is:ie)
    lclimate%iveg => climate%iveg(is:ie)
    lclimate%biome => climate%biome(is:ie)

    lclimate%dtemp => climate%dtemp(is:ie)
    lclimate%dmoist => climate%dmoist(is:ie)
    lclimate%mtemp => climate%mtemp(is:ie)
    lclimate%qtemp => climate%qtemp(is:ie)
    lclimate%mmoist => climate%mmoist(is:ie)
    lclimate%mtemp_min => climate%mtemp_min(is:ie)
    lclimate%mtemp_max => climate%mtemp_max(is:ie)
    lclimate%qtemp_max => climate%qtemp_max(is:ie)
    lclimate%qtemp_max_last_year => climate%qtemp_max_last_year(is:ie)
    lclimate%mtemp_min20 => climate%mtemp_min20(is:ie)
    lclimate%mtemp_max20 => climate%mtemp_max20(is:ie)
    lclimate%atemp_mean => climate%atemp_mean(is:ie)
    lclimate%dmoist_min => climate%dmoist_min(is:ie)
    lclimate%dmoist_max => climate%dmoist_max(is:ie)
    lclimate%dmoist_min20 => climate%dmoist_min20(is:ie)
    lclimate%dmoist_max20 => climate%dmoist_max20(is:ie)
    lclimate%AGDD5 => climate%AGDD5(is:ie)
    lclimate%GDD5 => climate%GDD5(is:ie)
    lclimate%AGDD0 => climate%AGDD0(is:ie)
    lclimate%GDD0 => climate%GDD0(is:ie)
    lclimate%gdd0_rec => climate%gdd0_rec(is:ie)
    lclimate%alpha_PT => climate%alpha_PT(is:ie)
    lclimate%evap_PT => climate%evap_PT(is:ie)
    lclimate%aevap => climate%aevap(is:ie)
    lclimate%alpha_PT20 => climate%alpha_PT20(is:ie)
    lclimate%dtemp_min => climate%dtemp_min(is:ie)
    lclimate%fdorm => climate%fdorm(is:ie)
    lclimate%frec => climate%frec(is:ie)
    lclimate%gmd => climate%gmd(is:ie)
    lclimate%fapar_ann_max => climate%fapar_ann_max(is:ie)
    lclimate%fapar_ann_max_last_year => climate%fapar_ann_max_last_year(is:ie)
  
    lclimate%mtemp_min_20 => climate%mtemp_min_20(is:ie,:)
    lclimate%mtemp_max_20 => climate%mtemp_max_20(is:ie,:)
    lclimate%dmoist_min_20 => climate%dmoist_min_20(is:ie,:)
    lclimate%dmoist_max_20 => climate%dmoist_max_20(is:ie,:)
    lclimate%dtemp_31 => climate%dtemp_31(is:ie,:)
    lclimate%dmoist_31 => climate%dmoist_31(is:ie,:)
    lclimate%alpha_PT_20 => climate%alpha_PT_20(is:ie,:)
    lclimate%dtemp_91 => climate%dtemp_91(is:ie,:)
  
    lclimate%APAR_leaf_sun => climate%APAR_leaf_sun(is:ie,:)
    lclimate%APAR_leaf_shade => climate%APAR_leaf_shade(is:ie,:)
    lclimate%Dleaf_sun => climate%Dleaf_sun(is:ie,:)
    lclimate%fwsoil => climate%fwsoil(is:ie,:)
    lclimate%Dleaf_shade => climate%Dleaf_shade(is:ie,:)
    lclimate%Tleaf_sun => climate%Tleaf_sun(is:ie,:)
    lclimate%Tleaf_shade => climate%Tleaf_shade(is:ie,:)
    lclimate%cs_sun => climate%cs_sun(is:ie,:)
    lclimate%cs_shade => climate%cs_shade(is:ie,:)
    lclimate%scalex_sun => climate%scalex_sun(is:ie,:)
    lclimate%scalex_shade => climate%scalex_shade(is:ie,:)
     
  end if  
  
end subroutine setp_climate

subroutine setp_met(met,lmet,tile)
  implicit none

  type(met_type), intent(in) :: met
  type(met_type), intent(inout) :: lmet
  integer, intent(in) :: tile
  integer :: is, ie

  is = tdata(tile)%toffset + 1
  ie = tdata(tile)%toffset + tdata(tile)%mp

  if ( is<=ie ) then
  
    lmet%year => met%year(is:ie)
    lmet%moy => met%moy(is:ie)

    lmet%ca => met%ca(is:ie)
    lmet%doy => met%doy(is:ie)
    lmet%hod => met%hod(is:ie)
    lmet%ofsd => met%ofsd(is:ie)
    lmet%fld => met%fld(is:ie)
    lmet%precip => met%precip(is:ie)
    lmet%precip_sn => met%precip_sn(is:ie)
    lmet%tk => met%tk(is:ie)
    lmet%tvair => met%tvair(is:ie)
    lmet%tvrad => met%tvrad(is:ie)
    lmet%pmb => met%pmb(is:ie)
    lmet%ua => met%ua(is:ie)
    lmet%qv => met%qv(is:ie)
    lmet%qvair => met%qvair(is:ie)
    lmet%da => met%da(is:ie)
    lmet%dva => met%dva(is:ie)
    lmet%coszen => met%coszen(is:ie)
    lmet%Ndep => met%Ndep(is:ie)

    lmet%fsd => met%fsd(is:ie,:)
    
  end if  

end subroutine setp_met

subroutine setp_phen(phen,lphen,tile)
  implicit none

  type(phen_variable), intent(in) :: phen
  type(phen_variable), intent(inout) :: lphen
  integer, intent(in) :: tile
  integer :: is, ie

  is = tdata(tile)%toffset + 1
  ie = tdata(tile)%toffset + tdata(tile)%mp

  if ( is<=ie ) then
  
    lphen%phase => phen%phase(is:ie)
    lphen%TKshed => phen%TKshed
    lphen%doyphase => phen%doyphase(is:ie,:)
    lphen%phen => phen%phen(is:ie)
    lphen%aphen => phen%aphen(is:ie)
    lphen%phasespin => phen%phasespin(is:ie,:)
    lphen%doyphasespin_1 => phen%doyphasespin_1(is:ie,:)
    lphen%doyphasespin_2 => phen%doyphasespin_2(is:ie,:)
    lphen%doyphasespin_3 => phen%doyphasespin_3(is:ie,:)
    lphen%doyphasespin_4 => phen%doyphasespin_4(is:ie,:)
    
  end if  

end subroutine setp_phen

subroutine setp_pop(pop,lpop,tile)
  implicit none

  type(pop_type), intent(in) :: pop
  type(pop_type), intent(inout) :: lpop
  integer, intent(in) :: tile
  integer :: is, ie

  is = tdata(tile)%poffset + 1
  ie = tdata(tile)%poffset + tdata(tile)%np

  if ( is<=ie ) then
    !lpop%pop_grid => pop%pop_grid(is:ie)
    !lpop%it_pop => pop%it_pop(is:ie)
    lpop%np = tdata(tile)%np
    !lpop%Iwood => pop%Iwood(is:ie)
  end if  

end subroutine setp_pop

subroutine setp_rad(rad,lrad,tile)
  implicit none

  type(radiation_type), intent(in) :: rad
  type(radiation_type), intent(inout) :: lrad
  integer, intent(in) :: tile
  integer :: is, ie

  is = tdata(tile)%toffset + 1
  ie = tdata(tile)%toffset + tdata(tile)%mp

  if ( is<=ie ) then
  
    lrad%transb => rad%transb(is:ie)
    lrad%albedo_T => rad%albedo_T(is:ie)
    lrad%longitude => rad%longitude(is:ie)
    lrad%workp1 => rad%workp1(is:ie)
    lrad%workp2 => rad%workp2(is:ie)
    lrad%workp3 => rad%workp3(is:ie)
    lrad%extkb => rad%extkb(is:ie)
    lrad%extkd2 => rad%extkd2(is:ie)
    lrad%extkd => rad%extkd(is:ie)
    lrad%flws => rad%flws(is:ie)
    lrad%latitude => rad%latitude(is:ie)
    lrad%lwabv => rad%lwabv(is:ie)
    lrad%qssabs => rad%qssabs(is:ie)
    lrad%transd => rad%transd(is:ie)
    lrad%trad => rad%trad(is:ie)

    lrad%fvlai => rad%fvlai(is:ie,:)
    lrad%rhocdf => rad%rhocdf(is:ie,:)
    lrad%rniso => rad%rniso(is:ie,:)
    lrad%scalex => rad%scalex(is:ie,:)
    lrad%albedo => rad%albedo(is:ie,:)
    lrad%reffdf => rad%reffdf(is:ie,:)
    lrad%reffbm => rad%reffbm(is:ie,:)
    lrad%extkbm => rad%extkbm(is:ie,:)
    lrad%extkdm => rad%extkdm(is:ie,:)
    lrad%fbeam => rad%fbeam(is:ie,:)
    lrad%cexpkbm => rad%cexpkbm(is:ie,:)
    lrad%cexpkdm => rad%cexpkdm(is:ie,:)
    lrad%rhocbm => rad%rhocbm(is:ie,:)
    lrad%gradis => rad%gradis(is:ie,:)

    lrad%qcan => rad%qcan(is:ie,:,:)
    
  end if  

end subroutine setp_rad

subroutine setp_rough(rough,lrough,tile)
  implicit none

  type(roughness_type), intent(in) :: rough
  type(roughness_type), intent(inout) :: lrough
  integer, intent(in) :: tile
  integer :: is, ie

  is = tdata(tile)%toffset + 1
  ie = tdata(tile)%toffset + tdata(tile)%mp

  if ( is<=ie ) then
  
    lrough%disp => rough%disp(is:ie)
    lrough%hruff => rough%hruff(is:ie)
    lrough%hruff_grmx => rough%hruff_grmx(is:ie)
    lrough%rt0us => rough%rt0us(is:ie)
    lrough%rt1usa => rough%rt1usa(is:ie)
    lrough%rt1usb => rough%rt1usb(is:ie)
    lrough%rt1 => rough%rt1(is:ie)
    lrough%za_uv => rough%za_uv(is:ie)
    lrough%za_tq => rough%za_tq(is:ie)
    lrough%z0m => rough%z0m(is:ie)
    lrough%zref_uv => rough%zref_uv(is:ie)
    lrough%zref_tq => rough%zref_tq(is:ie)
    lrough%zruffs => rough%zruffs(is:ie)
    lrough%z0soilsn => rough%z0soilsn(is:ie)
    lrough%z0soil => rough%z0soil(is:ie)

    lrough%coexp => rough%coexp(is:ie)

    lrough%usuh => rough%usuh(is:ie)

    lrough%term2 => rough%term2(is:ie)
    lrough%term3 => rough%term3(is:ie)
    lrough%term5 => rough%term5(is:ie)
    lrough%term6 => rough%term6(is:ie)
    lrough%term6a => rough%term6a(is:ie)
    
  end if  

end subroutine setp_rough

subroutine setp_soil(soil,lsoil,tile)
  implicit none

  type(soil_parameter_type), intent(in) :: soil
  type(soil_parameter_type), intent(inout) :: lsoil
  integer, intent(in) :: tile
  integer :: is, ie

  is = tdata(tile)%toffset + 1
  ie = tdata(tile)%toffset + tdata(tile)%mp

  if ( is<=ie ) then
  
    lsoil%isoilm => soil%isoilm(is:ie)

    lsoil%bch => soil%bch(is:ie)
    lsoil%c3 => soil%c3(is:ie)
    lsoil%clay => soil%clay(is:ie)
    lsoil%css => soil%css(is:ie)
    lsoil%hsbh => soil%hsbh(is:ie)
    lsoil%hyds => soil%hyds(is:ie)
    lsoil%i2bp3 => soil%i2bp3(is:ie)
    lsoil%ibp2 => soil%ibp2(is:ie)
    lsoil%rhosoil => soil%rhosoil(is:ie)
    lsoil%sand => soil%sand(is:ie)
    lsoil%sfc => soil%sfc(is:ie)
    lsoil%silt => soil%silt(is:ie)
    lsoil%ssat => soil%ssat(is:ie)
    lsoil%sucs => soil%sucs(is:ie)
    lsoil%swilt => soil%swilt(is:ie)
    lsoil%zse => soil%zse !ms
    lsoil%zshh => soil%zshh !ms
    lsoil%soilcol => soil%soilcol(is:ie)
    lsoil%albsoilf => soil%albsoilf(is:ie)

    lsoil%cnsd => soil%cnsd(is:ie)
    lsoil%pwb_min => soil%pwb_min(is:ie)

    lsoil%albsoil => soil%albsoil(is:ie,:)

    lsoil%nhorizons => soil%nhorizons(is:ie)
    lsoil%ishorizon => soil%ishorizon(is:ie,:)
    lsoil%clitt => soil%clitt(is:ie)
    lsoil%zeta => soil%zeta(is:ie)
    lsoil%fsatmax => soil%fsatmax(is:ie)
    lsoil%swilt_vec => soil%swilt_vec(is:ie,:)
    lsoil%ssat_vec => soil%ssat_vec(is:ie,:)
    lsoil%sfc_vec => soil%sfc_vec(is:ie,:)
    
    lsoil%heat_cap_lower_limit => soil%heat_cap_lower_limit(is:ie,:)
    
  end if  

end subroutine setp_soil

subroutine setp_ssnow(ssnow,lssnow,tile)
  implicit none

  type(soil_snow_type), intent(in) :: ssnow
  type(soil_snow_type), intent(inout) :: lssnow
  integer, intent(in) :: tile
  integer :: is, ie

  is = tdata(tile)%toffset + 1
  ie = tdata(tile)%toffset + tdata(tile)%mp

  if ( is<=ie ) then
  
    lssnow%isflag => ssnow%isflag(is:ie)

    lssnow%iantrct => ssnow%iantrct(is:ie)
    lssnow%pudsto => ssnow%pudsto(is:ie)
    lssnow%pudsmx => ssnow%pudsmx(is:ie)
    lssnow%cls => ssnow%cls(is:ie)
    lssnow%dfn_dtg => ssnow%dfn_dtg(is:ie)
    lssnow%dfh_dtg => ssnow%dfh_dtg(is:ie)
    lssnow%dfe_ddq => ssnow%dfe_ddq(is:ie)
    lssnow%ddq_dtg => ssnow%ddq_dtg(is:ie)
    lssnow%evapsn => ssnow%evapsn(is:ie)
    lssnow%fwtop => ssnow%fwtop(is:ie)
    lssnow%fwtop1 => ssnow%fwtop1(is:ie)
    lssnow%fwtop2 => ssnow%fwtop2(is:ie)
    lssnow%fwtop3 => ssnow%fwtop3(is:ie)
    lssnow%osnowd => ssnow%osnowd(is:ie)
    lssnow%potev => ssnow%potev(is:ie)
    lssnow%runoff => ssnow%runoff(is:ie)
    lssnow%rnof1 => ssnow%rnof1(is:ie)
    lssnow%rnof2 => ssnow%rnof2(is:ie)
    lssnow%rtsoil => ssnow%rtsoil(is:ie)
    lssnow%wbtot1 => ssnow%wbtot1(is:ie)
    lssnow%wbtot2 => ssnow%wbtot2(is:ie)
    lssnow%wb_lake => ssnow%wb_lake(is:ie)
    lssnow%sinfil => ssnow%sinfil(is:ie)
    lssnow%qstss => ssnow%qstss(is:ie)
    lssnow%wetfac => ssnow%wetfac(is:ie)
    lssnow%owetfac => ssnow%owetfac(is:ie)
    lssnow%t_snwlr => ssnow%t_snwlr(is:ie)
    lssnow%tggav => ssnow%tggav(is:ie)
    lssnow%otgg => ssnow%otgg(is:ie)
    lssnow%otss => ssnow%otss(is:ie)
    lssnow%otss_0 => ssnow%otss_0(is:ie)
    lssnow%tprecip => ssnow%tprecip(is:ie)
    lssnow%tevap => ssnow%tevap(is:ie)
    lssnow%trnoff => ssnow%trnoff(is:ie)
    lssnow%totenbal => ssnow%totenbal(is:ie)
    lssnow%totenbal2 => ssnow%totenbal2(is:ie)
    lssnow%fland => ssnow%fland(is:ie)
    lssnow%ifland => ssnow%ifland(is:ie)
    lssnow%qasrf => ssnow%qasrf(is:ie)
    lssnow%qfsrf => ssnow%qfsrf(is:ie)
    lssnow%qssrf => ssnow%qssrf(is:ie)
    lssnow%snage => ssnow%snage(is:ie)
    lssnow%snowd => ssnow%snowd(is:ie)
    lssnow%smelt => ssnow%smelt(is:ie)
    lssnow%ssdnn => ssnow%ssdnn(is:ie)
    lssnow%tss => ssnow%tss(is:ie)
    lssnow%tss_p => ssnow%tss_p(is:ie)
    lssnow%deltss => ssnow%deltss(is:ie)
    lssnow%owb1 => ssnow%owb1(is:ie)

    lssnow%sconds => ssnow%sconds(is:ie,:)
    lssnow%sdepth => ssnow%sdepth(is:ie,:)
    lssnow%smass => ssnow%smass(is:ie,:)
    lssnow%ssdn => ssnow%ssdn(is:ie,:)
    lssnow%tgg => ssnow%tgg(is:ie,:)
    lssnow%tggsn => ssnow%tggsn(is:ie,:)
    lssnow%dtmlt => ssnow%dtmlt(is:ie,:)
    lssnow%albsoilsn => ssnow%albsoilsn(is:ie,:)
    lssnow%evapfbl => ssnow%evapfbl(is:ie,:)
    lssnow%tilefrac => ssnow%tilefrac(is:ie,:)

    lssnow%wbtot => ssnow%wbtot(is:ie)

    lssnow%gammzz => ssnow%gammzz(is:ie,:)
    lssnow%wb => ssnow%wb(is:ie,:)
    lssnow%wbice => ssnow%wbice(is:ie,:)
    lssnow%wblf => ssnow%wblf(is:ie,:)
    lssnow%wbfice => ssnow%wbfice(is:ie,:)

    lssnow%S => ssnow%S(is:ie,:)
    lssnow%Tsoil => ssnow%Tsoil(is:ie,:)
    lssnow%SL => ssnow%SL(is:ie)
    lssnow%TL => ssnow%TL(is:ie)
    lssnow%h0 => ssnow%h0(is:ie)
    lssnow%rex => ssnow%rex(is:ie,:)
    lssnow%wflux => ssnow%wflux(is:ie,:)
    lssnow%delwcol => ssnow%delwcol(is:ie)
    lssnow%zdelta => ssnow%zdelta(is:ie)
    lssnow%kth => ssnow%kth(is:ie,:)
    lssnow%Tsurface => ssnow%Tsurface(is:ie)
    lssnow%lE => ssnow%lE(is:ie)
    lssnow%evap => ssnow%evap(is:ie)
    lssnow%ciso => ssnow%ciso(is:ie,:)
    lssnow%cisoL => ssnow%cisoL(is:ie)
    lssnow%rlitt => ssnow%rlitt(is:ie)
    lssnow%thetai => ssnow%thetai(is:ie,:)
    lssnow%snowliq => ssnow%snowliq(is:ie,:)
    lssnow%nsteps => ssnow%nsteps(is:ie)
    lssnow%TsurfaceFR => ssnow%TsurfaceFR(is:ie)
    lssnow%Ta_daily => ssnow%Ta_daily(is:ie,:)
    lssnow%nsnow => ssnow%nsnow(is:ie)
    lssnow%Qadv_daily => ssnow%Qadv_daily(is:ie)
    lssnow%G0_daily => ssnow%G0_daily(is:ie)
    lssnow%Qevap_daily => ssnow%Qevap_daily(is:ie)
    lssnow%Qprec_daily => ssnow%Qprec_daily(is:ie)
    lssnow%Qprec_snow_daily => ssnow%Qprec_snow_daily(is:ie)
    
    lssnow%GWwb => ssnow%GWwb(is:ie)
    lssnow%satfrac => ssnow%satfrac(is:ie)
    lssnow%rh_srf => ssnow%rh_srf(is:ie)
    lssnow%dfe_dtg => ssnow%dfe_dtg(is:ie)
    lssnow%wbliq => ssnow%wbliq(is:ie,:)
    
  end if  

end subroutine setp_ssnow

subroutine setp_sumflux(sum_flux,lsum_flux,tile)
  implicit none

  type(sum_flux_type), intent(in) :: sum_flux
  type(sum_flux_type), intent(inout) :: lsum_flux
  integer, intent(in) :: tile
  integer :: is, ie

  is = tdata(tile)%toffset + 1
  ie = tdata(tile)%toffset + tdata(tile)%mp

  if ( is<=ie ) then
  
    lsum_flux%sumpn => sum_flux%sumpn(is:ie)
    lsum_flux%sumrp => sum_flux%sumrp(is:ie)
    lsum_flux%sumrpw => sum_flux%sumrpw(is:ie)
    lsum_flux%sumrpr => sum_flux%sumrpr(is:ie)
    lsum_flux%sumrs => sum_flux%sumrs(is:ie)
    lsum_flux%sumrd => sum_flux%sumrd(is:ie)
    lsum_flux%dsumpn => sum_flux%dsumpn(is:ie)
    lsum_flux%dsumrp => sum_flux%dsumrp(is:ie)
    lsum_flux%dsumrs => sum_flux%dsumrs(is:ie)
    lsum_flux%dsumrd => sum_flux%dsumrd(is:ie)
    lsum_flux%sumxrp => sum_flux%sumxrp(is:ie)
    lsum_flux%sumxrs => sum_flux%sumxrs(is:ie)
    
  end if  

end subroutine setp_sumflux

subroutine setp_veg(veg,lveg,tile)
  implicit none

  type(veg_parameter_type), intent(in) :: veg
  type(veg_parameter_type), intent(inout) :: lveg
  integer, intent(in) :: tile
  integer :: is, ie

  is = tdata(tile)%toffset + 1
  ie = tdata(tile)%toffset + tdata(tile)%mp

  if ( is<=ie ) then
  
    lveg%iveg => veg%iveg(is:ie)
    lveg%iLU => veg%iLU(is:ie)
    lveg%canst1 => veg%canst1(is:ie)
    lveg%dleaf => veg%dleaf(is:ie)
    lveg%ejmax => veg%ejmax(is:ie)
    lveg%meth => veg%meth(is:ie)
    lveg%frac4 => veg%frac4(is:ie)
    lveg%hc => veg%hc(is:ie)
    lveg%vlai => veg%vlai(is:ie)
    lveg%xalbnir => veg%xalbnir(is:ie)
    lveg%rp20 => veg%rp20(is:ie)
    lveg%rpcoef => veg%rpcoef(is:ie)
    lveg%rs20 => veg%rs20(is:ie)
    lveg%shelrb => veg%shelrb(is:ie)
    lveg%vegcf => veg%vegcf(is:ie)
    lveg%tminvj => veg%tminvj(is:ie)
    lveg%toptvj => veg%toptvj(is:ie)
    lveg%tmaxvj => veg%tmaxvj(is:ie)
    lveg%vbeta => veg%vbeta(is:ie)
    lveg%vcmax => veg%vcmax(is:ie)
    lveg%xfang => veg%xfang(is:ie)
    lveg%extkn => veg%extkn(is:ie)
    lveg%vlaimax => veg%vlaimax(is:ie)
    lveg%wai => veg%wai(is:ie)
    lveg%a1gs => veg%a1gs(is:ie)
    lveg%d0gs => veg%d0gs(is:ie)
    lveg%alpha => veg%alpha(is:ie)
    lveg%convex => veg%convex(is:ie)
    lveg%cfrd => veg%cfrd(is:ie)
    lveg%gswmin => veg%gswmin(is:ie)
    lveg%conkc0 => veg%conkc0(is:ie)
    lveg%conko0 => veg%conko0(is:ie)
    lveg%ekc => veg%ekc(is:ie)
    lveg%eko => veg%eko(is:ie)
    lveg%g0 => veg%g0(is:ie)
    lveg%g1 => veg%g1(is:ie)

    lveg%deciduous => veg%deciduous(is:ie)

    lveg%refl => veg%refl(is:ie,:)
    lveg%taul => veg%taul(is:ie,:)
    lveg%froot => veg%froot(is:ie,:)

    lveg%rootbeta => veg%rootbeta(is:ie)
    lveg%gamma => veg%gamma(is:ie)
    lveg%ZR => veg%ZR(is:ie)
    lveg%F10 => veg%F10(is:ie)

    lveg%clitt => veg%clitt(is:ie)

    lveg%disturbance_interval => veg%disturbance_interval(is:ie,:)
    lveg%disturbance_intensity => veg%disturbance_intensity(is:ie,:)
  
    lveg%vcmax_shade => veg%vcmax_shade(is:ie)
    lveg%ejmax_shade => veg%ejmax_shade(is:ie)
    lveg%vcmax_sun => veg%vcmax_sun(is:ie)
    lveg%ejmax_sun => veg%ejmax_sun(is:ie)
    
  end if  

end subroutine setp_veg

subroutine setp_climate_save(climate_save,lclimate_save,tile)
  implicit none

  type(climate_save_type), intent(in) :: climate_save
  type(climate_save_type), intent(inout) :: lclimate_save
  integer, intent(in) :: tile
  integer :: is, ie

  is = tdata(tile)%toffset + 1
  ie = tdata(tile)%toffset + tdata(tile)%mp
  
  if ( is<=ie ) then
    
    lclimate_save%APAR_leaf_sun => climate_save%APAR_leaf_sun(is:ie)
    lclimate_save%APAR_leaf_shade => climate_save%APAR_leaf_shade(is:ie)
    lclimate_save%Dleaf_sun => climate_save%Dleaf_sun(is:ie)
    lclimate_save%fwsoil => climate_save%fwsoil(is:ie)
    lclimate_save%Dleaf_shade => climate_save%Dleaf_shade(is:ie)
    lclimate_save%Tleaf_sun => climate_save%Tleaf_sun(is:ie)
    lclimate_save%Tleaf_shade => climate_save%Tleaf_shade(is:ie)
    lclimate_save%cs_sun => climate_save%cs_sun(is:ie)
    lclimate_save%cs_shade => climate_save%cs_shade(is:ie)
    lclimate_save%scalex_sun => climate_save%scalex_sun(is:ie)
    lclimate_save%scalex_shade => climate_save%scalex_shade(is:ie)
    
  end if  
  
end subroutine setp_climate_save

end module cable_ccam4
