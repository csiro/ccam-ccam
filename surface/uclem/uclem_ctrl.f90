! UCLEM urban canopy model
    
! Copyright 2015-2024 Commonwealth Scientific Industrial Research Organisation (CSIRO)
    
! This file is part of the UCLEM urban canopy model
!
! UCLEM is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! UCLEM is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with UCLEM.  If not, see <http://www.gnu.org/licenses/>.

!------------------------------------------------------------------------------ 

! Usual practice is:
!   call uclem_init            ! to initalise state arrays, etc (use uclem_disable to disable calls to uclem subroutines)
!   call uclem_loadd           ! to load previous state arrays (from uclem_savem)
!   call uclem_deftype         ! to define urban type (or use uclem_fndef to define urban properties at each grid point)
!   ...
!   do t=1,tmax
!     ...
!     call uclem_newangle1     ! store current solar zenith and azimuthal angle (use uclem_ccangle for CCAM)
!     call uclem_alb1(split=1) ! returns urban albedo for direct component of shortwave radiation
!     call uclem_alb1(split=2) ! returns urban albedo for diffuse component of shortwave radiation
!     ...
!     call uclem_fbeam         ! store fraction of direct shortwave radiation (or use uclem_spitter to estimate fraction)
!     call uclem_alb1          ! returns net urban albedo (i.e., default split=0)
!     ...
!     call uclem_calc          ! calculates urban temperatures, fluxes, etc and blends with input
!     call uclem_cd            ! returns urban drag coefficient (or use uclem_zo for roughness length)
!     call uclem_scrnout       ! returns screen level diagnostics
!     ...
!   end do
!   ...
!   call uclem_saved           ! to save current state arrays (for use by tebloadm)
!   call uclem_end             ! to deallocate memory before quitting

! only uclem_init and uclem_calc are mandatory.  All other subroutine calls are optional.


! DEPENDICES:

! uclem_alb1(split=1)                   depends on     uclem_newangle1 (or uclem_ccangle)
! uclem_alb1(split=2)                   depends on     uclem_newangle1 (or uclem_ccangle)
! uclem_alb1 (i.e., default split=0)    depends on     uclem_newangle1 (or uclem_ccangle) and uclem_fbeam (or uclem_spitter)
! uclem_calc                            depends on     uclem_newangle1 (or uclem_ccangle) and uclem_fbeam (or uclem_spitter)  
! uclem_cd (or uclem_zo)                depends on     uclem_calc
! uclem_scrnout                         depends on     uclem_calc


module uclem_ctrl

#ifdef CCAM
use cc_omp, only : imax, ntiles
#endif
use uclem

implicit none

private
public uclem_init, uclem_end, uclem_calc, uclem_zo, uclem_type, uclem_alb1,            &
       uclem_newangle1, uclem_ccangle, uclem_disable, uclem_cd,                        &
       uclem_scrnout, uclem_fbeam, uclem_spitter, uclem_sigmau,                        &
       uclem_deftype, uclem_hydro, uclem_energy, uclem_misc,                           &
       uclem_deftype_export, uclem_avetemp, uclem_loadd_2, uclem_saved_2,              &
       uclem_loadd_3, uclem_saved_3
public urbtemp, refheight, soilunder

public upack_g, ufull_g, nl, nfrac
public f_roof,f_wall,f_road,f_slab,f_intm
public intm_g,rdhyd_g,rfhyd_g,rfveg_g
public road_g,roof_g,room_g,slab_g,walle_g,wallw_g,cnveg_g,intl_g
public f_g,p_g
public facetparams,facetdata,hydrodata,vegdata,intldata
public fparmdata,pdiagdata

! state arrays
integer, save :: ifull, nfrac
#ifndef CCAM
integer, save :: ntiles = 1     ! Emulate OMP
integer, save :: imax = 0       ! Emulate OMP
#endif
integer, dimension(:), allocatable, save :: ufull_g
logical, save :: uclem_active = .false.
logical, dimension(:,:), allocatable, save :: upack_g

type(facetdata), dimension(:,:), allocatable, save :: roof_g, road_g, walle_g, wallw_g, slab_g, intm_g, room_g
type(hydrodata), dimension(:,:), allocatable, save :: rfhyd_g, rdhyd_g
type(pdiagdata), dimension(:,:), allocatable, save :: p_g
type(facetparams), dimension(:), allocatable, save :: f_roof, f_road, f_wall, f_slab, f_intm
type(vegdata), dimension(:), allocatable,     save :: cnveg_g, rfveg_g
type(intldata), dimension(:), allocatable,    save :: intl_g
type(fparmdata), dimension(:), allocatable,   save :: f_g

integer, save      :: soilunder=1          ! Modify road heat capacity to extend under
                                           ! (0=road only, 1=canveg, 2=bld, 3=canveg & bld)

interface uclem_calc
  module procedure uclem_calc_standard, uclem_calc_thread
  module procedure uclem_calc2_standard, uclem_calc2_thread
end interface

interface uclem_deftype_export
  module procedure uclem_deftype_export_standard, uclem_deftype_export_thread
end interface

interface uclem_energy
  module procedure uclem_energy_standard, uclem_energy_thread
end interface

interface uclem_zo
  module procedure uclem_zo_standard, uclem_zo_thread
end interface

interface uclem_cd
  module procedure uclem_cd_standard, uclem_cd_thread
end interface

interface uclem_type
  module procedure uclem_type_standard, uclem_type_thread
end interface

interface uclem_hydro
  module procedure uclem_hydro_standard, uclem_hydro_thread
end interface

interface uclem_misc
  module procedure uclem_misc_standard, uclem_misc_thread
end interface

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine prepare the arrays used by the UCLEM scheme
! This is a compulsory subroutine that must be called during
! model initialisation

subroutine uclem_init(ifin,sigu,diag)

implicit none

integer, intent(in) :: ifin,diag
integer, dimension(:), allocatable, save :: utype
integer tile, is, ie, ifrac
real, dimension(ifin), intent(in) :: sigu

if ( diag>=1 ) write(6,*) "Initialising UCLEM"

uclem_active = .true.

ifull = ifin

if ( ntiles<1 ) then
  write(6,*) "ERROR: Invalid ntiles ",ntiles
  stop
end if

#ifndef CCAM
imax = ifull/ntiles
if ( mod(ifull,ntiles)/=0 ) then
  write(6,*) "ERROR: Invalid ntiles ",ntiles," for ifull ",ifull
  stop
end if
#endif

allocate( f_roof(ntiles), f_road(ntiles), f_wall(ntiles), f_slab(ntiles), f_intm(ntiles) )
allocate( cnveg_g(ntiles), rfveg_g(ntiles) )
allocate( f_g(ntiles) )
allocate( intl_g(ntiles) )
allocate( ufull_g(ntiles) )
allocate( upack_g(imax,ntiles) )

allocate( utype(imax) )

do tile = 1,ntiles
  is = (tile-1)*imax + 1
  ie = tile*imax

  upack_g(1:imax,tile) = sigu(is:ie)>0.
  ufull_g(tile) = count( upack_g(1:imax,tile) )

  allocate(f_g(tile)%rfvegdepth(ufull_g(tile)))
  allocate(f_g(tile)%ctime(ufull_g(tile)),f_g(tile)%bldairtemp(ufull_g(tile)),f_g(tile)%weekdayload(ufull_g(tile)))
  allocate(f_g(tile)%hangle(ufull_g(tile)),f_g(tile)%vangle(ufull_g(tile)),f_g(tile)%fbeam(ufull_g(tile)))
  allocate(f_g(tile)%hwratio(ufull_g(tile)))
  allocate(f_g(tile)%effhwratio(ufull_g(tile)),f_g(tile)%bldheight(ufull_g(tile)))
  allocate(f_g(tile)%sigmabld(ufull_g(tile)),f_g(tile)%industryfg(ufull_g(tile)))
  allocate(f_g(tile)%intgains_flr(ufull_g(tile)),f_g(tile)%trafficfg(ufull_g(tile)))
  allocate(f_g(tile)%swilt(ufull_g(tile)),f_g(tile)%sfc(ufull_g(tile)),f_g(tile)%ssat(ufull_g(tile)))
  allocate(f_g(tile)%heatprop(ufull_g(tile)),f_g(tile)%coolprop(ufull_g(tile)))
  allocate(f_g(tile)%intmassn(ufull_g(tile)),f_g(tile)%bldwidth(ufull_g(tile)))
  allocate(f_g(tile)%infilach(ufull_g(tile)),f_g(tile)%ventilach(ufull_g(tile)))
  allocate(f_g(tile)%sigmau(ufull_g(tile)))
  
  allocate(f_roof(tile)%depth(ufull_g(tile),nl),f_roof(tile)%lambda(ufull_g(tile),nl))
  allocate(f_roof(tile)%volcp(ufull_g(tile),nl))
  allocate(f_wall(tile)%depth(ufull_g(tile),nl),f_wall(tile)%lambda(ufull_g(tile),nl))
  allocate(f_wall(tile)%volcp(ufull_g(tile),nl))
  allocate(f_road(tile)%depth(ufull_g(tile),nl),f_road(tile)%lambda(ufull_g(tile),nl))
  allocate(f_road(tile)%volcp(ufull_g(tile),nl))
  allocate(f_slab(tile)%depth(ufull_g(tile),nl),f_slab(tile)%lambda(ufull_g(tile),nl))
  allocate(f_slab(tile)%volcp(ufull_g(tile),nl))
  allocate(f_intm(tile)%depth(ufull_g(tile),nl),f_intm(tile)%lambda(ufull_g(tile),nl))
  allocate(f_intm(tile)%volcp(ufull_g(tile),nl))
  allocate(f_roof(tile)%emiss(ufull_g(tile)),f_roof(tile)%alpha(ufull_g(tile)))
  allocate(f_wall(tile)%emiss(ufull_g(tile)),f_wall(tile)%alpha(ufull_g(tile)))
  allocate(f_road(tile)%emiss(ufull_g(tile)),f_road(tile)%alpha(ufull_g(tile)))
  allocate(f_slab(tile)%emiss(ufull_g(tile)))
 
  ! veg temperature is recalculated each time-step and hence does not need to be stored for each ifrac
  allocate(cnveg_g(tile)%emiss(ufull_g(tile)),cnveg_g(tile)%sigma(ufull_g(tile)))
  allocate(cnveg_g(tile)%alpha(ufull_g(tile)),cnveg_g(tile)%rsmin(ufull_g(tile)))
  allocate(cnveg_g(tile)%zo(ufull_g(tile)),cnveg_g(tile)%lai(ufull_g(tile)))
  allocate(rfveg_g(tile)%emiss(ufull_g(tile)),rfveg_g(tile)%sigma(ufull_g(tile)))
  allocate(rfveg_g(tile)%alpha(ufull_g(tile)),rfveg_g(tile)%rsmin(ufull_g(tile)))
  allocate(rfveg_g(tile)%zo(ufull_g(tile)),rfveg_g(tile)%lai(ufull_g(tile)))
  allocate(rfveg_g(tile)%temp(ufull_g(tile)),cnveg_g(tile)%temp(ufull_g(tile)))
  
  allocate(intl_g(tile)%viewf(ufull_g(tile),4,4),intl_g(tile)%psi(ufull_g(tile),4,4))
  
  if ( ufull_g(tile)>0 ) then
      
    ! define grid arrays
    f_g(tile)%sigmau = pack(sigu(is:ie),upack_g(1:imax,tile))
    f_g(tile)%rfvegdepth=0.1
    f_g(tile)%hwratio=1.
    f_g(tile)%sigmabld=0.5
    f_g(tile)%industryfg=0.
    f_g(tile)%intgains_flr=0.
    f_g(tile)%trafficfg=0.
    f_g(tile)%bldheight=10.
    f_g(tile)%bldairtemp=1. ! + urbtemp
    f_g(tile)%vangle=0.
    f_g(tile)%hangle=0.
    f_g(tile)%ctime=0.
    f_g(tile)%fbeam=1.
    f_g(tile)%swilt=0.
    f_g(tile)%sfc=0.5
    f_g(tile)%ssat=1.
    f_g(tile)%infilach=0.5
    f_g(tile)%ventilach=2.
    f_g(tile)%weekdayload=1.0
    
    f_roof(tile)%depth=0.1
    f_roof(tile)%volcp=2.E6
    f_roof(tile)%lambda=2.
    f_roof(tile)%alpha=0.2
    f_roof(tile)%emiss=0.97
    f_wall(tile)%depth=0.1
    f_wall(tile)%volcp=2.E6
    f_wall(tile)%lambda=2.
    f_wall(tile)%alpha=0.2
    f_wall(tile)%emiss=0.97
    f_road(tile)%depth=0.1
    f_road(tile)%volcp=2.E6
    f_road(tile)%lambda=2.
    f_road(tile)%alpha=0.2
    f_road(tile)%emiss=0.97
    f_slab(tile)%depth=0.1
    f_slab(tile)%volcp=2.E6
    f_slab(tile)%lambda=2.
    f_slab(tile)%emiss=0.97
    f_intm(tile)%depth=0.1
    f_intm(tile)%lambda=2.
    f_intm(tile)%volcp=2.E6
    
    cnveg_g(tile)%sigma=0.5
    cnveg_g(tile)%alpha=0.2
    cnveg_g(tile)%emiss=0.97
    cnveg_g(tile)%zo=0.1
    cnveg_g(tile)%lai=1.
    cnveg_g(tile)%rsmin=200.
    cnveg_g(tile)%temp=1. ! + urbtemp             ! updated in uclem_calc
    rfveg_g(tile)%sigma=0.
    rfveg_g(tile)%alpha=0.2
    rfveg_g(tile)%emiss=0.97
    rfveg_g(tile)%zo=0.1
    rfveg_g(tile)%lai=1.
    rfveg_g(tile)%rsmin=200.
    rfveg_g(tile)%temp=1. ! + urbtemp             ! updated in uclem_calc
    
    utype=1 ! default urban
    call uclem_type(utype,diag,f_g(tile),cnveg_g(tile),rfveg_g(tile),      &
                    f_roof(tile),f_road(tile),f_wall(tile),f_slab(tile),   &
                    f_intm(tile),intl_g(tile),upack_g(:,tile),             &
                    ufull_g(tile))

  end if
    
end do

deallocate( utype )


! moved here so uclem_type can read namelist
select case( intairtmeth )
  case(0,1)
    nfrac = 1
  case(2)
    nfrac = 3
  case default
    write(6,*) "ERROR: Unknown intairtmeth option ",intairtmeth
    stop
end select

  
allocate( roof_g(nfrac,ntiles), road_g(nfrac,ntiles), walle_g(nfrac,ntiles), wallw_g(nfrac,ntiles) )
allocate( slab_g(nfrac,ntiles), intm_g(nfrac,ntiles) )
allocate( room_g(nfrac,ntiles) )
allocate( rfhyd_g(nfrac,ntiles), rdhyd_g(nfrac,ntiles) )
allocate( p_g(nfrac,ntiles) )

do tile = 1,ntiles
  is = (tile-1)*imax + 1
  ie = tile*imax

  upack_g(1:imax,tile) = sigu(is:ie)>0.
  ufull_g(tile) = count( upack_g(1:imax,tile) )
    
  do ifrac = 1,nfrac
  
    allocate(roof_g(ifrac,tile)%nodetemp(ufull_g(tile),0:nl),road_g(ifrac,tile)%nodetemp(ufull_g(tile),0:nl))
    allocate(walle_g(ifrac,tile)%nodetemp(ufull_g(tile),0:nl),wallw_g(ifrac,tile)%nodetemp(ufull_g(tile),0:nl))
    allocate(slab_g(ifrac,tile)%nodetemp(ufull_g(tile),0:nl),intm_g(ifrac,tile)%nodetemp(ufull_g(tile),0:nl))
    allocate(room_g(ifrac,tile)%nodetemp(ufull_g(tile),1))
    allocate(road_g(ifrac,tile)%storage(ufull_g(tile),nl),roof_g(ifrac,tile)%storage(ufull_g(tile),nl))
    allocate(walle_g(ifrac,tile)%storage(ufull_g(tile),nl),wallw_g(ifrac,tile)%storage(ufull_g(tile),nl))
    allocate(slab_g(ifrac,tile)%storage(ufull_g(tile),nl),intm_g(ifrac,tile)%storage(ufull_g(tile),nl))
    allocate(room_g(ifrac,tile)%storage(ufull_g(tile),1)) 
  
    allocate(rfhyd_g(ifrac,tile)%surfwater(ufull_g(tile)),rfhyd_g(ifrac,tile)%snow(ufull_g(tile)))
    allocate(rfhyd_g(ifrac,tile)%den(ufull_g(tile)),rfhyd_g(ifrac,tile)%snowalpha(ufull_g(tile)))
    allocate(rfhyd_g(ifrac,tile)%leafwater(ufull_g(tile)),rfhyd_g(ifrac,tile)%soilwater(ufull_g(tile)))
    allocate(rdhyd_g(ifrac,tile)%surfwater(ufull_g(tile)),rdhyd_g(ifrac,tile)%snow(ufull_g(tile)))
    allocate(rdhyd_g(ifrac,tile)%den(ufull_g(tile)),rdhyd_g(ifrac,tile)%snowalpha(ufull_g(tile)))
    allocate(rdhyd_g(ifrac,tile)%leafwater(ufull_g(tile)),rdhyd_g(ifrac,tile)%soilwater(ufull_g(tile)))

    allocate(p_g(ifrac,tile)%lzom(ufull_g(tile)),p_g(ifrac,tile)%lzoh(ufull_g(tile)),p_g(ifrac,tile)%cndzmin(ufull_g(tile)))
    allocate(p_g(ifrac,tile)%cduv(ufull_g(tile)),p_g(ifrac,tile)%cdtq(ufull_g(tile)))
    allocate(p_g(ifrac,tile)%tscrn(ufull_g(tile)),p_g(ifrac,tile)%qscrn(ufull_g(tile)),p_g(ifrac,tile)%uscrn(ufull_g(tile)))
    allocate(p_g(ifrac,tile)%u10(ufull_g(tile)),p_g(ifrac,tile)%emiss(ufull_g(tile)))
    allocate(p_g(ifrac,tile)%bldheat(ufull_g(tile)),p_g(ifrac,tile)%bldcool(ufull_g(tile)),p_g(ifrac,tile)%traf(ufull_g(tile)))
    allocate(p_g(ifrac,tile)%surfrunoff(ufull_g(tile)),p_g(ifrac,tile)%irrig(ufull_g(tile)))
    allocate(p_g(ifrac,tile)%transveg(ufull_g(tile)), p_g(ifrac,tile)%acond_vegw(ufull_g(tile)))
    allocate(p_g(ifrac,tile)%soilmoist(ufull_g(tile)),p_g(ifrac,tile)%delsoilmoist(ufull_g(tile)))
    allocate(p_g(ifrac,tile)%rootmoistc(ufull_g(tile)), p_g(ifrac,tile)%ulai(ufull_g(tile)))
    allocate(p_g(ifrac,tile)%intgains_full(ufull_g(tile)),p_g(ifrac,tile)%storage_flux(ufull_g(tile)))
    allocate(p_g(ifrac,tile)%surferr(ufull_g(tile)),p_g(ifrac,tile)%atmoserr(ufull_g(tile)))
    allocate(p_g(ifrac,tile)%surferr_bias(ufull_g(tile)),p_g(ifrac,tile)%atmoserr_bias(ufull_g(tile)))
    allocate(p_g(ifrac,tile)%soilwetness(ufull_g(tile)),p_g(ifrac,tile)%soilwater(ufull_g(tile)))
    allocate(p_g(ifrac,tile)%snowmelt(ufull_g(tile)),p_g(ifrac,tile)%frac_sigma(ufull_g(tile)))
    allocate(p_g(ifrac,tile)%delintercept(ufull_g(tile)),p_g(ifrac,tile)%snowt(ufull_g(tile)))
    allocate(p_g(ifrac,tile)%vegt(ufull_g(tile)),p_g(ifrac,tile)%swe(ufull_g(tile)))
    allocate(p_g(ifrac,tile)%surfstor(ufull_g(tile)),p_g(ifrac,tile)%snowfrac(ufull_g(tile)))
    allocate(p_g(ifrac,tile)%salbedo(ufull_g(tile)),p_g(ifrac,tile)%calbedo(ufull_g(tile)))
    allocate(p_g(ifrac,tile)%taircanyon(ufull_g(tile)),p_g(ifrac,tile)%delswe(ufull_g(tile)))
    allocate(p_g(ifrac,tile)%roof_water_runoff(ufull_g(tile)),p_g(ifrac,tile)%roof_snow_runoff(ufull_g(tile)))
    allocate(p_g(ifrac,tile)%roof_soil_runoff(ufull_g(tile)),p_g(ifrac,tile)%road_water_runoff(ufull_g(tile)))
    allocate(p_g(ifrac,tile)%road_snow_runoff(ufull_g(tile)),p_g(ifrac,tile)%road_soil_runoff(ufull_g(tile)))

    if ( ufull_g(tile)>0 ) then

      ! Initialise state variables
      roof_g(ifrac,tile)%nodetemp=1._8  ! + urbtemp
      roof_g(ifrac,tile)%storage =0._8
      road_g(ifrac,tile)%nodetemp=1._8  ! + urbtemp
      road_g(ifrac,tile)%storage =0._8
      walle_g(ifrac,tile)%nodetemp=1._8 ! + urbtemp
      walle_g(ifrac,tile)%storage=0._8
      wallw_g(ifrac,tile)%nodetemp=1._8 ! + urbtemp
      wallw_g(ifrac,tile)%storage=0._8
      slab_g(ifrac,tile)%nodetemp=1._8 ! + urbtemp
      slab_g(ifrac,tile)%storage=0._8
      intm_g(ifrac,tile)%nodetemp=1._8 ! + urbtemp
      intm_g(ifrac,tile)%storage=0._8
      room_g(ifrac,tile)%nodetemp=1._8  ! + urbtemp
      room_g(ifrac,tile)%storage=0._8

      rfhyd_g(ifrac,tile)%surfwater=0.
      rfhyd_g(ifrac,tile)%snow=0.
      rfhyd_g(ifrac,tile)%den=minsnowden
      rfhyd_g(ifrac,tile)%snowalpha=maxsnowalpha
      rfhyd_g(ifrac,tile)%leafwater=0.
      rfhyd_g(ifrac,tile)%soilwater=0.
      rdhyd_g(ifrac,tile)%surfwater=0.
      rdhyd_g(ifrac,tile)%snow=0.
      rdhyd_g(ifrac,tile)%den=minsnowden
      rdhyd_g(ifrac,tile)%snowalpha=maxsnowalpha
      rdhyd_g(ifrac,tile)%leafwater=0.
      rdhyd_g(ifrac,tile)%soilwater=0.25
    
      p_g(ifrac,tile)%cndzmin=max(10.,0.1*f_g(tile)%bldheight+2.)                 ! updated in uclem_calc
      p_g(ifrac,tile)%lzom=log(p_g(ifrac,tile)%cndzmin/(0.1*f_g(tile)%bldheight)) ! updated in uclem_calc
      p_g(ifrac,tile)%lzoh=6.+p_g(ifrac,tile)%lzom ! (Kanda et al 2005)           ! updated in uclem_calc
      p_g(ifrac,tile)%cduv=(vkar/p_g(ifrac,tile)%lzom)**2                         ! updated in uclem_calc
      p_g(ifrac,tile)%cdtq=vkar**2/(p_g(ifrac,tile)%lzom*p_g(ifrac,tile)%lzoh)    ! updated in uclem_calc
      p_g(ifrac,tile)%tscrn=1.      ! + urbtemp                                   ! updated in uclem_calc
      p_g(ifrac,tile)%qscrn=0.                                                    ! updated in uclem_calc
      p_g(ifrac,tile)%uscrn=0.                                                    ! updated in uclem_calc
      p_g(ifrac,tile)%u10=0.                                                      ! updated in uclem_calc
      p_g(ifrac,tile)%emiss=0.97                                                  ! updated in uclem_calc
      p_g(ifrac,tile)%bldheat=0.
      p_g(ifrac,tile)%bldcool=0.
      p_g(ifrac,tile)%surfrunoff=0.
      p_g(ifrac,tile)%soilwetness=0.
      p_g(ifrac,tile)%soilwater=0.
      p_g(ifrac,tile)%transveg=0.
      p_g(ifrac,tile)%soilmoist=rdhyd_g(ifrac,tile)%soilwater*waterden*4.
      p_g(ifrac,tile)%rootmoistc=rdhyd_g(ifrac,tile)%soilwater*waterden*4.
      p_g(ifrac,tile)%delsoilmoist=0.
      p_g(ifrac,tile)%irrig=0.
      p_g(ifrac,tile)%acond_vegw=0.1
      p_g(ifrac,tile)%traf=0.
      p_g(ifrac,tile)%intgains_full=0.
      if ( ifrac==1 ) then 
        p_g(ifrac,tile)%frac_sigma=1.
      else
        p_g(ifrac,tile)%frac_sigma=0.  
      end if    
      p_g(ifrac,tile)%surferr=0._8
      p_g(ifrac,tile)%atmoserr=0._8
      p_g(ifrac,tile)%surferr_bias=0._8
      p_g(ifrac,tile)%atmoserr_bias=0._8
      p_g(ifrac,tile)%storage_flux=0._8
      p_g(ifrac,tile)%first_call=.true.
      p_g(ifrac,tile)%delintercept=0.
      p_g(ifrac,tile)%snowt=273.
      p_g(ifrac,tile)%vegt=273.
      p_g(ifrac,tile)%swe=0.
      p_g(ifrac,tile)%surfstor=0.
      p_g(ifrac,tile)%snowfrac=0.
      p_g(ifrac,tile)%salbedo=0.
      p_g(ifrac,tile)%calbedo=0.
      p_g(ifrac,tile)%taircanyon=273.
      p_g(ifrac,tile)%delswe=0.
      p_g(ifrac,tile)%ulai=1.
      p_g(ifrac,tile)%roof_water_runoff=0.
      p_g(ifrac,tile)%roof_snow_runoff=0.
      p_g(ifrac,tile)%roof_soil_runoff=0.
      p_g(ifrac,tile)%road_water_runoff=0.
      p_g(ifrac,tile)%road_snow_runoff=0.
      p_g(ifrac,tile)%road_soil_runoff=0.

    end if
    
  end do  
  
end do

return
end subroutine uclem_init

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine deallocates arrays used by the TEB scheme

subroutine uclem_end(diag)

implicit none

integer, intent(in) :: diag
integer tile, ifrac

if ( diag>=1 ) write(6,*) "Deallocating UCLEM arrays"

if ( uclem_active ) then

  do tile = 1,ntiles

    deallocate(f_roof(tile)%depth,f_wall(tile)%depth,f_road(tile)%depth,f_slab(tile)%depth,f_intm(tile)%depth)
    deallocate(f_roof(tile)%volcp,f_wall(tile)%volcp,f_road(tile)%volcp,f_slab(tile)%volcp,f_intm(tile)%volcp)
    deallocate(f_roof(tile)%lambda,f_wall(tile)%lambda,f_road(tile)%lambda,f_slab(tile)%lambda,f_intm(tile)%lambda)
    deallocate(f_roof(tile)%alpha,f_wall(tile)%alpha,f_road(tile)%alpha)
    deallocate(f_roof(tile)%emiss,f_wall(tile)%emiss,f_road(tile)%emiss)
    deallocate(f_slab(tile)%emiss)

    deallocate(cnveg_g(tile)%sigma,cnveg_g(tile)%alpha)
    deallocate(cnveg_g(tile)%emiss,rfveg_g(tile)%sigma,rfveg_g(tile)%alpha,rfveg_g(tile)%emiss)
    deallocate(cnveg_g(tile)%zo,cnveg_g(tile)%lai,cnveg_g(tile)%rsmin,rfveg_g(tile)%zo,rfveg_g(tile)%lai)
    deallocate(rfveg_g(tile)%rsmin,cnveg_g(tile)%temp,rfveg_g(tile)%temp)

    deallocate(intl_g(tile)%viewf,intl_g(tile)%psi)

    deallocate(f_g(tile)%sigmabld,f_g(tile)%hwratio,f_g(tile)%bldheight)
    deallocate(f_g(tile)%effhwratio)
    deallocate(f_g(tile)%industryfg,f_g(tile)%intgains_flr,f_g(tile)%trafficfg,f_g(tile)%vangle)
    deallocate(f_g(tile)%ctime,f_g(tile)%hangle,f_g(tile)%fbeam,f_g(tile)%weekdayload)
    deallocate(f_g(tile)%swilt,f_g(tile)%sfc,f_g(tile)%ssat)
    deallocate(f_g(tile)%bldairtemp,f_g(tile)%rfvegdepth)
    deallocate(f_g(tile)%intmassn,f_g(tile)%infilach,f_g(tile)%ventilach,f_g(tile)%heatprop,f_g(tile)%coolprop)
    deallocate(f_g(tile)%sigmau)
    
    do ifrac = 1,nfrac  
      
      deallocate(rfhyd_g(ifrac,tile)%surfwater,rfhyd_g(ifrac,tile)%snow,rfhyd_g(ifrac,tile)%den,rfhyd_g(ifrac,tile)%snowalpha)
      deallocate(rdhyd_g(ifrac,tile)%surfwater,rdhyd_g(ifrac,tile)%snow,rdhyd_g(ifrac,tile)%den,rdhyd_g(ifrac,tile)%snowalpha)
      deallocate(rdhyd_g(ifrac,tile)%leafwater,rdhyd_g(ifrac,tile)%soilwater,rfhyd_g(ifrac,tile)%leafwater)
      deallocate(rfhyd_g(ifrac,tile)%soilwater)

      deallocate(roof_g(ifrac,tile)%nodetemp,road_g(ifrac,tile)%nodetemp,walle_g(ifrac,tile)%nodetemp)
      deallocate(wallw_g(ifrac,tile)%nodetemp)
      deallocate(slab_g(ifrac,tile)%nodetemp,intm_g(ifrac,tile)%nodetemp,room_g(ifrac,tile)%nodetemp)
      deallocate(road_g(ifrac,tile)%storage,roof_g(ifrac,tile)%storage,walle_g(ifrac,tile)%storage,wallw_g(ifrac,tile)%storage)
      deallocate(slab_g(ifrac,tile)%storage,intm_g(ifrac,tile)%storage,room_g(ifrac,tile)%storage)
      
      deallocate(p_g(ifrac,tile)%lzom,p_g(ifrac,tile)%lzoh,p_g(ifrac,tile)%cndzmin,p_g(ifrac,tile)%cduv,p_g(ifrac,tile)%cdtq)
      deallocate(p_g(ifrac,tile)%tscrn,p_g(ifrac,tile)%qscrn,p_g(ifrac,tile)%uscrn,p_g(ifrac,tile)%u10,p_g(ifrac,tile)%emiss)
      deallocate(p_g(ifrac,tile)%surferr,p_g(ifrac,tile)%atmoserr,p_g(ifrac,tile)%surferr_bias,p_g(ifrac,tile)%atmoserr_bias)
      deallocate(p_g(ifrac,tile)%bldheat,p_g(ifrac,tile)%bldcool,p_g(ifrac,tile)%traf,p_g(ifrac,tile)%intgains_full)
      deallocate(p_g(ifrac,tile)%surfrunoff,p_g(ifrac,tile)%irrig,p_g(ifrac,tile)%acond_vegw)
      deallocate(p_g(ifrac,tile)%soilwetness,p_g(ifrac,tile)%soilwater)
      deallocate(p_g(ifrac,tile)%transveg,p_g(ifrac,tile)%soilmoist,p_g(ifrac,tile)%delsoilmoist)
      deallocate(p_g(ifrac,tile)%rootmoistc,p_g(ifrac,tile)%ulai)
      deallocate(p_g(ifrac,tile)%storage_flux,p_g(ifrac,tile)%snowmelt)
      deallocate(p_g(ifrac,tile)%delintercept,p_g(ifrac,tile)%snowt)
      deallocate(p_g(ifrac,tile)%vegt,p_g(ifrac,tile)%swe)
      deallocate(p_g(ifrac,tile)%surfstor,p_g(ifrac,tile)%snowfrac)
      deallocate(p_g(ifrac,tile)%salbedo,p_g(ifrac,tile)%calbedo)
      deallocate(p_g(ifrac,tile)%taircanyon,p_g(ifrac,tile)%delswe)
      
    end do  

  end do
  
  deallocate( roof_g, road_g, walle_g, wallw_g, slab_g, intm_g )
  deallocate( room_g )
  deallocate( f_roof, f_road, f_wall, f_slab, f_intm )
  deallocate( rfhyd_g, rdhyd_g )
  deallocate( cnveg_g, rfveg_g )
  deallocate( intl_g )
  deallocate( f_g )
  deallocate( p_g )
  deallocate( ufull_g )
  deallocate( upack_g )

end if

return
end subroutine uclem_end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! general version of tebload

subroutine uclem_loadd_3(urban,mode,ifrac,diag)

implicit none

integer, intent(in) :: ifrac, diag
integer ii, tile, is, ie
real(kind=8), dimension(ifull), intent(in) :: urban
character(len=*), intent(in) :: mode
character(len=10) :: teststr

if ( diag>=1 ) write(6,*) "Load UCLEM state array"
if (.not.uclem_active) return

do ii = 0,nl
  write(teststr,'("rooftemp",I1.1)') ii+1
  if ( trim(teststr)==trim(mode) ) then
    do tile = 1,ntiles
      is = (tile-1)*imax + 1
      ie = tile*imax
      if ( ufull_g(tile)>0 ) then
        roof_g(ifrac,tile)%nodetemp(:,ii)=pack(urban(is:ie),upack_g(:,tile))
        where ( roof_g(ifrac,tile)%nodetemp(:,ii)>150._8 )
          roof_g(ifrac,tile)%nodetemp(:,ii) = roof_g(ifrac,tile)%nodetemp(:,ii) - urbtemp  
        end where    
      end if
    end do
    return
  end if
  write(teststr,'("walletemp",I1.1)') ii+1
  if ( trim(teststr)==trim(mode) ) then
    do tile = 1,ntiles
      is = (tile-1)*imax + 1
      ie = tile*imax
      if ( ufull_g(tile)>0 ) then
        walle_g(ifrac,tile)%nodetemp(:,ii)=pack(urban(is:ie),upack_g(:,tile))
        where ( walle_g(ifrac,tile)%nodetemp(:,ii)>150._8 )
          walle_g(ifrac,tile)%nodetemp(:,ii) = walle_g(ifrac,tile)%nodetemp(:,ii) - urbtemp  
        end where  
      end if
    end do
    return
  end if
  write(teststr,'("wallwtemp",I1.1)') ii+1
  if ( trim(teststr)==trim(mode) ) then
    do tile = 1,ntiles
      is = (tile-1)*imax + 1
      ie = tile*imax
      if ( ufull_g(tile)>0 ) then
        wallw_g(ifrac,tile)%nodetemp(:,ii)=pack(urban(is:ie),upack_g(:,tile))
        where ( wallw_g(ifrac,tile)%nodetemp(:,ii)>150._8 )
          wallw_g(ifrac,tile)%nodetemp(:,ii) = wallw_g(ifrac,tile)%nodetemp(:,ii) - urbtemp  
        end where  
      end if
    end do
    return
  end if
  write(teststr,'("roadtemp",I1.1)') ii+1
  if ( trim(teststr)==trim(mode) ) then
    do tile = 1,ntiles
      is = (tile-1)*imax + 1
      ie = tile*imax
      if ( ufull_g(tile)>0 ) then
        road_g(ifrac,tile)%nodetemp(:,ii)=pack(urban(is:ie),upack_g(:,tile))
        where ( road_g(ifrac,tile)%nodetemp(:,ii)>150._8 )
          road_g(ifrac,tile)%nodetemp(:,ii) = road_g(ifrac,tile)%nodetemp(:,ii) - urbtemp  
        end where  
      end if
    end do
    return
  end if
  write(teststr,'("slabtemp",I1.1)') ii+1
  if ( trim(teststr)==trim(mode) ) then
    do tile = 1,ntiles
      is = (tile-1)*imax + 1
      ie = tile*imax
      if ( ufull_g(tile)>0 ) then
        slab_g(ifrac,tile)%nodetemp(:,ii)=pack(urban(is:ie),upack_g(:,tile))
        where ( slab_g(ifrac,tile)%nodetemp(:,ii)>150._8 )
          slab_g(ifrac,tile)%nodetemp(:,ii) = slab_g(ifrac,tile)%nodetemp(:,ii) - urbtemp  
        end where  
      end if
    end do
    return
  end if  
  write(teststr,'("intmtemp",I1.1)') ii+1
  if ( trim(teststr)==trim(mode) ) then
    do tile = 1,ntiles
      is = (tile-1)*imax + 1
      ie = tile*imax
      if ( ufull_g(tile)>0 ) then
        intm_g(ifrac,tile)%nodetemp(:,ii)=pack(urban(is:ie),upack_g(:,tile))
        where ( intm_g(ifrac,tile)%nodetemp(:,ii)>150._8 )
          intm_g(ifrac,tile)%nodetemp(:,ii) = intm_g(ifrac,tile)%nodetemp(:,ii) - urbtemp  
        end where  
      end if
    end do
    return
  end if   
end do  
  
select case(mode)
  case("roomtemp")
    do tile = 1,ntiles
      is = (tile-1)*imax + 1
      ie = tile*imax
      if ( ufull_g(tile)>0 ) then
        room_g(ifrac,tile)%nodetemp(:,1)=pack(urban(is:ie),upack_g(:,tile))  
        where ( room_g(ifrac,tile)%nodetemp(:,1)>150._8 )
          room_g(ifrac,tile)%nodetemp(:,1) = room_g(ifrac,tile)%nodetemp(:,1) - urbtemp  
        end where 
      end if
    end do
    return
end select
  
write(6,*) "ERROR: Unknown mode for uclem_loadd ",trim(mode)
stop
    
return
end subroutine uclem_loadd_3

subroutine uclem_loadd_2(urban,mode,ifrac,diag)

implicit none

integer, intent(in) :: ifrac, diag
integer tile, is, ie
real, dimension(ifull), intent(in) :: urban
character(len=*), intent(in) :: mode
character(len=10) :: teststr

if ( diag>=1 ) write(6,*) "Load UCLEM state array"
if (.not.uclem_active) return
  
select case(mode)
  case("canyonsoilmoisture")
    do tile = 1,ntiles
      is = (tile-1)*imax + 1
      ie = tile*imax
      if ( ufull_g(tile)>0 ) then
        rdhyd_g(ifrac,tile)%soilwater=pack(urban(is:ie),upack_g(:,tile))
      end if
    end do
    return
  case("roofsoilmoisture")
    do tile = 1,ntiles
      is = (tile-1)*imax + 1
      ie = tile*imax
      if ( ufull_g(tile)>0 ) then
        rfhyd_g(ifrac,tile)%soilwater=pack(urban(is:ie),upack_g(:,tile))
      end if
    end do
    return
  case("roadsurfacewater")
    do tile = 1,ntiles
      is = (tile-1)*imax + 1
      ie = tile*imax
      if ( ufull_g(tile)>0 ) then
        rdhyd_g(ifrac,tile)%surfwater=pack(urban(is:ie),upack_g(:,tile))  
      end if
    end do
    return
  case("roofsurfacewater")
    do tile = 1,ntiles
      is = (tile-1)*imax + 1
      ie = tile*imax
      if ( ufull_g(tile)>0 ) then
        rfhyd_g(ifrac,tile)%surfwater=pack(urban(is:ie),upack_g(:,tile))  
      end if
    end do
    return
  case("canyonleafwater")
    do tile = 1,ntiles
      is = (tile-1)*imax + 1
      ie = tile*imax
      if ( ufull_g(tile)>0 ) then
        rdhyd_g(ifrac,tile)%leafwater=pack(urban(is:ie),upack_g(:,tile))  
      end if
    end do
    return
  case("roofleafwater")
    do tile = 1,ntiles
      is = (tile-1)*imax + 1
      ie = tile*imax
      if ( ufull_g(tile)>0 ) then
        rfhyd_g(ifrac,tile)%leafwater=pack(urban(is:ie),upack_g(:,tile))  
      end if
    end do
    return
  case("roadsnowdepth")
    do tile = 1,ntiles
      is = (tile-1)*imax + 1
      ie = tile*imax
      if ( ufull_g(tile)>0 ) then
        rdhyd_g(ifrac,tile)%snow=pack(urban(is:ie),upack_g(:,tile))  
      end if
    end do
    return
  case("roofsnowdepth")
    do tile = 1,ntiles
      is = (tile-1)*imax + 1
      ie = tile*imax
      if ( ufull_g(tile)>0 ) then
        rfhyd_g(ifrac,tile)%snow=pack(urban(is:ie),upack_g(:,tile))  
      end if
    end do
    return
  case("roadsnowdensity")
    do tile = 1,ntiles
      is = (tile-1)*imax + 1
      ie = tile*imax
      if ( ufull_g(tile)>0 ) then
        rdhyd_g(ifrac,tile)%den=pack(urban(is:ie),upack_g(:,tile))  
      end if
    end do
    return
  case("roofsnowdensity")
    do tile = 1,ntiles
      is = (tile-1)*imax + 1
      ie = tile*imax
      if ( ufull_g(tile)>0 ) then
        rfhyd_g(ifrac,tile)%den=pack(urban(is:ie),upack_g(:,tile))  
      end if
    end do
    return
  case("roadsnowalbedo")
    do tile = 1,ntiles
      is = (tile-1)*imax + 1
      ie = tile*imax
      if ( ufull_g(tile)>0 ) then
        rdhyd_g(ifrac,tile)%snowalpha=pack(urban(is:ie),upack_g(:,tile))  
      end if
    end do
    return
  case("roofsnowalbedo")
    do tile = 1,ntiles
      is = (tile-1)*imax + 1
      ie = tile*imax
      if ( ufull_g(tile)>0 ) then
        rfhyd_g(ifrac,tile)%snowalpha=pack(urban(is:ie),upack_g(:,tile))  
      end if
    end do
    return
end select
  
write(6,*) "ERROR: Unknown mode for uclem_loadd ",trim(mode)
stop

return
end subroutine uclem_loadd_2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! this subroutine loads UCLEM type arrays (not compulsory)

! URBAN TYPES:
 
! 1 = Urban
! 2 = Urban (low)
! 3 = Urban (medium)
! 4 = Urban (high)
! 5 = Urban (cbd)
! 6 = Industrial (low)
! 7 = Industrial (medium)
! 8 = Industrial (high)

subroutine uclem_type_standard(itype,diag)

implicit none

integer, intent(in) :: diag
integer tile, is, ie
integer, dimension(ifull), intent(in) :: itype

if (diag>=1) write(6,*) "Load UCLEM building properties"
if (.not.uclem_active) return

do tile = 1,ntiles
  is = (tile-1)*imax + 1
  ie = tile*imax
  if ( ufull_g(tile)>0 ) then
    call uclem_type_thread(itype(is:ie),diag,f_g(tile),cnveg_g(tile),rfveg_g(tile), &
                           f_roof(tile),f_road(tile),f_wall(tile),f_slab(tile),     &
                           f_intm(tile),intl_g(tile),upack_g(:,tile),               &
                           ufull_g(tile))
  end if
end do

return
end subroutine uclem_type_standard

subroutine uclem_type_thread(itype,diag,fp,cnveg,rfveg,fp_roof,fp_road,fp_wall,fp_slab, &
                             fp_intm,intl,upack,ufull)

implicit none

integer, intent(in) :: diag, ufull
integer ii,j,ierr,nlp
integer, dimension(imax), intent(in) :: itype
integer, dimension(ufull) :: itmp
integer, parameter :: maxtype = 8
real x
real, dimension(ufull) :: tsigveg,tsigmabld,fp_coeffbldheight
! In-canyon vegetation fraction
real, dimension(maxtype) ::    csigvegc=(/ 0.38, 0.45, 0.38, 0.34, 0.05, 0.40, 0.30, 0.20 /)
! Green roof vegetation fraction
real, dimension(maxtype) ::    csigvegr=(/ 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00 /)
! Area fraction occupied by buildings
real, dimension(maxtype) ::   csigmabld=(/ 0.45, 0.40, 0.45, 0.46, 0.65, 0.40, 0.45, 0.50 /)
! Building height (m)
real, dimension(maxtype) ::  cbldheight=(/   6.,   4.,   6.,   8.,  18.,   4.,   8.,  12. /)
! Building height to width ratio
real, dimension(maxtype) ::    chwratio=(/  0.4,  0.2,  0.4,  0.6,   2.,  0.5,   1.,  1.5 /)
! Industrial sensible heat flux (W m^-2)
real, dimension(maxtype) :: cindustryfg=(/   0.,   0.,   0.,   0.,   0.,  10.,  20.,  30. /)
! Internal gains sensible heat flux [floor] (W m^-2)
real, dimension(maxtype) ::   cintgains=(/   5.,   5.,   5.,   5.,   5.,   5.,   5.,   5. /)
! Daily averaged traffic sensible heat flux (W m^-2)
real, dimension(maxtype) ::  ctrafficfg=(/  1.5,  1.5,  1.5,  1.5,  1.5,  1.5,  1.5,  1.5 /)
! Comfort temperature (K)
real, dimension(maxtype) :: cbldtemp=(/ 295.16, 295.16, 295.16, 295.16, 295.16, 295.16, 295.16, 295.16 /)
! Roof albedo
real, dimension(maxtype) ::  croofalpha=(/ 0.20, 0.20, 0.20, 0.20, 0.20, 0.20, 0.20, 0.20 /)    ! (Fortuniak 08) Masson = 0.15
! Wall albedo
real, dimension(maxtype) ::  cwallalpha=(/ 0.30, 0.30, 0.30, 0.30, 0.30, 0.30, 0.30, 0.30 /)    ! (Fortuniak 08) Masson = 0.25
! Road albedo
real, dimension(maxtype) ::  croadalpha=(/ 0.10, 0.10, 0.10, 0.10, 0.10, 0.10, 0.10, 0.10 /)    ! (Fortuniak 08) Masson = 0.08
! Canyon veg albedo
real, dimension(maxtype) ::  cvegalphac=(/ 0.20, 0.20, 0.20, 0.20, 0.20, 0.20, 0.20, 0.20 /)
! Roof veg albedo
real, dimension(maxtype) ::  cvegalphar=(/ 0.20, 0.20, 0.20, 0.20, 0.20, 0.20, 0.20, 0.20 /)
! Roof emissivity
real, dimension(maxtype) ::  croofemiss=(/ 0.90, 0.90, 0.90, 0.90, 0.90, 0.90, 0.90, 0.90 /)
! Wall emissivity
real, dimension(maxtype) ::  cwallemiss=(/ 0.85, 0.85, 0.85, 0.85, 0.85, 0.85, 0.85, 0.85 /) 
! Road emissivity
real, dimension(maxtype) ::  croademiss=(/ 0.94, 0.94, 0.94, 0.94, 0.94, 0.94, 0.94, 0.94 /)
! Slab emissivity
real, dimension(maxtype) ::  cslabemiss=(/ 0.90, 0.90, 0.90, 0.90, 0.90, 0.90, 0.90, 0.90 /) 
! Canyon veg emissivity
real, dimension(maxtype) ::  cvegemissc=(/ 0.96, 0.96, 0.96, 0.96, 0.96, 0.96, 0.96, 0.96 /)
! Roof veg emissivity
real, dimension(maxtype) ::  cvegemissr=(/ 0.96, 0.96, 0.96, 0.96, 0.96, 0.96, 0.96, 0.96 /)
! Green roof soil depth
real, dimension(maxtype) ::   cvegdeptr=(/ 0.10, 0.10, 0.10, 0.10, 0.10, 0.10, 0.10, 0.10 /)
! Roughness length of in-canyon vegetation (m)
real, dimension(maxtype) ::    czovegc=(/   0.1,   0.1,   0.1,   0.1,   0.1,   0.1,   0.1,   0.1 /)
! In-canyon vegetation LAI
real, dimension(maxtype) ::  cvegrlaic=(/   2.0,   2.0,   2.0,   2.0,   2.0,   2.0,   2.0,   2.0 /)
! Unconstrained canopy stomatal resistance
real, dimension(maxtype) :: cvegrsminc=(/ 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0 /)
! Roughness length of green roof vegetation (m)
real, dimension(maxtype) ::    czovegr=(/   0.1,   0.1,   0.1,   0.1,   0.1,   0.1,   0.1,   0.1 /)
! Green roof vegetation LAI
real, dimension(maxtype) ::  cvegrlair=(/   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0 /)
! Unconstrained canopy stomatal resistance
real, dimension(maxtype) :: cvegrsminr=(/ 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0 /)
! Soil wilting point (m^3 m^-3)
real, dimension(maxtype) ::     cswilt=(/  0.18,  0.18,  0.18,  0.18,  0.18,  0.18,  0.18,  0.18 /)
! Soil field capacity (m^3 m^-3)
real, dimension(maxtype) ::       csfc=(/  0.26,  0.26,  0.26,  0.26,  0.26,  0.26,  0.26,  0.26 /)
! Soil saturation point (m^3 m^-3)
real, dimension(maxtype) ::      cssat=(/  0.42,  0.42,  0.42,  0.42,  0.42,  0.42,  0.42,  0.42 /)
! Infiltration air volume changes per hour (m^3 m^-3)
real, dimension(maxtype) ::  cinfilach=(/  0.50,  0.50,  0.50,  0.50,  0.50,  0.50,  0.50,  0.50 /)
! Ventilation air volume changes per hour (m^3 m^-3)
real, dimension(maxtype) :: cventilach=(/  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00 /)
! Fraction of spaces with heating devices
real, dimension(maxtype) ::  cheatprop=(/  0.5,  0.5,  0.5,  0.5,  1.0,  0.00,  0.00,  0.00 /)
! Fraction of spaces with cooling devices
real, dimension(maxtype) ::  ccoolprop=(/  0.5,  0.5,  0.5,  0.5,  1.0,  0.00,  0.00,  0.00 /)

real, dimension(maxtype,nl) :: croofdepth
real, dimension(maxtype,nl) :: cwalldepth
real, dimension(maxtype,nl) :: croaddepth
real, dimension(maxtype,nl) :: cslabdepth
real, dimension(maxtype,nl) :: croofcp
real, dimension(maxtype,nl) :: cwallcp
real, dimension(maxtype,nl) :: croadcp
real, dimension(maxtype,nl) :: cslabcp
real, dimension(maxtype,nl) :: crooflambda
real, dimension(maxtype,nl) :: cwalllambda
real, dimension(maxtype,nl) :: croadlambda
real, dimension(maxtype,nl) :: cslablambda

logical, dimension(imax), intent(in) :: upack

type(fparmdata), intent(inout) :: fp
type(vegdata), intent(inout) :: cnveg, rfveg
type(facetparams), intent(inout) :: fp_roof, fp_road, fp_wall, fp_slab, fp_intm
type(intldata), intent(inout) :: intl


! facet array where: rows=maxtypes (landtypes) and columns=nl (material layers)
nlp=nl/4 ! number of layers in each material segment (over 4 material segments)
! depths (m)
croofdepth= reshape((/ ((0.01/real(nlp), ii=1,maxtype),j=1,nlp),    &
                       ((0.09/real(nlp), ii=1,maxtype),j=1,nlp),    & 
                       ((0.40/real(nlp), ii=1,maxtype),j=1,nlp),    &
                       ((0.10/real(nlp), ii=1,maxtype),j=1,nlp) /), &
                       (/maxtype,nl/))
cwalldepth= reshape((/ ((0.01/real(nlp), ii=1,maxtype),j=1,nlp),    &
                       ((0.04/real(nlp), ii=1,maxtype),j=1,nlp),    & 
                       ((0.10/real(nlp), ii=1,maxtype),j=1,nlp),    &
                       ((0.05/real(nlp), ii=1,maxtype),j=1,nlp) /), &
                       (/maxtype,nl/))
croaddepth= reshape((/ ((0.01/real(nlp), ii=1,maxtype),j=1,nlp),    &
                       ((0.04/real(nlp), ii=1,maxtype),j=1,nlp),    & 
                       ((0.45/real(nlp), ii=1,maxtype),j=1,nlp),    &
                       ((3.50/real(nlp), ii=1,maxtype),j=1,nlp) /), &
                       (/maxtype,nl/))
cslabdepth=reshape((/  ((0.05/real(nlp), ii=1,maxtype),j=1,nlp),    &
                       ((0.05/real(nlp), ii=1,maxtype),j=1,nlp),    & 
                       ((0.05/real(nlp), ii=1,maxtype),j=1,nlp),    &
                       ((0.05/real(nlp), ii=1,maxtype),j=1,nlp) /), &
                       (/maxtype,nl/))
! heat capacity (J m^-3 K^-1)
croofcp =   reshape((/ ((2.11E6, ii=1,maxtype),j=1,nlp),    & ! dense concrete (Oke 87)
                       ((2.11E6, ii=1,maxtype),j=1,nlp),    & ! dense concrete (Oke 87)
                       ((0.28E6, ii=1,maxtype),j=1,nlp),    & ! light concrete (Oke 87)
                       ((0.29E6, ii=1,maxtype),j=1,nlp) /), & ! insulation (Oke 87)
                       (/maxtype,nl/))
cwallcp =   reshape((/ ((1.55E6, ii=1,maxtype),j=1,nlp),    & ! concrete (Mills 93)
                       ((1.55E6, ii=1,maxtype),j=1,nlp),    & ! concrete (Mills 93)
                       ((1.55E6, ii=1,maxtype),j=1,nlp),    & ! concrete (Mills 93)
                       ((0.29E6, ii=1,maxtype),j=1,nlp) /), & ! insulation (Oke 87)
                       (/maxtype,nl/))
croadcp =   reshape((/ ((1.94E6, ii=1,maxtype),j=1,nlp),    & ! asphalt (Mills 93)
                       ((1.94E6, ii=1,maxtype),j=1,nlp),    & ! asphalt (Mills 93)
                       ((1.28E6, ii=1,maxtype),j=1,nlp),    & ! dry soil (Mills 93)
                       ((1.28E6, ii=1,maxtype),j=1,nlp) /), & ! dry soil (Mills 93)
                       (/maxtype,nl/))
cslabcp=   reshape((/  ((1.55E6, ii=1,maxtype),j=1,nlp),    & ! concrete (Mills 93)
                       ((1.55E6, ii=1,maxtype),j=1,nlp),    & ! concrete (Mills 93)
                       ((1.55E6, ii=1,maxtype),j=1,nlp),    & ! concrete (Mills 93)
                       ((1.55E6, ii=1,maxtype),j=1,nlp) /), & ! concrete (Mills 93)
                       (/maxtype,nl/))
! heat conductivity (W m^-1 K^-1)
crooflambda=reshape((/ ((1.5100, ii=1,maxtype),j=1,nlp),    & ! dense concrete (Oke 87)
                       ((1.5100, ii=1,maxtype),j=1,nlp),    & ! dense concrete (Oke 87)
                       ((0.0800, ii=1,maxtype),j=1,nlp),    & ! light concrete (Oke 87)
                       ((0.0500, ii=1,maxtype),j=1,nlp) /), & ! insulation (Oke 87)
                       (/maxtype,nl/))
cwalllambda=reshape((/ ((0.9338, ii=1,maxtype),j=1,nlp),    & ! concrete (Mills 93)
                       ((0.9338, ii=1,maxtype),j=1,nlp),    & ! concrete (Mills 93)
                       ((0.9338, ii=1,maxtype),j=1,nlp),    & ! concrete (Mills 93)
                       ((0.0500, ii=1,maxtype),j=1,nlp) /), & ! insulation (Oke 87)
                       (/maxtype,nl/))
croadlambda=reshape((/ ((0.7454, ii=1,maxtype),j=1,nlp),    & ! asphalt (Mills 93)
                       ((0.7454, ii=1,maxtype),j=1,nlp),    & ! asphalt (Mills 93)
                       ((0.2513, ii=1,maxtype),j=1,nlp),    & ! dry soil (Mills 93)
                       ((0.2513, ii=1,maxtype),j=1,nlp) /), & ! dry soil (Mills 93)
                       (/maxtype,nl/))
cslablambda=reshape((/ ((0.9338, ii=1,maxtype),j=1,nlp),    & ! concrete (Mills 93)
                       ((0.9338, ii=1,maxtype),j=1,nlp),    & ! concrete (Mills 93)
                       ((0.9338, ii=1,maxtype),j=1,nlp),    & ! concrete (Mills 93)
                       ((0.9338, ii=1,maxtype),j=1,nlp) /), & ! concrete (Mills 93)
                       (/maxtype,nl/))

itmp=pack(itype,upack)
where ( itmp>maxtype ) ! set default to new urban types as generic
  itmp = 1
end where
if ((minval(itmp)<1).or.(maxval(itmp)>maxtype)) then
  write(6,*) "ERROR: Urban type is out of range"
  stop
end if

tsigveg = csigvegc(itmp)
tsigmabld = csigmabld(itmp)
cnveg%sigma = max(min(tsigveg/(1.-tsigmabld),1.),0.)
rfveg%sigma = max(min(csigvegr(itmp),1.),0.)
fp%sigmabld = max(min(tsigmabld,1.),0.)
fp%hwratio = chwratio(itmp)

fp%industryfg = cindustryfg(itmp)
fp%intgains_flr = cintgains(itmp)
fp%trafficfg = ctrafficfg(itmp)
fp%bldheight = cbldheight(itmp)
fp_roof%alpha = croofalpha(itmp)
fp_wall%alpha = cwallalpha(itmp)
fp_road%alpha = croadalpha(itmp)
cnveg%alpha = cvegalphac(itmp)
rfveg%alpha = cvegalphar(itmp)
fp_roof%emiss = croofemiss(itmp)
fp_wall%emiss = cwallemiss(itmp)
fp_road%emiss = croademiss(itmp)
cnveg%emiss = cvegemissc(itmp)
rfveg%emiss = cvegemissr(itmp)
fp%bldairtemp = cbldtemp(itmp) - urbtemp
fp%rfvegdepth = cvegdeptr(itmp)
do ii = 1,nl
  fp_roof%depth(:,ii) = croofdepth(itmp,ii)
  fp_wall%depth(:,ii) = cwalldepth(itmp,ii)
  fp_road%depth(:,ii) = croaddepth(itmp,ii)
  fp_roof%lambda(:,ii) = crooflambda(itmp,ii)
  fp_wall%lambda(:,ii) = cwalllambda(itmp,ii)
  fp_road%lambda(:,ii) = croadlambda(itmp,ii)
  fp_roof%volcp(:,ii) = croofcp(itmp,ii)
  fp_wall%volcp(:,ii) = cwallcp(itmp,ii)
  select case(soilunder)
    case(0) ! storage under road only
      fp_road%volcp(:,ii) = croadcp(itmp,ii)
    case(1) ! storage under road and canveg
      fp_road%volcp(:,ii) = croadcp(itmp,ii)/(1.-cnveg%sigma)
    case(2) ! storage under road and bld
      fp_road%volcp(:,ii) = croadcp(itmp,ii)*(1./(1.-cnveg%sigma)*(1./(1.-fp%sigmabld)-1.) +1.)
    case(3) ! storage under road and canveg and bld (100% of grid point)
      fp_road%volcp(:,ii) = croadcp(itmp,ii)/(1.-cnveg%sigma)/(1.-fp%sigmabld)
    case DEFAULT
      write(6,*) "ERROR: Unknown soilunder mode ",soilunder
      stop
  end select
end do
cnveg%zo = czovegc(itmp)
cnveg%lai = cvegrlaic(itmp)
cnveg%rsmin = cvegrsminc(itmp)/max(cnveg%lai,1.E-8)
rfveg%zo = czovegr(itmp)
rfveg%lai = cvegrlair(itmp)
rfveg%rsmin = cvegrsminr(itmp)/max(rfveg%lai,1.E-8)
fp%swilt = cswilt(itmp)
fp%sfc = csfc(itmp)
fp%ssat = cssat(itmp)

! for varying internal temperature
fp_slab%emiss = cslabemiss(itmp)
fp%infilach = cinfilach(itmp)
fp%ventilach = cventilach(itmp)
fp%heatprop = cheatprop(itmp)
fp%coolprop = ccoolprop(itmp)
do ii = 1,nl
  fp_slab%depth(:,ii) = cslabdepth(itmp,ii)
  fp_intm%depth(:,ii) = cslabdepth(itmp,ii)  ! internal mass material same as slab
  fp_slab%lambda(:,ii) = cslablambda(itmp,ii)
  fp_intm%lambda(:,ii) = cslablambda(itmp,ii)
  fp_slab%volcp(:,ii) = cslabcp(itmp,ii)
  fp_intm%volcp(:,ii) = cslabcp(itmp,ii)
end do

! Here we modify the effective canyon geometry to account for in-canyon vegetation tall vegetation
fp_coeffbldheight = max(fp%bldheight-6.*cnveg%zo,0.2)/fp%bldheight
fp%effhwratio     = fp%hwratio*fp_coeffbldheight

call init_internal(fp)
call init_lwcoeff(fp,intl,ufull)

if ( diag>0 ) then
  write(6,*) 'hwratio, eff',fp%hwratio, fp%effhwratio
  write(6,*) 'bldheight, eff',fp%bldheight, fp_coeffbldheight
  write(6,*) 'sigmabld, sigmavegc', fp%sigmabld, cnveg%sigma
  write(6,*) 'roadcp multiple for soilunder:', soilunder,fp_road%volcp(itmp,1)/croadcp(itmp,1)
end if

return
end subroutine uclem_type_thread

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine uclem_deftype(paramdata,typedata,paramname,diag)

implicit none

integer, intent(in) :: diag
integer tile, is, ie, i, ifrac, maxtype
integer iq
integer, dimension(ifull), intent(in) :: typedata
integer, dimension(imax) :: itmp
real, dimension(:), intent(in) :: paramdata
real fp_coeffbldheight
logical found
character(len=*), intent(in) :: paramname
character(len=20) :: vname

if ( diag>=1 ) write(6,*) "Load UCLEM parameter: ",trim(paramname)
if ( .not.uclem_active ) return

maxtype = size(paramdata)

select case(paramname)
  case('bldheight')
    do tile = 1,ntiles
      is = (tile-1)*imax + 1
      ie = tile*imax  
      if ( ufull_g(tile)>0 ) then
        itmp(1:ufull_g(tile)) = pack(typedata(is:ie),upack_g(:,tile))
        if ( minval(itmp(1:ufull_g(tile)))<1 .or. maxval(itmp(1:ufull_g(tile)))>maxtype ) then
          write(6,*) "ERROR: Urban type is out of range"
          stop 
        end if
        f_g(tile)%bldheight = paramdata(itmp(1:ufull_g(tile)))
      end if
    end do  
  case('hwratio')
    do tile = 1,ntiles
      is = (tile-1)*imax + 1
      ie = tile*imax  
      if ( ufull_g(tile)>0 ) then
        itmp(1:ufull_g(tile)) = pack(typedata(is:ie),upack_g(:,tile))
        if ( minval(itmp(1:ufull_g(tile)))<1 .or. maxval(itmp(1:ufull_g(tile)))>maxtype ) then
          write(6,*) "ERROR: Urban type is out of range"
          stop 
        end if
        f_g(tile)%hwratio = paramdata(itmp(1:ufull_g(tile)))
      end if
    end do
  case('sigvegc')
    do tile = 1,ntiles
      is = (tile-1)*imax + 1
      ie = tile*imax  
      if ( ufull_g(tile)>0 ) then
        itmp(1:ufull_g(tile)) = pack(typedata(is:ie),upack_g(:,tile))
        if ( minval(itmp(1:ufull_g(tile)))<1 .or. maxval(itmp(1:ufull_g(tile)))>maxtype ) then
          write(6,*) "ERROR: Urban type is out of range"
          stop 
        end if
        cnveg_g(tile)%sigma = paramdata(itmp(1:ufull_g(tile)))/(1.-f_g(tile)%sigmabld)  
      end if
    end do  
  case('sigmabld')
    do tile = 1,ntiles
      is = (tile-1)*imax + 1
      ie = tile*imax  
      if ( ufull_g(tile)>0 ) then
        itmp(1:ufull_g(tile)) = pack(typedata(is:ie),upack_g(:,tile))
        if ( minval(itmp(1:ufull_g(tile)))<1 .or. maxval(itmp(1:ufull_g(tile)))>maxtype ) then
          write(6,*) "ERROR: Urban type is out of range"
          stop 
        end if
        cnveg_g(tile)%sigma = cnveg_g(tile)%sigma*(1.-f_g(tile)%sigmabld)/(1.-paramdata(itmp(1:ufull_g(tile))))
        f_g(tile)%sigmabld = paramdata(itmp(1:ufull_g(tile)))
      end if
    end do  
  case('industryfg')
    do tile = 1,ntiles
      is = (tile-1)*imax + 1
      ie = tile*imax  
      if ( ufull_g(tile)>0 ) then
        itmp(1:ufull_g(tile)) = pack(typedata(is:ie),upack_g(:,tile))
        if ( minval(itmp(1:ufull_g(tile)))<1 .or. maxval(itmp(1:ufull_g(tile)))>maxtype ) then
          write(6,*) "ERROR: Urban type is out of range"
          stop 
        end if
        f_g(tile)%industryfg = paramdata(itmp(1:ufull_g(tile)))
      end if
    end do  
  case('trafficfg')
    do tile = 1,ntiles
      is = (tile-1)*imax + 1
      ie = tile*imax  
      if ( ufull_g(tile)>0 ) then
        itmp(1:ufull_g(tile)) = pack(typedata(is:ie),upack_g(:,tile))
        if ( minval(itmp(1:ufull_g(tile)))<1 .or. maxval(itmp(1:ufull_g(tile)))>maxtype ) then
          write(6,*) "ERROR: Urban type is out of range"
          stop 
        end if
        f_g(tile)%trafficfg = paramdata(itmp(1:ufull_g(tile)))
      end if
    end do  
  case('roofalpha')
    do tile = 1,ntiles
      is = (tile-1)*imax + 1
      ie = tile*imax  
      if ( ufull_g(tile)>0 ) then
        itmp(1:ufull_g(tile)) = pack(typedata(is:ie),upack_g(:,tile))
        if ( minval(itmp(1:ufull_g(tile)))<1 .or. maxval(itmp(1:ufull_g(tile)))>maxtype ) then
          write(6,*) "ERROR: Urban type is out of range"
          stop 
        end if
        f_roof(tile)%alpha = paramdata(itmp(1:ufull_g(tile)))
      end if
    end do  
  case('wallalpha')
    do tile = 1,ntiles
      is = (tile-1)*imax + 1
      ie = tile*imax  
      if ( ufull_g(tile)>0 ) then
        itmp(1:ufull_g(tile)) = pack(typedata(is:ie),upack_g(:,tile))
        if ( minval(itmp(1:ufull_g(tile)))<1 .or. maxval(itmp(1:ufull_g(tile)))>maxtype ) then
          write(6,*) "ERROR: Urban type is out of range"
          stop 
        end if
        f_wall(tile)%alpha = paramdata(itmp(1:ufull_g(tile)))
      end if
    end do  
  case('roadalpha')
    do tile = 1,ntiles
      is = (tile-1)*imax + 1
      ie = tile*imax  
      if ( ufull_g(tile)>0 ) then
        itmp(1:ufull_g(tile)) = pack(typedata(is:ie),upack_g(:,tile))
        if ( minval(itmp(1:ufull_g(tile)))<1 .or. maxval(itmp(1:ufull_g(tile)))>maxtype ) then
          write(6,*) "ERROR: Urban type is out of range"
          stop 
        end if
        f_road(tile)%alpha = paramdata(itmp(1:ufull_g(tile)))
      end if
    end do  
  case('roofemiss')
    do tile = 1,ntiles
      is = (tile-1)*imax + 1
      ie = tile*imax  
      if ( ufull_g(tile)>0 ) then
        itmp(1:ufull_g(tile)) = pack(typedata(is:ie),upack_g(:,tile))
        if ( minval(itmp(1:ufull_g(tile)))<1 .or. maxval(itmp(1:ufull_g(tile)))>maxtype ) then
          write(6,*) "ERROR: Urban type is out of range"
          stop 
        end if
        f_roof(tile)%emiss = paramdata(itmp(1:ufull_g(tile)))
      end if
    end do  
  case('wallemiss')
    do tile = 1,ntiles
      is = (tile-1)*imax + 1
      ie = tile*imax  
      if ( ufull_g(tile)>0 ) then
        itmp(1:ufull_g(tile)) = pack(typedata(is:ie),upack_g(:,tile))
        if ( minval(itmp(1:ufull_g(tile)))<1 .or. maxval(itmp(1:ufull_g(tile)))>maxtype ) then
          write(6,*) "ERROR: Urban type is out of range"
          stop 
        end if
        f_wall(tile)%emiss = paramdata(itmp(1:ufull_g(tile)))
      end if
    end do  
  case('roademiss')
    do tile = 1,ntiles
      is = (tile-1)*imax + 1
      ie = tile*imax  
      if ( ufull_g(tile)>0 ) then
        itmp(1:ufull_g(tile)) = pack(typedata(is:ie),upack_g(:,tile))
        if ( minval(itmp(1:ufull_g(tile)))<1 .or. maxval(itmp(1:ufull_g(tile)))>maxtype ) then
          write(6,*) "ERROR: Urban type is out of range"
          stop 
        end if
        f_road(tile)%emiss = paramdata(itmp(1:ufull_g(tile)))
      end if
    end do   
  case('vegalphac')
    do tile = 1,ntiles
      is = (tile-1)*imax + 1
      ie = tile*imax  
      if ( ufull_g(tile)>0 ) then
        itmp(1:ufull_g(tile)) = pack(typedata(is:ie),upack_g(:,tile))
        if ( minval(itmp(1:ufull_g(tile)))<1 .or. maxval(itmp(1:ufull_g(tile)))>maxtype ) then
          write(6,*) "ERROR: Urban type is out of range"
          stop 
        end if
        cnveg_g(tile)%alpha = paramdata(itmp(1:ufull_g(tile)))
      end if
    end do  
  case('zovegc')
    do tile = 1,ntiles
      is = (tile-1)*imax + 1
      ie = tile*imax  
      if ( ufull_g(tile)>0 ) then
        itmp(1:ufull_g(tile)) = pack(typedata(is:ie),upack_g(:,tile))
        if ( minval(itmp(1:ufull_g(tile)))<1 .or. maxval(itmp(1:ufull_g(tile)))>maxtype ) then
          write(6,*) "ERROR: Urban type is out of range"
          stop 
        end if
        cnveg_g(tile)%zo = paramdata(itmp(1:ufull_g(tile)))
      end if
    end do  
  case('infilach')
    do tile = 1,ntiles
      is = (tile-1)*imax + 1
      ie = tile*imax  
      if ( ufull_g(tile)>0 ) then
        itmp(1:ufull_g(tile)) = pack(typedata(is:ie),upack_g(:,tile))
        if ( minval(itmp(1:ufull_g(tile)))<1 .or. maxval(itmp(1:ufull_g(tile)))>maxtype ) then
          write(6,*) "ERROR: Urban type is out of range"
          stop 
        end if
        f_g(tile)%infilach = paramdata(itmp(1:ufull_g(tile)))
      end if
    end do  
  case('intgains')
    do tile = 1,ntiles
      is = (tile-1)*imax + 1
      ie = tile*imax  
      if ( ufull_g(tile)>0 ) then
        itmp(1:ufull_g(tile)) = pack(typedata(is:ie),upack_g(:,tile))
        if ( minval(itmp(1:ufull_g(tile)))<1 .or. maxval(itmp(1:ufull_g(tile)))>maxtype ) then
          write(6,*) "ERROR: Urban type is out of range"
          stop 
        end if
        f_g(tile)%intgains_flr = paramdata(itmp(1:ufull_g(tile)))
      end if
    end do  
  case('bldairtemp')
    do tile = 1,ntiles
      is = (tile-1)*imax + 1
      ie = tile*imax  
      if ( ufull_g(tile)>0 ) then
        itmp(1:ufull_g(tile)) = pack(typedata(is:ie),upack_g(:,tile))
        if ( minval(itmp(1:ufull_g(tile)))<1 .or. maxval(itmp(1:ufull_g(tile)))>maxtype ) then
          write(6,*) "ERROR: Urban type is out of range"
          stop 
        end if
        f_g(tile)%bldairtemp = paramdata(itmp(1:ufull_g(tile))) - urbtemp
        do ifrac = 1,nfrac
          room_g(ifrac,tile)%nodetemp(:,1) = real(f_g(tile)%bldairtemp,8)
        end do  
      end if
    end do   
  case('heatprop')
    do tile = 1,ntiles
      is = (tile-1)*imax + 1
      ie = tile*imax  
      if ( ufull_g(tile)>0 ) then
        itmp(1:ufull_g(tile)) = pack(typedata(is:ie),upack_g(:,tile))
        if ( minval(itmp(1:ufull_g(tile)))<1 .or. maxval(itmp(1:ufull_g(tile)))>maxtype ) then
          write(6,*) "ERROR: Urban type is out of range"
          stop 
        end if
        f_g(tile)%heatprop = paramdata(itmp(1:ufull_g(tile)))
      end if
    end do   
  case('coolprop')
    do tile = 1,ntiles
      is = (tile-1)*imax + 1
      ie = tile*imax  
      if ( ufull_g(tile)>0 ) then
        itmp(1:ufull_g(tile)) = pack(typedata(is:ie),upack_g(:,tile))
        if ( minval(itmp(1:ufull_g(tile)))<1 .or. maxval(itmp(1:ufull_g(tile)))>maxtype ) then
          write(6,*) "ERROR: Urban type is out of range"
          stop 
        end if
        f_g(tile)%coolprop = paramdata(itmp(1:ufull_g(tile)))
      end if
    end do  
  case default
    found = .false.
    do i = 1,4
      write(vname,'("roofthick",(I1.1))') i
      if ( trim(paramname)==trim(vname) ) then
        do tile = 1,ntiles
          is = (tile-1)*imax + 1
          ie = tile*imax  
          if ( ufull_g(tile)>0 ) then
            itmp(1:ufull_g(tile)) = pack(typedata(is:ie),upack_g(:,tile))
            if ( minval(itmp(1:ufull_g(tile)))<1 .or. maxval(itmp(1:ufull_g(tile)))>maxtype ) then
              write(6,*) "ERROR: Urban type is out of range"
              stop 
            end if
            f_roof(tile)%depth(:,i) = paramdata(itmp(1:ufull_g(tile)))
          end if
        end do 
        found = .true.
        exit
      end if
      write(vname,'("roofcp",(I1.1))') i
      if ( trim(paramname)==trim(vname) ) then
        do tile = 1,ntiles
          is = (tile-1)*imax + 1
          ie = tile*imax  
          if ( ufull_g(tile)>0 ) then
            itmp(1:ufull_g(tile)) = pack(typedata(is:ie),upack_g(:,tile))
            if ( minval(itmp(1:ufull_g(tile)))<1 .or. maxval(itmp(1:ufull_g(tile)))>maxtype ) then
              write(6,*) "ERROR: Urban type is out of range"
              stop 
            end if
            f_roof(tile)%volcp(:,i) = paramdata(itmp(1:ufull_g(tile)))
          end if
        end do 
        found = .true.
        exit
      end if
      write(vname,'("roofcond",(I1.1))') i
      if ( trim(paramname)==trim(vname) ) then
        do tile = 1,ntiles
          is = (tile-1)*imax + 1
          ie = tile*imax  
          if ( ufull_g(tile)>0 ) then
            itmp(1:ufull_g(tile)) = pack(typedata(is:ie),upack_g(:,tile))
            if ( minval(itmp(1:ufull_g(tile)))<1 .or. maxval(itmp(1:ufull_g(tile)))>maxtype ) then
              write(6,*) "ERROR: Urban type is out of range"
              stop 
            end if
            f_roof(tile)%lambda(:,i) = paramdata(itmp(1:ufull_g(tile)))
          end if
        end do 
        found = .true.
        exit
      end if
      write(vname,'("wallthick",(I1.1))') i
      if ( trim(paramname)==trim(vname) ) then
        do tile = 1,ntiles
          is = (tile-1)*imax + 1
          ie = tile*imax  
          if ( ufull_g(tile)>0 ) then
            itmp(1:ufull_g(tile)) = pack(typedata(is:ie),upack_g(:,tile))
            if ( minval(itmp(1:ufull_g(tile)))<1 .or. maxval(itmp(1:ufull_g(tile)))>maxtype ) then
              write(6,*) "ERROR: Urban type is out of range"
              stop 
            end if
            f_wall(tile)%depth(:,i) = paramdata(itmp(1:ufull_g(tile)))
          end if
        end do 
        found = .true.
        exit
      end if
      write(vname,'("wallcp",(I1.1))') i
      if ( trim(paramname)==trim(vname) ) then
        do tile = 1,ntiles
          is = (tile-1)*imax + 1
          ie = tile*imax  
          if ( ufull_g(tile)>0 ) then
            itmp(1:ufull_g(tile)) = pack(typedata(is:ie),upack_g(:,tile))
            if ( minval(itmp(1:ufull_g(tile)))<1 .or. maxval(itmp(1:ufull_g(tile)))>maxtype ) then
              write(6,*) "ERROR: Urban type is out of range"
              stop 
            end if
            f_wall(tile)%volcp(:,i) = paramdata(itmp(1:ufull_g(tile)))
          end if
        end do 
        found = .true.
        exit
      end if
      write(vname,'("wallcond",(I1.1))') i
      if ( trim(paramname)==trim(vname) ) then
        do tile = 1,ntiles
          is = (tile-1)*imax + 1
          ie = tile*imax  
          if ( ufull_g(tile)>0 ) then
            itmp(1:ufull_g(tile)) = pack(typedata(is:ie),upack_g(:,tile))
            if ( minval(itmp(1:ufull_g(tile)))<1 .or. maxval(itmp(1:ufull_g(tile)))>maxtype ) then
              write(6,*) "ERROR: Urban type is out of range"
              stop 
            end if
            f_wall(tile)%lambda(:,i) = paramdata(itmp(1:ufull_g(tile)))
          end if
        end do 
        found = .true.
        exit
      end if
      write(vname,'("roadthick",(I1.1))') i
      if ( trim(paramname)==trim(vname) ) then
        do tile = 1,ntiles
          is = (tile-1)*imax + 1
          ie = tile*imax  
          if ( ufull_g(tile)>0 ) then
            itmp(1:ufull_g(tile)) = pack(typedata(is:ie),upack_g(:,tile))
            if ( minval(itmp(1:ufull_g(tile)))<1 .or. maxval(itmp(1:ufull_g(tile)))>maxtype ) then
              write(6,*) "ERROR: Urban type is out of range"
              stop 
            end if
            f_road(tile)%depth(:,i) = paramdata(itmp(1:ufull_g(tile)))
          end if
        end do 
        found = .true.
        exit
      end if
      write(vname,'("roadcp",(I1.1))') i
      if ( trim(paramname)==trim(vname) ) then
        do tile = 1,ntiles
          is = (tile-1)*imax + 1
          ie = tile*imax  
          if ( ufull_g(tile)>0 ) then
            itmp(1:ufull_g(tile)) = pack(typedata(is:ie),upack_g(:,tile))
            if ( minval(itmp(1:ufull_g(tile)))<1 .or. maxval(itmp(1:ufull_g(tile)))>maxtype ) then
              write(6,*) "ERROR: Urban type is out of range"
              stop 
            end if
            f_road(tile)%volcp(:,i) = paramdata(itmp(1:ufull_g(tile)))
          end if
        end do 
        found = .true.
        exit
      end if
      write(vname,'("roadcond",(I1.1))') i
      if ( trim(paramname)==trim(vname) ) then
        do tile = 1,ntiles
          is = (tile-1)*imax + 1
          ie = tile*imax  
          if ( ufull_g(tile)>0 ) then
            itmp(1:ufull_g(tile)) = pack(typedata(is:ie),upack_g(:,tile))
            if ( minval(itmp(1:ufull_g(tile)))<1 .or. maxval(itmp(1:ufull_g(tile)))>maxtype ) then
              write(6,*) "ERROR: Urban type is out of range"
              stop 
            end if
            f_road(tile)%lambda(:,i) = paramdata(itmp(1:ufull_g(tile)))
          end if
        end do 
        found = .true.
        exit
      end if
      write(vname,'("slabthick",(I1.1))') i
      if ( trim(paramname)==trim(vname) ) then
        do tile = 1,ntiles
          is = (tile-1)*imax + 1
          ie = tile*imax  
          if ( ufull_g(tile)>0 ) then
            itmp(1:ufull_g(tile)) = pack(typedata(is:ie),upack_g(:,tile))
            if ( minval(itmp(1:ufull_g(tile)))<1 .or. maxval(itmp(1:ufull_g(tile)))>maxtype ) then
              write(6,*) "ERROR: Urban type is out of range"
              stop 
            end if
            f_slab(tile)%depth(:,i) = paramdata(itmp(1:ufull_g(tile)))
            f_intm(tile)%depth(:,i) = paramdata(itmp(1:ufull_g(tile)))
          end if
        end do 
        found = .true.
        exit
      end if
      write(vname,'("slabcp",(I1.1))') i
      if ( trim(paramname)==trim(vname) ) then
        do tile = 1,ntiles
          is = (tile-1)*imax + 1
          ie = tile*imax  
          if ( ufull_g(tile)>0 ) then
            itmp(1:ufull_g(tile)) = pack(typedata(is:ie),upack_g(:,tile))
            if ( minval(itmp(1:ufull_g(tile)))<1 .or. maxval(itmp(1:ufull_g(tile)))>maxtype ) then
              write(6,*) "ERROR: Urban type is out of range"
              stop 
            end if
            f_slab(tile)%volcp(:,i) = paramdata(itmp(1:ufull_g(tile)))
            f_intm(tile)%volcp(:,i) = paramdata(itmp(1:ufull_g(tile)))
          end if
        end do 
        found = .true.
        exit
      end if
      write(vname,'("slabcond",(I1.1))') i
      if ( trim(paramname)==trim(vname) ) then
        do tile = 1,ntiles
          is = (tile-1)*imax + 1
          ie = tile*imax  
          if ( ufull_g(tile)>0 ) then
            itmp(1:ufull_g(tile)) = pack(typedata(is:ie),upack_g(:,tile))
            if ( minval(itmp(1:ufull_g(tile)))<1 .or. maxval(itmp(1:ufull_g(tile)))>maxtype ) then
              write(6,*) "ERROR: Urban type is out of range"
              stop 
            end if
            f_slab(tile)%lambda(:,i) = paramdata(itmp(1:ufull_g(tile)))
            f_intm(tile)%lambda(:,i) = paramdata(itmp(1:ufull_g(tile)))
          end if
        end do 
        found = .true.
        exit
      end if
    end do
    if ( .not.found ) then
      write(6,*) "ERROR: Unknown UCLEM parameter name ",trim(paramname)
      stop
    end if  
end select

  
! Here we modify the effective canyon geometry to account for in-canyon vegetation tall vegetation
do tile = 1,ntiles
  is = (tile-1)*imax + 1
  ie = tile*imax
  if ( ufull_g(tile)>0 ) then
    do iq = 1,ufull_g(tile)
      fp_coeffbldheight = max(f_g(tile)%bldheight(iq)-6.*cnveg_g(tile)%zo(iq),0.2)/f_g(tile)%bldheight(iq)
      f_g(tile)%effhwratio(iq) = f_g(tile)%hwratio(iq)*fp_coeffbldheight   
    end do  
    
    call init_internal(f_g(tile))
    call init_lwcoeff(f_g(tile),intl_g(tile),ufull_g(tile))
  end if
end do

return
end subroutine uclem_deftype

subroutine uclem_deftype_export_standard(paramdata,paramname,diag)

implicit none

integer, intent(in) :: diag
integer tile, is, ie
real, dimension(ifull), intent(inout) :: paramdata
character(len=*), intent(in) :: paramname
character(len=20) :: vname

if ( diag>=1 ) write(6,*) "Export UCLEM parameter: ",trim(paramname)
if ( .not.uclem_active ) return

do tile = 1,ntiles
  is = (tile-1)*imax + 1
  ie = tile*imax  
  if ( ufull_g(tile)>0 ) then
    call uclem_deftype_export_thread(paramdata(is:ie),paramname,f_g(tile),cnveg_g(tile),f_roof(tile),f_wall(tile), &
                                     f_road(tile),f_slab(tile),upack_g(:,tile),ufull_g(tile),diag)
  end if
end do  

return
end subroutine uclem_deftype_export_standard

subroutine uclem_deftype_export_thread(paramdata,paramname,fp,cnveg,fp_roof,fp_wall, &
                                       fp_road,fp_slab,upack,ufull,diag)

implicit none

integer, intent(in) :: ufull, diag
integer tile, is, ie, i
real, dimension(:), intent(inout) :: paramdata
logical, dimension(:), intent(in) :: upack
logical found
character(len=*), intent(in) :: paramname
character(len=20) :: vname
type(fparmdata), intent(in) :: fp
type(vegdata), intent(in) :: cnveg
type(facetparams), intent(in) :: fp_roof, fp_wall, fp_road, fp_slab

if ( diag>=2 ) write(6,*) "THREAD: Export UCLEM parameter: ",trim(paramname)
if ( ufull==0 ) return

select case(paramname)
  case('bldheight')
    paramdata=unpack(fp%bldheight,upack,paramdata)
  case('hwratio')
    paramdata=unpack(fp%hwratio,upack,paramdata)
  case('sigvegc')
    paramdata=unpack(cnveg%sigma*(1.-fp%sigmabld),upack,paramdata)
  case('sigmabld')
    paramdata=unpack(fp%sigmabld,upack,paramdata)
  case('industryfg')
    paramdata=unpack(fp%industryfg,upack,paramdata)
  case('trafficfg')
    paramdata=unpack(fp%trafficfg,upack,paramdata)
  case('roofalpha')
    paramdata=unpack(fp_roof%alpha,upack,paramdata)
  case('wallalpha')
    paramdata=unpack(fp_wall%alpha,upack,paramdata)
  case('roadalpha')
    paramdata=unpack(fp_road%alpha,upack,paramdata)
  case('roofemiss')
    paramdata=unpack(fp_roof%emiss,upack,paramdata)
  case('wallemiss')
    paramdata=unpack(fp_wall%emiss,upack,paramdata)
  case('roademiss')
    paramdata=unpack(fp_road%emiss,upack,paramdata)
  case('vegalphac')
    paramdata=unpack(cnveg%alpha,upack,paramdata)
  case('zovegc')
    paramdata=unpack(cnveg%zo,upack,paramdata)
  case('infilach')
    paramdata=unpack(fp%infilach,upack,paramdata)
  case('intgains')
    paramdata=unpack(fp%intgains_flr,upack,paramdata)
  case default
    found = .false.
    do i = 1,4
      write(vname,'("roofthick",(I1.1))') i
      if ( trim(paramname)==trim(vname) ) then
        paramdata=unpack(fp_roof%depth(:,i),upack,paramdata)
        found = .true.
        exit
      end if
      write(vname,'("roofcp",(I1.1))') i
      if ( trim(paramname)==trim(vname) ) then
        paramdata=unpack(fp_roof%volcp(:,i),upack,paramdata)
        found = .true.
        exit
      end if
      write(vname,'("roofcond",(I1.1))') i
      if ( trim(paramname)==trim(vname) ) then
        paramdata=unpack(fp_roof%lambda(:,i),upack,paramdata)
        found = .true.
        exit
      end if
      write(vname,'("wallthick",(I1.1))') i
      if ( trim(paramname)==trim(vname) ) then
        paramdata=unpack(fp_wall%depth(:,i),upack,paramdata)
        found = .true.
        exit
      end if
      write(vname,'("wallcp",(I1.1))') i
      if ( trim(paramname)==trim(vname) ) then
        paramdata=unpack(fp_wall%volcp(:,i),upack,paramdata)
        found = .true.
        exit
      end if
      write(vname,'("wallcond",(I1.1))') i
      if ( trim(paramname)==trim(vname) ) then
        paramdata=unpack(fp_wall%lambda(:,i),upack,paramdata)
        found = .true.
        exit
      end if
      write(vname,'("roadthick",(I1.1))') i
      if ( trim(paramname)==trim(vname) ) then
        paramdata=unpack(fp_road%depth(:,i),upack,paramdata)
        found = .true.
        exit
      end if
      write(vname,'("roadcp",(I1.1))') i
      if ( trim(paramname)==trim(vname) ) then
        paramdata=unpack(fp_road%volcp(:,i),upack,paramdata)
        found = .true.
        exit
      end if
      write(vname,'("roadcond",(I1.1))') i
      if ( trim(paramname)==trim(vname) ) then
        paramdata=unpack(fp_road%lambda(:,i),upack,paramdata)
        found = .true.
        exit
      end if
      write(vname,'("slabthick",(I1.1))') i
      if ( trim(paramname)==trim(vname) ) then
        paramdata=unpack(fp_slab%depth(:,i),upack,paramdata)
        found = .true.
        exit
      end if
      write(vname,'("slabcp",(I1.1))') i
      if ( trim(paramname)==trim(vname) ) then
        paramdata=unpack(fp_slab%volcp(:,i),upack,paramdata)
        found = .true.
        exit
      end if
      write(vname,'("slabcond",(I1.1))') i
      if ( trim(paramname)==trim(vname) ) then
        paramdata=unpack(fp_slab%lambda(:,i),upack,paramdata)
        found = .true.
        exit
      end if
    end do
    if ( .not.found ) then
      write(6,*) "ERROR: Unknown UCLEM parameter name ",trim(paramname)
      stop
    end if  
end select

return
end subroutine uclem_deftype_export_thread

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! general version of uclem_save

subroutine uclem_saved_3(urban,mode,ifrac,diag,rawtemp)

implicit none

integer, intent(in) :: ifrac, diag
integer ii, tile, is, ie
real(kind=8), dimension(ifull), intent(inout) :: urban
real(kind=8) urbtempadj
logical, intent(in), optional :: rawtemp
logical rawmode
character(len=*), intent(in) :: mode
character(len=10) :: teststr

if (diag>=1) write(6,*) "Save UCLEM state array"
if (.not.uclem_active) return

rawmode = .false.
if ( present(rawtemp) ) then
  rawmode = rawtemp
end if

if ( rawmode ) then
  urbtempadj = 0._8
else
  urbtempadj = real(urbtemp,8)
end if

do ii = 0,nl
  write(teststr,'("rooftemp",I1.1)') ii+1
  if ( trim(teststr)==trim(mode) ) then
    do tile = 1,ntiles
      is = (tile-1)*imax + 1
      ie = tile*imax
      if ( ufull_g(tile)>0 ) then
        urban(is:ie)=unpack(roof_g(ifrac,tile)%nodetemp(:,ii)+urbtempadj,upack_g(:,tile),urban(is:ie))
      end if
    end do
    return
  end if
  write(teststr,'("walletemp",I1.1)') ii+1
  if ( trim(teststr)==trim(mode) ) then
    do tile = 1,ntiles
      is = (tile-1)*imax + 1
      ie = tile*imax
      if ( ufull_g(tile)>0 ) then
        urban(is:ie)=unpack(walle_g(ifrac,tile)%nodetemp(:,ii)+urbtempadj,upack_g(:,tile),urban(is:ie))
      end if
    end do
    return
  end if
  write(teststr,'("wallwtemp",I1.1)') ii+1
  if ( trim(teststr)==trim(mode) ) then
    do tile = 1,ntiles
      is = (tile-1)*imax + 1
      ie = tile*imax
      if ( ufull_g(tile)>0 ) then
        urban(is:ie)=unpack(wallw_g(ifrac,tile)%nodetemp(:,ii)+urbtempadj,upack_g(:,tile),urban(is:ie))
      end if
    end do
    return
  end if
  write(teststr,'("roadtemp",I1.1)') ii+1
  if ( trim(teststr)==trim(mode) ) then
    do tile = 1,ntiles
      is = (tile-1)*imax + 1
      ie = tile*imax
      if ( ufull_g(tile)>0 ) then
        urban(is:ie)=unpack(road_g(ifrac,tile)%nodetemp(:,ii)+urbtempadj,upack_g(:,tile),urban(is:ie))
      end if
    end do
    return
  end if
  write(teststr,'("slabtemp",I1.1)') ii+1
  if ( trim(teststr)==trim(mode) ) then
    do tile = 1,ntiles
      is = (tile-1)*imax + 1
      ie = tile*imax
      if ( ufull_g(tile)>0 ) then
        urban(is:ie)=unpack(slab_g(ifrac,tile)%nodetemp(:,ii)+urbtempadj,upack_g(:,tile),urban(is:ie))
      end if
    end do
    return
  end if  
  write(teststr,'("intmtemp",I1.1)') ii+1
  if ( trim(teststr)==trim(mode) ) then
    do tile = 1,ntiles
      is = (tile-1)*imax + 1
      ie = tile*imax
      if ( ufull_g(tile)>0 ) then
        urban(is:ie)=unpack(intm_g(ifrac,tile)%nodetemp(:,ii)+urbtempadj,upack_g(:,tile),urban(is:ie))
      end if
    end do
    return
  end if   
end do  

select case(mode)
  case("roomtemp")
    do tile = 1,ntiles
      is = (tile-1)*imax + 1
      ie = tile*imax
      if ( ufull_g(tile)>0 ) then
        urban(is:ie)=unpack(room_g(ifrac,tile)%nodetemp(:,1)+urbtempadj,upack_g(:,tile),urban(is:ie))    
      end if
    end do
    return
end select

write(6,*) "ERROR: Unknown mode for uclem_saved ",trim(mode)
stop

return
end subroutine uclem_saved_3

subroutine uclem_saved_2(urban,mode,ifrac,diag)

implicit none

integer, intent(in) :: ifrac, diag
integer ii, tile, is, ie
real, dimension(ifull), intent(inout) :: urban
character(len=*), intent(in) :: mode
character(len=10) :: teststr

if (diag>=1) write(6,*) "Save UCLEM state array"
if (.not.uclem_active) return

select case(mode)
  case("canyonsoilmoisture")
    do tile = 1,ntiles
      is = (tile-1)*imax + 1
      ie = tile*imax
      if ( ufull_g(tile)>0 ) then
        urban(is:ie)=unpack(rdhyd_g(ifrac,tile)%soilwater,upack_g(:,tile),urban(is:ie))  
      end if
    end do
    return
  case("roofsoilmoisture")
    do tile = 1,ntiles
      is = (tile-1)*imax + 1
      ie = tile*imax
      if ( ufull_g(tile)>0 ) then
        urban(is:ie)=unpack(rfhyd_g(ifrac,tile)%soilwater,upack_g(:,tile),urban(is:ie))    
      end if
    end do
    return
  case("roadsurfacewater")
    do tile = 1,ntiles
      is = (tile-1)*imax + 1
      ie = tile*imax
      if ( ufull_g(tile)>0 ) then
        urban(is:ie)=unpack(rdhyd_g(ifrac,tile)%surfwater,upack_g(:,tile),urban(is:ie))    
      end if
    end do
    return
  case("roofsurfacewater")
    do tile = 1,ntiles
      is = (tile-1)*imax + 1
      ie = tile*imax
      if ( ufull_g(tile)>0 ) then
        urban(is:ie)=unpack(rfhyd_g(ifrac,tile)%surfwater,upack_g(:,tile),urban(is:ie))    
      end if
    end do
    return
  case("canyonleafwater")
    do tile = 1,ntiles
      is = (tile-1)*imax + 1
      ie = tile*imax
      if ( ufull_g(tile)>0 ) then
        urban(is:ie)=unpack(rdhyd_g(ifrac,tile)%leafwater,upack_g(:,tile),urban(is:ie))    
      end if
    end do
    return
  case("roofleafwater")
    do tile = 1,ntiles
      is = (tile-1)*imax + 1
      ie = tile*imax
      if ( ufull_g(tile)>0 ) then
        urban(is:ie)=unpack(rfhyd_g(ifrac,tile)%leafwater,upack_g(:,tile),urban(is:ie))    
      end if
    end do
    return
  case("roadsnowdepth")
    do tile = 1,ntiles
      is = (tile-1)*imax + 1
      ie = tile*imax
      if ( ufull_g(tile)>0 ) then
        urban(is:ie)=unpack(rdhyd_g(ifrac,tile)%snow,upack_g(:,tile),urban(is:ie))
      end if
    end do
    return
  case("roofsnowdepth")
    do tile = 1,ntiles
      is = (tile-1)*imax + 1
      ie = tile*imax
      if ( ufull_g(tile)>0 ) then
        urban(is:ie)=unpack(rfhyd_g(ifrac,tile)%snow,upack_g(:,tile),urban(is:ie))    
      end if
    end do
    return
  case("roadsnowdensity")
    do tile = 1,ntiles
      is = (tile-1)*imax + 1
      ie = tile*imax
      if ( ufull_g(tile)>0 ) then
        urban(is:ie)=unpack(rdhyd_g(ifrac,tile)%den,upack_g(:,tile),urban(is:ie))    
      end if
    end do
    return
  case("roofsnowdensity")
    do tile = 1,ntiles
      is = (tile-1)*imax + 1
      ie = tile*imax
      if ( ufull_g(tile)>0 ) then
        urban(is:ie)=unpack(rfhyd_g(ifrac,tile)%den,upack_g(:,tile),urban(is:ie))    
      end if
    end do
    return
  case("roadsnowalbedo")
    do tile = 1,ntiles
      is = (tile-1)*imax + 1
      ie = tile*imax
      if ( ufull_g(tile)>0 ) then
        urban(is:ie)=unpack(rdhyd_g(ifrac,tile)%snowalpha,upack_g(:,tile),urban(is:ie))    
      end if
    end do
    return
  case("roofsnowalbedo")
    do tile = 1,ntiles
      is = (tile-1)*imax + 1
      ie = tile*imax
      if ( ufull_g(tile)>0 ) then
        urban(is:ie)=unpack(rfhyd_g(ifrac,tile)%snowalpha,upack_g(:,tile),urban(is:ie))    
      end if
    end do
    return
end select

write(6,*) "ERROR: Unknown mode for uclem_saved ",trim(mode)
stop

return
end subroutine uclem_saved_2

subroutine uclem_energy_standard(o_data,mode,diag)

implicit none

integer, intent(in) :: diag
integer tile, is, ie
real, dimension(:), intent(inout) :: o_data
character(len=*), intent(in) :: mode

if ( diag>=1 ) write(6,*) "Extract energy output"
if (.not.uclem_active) return

do tile = 1,ntiles
  is = (tile-1)*imax + 1
  ie = tile*imax
  if ( ufull_g(tile)>0 ) then
    call uclem_energy_thread(o_data(is:ie),mode,diag,f_g(tile),p_g(:,tile),        & 
                             rdhyd_g(:,tile),rfhyd_g(:,tile),upack_g(:,tile),ufull_g(tile))
  end if
end do

return
end subroutine uclem_energy_standard

subroutine uclem_energy_thread(o_data,mode,diag,fp,pd,rdhyd,rfhyd,upack,ufull)

implicit none

integer, intent(in) :: ufull, diag
integer ifrac
real, dimension(:), intent(inout) :: o_data
real, dimension(ufull) :: ctmp, dtmp
character(len=*), intent(in) :: mode
logical, dimension(:), intent(in) :: upack
type(fparmdata), intent(in) :: fp
type(pdiagdata), dimension(nfrac), intent(in) :: pd
type(hydrodata), dimension(nfrac), intent(in) :: rdhyd, rfhyd

if ( diag>=2 ) write(6,*) "THREAD: Extract energy output"
if ( ufull==0 ) return

select case(mode)
  case("anthropogenic")
    dtmp = 0.
    do ifrac = 1,nfrac
      dtmp = dtmp + (pd(ifrac)%bldheat+pd(ifrac)%bldcool+pd(ifrac)%traf+fp%industryfg+pd(ifrac)%intgains_full) &
                    *pd(ifrac)%frac_sigma
    end do  
    ctmp = pack(o_data, upack)
    ctmp = (1.-fp%sigmau)*ctmp + fp%sigmau*dtmp
    o_data = unpack(ctmp, upack, o_data)
  case("elecgas")
    dtmp = 0.
    do ifrac = 1,nfrac   
      dtmp = dtmp + real(pd(ifrac)%bldheat+pd(ifrac)%bldcool+pd(ifrac)%intgains_full)*pd(ifrac)%frac_sigma
    end do  
    ctmp = pack(o_data, upack)
    ctmp = (1.-fp%sigmau)*ctmp + fp%sigmau*dtmp
    o_data = unpack(ctmp, upack, o_data)
  case("heating")
    dtmp = 0.
    do ifrac = 1,nfrac
      dtmp = dtmp + real(pd(ifrac)%bldheat)*pd(ifrac)%frac_sigma
    end do  
    ctmp = pack(o_data, upack)
    ctmp = (1.-fp%sigmau)*ctmp + fp%sigmau*dtmp
    o_data = unpack(ctmp, upack, o_data)
  case("cooling")
    dtmp = 0.
    do ifrac = 1,nfrac
      dtmp = dtmp + real(pd(ifrac)%bldcool)*pd(ifrac)%frac_sigma
    end do  
    ctmp = pack(o_data, upack)
    ctmp = (1.-fp%sigmau)*ctmp + fp%sigmau*dtmp
    o_data = unpack(ctmp, upack, o_data)
  case("storage")
    dtmp = 0.
    do ifrac = 1,nfrac
      dtmp = dtmp + real(pd(ifrac)%storage_flux)*pd(ifrac)%frac_sigma
    end do  
    ctmp = pack(o_data, upack)
    ctmp = (1.-fp%sigmau)*ctmp + fp%sigmau*dtmp
    o_data = unpack(ctmp, upack, o_data)
  case default
    write(6,*) "ERROR: Unknown uclem_energy mode ",trim(mode)
    stop
end select    

return
end subroutine uclem_energy_thread

subroutine uclem_misc_standard(o_data,mode,diag)

implicit none

integer, intent(in) :: diag
integer tile, is, ie
real, dimension(:), intent(inout) :: o_data
character(len=*), intent(in) :: mode

if ( diag>=1 ) write(6,*) "Extract energy output"
if (.not.uclem_active) return

do tile = 1,ntiles
  is = (tile-1)*imax + 1
  ie = tile*imax
  if ( ufull_g(tile)>0 ) then
    call uclem_misc_thread(o_data(is:ie),mode,diag,room_g(:,tile),f_g(tile),p_g(:,tile),upack_g(:,tile),ufull_g(tile))
  end if
end do

return
end subroutine uclem_misc_standard

subroutine uclem_misc_thread(o_data,mode,diag,room,fp,pd,upack,ufull)

implicit none

integer, intent(in) :: ufull, diag
integer ifrac
real, dimension(:), intent(inout) :: o_data
real, dimension(ufull) :: dtmp
character(len=*), intent(in) :: mode
logical, dimension(:), intent(in) :: upack
type(fparmdata), intent(in) :: fp
type(facetdata), dimension(nfrac), intent(in) :: room
type(pdiagdata), dimension(nfrac), intent(in) :: pd

if ( diag>=2 ) write(6,*) "THREAD: Extract energy output"
if ( ufull==0 ) return

select case(mode)
  case("emissivity")
    dtmp = 0.
    do ifrac = 1,nfrac
      dtmp = dtmp + pd(ifrac)%emiss*pd(ifrac)%frac_sigma
    end do  
    o_data = unpack(dtmp, upack, o_data)
  case("snowt")
    dtmp = 0.
    do ifrac = 1,nfrac
      dtmp = dtmp + pd(ifrac)%snowt*pd(ifrac)%frac_sigma  
    end do    
    o_data = unpack(dtmp, upack, o_data)
  case("vegt")
    dtmp = 0.
    do ifrac = 1,nfrac
      dtmp = dtmp + pd(ifrac)%vegt*pd(ifrac)%frac_sigma  
    end do    
    o_data = unpack(dtmp, upack, o_data)
  case("tairbuilding")
    dtmp = 0.
    do ifrac = 1,nfrac
      dtmp = dtmp + (real(room(ifrac)%nodetemp(:,1))+urbtemp)*pd(ifrac)%frac_sigma  
    end do    
    o_data = unpack(dtmp, upack, o_data)
  case("salbedo")
    dtmp = 0.
    do ifrac = 1,nfrac
      dtmp = dtmp + pd(ifrac)%salbedo*pd(ifrac)%frac_sigma  
    end do    
    o_data = unpack(dtmp, upack, o_data)
  case("calbedo")
    dtmp = 0.
    do ifrac = 1,nfrac
      dtmp = dtmp + pd(ifrac)%calbedo*pd(ifrac)%frac_sigma  
    end do    
    o_data = unpack(dtmp, upack, o_data)
  case("urbanlai")
    dtmp = 0.
    do ifrac = 1,nfrac
      dtmp = dtmp + pd(ifrac)%ulai*pd(ifrac)%frac_sigma  
    end do
    o_data = unpack(dtmp, upack, o_data)
  case("taircanyon")
    dtmp = 0.
    do ifrac = 1,nfrac
      dtmp = dtmp + (pd(ifrac)%taircanyon+urbtemp)*pd(ifrac)%frac_sigma  
    end do    
    o_data = unpack(dtmp, upack, o_data)      
  case default
    write(6,*) "ERROR: Unknown uclem_misc mode ",trim(mode)
    stop
end select    

return
end subroutine uclem_misc_thread

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine blends urban momentum and heat roughness lengths
! (This version neglects the displacement height (e.g., for CCAM))
!

subroutine uclem_zo_standard(zom,zoh,zoq,diag,raw)

implicit none

integer, intent(in) :: diag
integer tile, is, ie
real, dimension(:), intent(inout) :: zom, zoh, zoq
logical, intent(in), optional :: raw
logical mode

if ( diag>=1 ) write(6,*) "Calculate urban roughness lengths"
if (.not.uclem_active) return

mode=.false.
if (present(raw)) mode=raw

do tile = 1,ntiles
  is = (tile-1)*imax + 1
  ie = tile*imax
  if ( ufull_g(tile)>0 ) then
    call uclem_zo_thread(zom(is:ie),zoh(is:ie),zoq(is:ie),diag,p_g(:,tile),f_g(tile),upack_g(:,tile),ufull_g(tile),raw=mode)
  end if
end do

return
end subroutine uclem_zo_standard
                             
subroutine uclem_zo_thread(zom,zoh,zoq,diag,pd,fp,upack,ufull,raw)

implicit none

integer, intent(in) :: ufull, diag
integer ifrac
real, dimension(:), intent(inout) :: zom, zoh, zoq
real, dimension(ufull) :: workb,workc,workd,zmtmp,zhtmp,zqtmp
real, dimension(ufull) :: sumlzom, sumlzoh, sumlzoq, sumcndzmin
real, parameter :: zr=1.e-15 ! limits minimum roughness length for heat
logical, intent(in), optional :: raw
logical mode
logical, dimension(:), intent(in) :: upack
type(fparmdata), intent(in) :: fp
type(pdiagdata), dimension(nfrac), intent(in) :: pd

if ( diag>=2 ) write(6,*) "THREAD: Calculate urban roughness length"
if ( ufull==0 ) return

mode=.false.
if (present(raw)) mode=raw

sumlzom = 0.
sumlzoh = 0.
sumlzoq = 0.
sumcndzmin = 0.
do ifrac = 1,nfrac
  sumlzom = sumlzom + pd(ifrac)%frac_sigma/pd(ifrac)%lzom**2
  sumlzoh = sumlzoh + pd(ifrac)%frac_sigma/(pd(ifrac)%lzom*pd(ifrac)%lzoh)
  sumlzoq = sumlzoq + pd(ifrac)%frac_sigma/(pd(ifrac)%lzom*pd(ifrac)%lzoh)
  sumcndzmin = sumcndzmin + pd(ifrac)%frac_sigma*pd(ifrac)%cndzmin
end do

if (mode) then
  sumlzom = sqrt(sumlzom)
  sumlzoh = sumlzoh/sumlzom
  sumlzoq = sumlzoq/sumlzom
  zom=unpack(sumcndzmin*exp(-1./sumlzom),upack,zom)
  zoh=unpack(max(sumcndzmin*exp(-1./sumlzoh),zr),upack,zoh)
  zoq=unpack(max(sumcndzmin*exp(-1./sumlzoq),zr),upack,zoq)
else 
  ! evaluate at canyon displacement height (really the atmospheric model should provide a displacement height)
  zmtmp=pack(zom,upack)
  zhtmp=pack(zoh,upack)
  zqtmp=pack(zoq,upack)
  workb=sqrt((1.-fp%sigmau)/log(sumcndzmin/zmtmp)**2+fp%sigmau*sumlzom)
  workc=(1.-fp%sigmau)/(log(sumcndzmin/zmtmp)*log(sumcndzmin/zhtmp))+fp%sigmau*sumlzoh
  workc=workc/workb
  workd=(1.-fp%sigmau)/(log(sumcndzmin/zmtmp)*log(sumcndzmin/zqtmp))+fp%sigmau*sumlzoq
  workd=workd/workb
  workb=sumcndzmin*exp(-1./workb)
  workc=max(sumcndzmin*exp(-1./workc),zr)
  workd=max(sumcndzmin*exp(-1./workd),zr)
  zom=unpack(workb,upack,zom)
  zoh=unpack(workc,upack,zoh)
  zoq=unpack(workd,upack,zoq)
  if (minval(workc)<=zr) write(6,*) "WARN: minimum zoh reached"
end if

return
end subroutine uclem_zo_thread

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine blends the urban drag coeff
!

subroutine uclem_cd_standard(cduv,cdtq,diag,raw)
 
implicit none
 
integer, intent(in) :: diag
integer tile, is, ie
real, dimension(:), intent(inout) :: cduv, cdtq
logical, intent(in), optional :: raw
logical outmode

if ( diag>=1 ) write(6,*) "Calculate urban drag coeff"
if (.not.uclem_active) return

outmode=.false.
if (present(raw)) outmode=raw

do tile = 1,ntiles
  is = (tile-1)*imax + 1
  ie = tile*imax
  if ( ufull_g(tile)>0 ) then
    call uclem_cd_thread(cduv(is:ie),cdtq(is:ie),diag,p_g(:,tile),f_g(tile),upack_g(:,tile),ufull_g(tile),raw=outmode)
  end if
end do

return
end subroutine uclem_cd_standard

subroutine uclem_cd_thread(cduv,cdtq,diag,pd,fp,upack,ufull,raw)
 
implicit none
 
integer, intent(in) :: ufull, diag
integer ifrac
real, dimension(:), intent(inout) :: cduv, cdtq
real, dimension(ufull) :: ctmp, dtmp
logical, intent(in), optional :: raw
logical outmode
logical, dimension(:), intent(in) :: upack
type(fparmdata), intent(in) :: fp
type(pdiagdata), dimension(nfrac), intent(in) :: pd
 
if (diag>=2) write(6,*) "THREAD: Calculate urban drag coeff"
if ( ufull==0 ) return
 
outmode=.false.
if (present(raw)) outmode=raw
 

dtmp = 0.
do ifrac = 1,nfrac
  dtmp = dtmp + pd(ifrac)%cduv*pd(ifrac)%frac_sigma  
end do
if ( outmode ) then
  ctmp=dtmp
else
  ctmp=pack(cduv,upack)  
  ctmp=(1.-fp%sigmau)*ctmp+fp%sigmau*dtmp
end if
cduv=unpack(ctmp,upack,cduv)
 

dtmp = 0.
do ifrac = 1,nfrac
  dtmp = dtmp + pd(ifrac)%cdtq*pd(ifrac)%frac_sigma  
end do
if ( outmode ) then
  ctmp=dtmp
else
  ctmp=pack(cdtq,upack)  
  ctmp=(1.-fp%sigmau)*ctmp+fp%sigmau*dtmp
end if
cdtq=unpack(ctmp,upack,cdtq)
 
return
end subroutine uclem_cd_thread

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine is for hydrological outputs
!
 
subroutine uclem_hydro_standard(hydroout,mode,diag)

implicit none
 
integer, intent(in) :: diag
integer tile, is, ie
real, dimension(:), intent(inout) :: hydroout
character(len=*), intent(in) :: mode

if ( diag>=1 ) write(6,*) "Calculate hydrological outputs"
if (.not.uclem_active) return

do tile = 1,ntiles
  is = (tile-1)*imax + 1
  ie = tile*imax
  if ( ufull_g(tile)>0 ) then
    call uclem_hydro_thread(hydroout(is:ie),mode,diag,p_g(:,tile),f_g(tile),upack_g(:,tile),ufull_g(tile))
  end if
end do

return
end subroutine uclem_hydro_standard

subroutine uclem_hydro_thread(hydroout,mode,diag,pd,fp,upack,ufull)
 
implicit none
 
integer, intent(in) :: ufull, diag
integer ifrac
real, dimension(:), intent(inout) :: hydroout
real, dimension(ufull) :: ctmp, dtmp
character(len=*), intent(in) :: mode
logical, dimension(:), intent(in) :: upack
type(fparmdata), intent(in) :: fp
type(pdiagdata), dimension(nfrac), intent(in) :: pd
 
if ( diag>=2 ) write(6,*) "THREAD: Calculate hydrological outputs"
if ( ufull==0 ) return
 
select case(mode)
  case("snowmelt")
    dtmp = 0.
    do ifrac = 1,nfrac
      dtmp = dtmp + pd(ifrac)%snowmelt*pd(ifrac)%frac_sigma  
    end do    
    ctmp=pack(hydroout,upack)
    ctmp=(1.-fp%sigmau)*ctmp+fp%sigmau*dtmp
    hydroout=unpack(ctmp,upack,hydroout)
  case("soilmoisture")
    dtmp = 0.
    do ifrac = 1,nfrac
      dtmp = dtmp + pd(ifrac)%soilmoist*pd(ifrac)%frac_sigma
    end do  
    ctmp = pack(hydroout, upack)
    ctmp =(1.-fp%sigmau)*ctmp+fp%sigmau*dtmp
    hydroout=unpack(ctmp,upack,hydroout)
  case("soilwet")
    dtmp = 0.
    do ifrac = 1,nfrac
      dtmp = dtmp + pd(ifrac)%soilwetness*pd(ifrac)%frac_sigma
    end do  
    ctmp = pack(hydroout, upack)
    ctmp =(1.-fp%sigmau)*ctmp+fp%sigmau*dtmp
    hydroout=unpack(ctmp,upack,hydroout)
  case("delsoilmoist")
    dtmp = 0.
    do ifrac = 1,nfrac
      dtmp = dtmp + pd(ifrac)%delsoilmoist*pd(ifrac)%frac_sigma
    end do  
    ctmp = pack(hydroout, upack)
    ctmp = (1.-fp%sigmau)*ctmp + fp%sigmau*dtmp
    hydroout=unpack(ctmp,upack,hydroout)
  case("rootmoisture")
    dtmp = 0.
    do ifrac = 1,nfrac
      dtmp = dtmp + pd(ifrac)%rootmoistc*pd(ifrac)%frac_sigma
    end do  
    ctmp = pack(hydroout, upack)
    ctmp =(1.-fp%sigmau)*ctmp+fp%sigmau*dtmp
    hydroout=unpack(ctmp,upack,hydroout)
  case("irrigation")
    dtmp = 0.
    do ifrac = 1,nfrac
      dtmp = dtmp + pd(ifrac)%irrig*pd(ifrac)%frac_sigma
    end do  
    ctmp = pack(hydroout, upack)
    ctmp = (1.-fp%sigmau)*ctmp + fp%sigmau*dtmp
    hydroout=unpack(ctmp,upack,hydroout)
  case("transpirationveg")
    dtmp = 0.
    do ifrac = 1,nfrac
      dtmp = dtmp + pd(ifrac)%transveg*pd(ifrac)%frac_sigma
    end do  
    ctmp = pack(hydroout, upack)
    ctmp = (1.-fp%sigmau)*ctmp + fp%sigmau*dtmp
    hydroout = unpack(ctmp, upack, hydroout)
  case("aeroconductionveg")
    dtmp = 0.
    do ifrac = 1,nfrac
      dtmp = dtmp + pd(ifrac)%acond_vegw*pd(ifrac)%frac_sigma
    end do  
    ctmp = pack(hydroout, upack)
    ctmp = (1.-fp%sigmau)*ctmp+fp%sigmau*dtmp
    hydroout=unpack(ctmp,upack,hydroout)
  case("delintercept")
    dtmp = 0.
    do ifrac = 1,nfrac
      dtmp = dtmp + pd(ifrac)%delintercept*pd(ifrac)%frac_sigma  
    end do    
    ctmp=pack(hydroout,upack)
    ctmp=(1.-fp%sigmau)*ctmp+fp%sigmau*dtmp
    hydroout=unpack(ctmp,upack,hydroout)
  case("swe")
    dtmp = 0.
    do ifrac = 1,nfrac
      dtmp = dtmp + pd(ifrac)%swe*pd(ifrac)%frac_sigma  
    end do    
    ctmp=pack(hydroout,upack)
    ctmp=(1.-fp%sigmau)*ctmp+fp%sigmau*dtmp
    hydroout=unpack(ctmp,upack,hydroout)    
  case("surfstor")
    dtmp = 0.
    do ifrac = 1,nfrac
      dtmp = dtmp + pd(ifrac)%surfstor*pd(ifrac)%frac_sigma  
    end do    
    ctmp=pack(hydroout,upack)
    ctmp=(1.-fp%sigmau)*ctmp+fp%sigmau*dtmp
    hydroout=unpack(ctmp,upack,hydroout)   
  case("snowfrac")
    dtmp = 0.
    do ifrac = 1,nfrac
      dtmp = dtmp + pd(ifrac)%snowfrac*pd(ifrac)%frac_sigma  
    end do    
    ctmp=pack(hydroout,upack)
    ctmp=(1.-fp%sigmau)*ctmp+fp%sigmau*dtmp
    hydroout=unpack(ctmp,upack,hydroout) 
  case("delswe")
    dtmp = 0.
    do ifrac = 1,nfrac
      dtmp = dtmp + pd(ifrac)%delswe*pd(ifrac)%frac_sigma  
    end do    
    ctmp=pack(hydroout,upack)
    ctmp=(1.-fp%sigmau)*ctmp+fp%sigmau*dtmp
    hydroout=unpack(ctmp,upack,hydroout)
  case("surfrunoff")
    dtmp = 0.
    do ifrac = 1,nfrac
      dtmp = dtmp + pd(ifrac)%surfrunoff*pd(ifrac)%frac_sigma
    end do  
    ctmp = pack(hydroout, upack)
    ctmp = (1.-fp%sigmau)*ctmp + fp%sigmau*dtmp
    hydroout = unpack(ctmp, upack, hydroout)
  case("soilwaterfrac")
    dtmp = 0.
    do ifrac = 1,nfrac
      dtmp = dtmp + pd(ifrac)%soilwater*pd(ifrac)%frac_sigma
    end do  
    ctmp = pack(hydroout, upack)
    ctmp = (1.-fp%sigmau)*ctmp + fp%sigmau*dtmp
    hydroout = unpack(ctmp, upack, hydroout)
  case default
    write(6,*) "ERROR: Unknown uclem_hydro mode ",trim(mode)
    stop
end select
 
return
end subroutine uclem_hydro_thread

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine is for temperature outputs
!

subroutine uclem_avetemp(tempout,mode,diag)

implicit none

integer, intent(in) :: diag
integer tile, is, ie
real, dimension(:), intent(inout) :: tempout
character(len=*), intent(in) :: mode

if ( diag>=1 ) write(6,*) "Calculate temperature outputs"
if (.not.uclem_active) return

do tile = 1,ntiles
  is = (tile-1)*imax + 1
  ie = tile*imax
  if ( ufull_g(tile)>0 ) then
    call uclem_avetemp_thread(tempout(is:ie),mode,diag,roof_g(:,tile),walle_g(:,tile),wallw_g(:,tile),road_g(:,tile), &
                              slab_g(:,tile),intm_g(:,tile),p_g(:,tile),upack_g(:,tile),ufull_g(tile))
  end if
end do

return
end subroutine uclem_avetemp

subroutine uclem_avetemp_thread(tempout,mode,diag,roof,walle,wallw,road,slab, &
                                intm,pd,upack,ufull)

implicit none

integer, intent(in) :: ufull, diag
integer ifrac, ii
real, dimension(:), intent(inout) :: tempout
real, dimension(ufull) :: ctmp
logical, dimension(:), intent(in) :: upack
character(len=*), intent(in) :: mode
character(len=10) :: teststr
type(facetdata), dimension(nfrac), intent(in) :: roof, walle, wallw, road, slab, intm
type(pdiagdata), dimension(nfrac), intent(in) :: pd

if ( diag>=2 ) write(6,*) "THREAD: Calculate temperature outputs"
if ( ufull==0 ) return

do ii = 0,nl
  write(teststr,'("rooftemp",I1.1)') ii+1
  if ( trim(teststr)==trim(mode) ) then
    ctmp = 0.
    do ifrac = 1,nfrac
      ctmp = ctmp + real(roof(ifrac)%nodetemp(:,ii))*pd(ifrac)%frac_sigma
    end do
    tempout = unpack( ctmp+urbtemp, upack, tempout )
    return
  end if  
  write(teststr,'("walletemp",I1.1)') ii+1
  if ( trim(teststr)==trim(mode) ) then
    ctmp = 0.
    do ifrac = 1,nfrac
      ctmp = ctmp + real(walle(ifrac)%nodetemp(:,ii))*pd(ifrac)%frac_sigma
    end do
    tempout = unpack( ctmp+urbtemp, upack, tempout )
    return
  end if 
  write(teststr,'("wallwtemp",I1.1)') ii+1
  if ( trim(teststr)==trim(mode) ) then
    ctmp = 0.
    do ifrac = 1,nfrac
      ctmp = ctmp + real(wallw(ifrac)%nodetemp(:,ii))*pd(ifrac)%frac_sigma
    end do
    tempout = unpack( ctmp+urbtemp, upack, tempout )
    return
  end if
  write(teststr,'("roadtemp",I1.1)') ii+1
  if ( trim(teststr)==trim(mode) ) then
    ctmp = 0.
    do ifrac = 1,nfrac
      ctmp = ctmp + real(road(ifrac)%nodetemp(:,ii))*pd(ifrac)%frac_sigma
    end do
    tempout = unpack( ctmp+urbtemp, upack, tempout )
    return
  end if
  write(teststr,'("slabtemp",I1.1)') ii+1
  if ( trim(teststr)==trim(mode) ) then
    ctmp = 0.
    do ifrac = 1,nfrac
      ctmp = ctmp + real(slab(ifrac)%nodetemp(:,ii))*pd(ifrac)%frac_sigma
    end do
    tempout = unpack( ctmp+urbtemp, upack, tempout )
    return
  end if
  write(teststr,'("intmtemp",I1.1)') ii+1
  if ( trim(teststr)==trim(mode) ) then
    ctmp = 0.
    do ifrac = 1,nfrac
      ctmp = ctmp + real(intm(ifrac)%nodetemp(:,ii))*pd(ifrac)%frac_sigma
    end do
    tempout = unpack( ctmp+urbtemp, upack, tempout )
    return
  end if
end do

write(6,*) "ERROR: unknown option for uclem_avetmp_thread ",trim(mode)
stop
  
return
end subroutine uclem_avetemp_thread

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Store fraction of direct radiation
!

subroutine uclem_fbeam(is,ifin,fbeam,diag)

implicit none

integer, intent(in) :: is,ifin,diag
integer ifinish,ib,ie
integer tile, js, je, kstart, kfinish, jstart, jfinish
real, dimension(ifin), intent(in) :: fbeam

if ( diag>=1 ) write(6,*) "Assign urban direct beam ratio"
if ( .not.uclem_active ) return

ifinish = is + ifin - 1

do tile = 1,ntiles
  js = (tile-1)*imax + 1 ! js:je is the tile portion of 1:ifull
  je = tile*imax         ! js:je is the tile portion of 1:ifull
  if ( ufull_g(tile)>0 ) then
    kstart = max( is - js + 1, 1)          ! kstart:kfinish is the requested portion of 1:imax
    kfinish = min( ifinish - js + 1, imax) ! kstart:kfinish is the requested portion of 1:imax
    if ( kstart<=kfinish ) then
      jstart = kstart + js - is                 ! jstart:jfinish is the tile portion of 1:ifin
      jfinish = kfinish + js - is               ! jstart:jfinish is the tile portion of 1:ifin
      ib = count(upack_g(1:kstart-1,tile))+1
      ie = count(upack_g(kstart:kfinish,tile))+ib-1
      if ( ib<=ie ) then
        f_g(tile)%fbeam(ib:ie)=pack(fbeam(jstart:jfinish),upack_g(kstart:kfinish,tile))
      end if
    end if
  end if
end do

return
end subroutine uclem_fbeam

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Use Spitter et al (1986) method to estimate fraction of direct
! shortwave radiation (from CABLE v1.4)
!

subroutine uclem_spitter(is,ifin,fjd,sg,cosin,diag)

implicit none

integer, intent(in) :: is,ifin,diag
integer ib,ie,ifinish
integer tile, js, je, kstart, kfinish, jstart, jfinish
real, dimension(ifin), intent(in) :: sg,cosin
real, dimension(imax) :: tmpr,tmpk,tmprat
real, dimension(imax) :: lsg,lcosin
real, intent(in) :: fjd
real, parameter :: solcon = 1370.

if ( diag>=1 ) write(6,*) "Diagnose urban direct beam ratio"
if ( .not.uclem_active ) return

ifinish = is + ifin - 1

do tile = 1,ntiles
  js = (tile-1)*imax + 1 ! js:je is the tile portion of 1:ifull
  je = tile*imax         ! js:je is the tile portion of 1:ifull
  if ( ufull_g(tile)>0 ) then
      
    kstart = max( is - js + 1, 1)          ! kstart:kfinish is the requested portion of 1:imax
    kfinish = min( ifinish - js + 1, imax) ! kstart:kfinish is the requested portion of 1:imax
    if ( kstart<=kfinish ) then
      jstart = kstart + js - is             ! jstart:jfinish is the tile portion of 1:ifin
      jfinish = kfinish + js - is           ! jstart:jfinish is the tile portion of 1:ifin
      ib = count(upack_g(1:kstart-1,tile))+1
      ie = count(upack_g(kstart:kfinish,tile))+ib-1
      if ( ib<=ie ) then

        lsg(ib:ie)   =pack(sg(jstart:jfinish),upack_g(kstart:kfinish,tile))
        lcosin(ib:ie)=pack(cosin(jstart:jfinish),upack_g(kstart:kfinish,tile))

        tmpr(ib:ie)=0.847+lcosin(ib:ie)*(1.04*lcosin(ib:ie)-1.61)
        tmpk(ib:ie)=(1.47-tmpr(ib:ie))/1.66
        where (lcosin(ib:ie)>1.0e-10 .and. lsg(ib:ie)>10.)
          tmprat(ib:ie)=lsg(ib:ie)/(solcon*(1.+0.033*cos(2.*pi*(fjd-10.)/365.))*lcosin(ib:ie))
        elsewhere
          tmprat(ib:ie)=0.
        end where
        where (tmprat(ib:ie)>tmpk(ib:ie))
          f_g(tile)%fbeam(ib:ie)=max(1.-tmpr(ib:ie),0.)
        elsewhere (tmprat(ib:ie)>0.35)
          f_g(tile)%fbeam(ib:ie)=min(1.66*tmprat(ib:ie)-0.4728,1.)
        elsewhere (tmprat(ib:ie)>0.22)
          f_g(tile)%fbeam(ib:ie)=6.4*(tmprat(ib:ie)-0.22)**2
        elsewhere
          f_g(tile)%fbeam(ib:ie)=0.
        end where
        
      end if
    end if
  end if
end do

return
end subroutine uclem_spitter

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine calculates the urban contribution to albedo.
! (selected grid points only)

! raw   (.false.=blend, .true.=output only)
! split (0=net albedo, 1=direct albedo, 2=diffuse albedo)

subroutine uclem_alb1(is,ifin,alb,diag,raw,split)

implicit none

integer, intent(in) :: is,ifin,diag
integer ucount,ib,ie,ifinish,albmode
integer tile, js, je, kstart, kfinish, jstart, jfinish
integer ifrac
integer, intent(in), optional :: split
real, dimension(ifin), intent(inout) :: alb
real, dimension(imax) :: ualb,utmp
real, dimension(imax) :: dumfbeam,snowdeltac,snowdeltar
real, dimension(imax) :: tmpalb,wallpsi,roadpsi
real, dimension(imax) :: sg_roof,sg_vegr,sg_road,sg_walle,sg_wallw,sg_vegc,sg_rfsn,sg_rdsn
logical, intent(in), optional :: raw
logical outmode

if ( diag>=1 ) write(6,*) "Calculate urban albedo (broad)"
if ( .not.uclem_active ) return

outmode=.false.
if (present(raw)) outmode=raw

albmode=0 ! net albedo
if (present(split)) albmode=split

ifinish = is + ifin - 1

do tile = 1,ntiles
  js = (tile-1)*imax + 1 ! js:je is the tile portion of 1:ifull
  je = tile*imax         ! js:je is the tile portion of 1:ifull
  if ( ufull_g(tile)>0 ) then
      
    kstart = max( is - js + 1, 1)          ! kstart:kfinish is the requested portion of 1:imax
    kfinish = min( ifinish - js + 1, imax) ! kstart:kfinish is the requested portion of 1:imax
    if ( kstart<=kfinish ) then
      jstart = kstart + js - is             ! jstart:jfinish is the tile portion of 1:ifin
      jfinish = kfinish + js - is           ! jstart:jfinish is the tile portion of 1:ifin
      ib = count(upack_g(1:kstart-1,tile))+1
      ie = count(upack_g(kstart:kfinish,tile))+ib-1
      if ( ib<=ie ) then

        ucount = ie - ib + 1  

        select case(albmode)
          case default ! net albedo
            dumfbeam(ib:ie)=f_g(tile)%fbeam(ib:ie)
          case(1)      ! direct albedo
            dumfbeam(ib:ie)=1.
          case(2)      ! diffuse albedo
            dumfbeam(ib:ie)=0.
        end select

        ualb(ib:ie) = 0.  
        do ifrac = 1,nfrac  

          call uclem_calc_alb(tmpalb(ib:ie),f_g(tile)%hwratio(ib:ie),f_g(tile)%effhwratio(ib:ie), &
                              f_g(tile)%vangle(ib:ie),f_g(tile)%hangle(ib:ie),dumfbeam(ib:ie),    &
                              cnveg_g(tile)%sigma(ib:ie),rfveg_g(tile)%sigma(ib:ie),              &
                              f_g(tile)%sigmabld(ib:ie),f_road(tile)%alpha(ib:ie),                &
                              f_wall(tile)%alpha(ib:ie),f_roof(tile)%alpha(ib:ie),                &
                              cnveg_g(tile)%alpha(ib:ie),rfveg_g(tile)%alpha(ib:ie),              &
                              rdhyd_g(ifrac,tile)%snowalpha(ib:ie),                               &
                              rfhyd_g(ifrac,tile)%snowalpha(ib:ie),                               &
                              rdhyd_g(ifrac,tile)%snow(ib:ie),                                    &
                              rfhyd_g(ifrac,tile)%snow(ib:ie))

          ualb(ib:ie) = ualb(ib:ie) + tmpalb(ib:ie)*p_g(ifrac,tile)%frac_sigma(ib:ie)
        end do  

        if ( outmode ) then
          alb(jstart:jfinish) = unpack(ualb(ib:ie),upack_g(kstart:kfinish,tile),alb(jstart:jfinish))
        else
          utmp(ib:ie) = pack(alb(jstart:jfinish),upack_g(kstart:kfinish,tile))
          utmp(ib:ie) = (1.-f_g(tile)%sigmau(ib:ie))*utmp(ib:ie) + f_g(tile)%sigmau(ib:ie)*ualb(ib:ie)
          alb(jstart:jfinish) = unpack(utmp(ib:ie),upack_g(kstart:kfinish,tile),alb(jstart:jfinish))
        end if
        
      end if
    end if
  end if
end do  

return
end subroutine uclem_alb1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine stores the zenith angle and the solar azimuth angle
! (single grid point)

subroutine uclem_newangle1(is,ifin,cosin,azimuthin,ctimein)

implicit none

integer, intent(in) :: is,ifin
integer ifinish,ib,ie
integer tile, js, je, kstart, kfinish, jstart, jfinish
real, dimension(ifin), intent(in) :: cosin     ! cosine of zenith angle
real, dimension(ifin), intent(in) :: azimuthin ! azimuthal angle
real, dimension(ifin), intent(in) :: ctimein   ! local hour (0<=ctime<=1)

if (.not.uclem_active) return

ifinish = is + ifin - 1

do tile = 1,ntiles
  js = (tile-1)*imax + 1 ! js:je is the tile portion of 1:ifull
  je = tile*imax         ! js:je is the tile portion of 1:ifull
  if ( ufull_g(tile)>0 ) then
    kstart = max( is - js + 1, 1)          ! kstart:kfinish is the requested portion of 1:imax
    kfinish = min( ifinish - js + 1, imax) ! kstart:kfinish is the requested portion of 1:imax
    if ( kstart<=kfinish ) then
      jstart = kstart + js - is             ! jstart:jfinish is the tile portion of 1:ifin
      jfinish = kfinish + js - is           ! jstart:jfinish is the tile portion of 1:ifin
      ib = count(upack_g(1:kstart-1,tile))+1
      ie = count(upack_g(kstart:kfinish,tile))+ib-1
      if ( ib<=ie ) then
        f_g(tile)%hangle(ib:ie)=0.5*pi-pack(azimuthin(jstart:jfinish),upack_g(kstart:kfinish,tile))
        f_g(tile)%vangle(ib:ie)=acos(pack(cosin(jstart:jfinish),upack_g(kstart:kfinish,tile)))
        f_g(tile)%ctime(ib:ie)=pack(ctimein(jstart:jfinish),upack_g(kstart:kfinish,tile))
      end if
    end if
  end if
end do  

return
end subroutine uclem_newangle1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This version of tebnewangle is for CCAM and TAPM
!

subroutine uclem_ccangle(is,ifin,cosin,rlon,rlat,fjd,slag,dt,sdlt)

implicit none

integer, intent(in) :: is,ifin
integer ifinish, ib, ie, iqu
integer tile, js, je, kstart, kfinish, jstart, jfinish
real, intent(in) :: fjd,slag,dt,sdlt
real cdlt, x, y
real, dimension(ifin), intent(in) :: cosin,rlon,rlat
real, dimension(imax) :: hloc,lattmp

! cosin = cosine of zenith angle
! rlon = longitude
! rlat = latitude
! fjd = day of year
! slag = sun lag angle
! sdlt = sin declination of sun

if (.not.uclem_active) return

ifinish = is + ifin - 1

cdlt=sqrt(min(max(1.-sdlt**2,0.),1.))

do tile = 1,ntiles
  js = (tile-1)*imax + 1 ! js:je is the tile portion of 1:ifull
  je = tile*imax         ! js:je is the tile portion of 1:ifull
  if ( ufull_g(tile)>0 ) then
      
    kstart = max( is - js + 1, 1)          ! kstart:kfinish is the requested portion of 1:imax
    kfinish = min( ifinish - js + 1, imax) ! kstart:kfinish is the requested portion of 1:imax
    if ( kstart<=kfinish ) then
      jstart = kstart + js - is             ! jstart:jfinish is the tile portion of 1:ifin
      jfinish = kfinish + js - is           ! jstart:jfinish is the tile portion of 1:ifin
      ib = count(upack_g(1:kstart-1,tile))+1
      ie = count(upack_g(kstart:kfinish,tile))+ib-1
      if ( ib<=ie ) then

        lattmp(ib:ie)=pack(rlat(jstart:jfinish),upack_g(kstart:kfinish,tile))

        ! from CCAM zenith.f90
        hloc(ib:ie)=2.*pi*fjd+slag+pi+pack(rlon(jstart:jfinish),upack_g(kstart:kfinish,tile))+dt*pi/86400.
        
        ! estimate azimuth angle
        do iqu = ib,ie
          x=sin(-hloc(iqu))*cdlt
          y=-cos(-hloc(iqu))*cdlt*sin(lattmp(iqu))+cos(lattmp(iqu))*sdlt
          !azimuth=atan2(x,y)
          f_g(tile)%hangle(iqu)=0.5*pi-atan2(x,y)
        end do  
        
        f_g(tile)%vangle(ib:ie)=acos(pack(cosin(jstart:jfinish),upack_g(kstart:kfinish,tile)))
        f_g(tile)%ctime(ib:ie)=min(max(mod(0.5*hloc(ib:ie)/pi-0.5,1.),0.),1.)
        ! calculate weekdayload loading
        f_g(tile)%weekdayload(ib:ie)=1.0
        if ( (mod(floor(fjd),7)==3) .or. (mod(floor(fjd),7)==4) ) then
          f_g(tile)%weekdayload(ib:ie)=0.9  ! weekend = 90%
        end if  
      end if
    end if
  end if
end do  

return
end subroutine uclem_ccangle

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine calculates screen level diagnostics
!

subroutine uclem_scrnout(tscrn,qscrn,uscrn,u10,diag,raw)

implicit none

integer, intent(in) :: diag
integer tile, is, ie
real, dimension(ifull), intent(inout) :: tscrn,qscrn,uscrn,u10
logical, intent(in), optional :: raw
logical mode

if (diag>=1) write(6,*) "Calculate urban 2m diagnostics"
if (.not.uclem_active) return

mode=.false.
if (present(raw)) mode=raw

do tile = 1,ntiles
  is = (tile-1)*imax + 1
  ie = tile*imax
  if ( ufull_g(tile)>0 ) then
    call uclem_scrnout_thread(tscrn(is:ie),qscrn(is:ie),uscrn(is:ie),u10(is:ie),f_g(tile),p_g(:,tile), &
                              upack_g(:,tile),ufull_g(tile),diag,mode)
  end if
end do  

return
end subroutine uclem_scrnout

subroutine uclem_scrnout_thread(tscrn,qscrn,uscrn,u10,fp,pd,upack,ufull,diag,mode)

implicit none

integer, intent(in) :: ufull, diag
integer ifrac
real, dimension(:), intent(inout) :: tscrn,qscrn,uscrn,u10
real, dimension(ufull) :: tmp
real, dimension(ufull) :: tmp_tscrn, tmp_qscrn, tmp_uscrn, tmp_u10
logical, intent(in) :: mode
logical, dimension(:), intent(in) :: upack
type(fparmdata), intent(in) :: fp
type(pdiagdata), dimension(nfrac), intent(in) :: pd

if (diag>=1) write(6,*) "THREAD: Calculate urban 2m diagnostics"
if (.not.uclem_active) return

tmp_tscrn = 0.
tmp_qscrn = 0.
tmp_uscrn = 0.
tmp_u10 = 0.
do ifrac = 1,nfrac
  tmp_tscrn = tmp_tscrn + (pd(ifrac)%tscrn+urbtemp)*pd(ifrac)%frac_sigma
  tmp_qscrn = tmp_qscrn + pd(ifrac)%qscrn*pd(ifrac)%frac_sigma
  tmp_uscrn = tmp_uscrn + pd(ifrac)%uscrn*pd(ifrac)%frac_sigma
  tmp_u10 = tmp_u10 + pd(ifrac)%u10*pd(ifrac)%frac_sigma
end do    

if (mode) then
  tscrn(:)=unpack(tmp_tscrn,upack,tscrn(:))
  qscrn(:)=unpack(tmp_qscrn,upack,qscrn(:))
  uscrn(:)=unpack(tmp_uscrn,upack,uscrn(:))
  u10(:)  =unpack(tmp_u10,  upack,u10(:)  )
else
  tmp(:)=pack(tscrn(:),upack)
  tmp(:)=fp%sigmau*tmp_tscrn + (1.-fp%sigmau)*tmp(:)
  tscrn(:)=unpack(tmp(:),upack,tscrn(:))
  tmp(:)=pack(qscrn(:),upack)
  tmp(:)=fp%sigmau*tmp_qscrn + (1.-fp%sigmau)*tmp(:)
  qscrn(:)=unpack(tmp(:),upack,qscrn(:))
  tmp(:)=pack(uscrn(:),upack)
  tmp(:)=fp%sigmau*tmp_uscrn + (1.-fp%sigmau)*tmp(:)
  uscrn(:)=unpack(tmp(:),upack,uscrn(:))
  tmp(:)=pack(u10(:),upack)
  tmp(:)=fp%sigmau*tmp_u10 + (1.-fp%sigmau)*tmp(:)
  u10(:)=unpack(tmp(:),upack,u10(:))
end if

return
end subroutine uclem_scrnout_thread

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Extract urban fraction
subroutine uclem_sigmau(sigu,diag)

implicit none

integer, intent(in) :: diag
integer tile, is, ie
real, dimension(ifull), intent(out) :: sigu

if (diag>=1) write(6,*) "Calculate urban cover fraction"
sigu=0.
if (.not.uclem_active) return
do tile = 1,ntiles
  is = (tile-1)*imax + 1
  ie = tile*imax
  if ( ufull_g(tile)>0 ) then
    sigu(is:ie)=unpack(f_g(tile)%sigmau,upack_g(:,tile),0.)
  end if
end do

return
end subroutine uclem_sigmau

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Main routine for calculating urban flux contribution

! ifull = number of horizontal grid points
! dt = model time step (sec)
! zmin = first model level height (m)
! sg = incoming short wave radiation (W/m^2)
! rg = incoming long wave radiation (W/m^2)
! rnd = incoming rainfall/snowfall rate (kg/(m^2 s))
! rho = atmospheric density at first model level (kg/m^3)
! temp = atmospheric temperature at first model level (K)
! mixr = atmospheric mixing ratio at first model level (kg/kg)
! ps = surface pressure (Pa)
! pa = pressure at first model level (Pa)
! uu = U component of wind speed at first model level (m/s)
! vv = V component of wind speed at first model level (m/s)
! umin = minimum wind speed (m/s)
! ofg = Input/Output sensible heat flux (W/m^2)
! oeg = Input/Output latent heat flux (W/m^2)
! ots = Input/Output radiative/skin temperature (K)
! owf = Input/Output wetness fraction/surface water (%)
! orn = runoff
! oevspsbl = evapotranspiration
! osbl = sublimation
! diag = diagnostic message mode (0=off, 1=basic messages, 2=more detailed messages, etc)

subroutine uclem_calc_standard(ofg,oeg,ots,owf,orn,oevspsbl,osbl,dt,zmin,sg,rg,rnd,snd,rho, &
                               temp,mixr,ps,uu,vv,umin,diag,raw)

implicit none

integer, intent(in) :: diag
integer tile, is, ie
real, intent(in) :: dt,umin
real, dimension(ifull), intent(in) :: sg,rg,rnd,snd,rho,temp,mixr,ps,uu,vv,zmin
real, dimension(ifull), intent(inout) :: ofg,oeg,ots,owf,orn,oevspsbl,osbl
logical, intent(in), optional :: raw
logical mode

! mode = .false. implies weight output with urban area cover fraction
! mode = .true. implies no weighting of output with urban area cover fraction (assumes 100% cover)
mode=.false.
if (present(raw)) mode=raw

call uclem_calc2_standard(ofg,oeg,ots,owf,orn,oevspsbl,osbl,dt,zmin,zmin,sg,rg,rnd,snd,rho, &
                          temp,temp,mixr,mixr,ps,uu,uu,vv,vv,umin,diag,raw=mode)

return
end subroutine uclem_calc_standard

subroutine uclem_calc2_standard(ofg,oeg,ots,owf,orn,oevspsbl,osbl,dt,zmin,zroof,sg,rg,rnd,snd,rho, &
                                temp,temproof,mixr,mixrroof,ps,uu,uuroof,vv,vvroof,umin,diag,raw)

implicit none

integer, intent(in) :: diag
integer tile, is, ie
real, intent(in) :: dt,umin
real, dimension(ifull), intent(in) :: sg,rg,rnd,snd,rho,temp,mixr,ps,uu,vv,zmin
real, dimension(ifull), intent(in) :: temproof,mixrroof,uuroof,vvroof,zroof
real, dimension(ifull), intent(inout) :: ofg,oeg,ots,owf,orn,oevspsbl,osbl
logical, intent(in), optional :: raw
logical mode

! mode = .false. implies weight output with urban area cover fraction
! mode = .true. implies no weighting of output with urban area cover fraction (assumes 100% cover)
mode=.false.
if (present(raw)) mode=raw

if ( .not.uclem_active ) return

do tile = 1,ntiles
  is = (tile-1)*imax + 1
  ie = tile*imax
  if ( ufull_g(tile)>0 ) then
    call uclem_calc2_thread(ofg(is:ie),oeg(is:ie),ots(is:ie),owf(is:ie),orn(is:ie),oevspsbl(is:ie),       &
                            osbl(is:ie),dt,zmin(is:ie),zroof(is:ie),sg(is:ie),rg(is:ie),rnd(is:ie),       &
                            snd(is:ie),rho(is:ie),temp(is:ie),temproof(is:ie),mixr(is:ie),                &
                            mixrroof(is:ie),ps(is:ie),uu(is:ie),uuroof(is:ie),vv(is:ie),vvroof(is:ie),    &
                            umin,f_g(tile),f_intm(tile),f_road(tile),f_roof(tile),f_slab(tile),           &
                            f_wall(tile),intm_g(:,tile),p_g(:,tile),rdhyd_g(:,tile),rfhyd_g(:,tile),      &
                            rfveg_g(tile),road_g(:,tile),roof_g(:,tile),room_g(:,tile),slab_g(:,tile),    &
                            walle_g(:,tile),wallw_g(:,tile),cnveg_g(tile),intl_g(tile),upack_g(:,tile),   &
                            ufull_g(tile),diag,raw=mode)
  end if
end do

return
end subroutine uclem_calc2_standard
                             
subroutine uclem_calc_thread(ofg,oeg,ots,owf,orn,oevspsbl,osbl,dt,zmin,sg,rg,rnd,snd,rho,temp,    &
                             mixr,ps,uu,vv,umin,fp,fp_intm,fp_road,fp_roof,fp_slab,fp_wall,intm,pd,      &
                             rdhyd,rfhyd,rfveg,road,roof,room,slab,walle,wallw,cnveg,intl,               &
                             upack,ufull,diag,raw)

implicit none

integer, intent(in) :: ufull, diag
integer ifrac
real, intent(in) :: dt, umin
real, dimension(:), intent(in) :: sg,rg,rnd,snd,rho,temp,mixr,ps,uu,vv,zmin
real, dimension(:), intent(inout) :: ofg,oeg,ots,owf,orn,oevspsbl,osbl
real, dimension(ufull) :: tmp
real, dimension(ufull) :: a_sg,a_rg,a_rho,a_temp,a_mixr,a_ps,a_umag,a_udir,a_rnd,a_snd,a_zmin
real, dimension(ufull) :: u_fg,u_eg,u_ts,u_wf,u_rn,u_evspsbl,u_sbl
real, dimension(ufull) :: tmp_fg, tmp_eg, tmp_ts, tmp_wf, tmp_rn, tmp_evspsbl, tmp_sbl
logical, intent(in), optional :: raw
logical mode
logical, dimension(:), intent(in) :: upack
type(hydrodata), dimension(nfrac), intent(inout) :: rdhyd, rfhyd
type(facetdata), dimension(nfrac), intent(inout) :: road, roof, room, slab, walle, wallw, intm
type(pdiagdata), dimension(nfrac), intent(inout) :: pd
type(facetparams), intent(in) :: fp_intm, fp_road, fp_roof, fp_slab, fp_wall
type(vegdata), intent(inout) :: rfveg
type(vegdata), intent(inout) :: cnveg
type(intldata), intent(in) :: intl
type(fparmdata), intent(in) :: fp

call uclem_calc2_thread(ofg,oeg,ots,owf,orn,oevspsbl,osbl,dt,zmin,zmin,sg,rg,rnd,snd,rho,           &
                        temp,temp,mixr,mixr,ps,uu,uu,vv,vv,                                         &
                        umin,fp,fp_intm,fp_road,fp_roof,fp_slab,fp_wall,intm,pd,                    &
                        rdhyd,rfhyd,rfveg,road,roof,room,slab,walle,wallw,cnveg,intl,               &
                        upack,ufull,diag,raw)

return
end subroutine uclem_calc_thread

subroutine uclem_calc2_thread(ofg,oeg,ots,owf,orn,oevspsbl,osbl,dt,zmin,zroof,sg,rg,rnd,snd,rho,  &
                              temp,temproof,mixr,mixrroof,ps,uu,uuroof,vv,vvroof,                         &
                              umin,fp,fp_intm,fp_road,fp_roof,fp_slab,fp_wall,intm,pd,                    &
                              rdhyd,rfhyd,rfveg,road,roof,room,slab,walle,wallw,cnveg,intl,               &
                              upack,ufull,diag,raw)

implicit none

integer, intent(in) :: ufull, diag
integer ifrac
real, intent(in) :: dt, umin
real, dimension(:), intent(in) :: sg,rg,rnd,snd,rho,temp,mixr,ps,uu,vv,zmin
real, dimension(:), intent(in) :: temproof,mixrroof,uuroof,vvroof,zroof
real, dimension(:), intent(inout) :: ofg,oeg,ots,owf,orn,oevspsbl,osbl
real, dimension(ufull) :: tmp
real, dimension(ufull) :: a_sg,a_rg,a_rho,a_temp,a_mixr,a_ps,a_umag,a_udir,a_rnd,a_snd,a_zmin
real, dimension(ufull) :: a_temproof, a_mixrroof, a_umagroof, a_zroof
real, dimension(ufull) :: u_fg,u_eg,u_ts,u_wf,u_rn,u_evspsbl,u_sbl
real, dimension(ufull) :: tmp_fg, tmp_eg, tmp_ts, tmp_wf, tmp_rn, tmp_evspsbl, tmp_sbl
logical, intent(in), optional :: raw
logical mode
logical, dimension(:), intent(in) :: upack
type(hydrodata), dimension(nfrac), intent(inout) :: rdhyd, rfhyd
type(facetdata), dimension(nfrac), intent(inout) :: road, roof, room, slab, walle, wallw, intm
type(pdiagdata), dimension(nfrac), intent(inout) :: pd
type(facetparams), intent(in) :: fp_intm, fp_road, fp_roof, fp_slab, fp_wall
type(vegdata), intent(inout) :: rfveg
type(vegdata), intent(inout) :: cnveg
type(intldata), intent(in) :: intl
type(fparmdata), intent(in) :: fp

if ( ufull==0 ) return ! no urban grid points

! mode = .false. implies weight output with urban area cover fraction
! mode = .true. implies no weighting of output with urban area cover fraction (assumes 100% cover)
mode=.false.
if (present(raw)) mode=raw

! Host model meteorological data
a_zmin    =pack(zmin,                         upack)
a_zroof   =pack(zroof,                        upack)
a_sg      =pack(sg,                           upack)
a_rg      =pack(rg,                           upack)
a_rho     =pack(rho,                          upack)
a_temp    =pack(temp-urbtemp,                 upack)
a_temproof=pack(temproof-urbtemp,             upack)
a_mixr    =pack(mixr,                         upack)
a_mixrroof=pack(mixrroof,                     upack)
a_ps      =pack(ps,                           upack)
a_umag    =max(pack(sqrt(uu**2+vv**2),        upack),umin)
a_umagroof=max(pack(sqrt(uuroof**2+vvroof**2),upack),umin)
a_udir    =pack(atan2(vv,uu),                 upack)
a_rnd     =pack(rnd-snd,                      upack)
a_snd     =pack(snd,                          upack)

u_fg = 0.
u_eg = 0.
u_ts = 0.
u_wf = 0.
u_rn = 0.
u_evspsbl = 0.
u_sbl = 0.

do ifrac = 1,nfrac

  ! Update urban prognostic variables
  call uclem_eval(tmp_fg,tmp_eg,tmp_ts,tmp_wf,tmp_rn,tmp_evspsbl,tmp_sbl,dt,a_sg,a_rg,a_rho,a_temp,a_temproof,  &
                  a_mixr,a_mixrroof,a_ps,a_umag,a_umagroof,a_udir,a_rnd,a_snd,a_zmin,a_zroof,fp,fp_intm,        &
                  fp_road,fp_roof,fp_slab,fp_wall,intm(ifrac),pd(ifrac),rdhyd(ifrac),rfhyd(ifrac),rfveg,        &
                  road(ifrac),roof(ifrac),room(ifrac),slab(ifrac),walle(ifrac),wallw(ifrac),cnveg,intl,ufull,   &
                  ifrac,0,diag)

  u_fg = u_fg + tmp_fg*pd(ifrac)%frac_sigma
  u_eg = u_eg + tmp_eg*pd(ifrac)%frac_sigma
  u_ts = u_ts + (tmp_ts+urbtemp)**4*pd(ifrac)%frac_sigma ! combine as longwave radiation
  u_wf = u_wf + tmp_wf*pd(ifrac)%frac_sigma
  u_rn = u_rn + tmp_rn*pd(ifrac)%frac_sigma
  u_evspsbl = u_evspsbl + tmp_evspsbl*pd(ifrac)%frac_sigma
  u_sbl = u_sbl + tmp_sbl*pd(ifrac)%frac_sigma
  
end do

u_ts = u_ts**0.25 ! convert back to temperature
  
! export urban fluxes on host grid
if (mode) then
  ofg = unpack(u_fg,upack,ofg)
  oeg = unpack(u_eg,upack,oeg)
  ots = unpack(u_ts,upack,ots)
  owf = unpack(u_wf,upack,owf)
  orn = unpack(u_rn,upack,orn)
  oevspsbl = unpack(u_evspsbl,upack,oevspsbl)
  osbl = unpack(u_sbl,upack,osbl)
else
  tmp = pack(ofg,upack)
  tmp = (1.-fp%sigmau)*tmp + fp%sigmau*u_fg
  ofg = unpack(tmp,upack,ofg)
  tmp = pack(oeg,upack)
  tmp = (1.-fp%sigmau)*tmp + fp%sigmau*u_eg
  oeg = unpack(tmp,upack,oeg)
  tmp = pack(ots,upack)
  tmp = ((1.-fp%sigmau)*tmp**4+fp%sigmau*u_ts**4)**0.25
  ots = unpack(tmp,upack,ots)
  tmp = pack(owf,upack)
  tmp = (1.-fp%sigmau)*tmp + fp%sigmau*u_wf
  owf = unpack(tmp,upack,owf)
  tmp = pack(orn,upack)
  tmp = (1.-fp%sigmau)*tmp + fp%sigmau*u_rn
  orn = unpack(tmp,upack,orn)
  tmp = pack(oevspsbl,upack)
  tmp = (1.-fp%sigmau)*tmp + fp%sigmau*u_evspsbl
  oevspsbl = unpack(tmp,upack,oevspsbl)
  tmp = pack(osbl,upack)
  tmp = (1.-fp%sigmau)*tmp + fp%sigmau*u_sbl
  osbl = unpack(tmp,upack,osbl)
end if

return
end subroutine uclem_calc2_thread
                    

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Disables UCLEM so subroutine calls have no effect

subroutine uclem_disable(diag)

implicit none

integer, intent(in) :: diag

if ( diag>=1 ) write(6,*) "Disable UCLEM"

uclem_active = .false.

return
end subroutine uclem_disable

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end module uclem_ctrl
