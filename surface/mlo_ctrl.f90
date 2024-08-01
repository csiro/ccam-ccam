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
    
! This is a 1D, mixed layer ocean model based on Large, et al (1994), for ensemble regional climate simulations.
! In CCAM this module interfaces withmlodynamics.f90 for river routing, diffusion and advection routines.

! This version supports a thermodynamic model of sea ice based on O'Farrell from Mk3.5.  We have included a
! free surface so that lakes can change depth, etc.

! This version can assimilate SSTs from GCMs, using a convolution based digital filter (see nesting.f90),
! which avoids problems with complex land-sea boundary conditions.

! Order of calls:
!    call mloinit
!    call mloload
!    ...
!    main loop
!       ...
!       call mloalb4
!       ...
!       call mloeval
!       ...
!       call mloscrnout
!       ...
!    end loop
!    ...
!    call mlosave
!    call mloend

module mlo_ctrl

#ifdef CCAM
use newmpar_m, only : imax, ntiles
#endif
use mlo

implicit none

private
public mloinit,mloend,mloeval,mloimport,mloexport,mloload,mlosave,mloregrid,mlodiag,mloalb2,mloalb4, &
       mloscrnout,mloextra,mloimpice,mloexpice,mloexpdep,mloexpdensity,mloexpmelt,mloexpgamm,        &
       mloimport3d,mloexport3d,mlovlevels,mlocheck,mlodiagice,mlo_updatediag,mlo_updatekm,           &
       mlosurf,mlonewice,mlo_ema
public minsfc,minsal,maxsal,icemax
public micdwn
public wlev,zomode,wrtemp,wrtrho,mxd,mindep,minwater,zoseaice,factchseaice,otaumode,mlosigma
public oclosure,pdl,pdu,usepice,minicemass,cdbot,cp0,ominl,omaxl,mlo_adjeta,mlo_limitsal
public mlo_timeave_length,kemaxdt,mlo_step,mlo_uvcoupl,fluxwgt,delwater

#ifdef CCAM
public water_g,ice_g
public alphavis_seaice, alphanir_seaice
public alphavis_seasnw, alphanir_seasnw
public dgwater_g,dgice_g,dgscrn_g
public depth_g
public waterdata,icedata
public dgwaterdata,dgicedata,dgscrndata,depthdata
public turb_g
public turbdata
public omink,omineps
public k_mode,eps_mode,limitL,fixedce3
public nops,nopb,fixedstabfunc
#endif

! parameters
integer, save :: ifull                                  ! Grid size
#ifndef CCAM
integer, save :: ntiles = 1                             ! Emulate OMP
integer, save :: imax = 0                               ! Emulate OMP
#endif
! model arrays
logical, save :: mlo_active = .false.                   ! Flag if MLO has been initialised
real, dimension(:,:), allocatable, save :: micdwn       ! This variable is for CCAM onthefly.f

type(waterdata), dimension(:), allocatable, save :: water_g
type(icedata), dimension(:), allocatable, save :: ice_g
type(dgwaterdata), dimension(:), allocatable, save :: dgwater_g
type(dgicedata), dimension(:), allocatable, save :: dgice_g
type(dgscrndata), dimension(:), allocatable, save :: dgscrn_g
type(depthdata), dimension(:), allocatable, save :: depth_g
type(turbdata), dimension(:), allocatable, save :: turb_g
  
interface mloeval
  module procedure mloeval_standard, mloeval_thread
end interface

interface mloimport
  module procedure mlo_import_ifull, mlo_import_imax
end interface

interface mloexport
  module procedure mlo_export_ifull, mlo_export_imax
end interface

interface mloimpice
   module procedure mlo_impice_ifull, mlo_impice_imax
end interface

interface mloexpice
  module procedure mlo_expice_ifull, mlo_expice_imax
end interface

interface mloextra
  module procedure mloextra_ifull, mloextra_imax
end interface

interface mloexpdep
  module procedure mloexpdep_old
  module procedure mlo_export_depth_ifull, mlo_export_depth_imax
end interface

interface mlosurf
  module procedure mlo_surf_ifull, mlo_surf_imax
end interface

interface mlo_updatekm
  module procedure mlo_updatekm_imax
end interface

interface mlodiag
  module procedure mlodiag_old
  module procedure mlo_diag_ifull, mlo_diag_imax
end interface

interface mlo_updatediag
  module procedure mlo_update_diag_imax
end interface

interface mlodiagice
  module procedure mlo_diagice_ifull, mlo_diagice_imax
end interface

contains


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  Initialise MLO

subroutine mloinit(ifin,depin,f,diag)

implicit none

integer, intent(in) :: ifin, diag
integer iqw, ii, is, ie, tile
real, dimension(ifin), intent(in) :: depin, f
real, dimension(ifin) :: deptmp
real, dimension(wlev) :: dumdf
real, dimension(wlev+1) :: dumdh
logical, dimension(ifin) :: lbottom

if (diag>=1) write(6,*) "Initialising MLO"

mlo_active = .false.

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

if ( minsfc>minwater ) then
  write(6,*) "ERROR: MLO parameters are invalid.  minsfc>minwater"
  stop
end if

allocate( water_g(ntiles) )
allocate( ice_g(ntiles) )
allocate( dgwater_g(ntiles) )
allocate( dgice_g(ntiles) )
allocate( dgscrn_g(ntiles) )
allocate( depth_g(ntiles) )
allocate( turb_g(ntiles) )

do tile = 1,ntiles
  is = (tile-1)*imax + 1
  ie = tile*imax
  
  depth_g(tile)%data_allocated = .false.
  
  if ( any( depin(is:ie)>minwater ) ) then
 
    allocate(water_g(tile)%temp(imax,wlev),water_g(tile)%sal(imax,wlev))
    allocate(water_g(tile)%u(imax,wlev),water_g(tile)%v(imax,wlev))
    allocate(water_g(tile)%w(imax,wlev))
    allocate(water_g(tile)%eta(imax))
    allocate(water_g(tile)%ubot(imax),water_g(tile)%vbot(imax))
    allocate(water_g(tile)%utop(imax),water_g(tile)%vtop(imax))
    allocate(water_g(tile)%u_ema(imax,wlev),water_g(tile)%v_ema(imax,wlev))
    allocate(water_g(tile)%w_ema(imax,wlev))
    allocate(water_g(tile)%temp_ema(imax,wlev),water_g(tile)%sal_ema(imax,wlev))
    allocate(water_g(tile)%dudz(imax,wlev),water_g(tile)%dvdz(imax,wlev))
    allocate(water_g(tile)%dwdx(imax,wlev),water_g(tile)%dwdy(imax,wlev))
    allocate(ice_g(tile)%temp(imax,0:2),ice_g(tile)%thick(imax),ice_g(tile)%snowd(imax))
    allocate(ice_g(tile)%fracice(imax),ice_g(tile)%tsurf(imax),ice_g(tile)%store(imax))
    allocate(ice_g(tile)%u(imax),ice_g(tile)%v(imax))
    allocate(dgwater_g(tile)%visdiralb(imax),dgwater_g(tile)%visdifalb(imax))
    allocate(dgwater_g(tile)%nirdiralb(imax),dgwater_g(tile)%nirdifalb(imax))
    allocate(dgwater_g(tile)%zo(imax),dgwater_g(tile)%zoh(imax),dgwater_g(tile)%zoq(imax))
    allocate(dgwater_g(tile)%cd(imax),dgwater_g(tile)%cdh(imax),dgwater_g(tile)%cdq(imax))
    allocate(dgwater_g(tile)%umod(imax))
    allocate(dgwater_g(tile)%fg(imax),dgwater_g(tile)%eg(imax))
    allocate(dgwater_g(tile)%mixdepth(imax),dgwater_g(tile)%bf(imax))
    allocate(dgwater_g(tile)%taux(imax),dgwater_g(tile)%tauy(imax))
    allocate(dgwater_g(tile)%rho(imax,wlev),dgwater_g(tile)%alpha(imax,wlev))
    allocate(dgwater_g(tile)%beta(imax,wlev),dgwater_g(tile)%rad(imax,wlev))
    allocate(dgwater_g(tile)%wt0(imax),dgwater_g(tile)%wt0_rad(imax))
    allocate(dgwater_g(tile)%wt0_melt(imax),dgwater_g(tile)%wt0_eg(imax))
    allocate(dgwater_g(tile)%ws0(imax),dgwater_g(tile)%wu0(imax))
    allocate(dgwater_g(tile)%ws0_subsurf(imax))
    allocate(dgwater_g(tile)%wv0(imax),dgwater_g(tile)%wt0_fb(imax))
    allocate(dgwater_g(tile)%cd_bot(imax))
    allocate(dgice_g(tile)%visdiralb(imax),dgice_g(tile)%visdifalb(imax))
    allocate(dgice_g(tile)%nirdiralb(imax),dgice_g(tile)%nirdifalb(imax))
    allocate(dgice_g(tile)%cd(imax),dgice_g(tile)%cdh(imax),dgice_g(tile)%cdq(imax))
    allocate(dgice_g(tile)%umod(imax))
    allocate(dgice_g(tile)%fg(imax),dgice_g(tile)%eg(imax))
    allocate(dgice_g(tile)%wetfrac(imax),dgwater_g(tile)%mixind(imax))
    allocate(dgice_g(tile)%tauxica(imax),dgice_g(tile)%tauyica(imax))
    allocate(dgice_g(tile)%tauxicw(imax),dgice_g(tile)%tauyicw(imax))
    allocate(dgice_g(tile)%imass(imax))
    allocate(dgice_g(tile)%cd_bot(imax))
    allocate(dgscrn_g(tile)%temp(imax),dgscrn_g(tile)%u2(imax),dgscrn_g(tile)%qg(imax))
    allocate(dgscrn_g(tile)%u10(imax))
    allocate(depth_g(tile)%depth(imax,wlev),depth_g(tile)%dz(imax,wlev))
    allocate(depth_g(tile)%depth_hl(imax,wlev+1),depth_g(tile)%dz_hl(imax,2:wlev))
    allocate(depth_g(tile)%f(imax),depth_g(tile)%ibot(imax))
    allocate(turb_g(tile)%k(imax,wlev),turb_g(tile)%eps(imax,wlev))

    water_g(tile)%temp=288.-wrtemp    ! K
    water_g(tile)%sal=34.72           ! PSU
    water_g(tile)%u=0.                ! m/s
    water_g(tile)%v=0.                ! m/s
    water_g(tile)%w=0.                ! m/s
    water_g(tile)%eta=0.              ! m
    water_g(tile)%ubot=0.             ! m/s
    water_g(tile)%vbot=0.             ! m/s
    water_g(tile)%utop=0.             ! m/s
    water_g(tile)%vtop=0.             ! m/s
    water_g(tile)%u_ema=0.            ! m/s
    water_g(tile)%v_ema=0.            ! m/s
    water_g(tile)%w_ema=0.            ! m/s
    water_g(tile)%temp_ema=0.         ! K
    water_g(tile)%sal_ema=0.          ! PSU
    water_g(tile)%dudz=0.
    water_g(tile)%dvdz=0.
    water_g(tile)%dwdx=0.
    water_g(tile)%dwdy=0.
  
    ice_g(tile)%thick=0.              ! m
    ice_g(tile)%snowd=0.              ! m
    ice_g(tile)%fracice=0.            ! %
    ice_g(tile)%tsurf=271.2           ! K
    ice_g(tile)%temp=271.2            ! K
    ice_g(tile)%store=0.
    ice_g(tile)%u=0.                  ! m/s
    ice_g(tile)%v=0.                  ! m/s

    dgwater_g(tile)%mixdepth=-1.      ! m
    dgwater_g(tile)%bf=0.
    dgwater_g(tile)%visdiralb=0.06
    dgwater_g(tile)%visdifalb=0.06
    dgwater_g(tile)%nirdiralb=0.06
    dgwater_g(tile)%nirdifalb=0.06
    dgwater_g(tile)%zo=0.001
    dgwater_g(tile)%zoh=0.001
    dgwater_g(tile)%zoq=0.001
    dgwater_g(tile)%cd=0.
    dgwater_g(tile)%cdh=0.
    dgwater_g(tile)%cdq=0.
    dgwater_g(tile)%cd_bot=0.
    dgwater_g(tile)%umod=0.
    dgwater_g(tile)%fg=0.
    dgwater_g(tile)%eg=0.
    dgwater_g(tile)%mixind=wlev-1
    dgwater_g(tile)%taux=0.
    dgwater_g(tile)%tauy=0.
    dgwater_g(tile)%rho=wrtrho
    dgwater_g(tile)%alpha=0.
    dgwater_g(tile)%beta=0.
    dgwater_g(tile)%wt0=0.
    dgwater_g(tile)%wt0_rad=0.
    dgwater_g(tile)%wt0_melt=0.
    dgwater_g(tile)%wt0_eg=0.
    dgwater_g(tile)%wt0_fb=0.
    dgwater_g(tile)%ws0=0.
    dgwater_g(tile)%ws0_subsurf=0.
    dgwater_g(tile)%wu0=0.
    dgwater_g(tile)%wv0=0.
    
    dgice_g(tile)%visdiralb=0.65
    dgice_g(tile)%visdifalb=0.65
    dgice_g(tile)%nirdiralb=0.65
    dgice_g(tile)%nirdifalb=0.65
    dgice_g(tile)%cd=0.
    dgice_g(tile)%cdh=0.
    dgice_g(tile)%cdq=0.
    dgice_g(tile)%cd_bot=0.
    dgice_g(tile)%umod=0.
    dgice_g(tile)%fg=0.
    dgice_g(tile)%eg=0.
    dgice_g(tile)%wetfrac=1.
    dgice_g(tile)%tauxica=0.
    dgice_g(tile)%tauyica=0.
    dgice_g(tile)%tauxicw=0.
    dgice_g(tile)%tauyicw=0.
    dgice_g(tile)%imass=0.
    
    dgscrn_g(tile)%temp=273.2
    dgscrn_g(tile)%qg=0.
    dgscrn_g(tile)%u2=0.
    dgscrn_g(tile)%u10=0.
    
    depth_g(tile)%ibot=wlev

    ! MLO - 30 level ( mlosigma=6 )
    !depth = (/   4.5,  14.2,  25.9,  39.8,  56.5,  76.3, 100.0, 128.1, 161.7, 201.5, &
    !           248.7, 304.7, 370.9, 449.0, 540.8, 648.4, 774.0, 920.0,1088.9,1283.0, &
    !          1504.5,1755.3,2036.8,2349.6,2693.2,3066.2,3456.6,3887.4,4325.2,4775.7 /)
    ! Mk3.5 - 31 level
    !depth = (/   5.0,  15.0,  28.2,  42.0,  59.7,  78.5, 102.1, 127.9, 159.5, 194.6, &
    !           237.0, 284.7, 341.7, 406.4, 483.2, 570.9, 674.9, 793.8, 934.1,1095.2, &
    !          1284.7,1502.9,1758.8,2054.3,2400.0,2800.0,3200.0,3600.0,4000.0,4400.0, &
    !          4800.0 /)
    !ACCESS 1.3 - 50 level
    !depth = (/   5.0,  15.0,  25.0,  35.0,  45.0,  55.0,  65.0,  75.0,  85.0,  95.0, &
    !           105.0, 115.0, 125.0, 135.0, 145.0, 155.0, 165.0, 175.0, 185.0, 195.0, &
    !           205.0, 216.8, 241.4, 280.8, 343.3, 427.3, 536.7, 665.4, 812.8, 969.1, &
    !          1130.9,1289.6,1455.8,1622.9,1801.6,1984.9,2182.9,2388.4,2610.9,2842.6, &
    !          3092.2,3351.3,3628.1,3913.3,4214.5,4521.9,4842.6,5166.1,5499.2,5831.3 /)

    deptmp(1:imax) = depin(is:ie)

    do iqw = 1,imax
      call vgrid(wlev,deptmp(iqw),dumdf,dumdh)
      depth_g(tile)%depth(iqw,:) = dumdf
      depth_g(tile)%depth_hl(iqw,:) = dumdh
    end do
    do ii = 1,wlev
      depth_g(tile)%dz(:,ii) = depth_g(tile)%depth_hl(:,ii+1) - depth_g(tile)%depth_hl(:,ii)
    end do
    ! dz_hl just used for mixing between full levels.  Hence ii=1 and ii=wlev+1 are not defined.
    do ii = 2,wlev
      where ( depth_g(tile)%dz(:,ii)*depth_g(tile)%dz(:,ii-1)>1.e-4 )
        depth_g(tile)%dz_hl(:,ii) = depth_g(tile)%depth(:,ii) - depth_g(tile)%depth(:,ii-1)
      elsewhere
        depth_g(tile)%dz_hl(:,ii) = 0.  
      end where    
    end do
    
    depth_g(tile)%f(:) = f(is:ie)
    depth_g(tile)%data_allocated = .true.
    mlo_active = .true.
    
    ! find bottom index
    do ii = 1,wlev
      lbottom(1:imax) = depth_g(tile)%depth_hl(:,ii+1)>=depth_g(tile)%depth_hl(:,wlev+1) .and.  &
                        depth_g(tile)%depth_hl(:,wlev+1)<mxd .and. depth_g(tile)%dz(:,ii)>1.e-4  
      where ( lbottom(1:imax) )
        depth_g(tile)%ibot(:) = ii
      end where
    end do

    turb_g(tile)%k = omink
    turb_g(tile)%eps = omineps
    
  end if  

end do    

return
end subroutine mloinit

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Deallocate MLO arrays

subroutine mloend

implicit none

integer tile

if ( mlo_active ) then

  do tile = 1,ntiles
      
    if ( depth_g(tile)%data_allocated ) then  

      deallocate(water_g(tile)%temp,water_g(tile)%sal,water_g(tile)%u,water_g(tile)%v)
      deallocate(water_g(tile)%w)
      deallocate(water_g(tile)%eta)
      deallocate(water_g(tile)%ubot,water_g(tile)%vbot)
      deallocate(water_g(tile)%utop,water_g(tile)%vtop)
      deallocate(water_g(tile)%u_ema,water_g(tile)%v_ema)
      deallocate(water_g(tile)%w_ema)
      deallocate(water_g(tile)%temp_ema,water_g(tile)%sal_ema)
      deallocate(water_g(tile)%dudz,water_g(tile)%dvdz)    
      deallocate(water_g(tile)%dwdx,water_g(tile)%dwdy)    
      deallocate(ice_g(tile)%temp,ice_g(tile)%thick,ice_g(tile)%snowd)
      deallocate(ice_g(tile)%fracice,ice_g(tile)%tsurf,ice_g(tile)%store)
      deallocate(ice_g(tile)%u,ice_g(tile)%v)
      deallocate(dgwater_g(tile)%mixdepth,dgwater_g(tile)%bf)
      deallocate(dgwater_g(tile)%visdiralb,dgwater_g(tile)%visdifalb)
      deallocate(dgwater_g(tile)%nirdiralb,dgwater_g(tile)%nirdifalb)
      deallocate(dgwater_g(tile)%zo,dgwater_g(tile)%zoh,dgwater_g(tile)%zoq)
      deallocate(dgwater_g(tile)%cd)
      deallocate(dgwater_g(tile)%cdq)
      deallocate(dgwater_g(tile)%cdh,dgwater_g(tile)%fg,dgwater_g(tile)%eg)
      deallocate(dgwater_g(tile)%umod)
      deallocate(dgwater_g(tile)%taux,dgwater_g(tile)%tauy)
      deallocate(dgwater_g(tile)%rho,dgwater_g(tile)%alpha)
      deallocate(dgwater_g(tile)%beta,dgwater_g(tile)%rad)
      deallocate(dgwater_g(tile)%wt0,dgwater_g(tile)%wt0_rad)
      deallocate(dgwater_g(tile)%wt0_melt,dgwater_g(tile)%wt0_eg)
      deallocate(dgwater_g(tile)%ws0,dgwater_g(tile)%wu0)
      deallocate(dgwater_g(tile)%ws0_subsurf)
      deallocate(dgwater_g(tile)%wv0,dgwater_g(tile)%wt0_fb)
      deallocate(dgwater_g(tile)%cd_bot)
      deallocate(dgice_g(tile)%visdiralb,dgice_g(tile)%visdifalb)
      deallocate(dgice_g(tile)%nirdiralb,dgice_g(tile)%nirdifalb)
      deallocate(dgice_g(tile)%cd)
      deallocate(dgice_g(tile)%cdq)
      deallocate(dgice_g(tile)%cdh,dgice_g(tile)%fg,dgice_g(tile)%eg)
      deallocate(dgice_g(tile)%umod)
      deallocate(dgice_g(tile)%wetfrac,dgwater_g(tile)%mixind)
      deallocate(dgice_g(tile)%tauxica,dgice_g(tile)%tauyica)
      deallocate(dgice_g(tile)%tauxicw,dgice_g(tile)%tauyicw)
      deallocate(dgice_g(tile)%imass)
      deallocate(dgice_g(tile)%cd_bot)
      deallocate(dgscrn_g(tile)%temp,dgscrn_g(tile)%u2,dgscrn_g(tile)%qg,dgscrn_g(tile)%u10)
      deallocate(depth_g(tile)%depth,depth_g(tile)%dz,depth_g(tile)%depth_hl,depth_g(tile)%dz_hl)
      deallocate(depth_g(tile)%f,depth_g(tile)%ibot)
      deallocate(turb_g(tile)%k,turb_g(tile)%eps)
      
      depth_g(tile)%data_allocated = .false.
      
    end if  

  end do

  deallocate( water_g )
  deallocate( ice_g )
  deallocate( dgwater_g )
  deallocate( dgice_g )
  deallocate( dgscrn_g )
  deallocate( depth_g )
  deallocate( turb_g )

end if

mlo_active = .false.

return
end subroutine mloend

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Load MLO data

subroutine mloload(datain,shin,icein,diag)

implicit none

integer, intent(in) :: diag
integer ii, tile, is, ie
real, dimension(ifull,wlev,6), intent(in) :: datain
real, dimension(ifull,10), intent(in) :: icein
real, dimension(ifull), intent(in) :: shin

if (diag>=1) write(6,*) "Load MLO data"
if (.not.mlo_active) return

do tile = 1,ntiles
  is = (tile-1)*imax + 1
  ie = tile*imax
  
  if ( depth_g(tile)%data_allocated ) then
  
    do ii=1,wlev
      where ( depth_g(tile)%dz(:,ii)>=1.e-4 )  
        water_g(tile)%temp(:,ii)=datain(is:ie,ii,1)
        water_g(tile)%sal(:,ii) =datain(is:ie,ii,2)
        water_g(tile)%u(:,ii)   =datain(is:ie,ii,3)
        water_g(tile)%v(:,ii)   =datain(is:ie,ii,4)
        turb_g(tile)%k(:,ii)   =datain(is:ie,ii,5)
        turb_g(tile)%eps(:,ii) =datain(is:ie,ii,6)
      elsewhere
        water_g(tile)%temp(:,ii) = 288.-wrtemp
        water_g(tile)%sal(:,ii) = 34.72
        water_g(tile)%u(:,ii) = 0.
        water_g(tile)%v(:,ii) = 0.
        turb_g(tile)%k(:,ii)   = omink
        turb_g(tile)%eps(:,ii) = omineps
      end where  
    end do
    where ( depth_g(tile)%dz(:,1)>=1.e-4 )
      water_g(tile)%eta(:)   =shin(is:ie)
      ice_g(tile)%tsurf(:)   =icein(is:ie,1)
      ice_g(tile)%temp(:,0)  =icein(is:ie,2)
      ice_g(tile)%temp(:,1)  =icein(is:ie,3)
      ice_g(tile)%temp(:,2)  =icein(is:ie,4)
      ice_g(tile)%fracice(:) =icein(is:ie,5)
      ice_g(tile)%thick(:)   =icein(is:ie,6)
      ice_g(tile)%snowd(:)   =icein(is:ie,7)
      ice_g(tile)%store(:)   =icein(is:ie,8)
      ice_g(tile)%u(:)       =icein(is:ie,9)
      ice_g(tile)%v(:)       =icein(is:ie,10)
    end where
    
    ! fix for old restart files
    !if ( any( water_g(tile)%sal<-1.e-4 .or. water_g(tile)%sal>100. ) ) then
    !  write(6,*) "WARN: Salinity is out-of-range in MLO-load"
    !  write(6,*) "minval,maxval ",minval(water_g(tile)%sal),maxval(water_g(tile)%sal)
    !  water_g(tile)%sal = min( max( water_g(tile)%sal, 0. ), maxsal )
    !end if

    call mlocheck("MLO-load",water_temp=water_g(tile)%temp,water_sal=water_g(tile)%sal, &
                  water_u=water_g(tile)%u,water_v=water_g(tile)%v,                      &
                  ice_tsurf=ice_g(tile)%tsurf,ice_temp=ice_g(tile)%temp)
    
  end if  

end do ! tile loop

return
end subroutine mloload

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Save MLO data

subroutine mlosave(dataout,depout,shout,iceout,diag)

implicit none

integer, intent(in) :: diag
integer ii, tile, is, ie
real, dimension(ifull,wlev,6), intent(inout) :: dataout
real, dimension(ifull,10), intent(inout) :: iceout
real, dimension(ifull), intent(inout) :: depout,shout

if (diag>=1) write(6,*) "Save MLO data"

iceout(:,8)=0.
depout=0.
shout=0.

if (.not.mlo_active) return

do tile = 1,ntiles
  is = (tile-1)*imax + 1
  ie = tile*imax
  
  if ( depth_g(tile)%data_allocated ) then
  
    do ii=1,wlev
      where ( depth_g(tile)%dz(:,ii)>=1.e-4 )
        dataout(is:ie,ii,1)=water_g(tile)%temp(:,ii)
        dataout(is:ie,ii,2)=water_g(tile)%sal(:,ii)
        dataout(is:ie,ii,3)=water_g(tile)%u(:,ii)
        dataout(is:ie,ii,4)=water_g(tile)%v(:,ii)
        dataout(is:ie,ii,5)=turb_g(tile)%k(:,ii)
        dataout(is:ie,ii,6)=turb_g(tile)%eps(:,ii)
      end where
    end do
    where ( depth_g(tile)%dz(:,1)>=1.e-4 )
      shout(is:ie)      =water_g(tile)%eta(:)
      iceout(is:ie,1)   =ice_g(tile)%tsurf(:)
      iceout(is:ie,2)   =ice_g(tile)%temp(:,0)
      iceout(is:ie,3)   =ice_g(tile)%temp(:,1)
      iceout(is:ie,4)   =ice_g(tile)%temp(:,2)
      iceout(is:ie,5)   =ice_g(tile)%fracice(:)
      iceout(is:ie,6)   =ice_g(tile)%thick(:)
      iceout(is:ie,7)   =ice_g(tile)%snowd(:)
      iceout(is:ie,8)   =ice_g(tile)%store(:)
      iceout(is:ie,9)   =ice_g(tile)%u(:)
      iceout(is:ie,10)  =ice_g(tile)%v(:)
      depout(is:ie)     =depth_g(tile)%depth_hl(:,wlev+1)
    end where  

    call mlocheck("MLO-save",water_temp=water_g(tile)%temp,water_sal=water_g(tile)%sal, &
                  water_u=water_g(tile)%u,water_v=water_g(tile)%v,                      &
                  ice_tsurf=ice_g(tile)%tsurf,ice_temp=ice_g(tile)%temp)
    
  end if  

end do ! tile loop

return
end subroutine mlosave

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Import sst for nudging

subroutine mlo_import_ifull(mode,sst,ilev,diag)

implicit none

integer, intent(in) :: ilev,diag
integer tile, is, ie
real, dimension(:), intent(in) :: sst
character(len=*), intent(in) :: mode

if (diag>=1) write(6,*) "Import MLO data"
if (.not.mlo_active) return

do tile = 1,ntiles
  is = (tile-1)*imax + 1
  ie = tile*imax
  call mlo_import_imax(mode,sst(is:ie),ilev,diag,water_g(tile),depth_g(tile))
end do
      
return
end subroutine mlo_import_ifull

subroutine mlo_import_imax(mode,sst,ilev,diag,water,depth)

implicit none

integer, intent(in) :: ilev, diag
real, dimension(imax), intent(in) :: sst
type(waterdata), intent(inout) :: water
type(depthdata), intent(in) :: depth
character(len=*), intent(in) :: mode

if (diag>=2) write(6,*) "THREAD: Import MLO data"
if ( .not.mlo_active ) return
if ( .not.depth%data_allocated ) return

select case(mode)
  case("temp")
    where ( depth%dz(:,ilev)>=1.e-4 )
      water%temp(:,ilev)=sst
    elsewhere
      water%temp(:,ilev) = 288.-wrtemp
    end where
  case("sal")
    where( depth%dz(:,ilev)>=1.e-4 )
      water%sal(:,ilev)=sst
    elsewhere
      water%sal(:,ilev) = 34.72
    end where
  case("u")
    where ( depth%dz(:,ilev)>=1.e-4 )
      water%u(:,ilev)=sst
    elsewhere
      water%u(:,ilev) = 0.
    end where
  case("v")
    where ( depth%dz(:,ilev)>=1.e-4 )
      water%v(:,ilev)=sst
    elsewhere
      water%v(:,ilev) = 0.
    end where
  case("w")
    where ( depth%dz(:,ilev)>=1.e-4 )
      water%w(:,ilev)=sst
    elsewhere
      water%w(:,ilev) = 0.
    end where    
  case("eta")
    where ( depth%dz(:,1)>=1.e-4 )  
      water%eta=sst
    end where  
  case("utop")
    where ( depth%dz(:,1)>=1.e-4 )  
      water%utop=sst
    end where  
  case("vtop")
    where ( depth%dz(:,1)>=1.e-4 )  
      water%vtop=sst
    end where  
  case("ubot")
    where ( depth%dz(:,1)>=1.e-4 )  
      water%ubot=sst
    end where  
  case("vbot")
    where ( depth%dz(:,1)>=1.e-4 )  
      water%vbot=sst
    end where  
  case("u_ema")
    where ( depth%dz(:,ilev)>=1.e-4 )
      water%u_ema(:,ilev)=sst
    elsewhere
      water%u_ema(:,ilev) = 0.
    end where
  case("v_ema")
    where ( depth%dz(:,ilev)>=1.e-4 )
      water%v_ema(:,ilev)=sst
    elsewhere
      water%v_ema(:,ilev) = 0.
    end where
  case("w_ema")
    where ( depth%dz(:,ilev)>=1.e-4 )
      water%w_ema(:,ilev)=sst
    elsewhere
      water%w_ema(:,ilev) = 0.
    end where
  case("temp_ema")
    where ( depth%dz(:,ilev)>=1.e-4 )
      water%temp_ema(:,ilev)=sst
    elsewhere
      water%temp_ema(:,ilev) = 0.
    end where    
  case("sal_ema")
    where ( depth%dz(:,ilev)>=1.e-4 )
      water%sal_ema(:,ilev)=sst
    elsewhere
      water%sal_ema(:,ilev) = 0.
    end where 
  case("dwdx")
    where ( depth%dz(:,ilev)>=1.e-4 )
      water%dwdx(:,ilev)=sst
    elsewhere
      water%dwdx(:,ilev) = 0.
    end where   
  case("dwdy")
    where ( depth%dz(:,ilev)>=1.e-4 )
      water%dwdy(:,ilev)=sst
    elsewhere
      water%dwdy(:,ilev) = 0.
    end where
  case default
    write(6,*) "ERROR: Invalid mode for mloimport with mode = ",trim(mode)
    stop
end select

return
end subroutine mlo_import_imax

subroutine mloimport3d(mode,sst,diag)

implicit none

integer, intent(in) :: mode,diag
integer ii, tile, is, ie
real, dimension(:,:), intent(in) :: sst

if (diag>=1) write(6,*) "Import 3D MLO data"
if (.not.mlo_active) return

select case(mode)
  case(0)
    do tile = 1,ntiles
      is = (tile-1)*imax + 1
      ie = tile*imax  
      if ( depth_g(tile)%data_allocated ) then
        do ii = 1,wlev
          where ( depth_g(tile)%dz(:,ii)>=1.e-4 )
            water_g(tile)%temp(:,ii) = sst(is:ie,ii)
          elsewhere
            water_g(tile)%temp(:,ii) = 288.-wrtemp
          end where
        end do  
      end if  
    end do
  case(1)
    do tile = 1,ntiles
      is = (tile-1)*imax + 1
      ie = tile*imax
      if ( depth_g(tile)%data_allocated ) then
        do ii = 1,wlev
          where ( depth_g(tile)%dz(:,ii)>=1.e-4 )
            water_g(tile)%sal(:,ii) = sst(is:ie,ii)
          elsewhere
            water_g(tile)%sal(:,ii) = 34.72
          end where
        end do  
      end if  
    end do
  case(2)
    do tile = 1,ntiles
      is = (tile-1)*imax + 1
      ie = tile*imax
      if ( depth_g(tile)%data_allocated ) then
        do ii = 1,wlev
          where ( depth_g(tile)%dz(:,ii)>=1.e-4 )
            water_g(tile)%u(:,ii) = sst(is:ie,ii)
          elsewhere
            water_g(tile)%u(:,ii) = 0.
          end where
        end do  
      end if  
    end do
  case(3)
    do tile = 1,ntiles
      is = (tile-1)*imax + 1
      ie = tile*imax
      if ( depth_g(tile)%data_allocated ) then
        do ii = 1,wlev
          where ( depth_g(tile)%dz(:,ii)>=1.e-4 )
            water_g(tile)%v(:,ii) = sst(is:ie,ii)
          elsewhere
            water_g(tile)%v(:,ii) = 0.
          end where
        end do  
      end if  
    end do
end select

return
end subroutine mloimport3d


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Import ice temp

subroutine mlo_impice_ifull(mode,tsn,diag)

implicit none

integer, intent(in) :: diag
integer tile, is, ie
real, dimension(:), intent(in) :: tsn
character(len=*), intent(in) :: mode

if (diag>=1) write(6,*) "Import MLO ice data"
if (.not.mlo_active) return

do tile = 1,ntiles
  is = (tile-1)*imax + 1
  ie = tile*imax
  call mlo_impice_imax(mode,tsn(is:ie),diag,ice_g(tile),depth_g(tile))
end do

return
end subroutine mlo_impice_ifull

subroutine mlo_impice_imax(mode,tsn,diag,ice,depth)

implicit none

integer, intent(in) :: diag
real, dimension(imax), intent(in) :: tsn
type(icedata), intent(inout) :: ice
type(depthdata), intent(in) :: depth
character(len=*), intent(in) :: mode

if (diag>=2) write(6,*) "THREAD: Import MLO ice data"
if ( .not.mlo_active ) return
if ( .not.depth%data_allocated ) return

select case(mode)
  case("tsurf")
    where ( depth%dz(:,1)>=1.e-4 )
      ice%tsurf=tsn
    end where  
  case("temp0")
    where ( depth%dz(:,1)>=1.e-4 )  
      ice%temp(:,0)=tsn
    end where  
  case("temp1")
    where ( depth%dz(:,1)>=1.e-4 )  
      ice%temp(:,1)=tsn
    end where  
  case("temp2")
    where ( depth%dz(:,1)>=1.e-4 )  
      ice%temp(:,2)=tsn
    end where  
  case("fracice")
    where ( depth%dz(:,1)>=1.e-4 )  
      ice%fracice=tsn
    end where  
  case("thick")
    where ( depth%dz(:,1)>=1.e-4 )  
      ice%thick=tsn
    end where  
  case("snowd")
    where ( depth%dz(:,1)>=1.e-4 )  
      ice%snowd=tsn
    end where  
  case("store")
    where ( depth%dz(:,1)>=1.e-4 )  
      ice%store=tsn
    end where  
  case("u")
    where ( depth%dz(:,1)>=1.e-4 )  
      ice%u=tsn
    end where  
  case("v")
    where ( depth%dz(:,1)>=1.e-4 )  
      ice%v=tsn
    end where  
  case DEFAULT
    write(6,*) "ERROR: Invalid mode ",trim(mode)
    stop
end select

return
end subroutine mlo_impice_imax

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Export sst for nudging

subroutine mlo_export_ifull(mode,sst,ilev,diag)

implicit none

integer, intent(in) :: ilev,diag
integer tile, is, ie
real, dimension(:), intent(inout) :: sst
character(len=*), intent(in) :: mode

if (diag>=1) write(6,*) "Export MLO SST data"
if (.not.mlo_active) return

do tile = 1,ntiles
  is = (tile-1)*imax + 1
  ie = tile*imax
  call mlo_export_imax(mode,sst(is:ie),ilev,diag,water_g(tile),depth_g(tile))
end do

return
end subroutine mlo_export_ifull

subroutine mlo_export_imax(mode,sst,ilev,diag,water,depth)

implicit none

integer, intent(in) :: ilev,diag
real, dimension(imax), intent(inout) :: sst
type(waterdata), intent(in) :: water
type(depthdata), intent(in) :: depth
character(len=*), intent(in) :: mode

if (diag>=2) write(6,*) "THREAD: Export MLO SST data"
if ( .not.mlo_active ) return
if ( .not.depth%data_allocated ) return

select case(mode)
  case("temp")
    where ( depth%dz(:,1)>=1.e-4 )
      sst=water%temp(:,ilev)
    end where
  case("sal")
    where ( depth%dz(:,1)>=1.e-4 )
      sst=water%sal(:,ilev)
    end where  
  case("u")
    where ( depth%dz(:,1)>=1.e-4 )
      sst=water%u(:,ilev)
    end where
  case("v")
    where ( depth%dz(:,1)>=1.e-4 )
      sst=water%v(:,ilev)
    end where 
  case("w")
    where ( depth%dz(:,1)>=1.e-4 )
      sst=water%w(:,ilev)
    end where    
  case("eta")
    where ( depth%dz(:,1)>=1.e-4 )  
      sst=water%eta
    end where   
  case("utop")
    where ( depth%dz(:,1)>=1.e-4 )
      sst=water%utop
    end where  
  case("vtop")
    where ( depth%dz(:,1)>=1.e-4 )  
      sst=water%vtop
    end where
  case("ubot")
    where ( depth%dz(:,1)>=1.e-4 )
      sst=water%ubot
    end where  
  case("vbot")
    where ( depth%dz(:,1)>=1.e-4 )  
      sst=water%vbot
    end where  
  case("ibot")
    where ( depth%dz(:,1)>=1.e-4 )      
      sst=real(depth%ibot)
    end where  
  case("u_ema")
    where ( depth%dz(:,1)>=1.e-4 )
      sst=water%u_ema(:,ilev)
    end where
  case("v_ema")
    where ( depth%dz(:,1)>=1.e-4 )
      sst=water%v_ema(:,ilev)
    end where  
  case("w_ema")
    where ( depth%dz(:,1)>=1.e-4 )
      sst=water%w_ema(:,ilev)
    end where
  case("temp_ema")
    where ( depth%dz(:,1)>=1.e-4 )
      sst=water%temp_ema(:,ilev)
    end where
  case("sal_ema")
    where ( depth%dz(:,1)>=1.e-4 )
      sst=water%sal_ema(:,ilev)
    end where
  case DEFAULT
    write(6,*) "ERROR: Invalid mode ",trim(mode)
    stop    
end select

return
end subroutine mlo_export_imax

subroutine mlo_surf_ifull(mode,sst,diag)

implicit none

integer, intent(in) :: diag
integer tile, is, ie
real, dimension(:), intent(inout) :: sst
character(len=*), intent(in) :: mode

if (diag>=1) write(6,*) "Export MLO surf data"
if (.not.mlo_active) return

do tile = 1,ntiles
  is = (tile-1)*imax + 1
  ie = tile*imax
  call mlo_surf_imax(mode,sst(is:ie),diag,water_g(tile),ice_g(tile),depth_g(tile))
end do

return
end subroutine mlo_surf_ifull

subroutine mlo_surf_imax(mode,sst,diag,water,ice,depth)

implicit none

integer, intent(in) :: diag
real, dimension(imax), intent(inout) :: sst
type(waterdata), intent(in) :: water
type(icedata), intent(in) :: ice
type(depthdata), intent(in) :: depth
real, dimension(imax) :: workb
character(len=*), intent(in) :: mode

if (diag>=2) write(6,*) "THREAD: Export MLO surf data"
if ( .not.mlo_active ) return
if ( .not.depth%data_allocated ) return

select case(mode)
  case("sst")
    workb = emisice**0.25*ice%tsurf
    where ( depth%dz(:,1)>1.e-4 )
      sst   = (1.-ice%fracice)*(water%temp(:,1)+wrtemp)+ice%fracice*workb
    end where  
end select

return
end subroutine mlo_surf_imax

subroutine mloexport3d(mode,sst,diag)

implicit none

integer, intent(in) :: mode,diag
integer ii, tile, is, ie
real, dimension(:,:), intent(inout) :: sst

if (diag>=1) write(6,*) "Export 3D MLO data"
if (.not.mlo_active) return

select case(mode)
  case(0)
    do tile = 1,ntiles
      is = (tile-1)*imax + 1
      ie = tile*imax
      if ( depth_g(tile)%data_allocated ) then
        do ii = 1,wlev 
          where ( depth_g(tile)%dz(:,1)>=1.e-4 ) 
            sst(is:ie,ii)=water_g(tile)%temp(:,ii)
          end where
        end do  
      end if  
    end do
  case(1)
    do tile = 1,ntiles
      is = (tile-1)*imax + 1
      ie = tile*imax
      if ( depth_g(tile)%data_allocated ) then
        do ii = 1,wlev 
          where ( depth_g(tile)%dz(:,1)>=1.e-4 )
            sst(is:ie,ii)=water_g(tile)%sal(:,ii)
          end where
        end do
      end if  
    end do
  case(2)
    do tile = 1,ntiles
      is = (tile-1)*imax + 1
      ie = tile*imax
      if ( depth_g(tile)%data_allocated ) then
        do ii = 1,wlev 
          where ( depth_g(tile)%dz(:,1)>=1.e-4 )
            sst(is:ie,ii)=water_g(tile)%u(:,ii)
          end where
        end do
      end if  
    end do
  case(3)
    do tile = 1,ntiles
      is = (tile-1)*imax + 1
      ie = tile*imax
      if ( depth_g(tile)%data_allocated ) then
        do ii = 1,wlev  
          where ( depth_g(tile)%dz(:,1)>=1.e-4 )
            sst(is:ie,ii)=water_g(tile)%v(:,ii)
          end where
        end do
      end if  
    end do
end select

return
end subroutine mloexport3d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Export ice temp

subroutine mlo_expice_ifull(mode,tsn,diag)

implicit none

integer, intent(in) :: diag
integer tile, is, ie
real, dimension(:), intent(inout) :: tsn
character(len=*), intent(in) :: mode

if (diag>=1) write(6,*) "Export MLO ice data"
if (.not.mlo_active) return

do tile = 1,ntiles
  is = (tile-1)*imax + 1
  ie = tile*imax
  call mlo_expice_imax(mode,tsn(is:ie),diag,ice_g(tile),depth_g(tile))
end do

return
end subroutine mlo_expice_ifull

subroutine mlo_expice_imax(mode,tsn,diag,ice,depth)

implicit none

integer, intent(in) :: diag
real, dimension(imax), intent(inout) :: tsn
type(icedata), intent(in) :: ice
type(depthdata), intent(in) :: depth
character(len=*), intent(in) :: mode

if (diag>=2) write(6,*) "THREAD: Export MLO ice data"
if ( .not.mlo_active ) return
if ( .not.depth%data_allocated ) return

select case(mode)
  case("tsurf")
    where ( depth%dz(:,1)>1.e-4 )
      tsn=ice%tsurf
    end where  
  case("temp0")
    where ( depth%dz(:,1)>1.e-4 )  
      tsn=ice%temp(:,0)
    end where  
  case("temp1")
    where ( depth%dz(:,1)>1.e-4 )  
      tsn=ice%temp(:,1)
    end where  
  case("temp2")
    where ( depth%dz(:,1)>1.e-4 )  
      tsn=ice%temp(:,2)
    end where  
  case("fracice")
    where ( depth%dz(:,1)>1.e-4 )  
      tsn=ice%fracice
    end where
  case("thick")
    where ( depth%dz(:,1)>1.e-4 )  
      tsn=ice%thick
    end where
  case("snowd")
    where ( depth%dz(:,1)>1.e-4 )  
      tsn=ice%snowd
    end where
  case("store")
    where ( depth%dz(:,1)>1.e-4 )  
      tsn=ice%store
    end where
  case("u")
    where ( depth%dz(:,1)>1.e-4 )  
      tsn=ice%u
    end where
  case("v")
    where ( depth%dz(:,1)>1.e-4 )  
      tsn=ice%v
    end where
  case DEFAULT
    write(6,*) "ERROR: Invalid mode ",mode
    stop
end select

return
end subroutine mlo_expice_imax

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Return mixed layer depth

subroutine mlodiag_old(mld,diag)

implicit none

integer, intent(in) :: diag
integer tile, is, ie
real, dimension(:), intent(out) :: mld

if (diag>=1) write(6,*) "Export MLO mixed layer depth"
mld=0.
if (.not.mlo_active) return

do tile = 1,ntiles
  is = (tile-1)*imax + 1
  ie = tile*imax
  if ( depth_g(tile)%data_allocated ) then
    where ( depth_g(tile)%dz(:,1)>1.e-4 )
      mld(is:ie)=dgwater_g(tile)%mixdepth
    end where    
  end if  
end do

return
end subroutine mlodiag_old

subroutine mlo_diag_ifull(mode,mld,ilev,diag)

implicit none

integer, intent(in) :: ilev,diag
integer tile, is, ie
real, dimension(:), intent(inout) :: mld
character(len=*), intent(in) :: mode

if (diag>=1) write(6,*) "Export MLO diagnostic data"
if (.not.mlo_active) return

do tile = 1,ntiles
  is = (tile-1)*imax + 1
  ie = tile*imax
  call mlo_diag_imax(mode,mld(is:ie),ilev,diag,dgwater_g(tile),depth_g(tile))
end do

return
end subroutine mlo_diag_ifull

subroutine mlo_diag_imax(mode,mld,ilev,diag,dgwater,depth)

implicit none

integer, intent(in) :: ilev, diag
integer :: iqw
real, dimension(imax), intent(inout) :: mld
real, dimension(imax) :: work
type(dgwaterdata), intent(in) :: dgwater
type(depthdata), intent(in) :: depth
character(len=*), intent(in) :: mode

if (diag>=2) write(6,*) "THREAD: Export MLO ice data"
mld=0.
if ( .not.mlo_active ) return
if ( .not.depth%data_allocated ) return

select case(mode)
  case("mixdepth")
    where ( depth%dz(:,1)>1.e-4 )  
      mld=dgwater%mixdepth
    end where  
  case("rad")
    where ( depth%dz(:,1)>=1.e-4 )
      mld=dgwater%rad(:,ilev)
    end where
  case("t0")
    if ( incradgam>0 ) then
      ! include radiation in counter-gradient term
      do iqw = 1,imax
        work(iqw) = dgwater%wt0(iqw) + sum(dgwater%rad(iqw,1:dgwater%mixind(iqw)))
      end do
      where ( depth%dz(:,1)>=1.e-4 )
        mld = work
      end where  
    else
      where ( depth%dz(:,1)>=1.e-4 )  
        mld = dgwater%wt0
      end where  
    end if
  case("cdh")
    where ( depth%dz(:,1)>=1.e-4 )  
      mld = dgwater%cdh
    end where  
  case("cd")
    where ( depth%dz(:,1)>=1.e-4 )  
      mld = dgwater%cd
    end where  
  case("cd_bot")
    where ( depth%dz(:,1)>=1.e-4 )  
      mld = dgwater%cd_bot
    end where  
  case("wt0")
    where ( depth%dz(:,1)>=1.e-4 )  
      mld = dgwater%wt0
    end where  
  case("wt0_rad")
    where ( depth%dz(:,1)>=1.e-4 )  
      mld = dgwater%wt0_rad
    end where  
  case("wt0_melt")
    where ( depth%dz(:,1)>=1.e-4 )  
      mld = dgwater%wt0_melt
    end where  
  case("wt0_eg")
    where ( depth%dz(:,1)>=1.e-4 )  
      mld = dgwater%wt0_eg
    end where  
  case("wt0_fb")
    where ( depth%dz(:,1)>=1.e-4 )  
      mld = dgwater%wt0_fb
    end where  
  case("ws0")
    where ( depth%dz(:,1)>=1.e-4 )  
      mld = dgwater%ws0
    end where  
  case("ws0_subsurf")
    where ( depth%dz(:,1)>=1.e-4 )  
      mld = dgwater%ws0_subsurf
    end where  
  case default
    write(6,*) "ERROR: Invalid mode for mlodiag with mode = ",trim(mode)
    stop
end select

return
end subroutine mlo_diag_imax

subroutine mlo_update_diag_imax(mode,dat,diag,dgwater,depth)

implicit none

integer, intent(in) :: diag
real, dimension(imax), intent(in) :: dat
type(dgwaterdata), intent(inout) :: dgwater
type(depthdata), intent(in) :: depth
character(len=*), intent(in) :: mode

if (diag>=2) write(6,*) "THREAD: Import MLO data"
if ( .not.mlo_active ) return
if ( .not.depth%data_allocated ) return

select case(mode)
  case("wt0")
    where ( depth%dz(:,1)>=1.e-4 )  
      dgwater%wt0=dat
    end where  
  case("wu0")
    where ( depth%dz(:,1)>=1.e-4 )  
      dgwater%wu0=dat
    end where  
  case("wv0")
    where ( depth%dz(:,1)>=1.e-4 )  
      dgwater%wv0=dat
    end where  
end select

return
end subroutine mlo_update_diag_imax

subroutine mlo_diagice_ifull(mode,mld,diag)

implicit none

integer, intent(in) :: diag
integer tile, is, ie
real, dimension(:), intent(inout) :: mld
character(len=*), intent(in) :: mode

if (diag>=1) write(6,*) "Export MLO diagnostic data"
if (.not.mlo_active) return

do tile = 1,ntiles
  is = (tile-1)*imax + 1
  ie = tile*imax
  call mlo_diagice_imax(mode,mld(is:ie),diag,dgice_g(tile),depth_g(tile))
end do

return
end subroutine mlo_diagice_ifull

subroutine mlo_diagice_imax(mode,mld,diag,dgice,depth)

implicit none

integer, intent(in) :: diag
real, dimension(imax), intent(inout) :: mld
type(dgicedata), intent(in) :: dgice
type(depthdata), intent(in) :: depth
character(len=*), intent(in) :: mode

if (diag>=2) write(6,*) "THREAD: Export MLO ice data"
mld = 0.
if ( .not.mlo_active) return
if ( .not.depth%data_allocated ) return

select case(mode)
  case("fg")
    where ( depth%dz(:,1)>=1.e-4 )
      mld = dgice%fg
    end where  
  case("mass")
    where ( depth%dz(:,1)>=1.e-4 )  
      mld = dgice%imass
    end where  
  case("cd")
    where ( depth%dz(:,1)>=1.e-4 )  
      mld = dgice%cd
    end where  
  case("cd_bot")
    where ( depth%dz(:,1)>=1.e-4 )  
      mld = dgice%cd_bot
    end where  
  case default
    write(6,*) "ERROR: Invalid mode for mlodiagice with mode = ",trim(mode)
    stop
end select

return
end subroutine mlo_diagice_imax

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Return roughness length for heat

subroutine mloextra_ifull(mode,zoh,zmin,diag)

implicit none

integer, intent(in) :: mode,diag
integer tile, is, ie
real, dimension(:), intent(out) :: zoh
real, dimension(:), intent(in) :: zmin

if (diag>=1) write(6,*) "Export additional MLO data"
zoh=0.
if (.not.mlo_active) return

do tile = 1,ntiles
  is = (tile-1)*imax + 1
  ie = tile*imax
  call mloextra_imax(mode,zoh(is:ie),zmin(is:ie),diag,          &
                     dgwater_g(tile),dgice_g(tile),ice_g(tile), &
                     depth_g(tile))
end do

return
end subroutine mloextra_ifull

subroutine mloextra_imax(mode,zoh,zmin,diag,dgwater,dgice,ice,depth)

implicit none

integer, intent(in) :: mode,diag
real, dimension(imax), intent(out) :: zoh
real, dimension(imax), intent(in) :: zmin
type(dgwaterdata), intent(in) :: dgwater
type(dgicedata), intent(in) :: dgice
type(icedata), intent(in) :: ice
type(depthdata), intent(in) :: depth
real, dimension(imax) :: workb,workc
real, dimension(imax) :: dumazmin
real zohseaice, zoqseaice

if (diag>=2) write(6,*) "THREAD: Export additional MLO data"
zoh=0.
if ( .not.mlo_active ) return
if ( .not.depth%data_allocated ) return

select case(mode)
  case(0) ! zoh
    zohseaice=zoseaice/(factchseaice*factchseaice)  
    dumazmin=max(zmin,dgwater%zo+0.2,zoseaice+0.2)
    workb=(1.-ice%fracice)/(log(dumazmin/dgwater%zo)*log(dumazmin/dgwater%zoh)) &
         +ice%fracice/(log(dumazmin/zoseaice)*log(dumazmin/zohseaice))
    workc=(1.-ice%fracice)/log(dumazmin/dgwater%zo)**2+ice%fracice/log(dumazmin/zoseaice)**2
    workc=sqrt(workc)
    where ( depth%dz(:,1)>=1.e-4 )
      zoh=dumazmin*exp(-workc/workb)
    end where  
  case(1) ! taux
    workb=(1.-ice%fracice)*dgwater%taux+ice%fracice*dgice%tauxica
    where ( depth%dz(:,1)>=1.e-4 )
      zoh=workb
    end where  
  case(2) ! tauy
    workb=(1.-ice%fracice)*dgwater%tauy+ice%fracice*dgice%tauyica
    where ( depth%dz(:,1)>=1.e-4 )
      zoh=workb
    end where  
  case(3) ! zoq
    zoqseaice=zoseaice/(factchseaice*factchseaice)  
    dumazmin=max(zmin,dgwater%zo+0.2,zoseaice+0.2)
    workb=(1.-ice%fracice)/(log(dumazmin/dgwater%zo)*log(dumazmin/dgwater%zoq)) &
         +ice%fracice/(log(dumazmin/zoseaice)*log(dumazmin/zoqseaice))
    workc=(1.-ice%fracice)/log(dumazmin/dgwater%zo)**2+ice%fracice/log(dumazmin/zoseaice)**2
    workc=sqrt(workc)
    where ( depth%dz(:,1)>=1.e-4 )
      zoh=dumazmin*exp(-workc/workb)
    end where  
  case default
    write(6,*) "ERROR: Invalid mode ",mode
    stop
end select

return
end subroutine mloextra_imax

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Calculate screen diagnostics

subroutine mloscrnout(tscrn,qgscrn,uscrn,u10,diag)

implicit none

integer, intent(in) :: diag
integer tile, is, ie
real, dimension(:), intent(inout) :: tscrn,qgscrn,uscrn,u10

if (diag>=1) write(6,*) "Export MLO 2m diagnostics"
if (.not.mlo_active) return

do tile = 1,ntiles
  is = (tile-1)*imax + 1
  ie = tile*imax
  if ( depth_g(tile)%data_allocated ) then
    where ( depth_g(tile)%dz(:,1)>=1.e-4 )
      tscrn(is:ie) =dgscrn_g(tile)%temp
      qgscrn(is:ie)=dgscrn_g(tile)%qg
      uscrn(is:ie) =dgscrn_g(tile)%u2
      u10(is:ie)   =dgscrn_g(tile)%u10
    end where  
  end if  
end do

return
end subroutine mloscrnout

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Calculate albedo data (VIS + NIR)

subroutine mloalb2(istart,ifin,coszro,ovisalb,oniralb,diag)

implicit none

integer, intent(in) :: istart,ifin,diag
integer ifinish,ib,ie
integer tile, js, je, kstart, kfinish, jstart, jfinish
real, dimension(ifin), intent(in) :: coszro
real, dimension(ifin), intent(inout) :: ovisalb,oniralb
real, dimension(imax) :: watervis,waternir,icevis,icenir
real, dimension(imax) :: snow

if (diag>=1) write(6,*) "Export MLO albedo data vis/nir"
if (.not.mlo_active) return

ifinish = istart + ifin - 1

do tile = 1,ntiles
  js = (tile-1)*imax + 1 ! js:je is the tile portion of 1:ifull
  je = tile*imax         ! js:je is the tile portion of 1:ifull
  
  if ( depth_g(tile)%data_allocated ) then
      
    kstart = max( istart - js + 1, 1)      ! kstart:kfinish is the requested portion of 1:imax
    kfinish = min( ifinish - js + 1, imax) ! kstart:kfinish is the requested portion of 1:imax
    if ( kstart<=kfinish ) then
      jstart = kstart + js - istart        ! jstart:jfinish is the tile portion of 1:ifin
      jfinish = kfinish + js - istart      ! jstart:jfinish is the tile portion of 1:ifin
      ib = kstart
      ie = kfinish
        
      !pond(ib:ie)=max(1.+.008*min(ice(tile)%tsurf(ib:ie)-273.16,0.),0.)
      !snow(ib:ie)=min(max(ice(tile)%snowd(ib:ie)/0.05,0.),1.)
      snow(ib:ie)=0.

      watervis(ib:ie)=.05/(coszro(jstart:jfinish)+0.15)
      waternir(ib:ie)=.05/(coszro(jstart:jfinish)+0.15)
      ! need to factor snow age into albedo
      !icevis(ib:ie)=(alphavis_seaice*(1.-pond(ib:ie))+0.75*pond(ib:ie))*(1.-snow(ib:ie)) &
      !    +(alphavis_seasnw*(1.-pond(ib:ie))+0.85*pond(ib:ie))*snow(ib:ie)
      !icenir(ib:ie)=(alphanir_seaice*(1.-pond(ib:ie))+0.35*pond(ib:ie))*(1.-snow(ib:ie)) &
      !    +(alphanir_seasnw*(1.-pond(ib:ie))+0.55*pond(ib:ie))*snow(ib:ie)
      icevis(ib:ie)=alphavis_seaice*(1.-snow(ib:ie))+alphavis_seasnw*snow(ib:ie)
      icenir(ib:ie)=alphanir_seaice*(1.-snow(ib:ie))+alphanir_seasnw*snow(ib:ie)

      where ( depth_g(tile)%dz(ib:ie,1)>=1.e-4 )
        ovisalb(jstart:jfinish)=icevis(ib:ie)*ice_g(tile)%fracice(ib:ie)               &
                               +(1.-ice_g(tile)%fracice(ib:ie))*watervis(ib:ie)
        oniralb(jstart:jfinish)=icenir(ib:ie)*ice_g(tile)%fracice(ib:ie)               &
                               +(1.-ice_g(tile)%fracice(ib:ie))*waternir(ib:ie)

        dgwater_g(tile)%visdiralb(ib:ie)=watervis(ib:ie)
        dgwater_g(tile)%visdifalb(ib:ie)=watervis(ib:ie)
        dgwater_g(tile)%nirdiralb(ib:ie)=waternir(ib:ie)
        dgwater_g(tile)%nirdifalb(ib:ie)=waternir(ib:ie)
        dgice_g(tile)%visdiralb(ib:ie)=icevis(ib:ie)
        dgice_g(tile)%visdifalb(ib:ie)=icevis(ib:ie)
        dgice_g(tile)%nirdiralb(ib:ie)=icenir(ib:ie)
        dgice_g(tile)%nirdifalb(ib:ie)=icenir(ib:ie)
      end where  

    end if  
    
  end if
  
end do ! tile loop

return
end subroutine mloalb2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Calculate albedo data (VIS-DIR, VIS-DIF, NIR-DIR & NIR-DIF)

subroutine mloalb4(istart,ifin,coszro,ovisdir,ovisdif,onirdir,onirdif,diag)

implicit none

integer, intent(in) :: istart,ifin,diag
integer ifinish,ib,ie
integer tile, js, je, kstart, kfinish, jstart, jfinish
real, dimension(ifin), intent(in) :: coszro
real, dimension(ifin), intent(inout) :: ovisdir,ovisdif,onirdir,onirdif
real, dimension(imax) :: snow

if (diag>=1) write(6,*) "Export MLO albedo vis/nir/dir/dif"
if (.not.mlo_active) return

ifinish = istart + ifin - 1 ! istart:ifinish is the requested portion of 1:ifull

do tile = 1,ntiles
  js = (tile-1)*imax + 1 ! js:je is the tile portion of 1:ifull
  je = tile*imax         ! js:je is the tile portion of 1:ifull
    
  if ( depth_g(tile)%data_allocated ) then
  
    kstart = max( istart - js + 1, 1)      ! kstart:kfinish is the requested portion of 1:imax
    kfinish = min( ifinish - js + 1, imax) ! kstart:kfinish is the requested portion of 1:imax
    if ( kstart<=kfinish ) then
      jstart = kstart + js - istart        ! jstart:jfinish is the tile portion of 1:ifin
      jfinish = kfinish + js - istart      ! jstart:jfinish is the tile portion of 1:ifin
      ib = kstart
      ie = kfinish

      !pond(ib:ie)=max(1.+.008*min(ice(tile)%tsurf(ib:ie)-273.16,0.),0.)
      !snow(ib:ie)=min(max(ice(tile)%snowd(ib:ie)/0.05,0.),1.)
      snow(ib:ie)=0.

      where ( coszro(jstart:jfinish)>1.e-8 .and. depth_g(tile)%dz(ib:ie,1)>=1.e-4 )
        dgwater_g(tile)%visdiralb(ib:ie)=0.026/(coszro(jstart:jfinish)**1.7+0.065)+0.15*(coszro(jstart:jfinish)-0.1)* &
                        (coszro(jstart:jfinish)-0.5)*(coszro(jstart:jfinish)-1.)
      elsewhere ( depth_g(tile)%dz(ib:ie,1)>=1.e-4 )
        dgwater_g(tile)%visdiralb(ib:ie)=0.3925
      end where
    
      where ( depth_g(tile)%dz(ib:ie,1)>=1.e-4 )
        dgwater_g(tile)%visdifalb(ib:ie)=0.06
        dgwater_g(tile)%nirdiralb(ib:ie)=dgwater_g(tile)%visdiralb(ib:ie)
        dgwater_g(tile)%nirdifalb(ib:ie)=0.06
        ! need to factor snow age into albedo
        !dgice_g(tile)%visdiralb(ib:ie)=(alphavis_seaice*(1.-pond(ib:ie))+0.75*pond(ib:ie))*(1.-snow(ib:ie)) &
        !              +(alphavis_seasnw*(1.-pond(ib:ie))+0.85*pond(ib:ie))*snow(ib:ie)
        dgice_g(tile)%visdiralb(ib:ie)=alphavis_seaice*(1.-snow(ib:ie))+alphavis_seasnw*snow(ib:ie)
        dgice_g(tile)%visdifalb(ib:ie)=dgice_g(tile)%visdiralb(ib:ie)
        !dgice_g(tile)%nirdiralb(ib:ie)=(alphanir_seaice*(1.-pond(ib:ie))+0.35*pond(ib:ie))*(1.-snow(ib:ie)) &
        !              +(alphanir_seasnw*(1.-pond(ib:ie))+0.55*pond(ib:ie))*snow(ib:ie)
        dgice_g(tile)%nirdiralb(ib:ie)=alphanir_seaice*(1.-snow(ib:ie))+alphanir_seasnw*snow(ib:ie)
        dgice_g(tile)%nirdifalb(ib:ie)=dgice_g(tile)%nirdiralb(ib:ie)
      
        ovisdir(jstart:jfinish)=ice_g(tile)%fracice(ib:ie)*dgice_g(tile)%visdiralb(ib:ie) &
            +(1.-ice_g(tile)%fracice(ib:ie))*dgwater_g(tile)%visdiralb(ib:ie)
        ovisdif(jstart:jfinish)=ice_g(tile)%fracice(ib:ie)*dgice_g(tile)%visdifalb(ib:ie)  &
            +(1.-ice_g(tile)%fracice(ib:ie))*dgwater_g(tile)%visdifalb(ib:ie)
        onirdir(jstart:jfinish)=ice_g(tile)%fracice(ib:ie)*dgice_g(tile)%nirdiralb(ib:ie)  &
            +(1.-ice_g(tile)%fracice(ib:ie))*dgwater_g(tile)%nirdiralb(ib:ie)
        onirdif(jstart:jfinish)=ice_g(tile)%fracice(ib:ie)*dgice_g(tile)%nirdifalb(ib:ie)  &
            +(1.-ice_g(tile)%fracice(ib:ie))*dgwater_g(tile)%nirdifalb(ib:ie)
      end where
    
    end if  
    
  end if  
  
end do

return
end subroutine mloalb4

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Regrid MLO data

subroutine mloregrid(wlin,sigin,depin,mloin,mlodat,mode)

implicit none

integer, intent(in) :: wlin,mode
integer tile, is, ie
real, dimension(:,:), intent(in) :: sigin
real, dimension(:), intent(in) :: depin
real, dimension(:,:), intent(in) :: mloin
real, dimension(:,:), intent(inout) :: mlodat
real, dimension(imax,wlin) :: mloin_tmp, sigin_tmp
real, dimension(imax,wlev) :: mlodat_tmp

if ( .not.mlo_active ) return

do tile = 1,ntiles
  is = (tile-1)*imax + 1
  ie = tile*imax
  sigin_tmp(1:imax,1:wlin) = sigin(is:ie,1:wlin)
  mloin_tmp(1:imax,1:wlin) = mloin(is:ie,1:wlin)
  mlodat_tmp(1:imax,1:wlev) = mlodat(is:ie,1:wlev)
  call mloregrid_work(wlin,sigin_tmp,depin(is:ie),mloin_tmp,mlodat_tmp,mode, &
                      depth_g(tile))
  mlodat(is:ie,1:wlev) = mlodat_tmp(1:imax,1:wlev)
end do

return
end subroutine mloregrid

subroutine mloregrid_work(wlin,sig_tmp,depin,mloin,mlodat,mode, &
                          depth)

implicit none

integer, intent(in) :: wlin,mode
integer iqw,ii,jj,jj_found,pos(1)
real, dimension(imax,wlin), intent(in) :: sig_tmp
real, dimension(imax), intent(in) :: depin
real, dimension(imax,wlin), intent(in) :: mloin
real, dimension(imax,wlev), intent(inout) :: mlodat
real, dimension(wlin) :: dpin, sigin_tmp
real, dimension(wlev) :: sig
real, dimension(imax,wlin) :: sigin
real x
type(depthdata), intent(in) :: depth

if ( .not.mlo_active ) return
if ( .not.depth%data_allocated ) return

do iqw = 1,imax
  if ( sig_tmp(iqw,1)>sig_tmp(iqw,wlin) ) then
    write(6,*) "ERROR: Input levels for MLO are in reverse order"
    write(6,*) "sig_tmp ",sig_tmp(iqw,:)
    stop
  end if
end do  

if ( any(sig_tmp>1.) ) then
  ! found z* levels
  do ii = 1,wlin
    do iqw = 1,imax  
      sigin(iqw,ii) = sig_tmp(iqw,ii)/max(depin(iqw),1.e-8)
    end do  
  end do  
else
  write(6,*) "ERROR: Simga levels for MLO are no longer supported"
  stop
end if

select case(mode)
  case(0,1) ! interpolate to depth
    do iqw = 1,imax
      dpin(1:wlin) = min( sigin(iqw,1:wlin)*depin(iqw), depin(iqw) )  
      if ( wlev==wlin ) then
        if ( all( abs(depth%depth(iqw,1:wlev)-dpin(1:wlev))/depth%depth(iqw,1:wlev)<1.e-6 ) ) then
          mlodat(iqw,1:wlev) = mloin(iqw,1:wlev)
          cycle
        end if
      end if
      do ii = 1,wlev
        if ( depth%depth(iqw,ii)<=dpin(1) ) then
          mlodat(iqw,ii) = mloin(iqw,1)
        else
          ! search down column.  May have multiple levels with same depth, so
          ! we want the first level of the same depth.
          jj_found = wlin  
          do jj = 2,wlin
            if ( depth%depth(iqw,ii)<dpin(jj) ) then
              jj_found = jj  
              exit
            end if
          end do  
          if ( depth%depth(iqw,ii)<dpin(jj_found) ) then
            x = (depth%depth(iqw,ii)-dpin(jj_found-1))/max(dpin(jj_found)-dpin(jj_found-1),1.e-20)
            mlodat(iqw,ii) = mloin(iqw,jj_found)*x + mloin(iqw,jj_found-1)*(1.-x)
          else
            mlodat(iqw,ii) = mloin(iqw,wlin)
          end if
        end if
      end do
    end do
  case(2,3) ! interpolate to sigma level
    do iqw = 1,imax
      sig = depth%depth(iqw,:)/max(depth%depth_hl(iqw,wlev),1.e-20)
      if ( wlev==wlin ) then
        if ( all( abs(sig(1:wlev)-sigin(iqw,1:wlev))<1.e-6 ) ) then
          mlodat(iqw,1:wlev) = mloin(iqw,1:wlev)
          cycle
        end if
      end if
      do ii = 1,wlev
        if ( sig(ii)>=sigin(iqw,wlin) ) then
          mlodat(iqw,ii) = mloin(iqw,wlin)
        else if ( sig(ii)<=sigin(iqw,1) ) then
          mlodat(iqw,ii) = mloin(iqw,1)
        else
          sigin_tmp = sigin(iqw,:)  
          pos = maxloc(sigin_tmp,sigin_tmp<sig(ii))
          pos(1) = max(1,min(wlin-1,pos(1)))
          x = (sig(ii)-sigin(iqw,pos(1)))/max(sigin(iqw,pos(1)+1)-sigin(iqw,pos(1)),1.e-20)
          x = max(0.,min(1.,x))
          mlodat(iqw,ii) = mloin(iqw,pos(1)+1)*x+mloin(iqw,pos(1))*(1.-x)
        end if
      end do
    end do
  case default
    write(6,*) "ERROR: Unknown mode for mloregrid ",mode
    stop
end select

return
end subroutine mloregrid_work    
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Export ocean depth data

subroutine mloexpdep_old(mode,odep,ii,diag)

implicit none

integer, intent(in) :: mode,ii,diag
integer tile, is, ie
real, dimension(:), intent(out) :: odep

if (diag>=1) write(6,*) "Export MLO ocean depth data"
odep=0.
if (.not.mlo_active) return

select case(mode)
  case(0)
    do tile = 1,ntiles
      is = (tile-1)*imax + 1
      ie = tile*imax
      if ( depth_g(tile)%data_allocated ) then
        where ( depth_g(tile)%dz(:,ii)>=1.e-4 )
          odep(is:ie)=depth_g(tile)%depth(:,ii)
        end where
      end if  
    end do
  case(1)
    do tile = 1,ntiles
      is = (tile-1)*imax + 1
      ie = tile*imax
      if ( depth_g(tile)%data_allocated ) then
        where ( depth_g(tile)%dz(:,ii)>=1.e-4 )
          odep(is:ie)=depth_g(tile)%dz(:,ii)
        end where
      end if  
    end do
  case default
    write(6,*) "ERROR: Unknown mloexpdep mode ",mode
    stop
end select
  
return
end subroutine mloexpdep_old

subroutine mlo_export_depth_ifull(mode,odep,ilev,diag)

implicit none

integer, intent(in) :: ilev,diag
integer tile, is, ie
real, dimension(:), intent(out) :: odep
character(len=*), intent(in) :: mode

if (diag>=1) write(6,*) "Export MLO depth data"
if (.not.mlo_active) return

do tile = 1,ntiles
  is = (tile-1)*imax + 1
  ie = tile*imax
  call mlo_export_depth_imax(mode,odep(is:ie),ilev,diag,depth_g(tile))
end do
      
return
end subroutine mlo_export_depth_ifull

subroutine mlo_export_depth_imax(mode,odep,ilev,diag,depth)

implicit none

integer, intent(in) :: ilev,diag
integer iqw, ii
real, dimension(imax), intent(out) :: odep
type(depthdata), intent(in) :: depth
character(len=*), intent(in) :: mode

if (diag>=2) write(6,*) "THREAD: Export MLO depth data"
odep=0.
if ( .not.mlo_active ) return
if ( .not.depth%data_allocated ) return

select case(mode)
  case("depth_dz_hl")
    where ( depth%dz(:,1)>=1.e-4 )
      odep = depth%dz_hl(:,ilev)
    end where  
  case("depth_dz_fl")  
    where ( depth%dz(:,1)>=1.e-4 )  
      odep = depth%dz(:,ilev)
    end where  
  case("depth_hl")
    where ( depth%dz(:,1)>=1.e-4 )  
      odep = depth%depth_hl(:,ilev)
    end where  
  case("depth_fl")
    where ( depth%dz(:,1)>=1.e-4 )  
      odep = depth%depth(:,ilev)
    end where 
  case("depth_p")
    do iqw = 1,imax
      ii = depth%ibot(iqw)
      if ( depth%dz(iqw,ii)>=1.e-4 ) then
        odep(iqw) = depth%depth(iqw,ii)
      end if
    end do
  case default
    write(6,*) "ERROR: Unknown option for mlo_export_depth with mode = ",trim(mode)
    stop
end select
  
return
end subroutine mlo_export_depth_imax

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Extract density

subroutine mloexpdensity(odensity,alpha,beta,tt,ss,ddz,pxtr,diag,rawrho)

implicit none

integer, intent(in) :: diag
real, dimension(:,:), intent(in) :: tt
real, dimension(:), intent(in) :: pxtr
real, dimension(:,:), intent(in) :: ss,ddz
real, dimension(:,:), intent(out) :: odensity,alpha,beta
logical, intent(in), optional :: rawrho
logical rawmode

if (diag>=1) write(6,*) "Calculate MLO density"

rawmode = .false.
if ( present( rawrho ) ) then
  rawmode = rawrho
end if

call calcdensity(odensity,alpha,beta,tt,ss,ddz,pxtr)

if ( .not.rawmode ) then
  odensity = odensity + wrtrho
end if

return
end subroutine mloexpdensity

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Extract melting temperature

subroutine mloexpmelt(omelt)

implicit none

integer tile, is, ie
real, dimension(:), intent(out) :: omelt

omelt=273.16
if (.not.mlo_active) return

do tile = 1,ntiles
  is = (tile-1)*imax + 1
  ie = tile*imax
  call mloexpmelt_work(omelt(is:ie),water_g(tile),depth_g(tile))
end do

return
end subroutine mloexpmelt
    
subroutine mloexpmelt_work(omelt,water,depth)

implicit none

real, dimension(imax), intent(inout) :: omelt
real, dimension(imax) :: tmelt
type(waterdata), intent(in) :: water
type(depthdata), intent(in) :: depth

if ( .not.mlo_active ) return
if ( .not.depth%data_allocated ) return

call calcmelt(tmelt,water)
where ( depth%dz(:,1)>=1.e-4 ) 
  omelt = tmelt
end where  

return
end subroutine mloexpmelt_work

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Extract factors for energy conservation

pure subroutine mloexpgamm(gamm,ip_dic,ip_dsn,diag)

implicit none

integer, intent(in) :: diag
real, dimension(:), intent(in) :: ip_dic, ip_dsn
real, dimension(:,:), intent(out) :: gamm

gamm(1:ifull,1)=gammi
gamm(1:ifull,2)=max(ip_dsn(1:ifull),0.)*cps
gamm(1:ifull,3)=max(ip_dic(1:ifull),0.)*0.5*cpi

return
end subroutine mloexpgamm

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Export vertical levels

subroutine mlovlevels(ans,sigma)

implicit none

real, dimension(wlev), intent(out) :: ans
real, dimension(wlev+1) :: ans_hl
logical, intent(in), optional :: sigma
logical usesigma

usesigma = .false.
if ( present(sigma) ) then
  usesigma = sigma
end if

select case(mlosigma)
  case(4,5,6,7)
    call vgrid(wlev,mxd,ans,ans_hl)
    if ( usesigma ) then
      ans = ans/mxd  
    end if    
  case default
    write(6,*) "ERROR: Unknown option in mlovlevels for mlosigma ",mlosigma
    stop
end select
  
end subroutine mlovlevels  


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Pack atmospheric data for MLO eval

! calcprog = 0 update sea-ice and water
! calcprog = 1 do not update sea-ice or water
! calcprog = 2 update sea-ice only
! calcprog = 3 update sea-ice thermodynamics only (not momentum)

subroutine mloeval_standard(sst,zo,cd,cds,umod,fg,eg,evspsbl,sbl,wetfac,epot,epan,fracice, &
                   siced,snowd,dt,zmin,zmins,sg,rg,precp,precs,uatm,vatm,temp,qg,ps,       &
                   visnirratio,fbvis,fbnir,inflow,diag,calcprog)                   

implicit none

integer, intent(in) :: diag, calcprog
integer tile, is, ie
real, intent(in) :: dt
real, dimension(:), intent(in) :: sg,rg,precp,precs,uatm,vatm,temp,qg,ps,visnirratio,fbvis,fbnir,inflow,zmin,zmins
real, dimension(:), intent(inout) :: sst,zo,cd,cds,umod,fg,eg,evspsbl,sbl,wetfac,fracice,siced,epot,epan,snowd

do tile = 1,ntiles
  is = (tile-1)*imax + 1
  ie = tile*imax
  call mloeval_thread(sst(is:ie),zo(is:ie),cd(is:ie),cds(is:ie),umod(is:ie),fg(is:ie),eg(is:ie),    &
                     evspsbl(is:ie),sbl(is:ie),wetfac(is:ie),epot(is:ie),epan(is:ie),               &
                     fracice(is:ie),siced(is:ie),snowd(is:ie),dt,zmin(is:ie),zmins(is:ie),          &
                     sg(is:ie),rg(is:ie),precp(is:ie),precs(is:ie),uatm(is:ie),vatm(is:ie),         &
                     temp(is:ie),qg(is:ie),ps(is:ie),visnirratio(is:ie),fbvis(is:ie),               &
                     fbnir(is:ie),inflow(is:ie),diag,calcprog,                                      &
                     depth_g(tile),dgice_g(tile),dgscrn_g(tile),dgwater_g(tile),                    &
                     ice_g(tile),water_g(tile),turb_g(tile))
end do

return
end subroutine mloeval_standard
                   
subroutine mloeval_thread(sst,zo,cd,cds,umod,fg,eg,evspsbl,sbl,wetfac,epot,epan,fracice,   &
                   siced,snowd,dt,zmin,zmins,sg,rg,precp,precs,uatm,vatm,temp,qg,ps,       &
                   visnirratio,fbvis,fbnir,inflow,diag,calcprog,                           &
                   depth,dgice,dgscrn,dgwater,ice,water,                                   &
                   turb)                   

implicit none

integer, intent(in) :: diag, calcprog
real, intent(in) :: dt
real, dimension(imax), intent(in) :: sg,rg,precp,precs,uatm,vatm,temp,qg,ps,visnirratio,fbvis,fbnir,inflow,zmin,zmins
real, dimension(imax), intent(inout) :: sst,zo,cd,cds,umod,fg,eg,evspsbl,sbl,wetfac,fracice,siced,epot,epan,snowd
real, dimension(imax) :: atm_rnd,atm_snd
real, dimension(imax) :: workb,workc,dumazmin
type(dgicedata), intent(inout) :: dgice
type(dgscrndata), intent(inout) :: dgscrn
type(dgwaterdata), intent(inout) :: dgwater
type(icedata), intent(inout) :: ice
type(waterdata), intent(inout) :: water
type(depthdata), intent(in) :: depth
type(turbdata), intent(inout) :: turb

if ( .not.mlo_active .or. .not.depth%data_allocated ) return

atm_rnd = precp - precs
atm_snd = precs

call mloeval_work(dt,zmin,zmins,sg,rg,atm_rnd,atm_snd,uatm,vatm,          &
                   temp,qg,ps,                                            &
                   visnirratio,fbvis,fbnir,inflow,diag,                   &
                   calcprog,depth,dgice,dgscrn,dgwater,ice,water,turb)

workb  =emisice**0.25*ice%tsurf
dumazmin=max(zmin,dgwater%zo+0.2,zoseaice+0.2)
workc  =(1.-ice%fracice)/log(dumazmin/dgwater%zo)**2+ice%fracice/log(dumazmin/zoseaice)**2
where ( depth%dz(:,1)>=1.e-4 )
  sst    =(1.-ice%fracice)*(water%temp(:,1)+wrtemp)+ice%fracice*workb
  zo     =dumazmin*exp(-1./sqrt(workc))
  cd     =(1.-ice%fracice)*dgwater%cd  +ice%fracice*dgice%cd  ! includes umod
  cds    =(1.-ice%fracice)*dgwater%cdh +ice%fracice*dgice%cdh ! includes umod
  umod   =(1.-ice%fracice)*dgwater%umod+ice%fracice*dgice%umod
  fg     =(1.-ice%fracice)*dgwater%fg  +ice%fracice*dgice%fg
  eg     =(1.-ice%fracice)*dgwater%eg  +ice%fracice*dgice%eg
  wetfac =(1.-ice%fracice)             +ice%fracice*dgice%wetfrac
  epan   =dgwater%eg
  epot   =(1.-ice%fracice)*dgwater%eg  +ice%fracice*dgice%eg/max(dgice%wetfrac,1.e-20)
  fracice=ice%fracice
  siced  =ice%thick
  snowd  =ice%snowd
  evspsbl=(1.-ice%fracice)*dgwater%eg/lv+ice%fracice*dgice%eg/ls
  sbl    =ice%fracice*dgice%eg/ls
end where  

return
end subroutine mloeval_thread

subroutine mlo_updatekm_imax(km_o,ks_o,gammas_o,dt,diag, &
                             depth,ice,dgwater,water,turb)                   

implicit none

integer, intent(in) :: diag
integer ii, iqw
real, intent(in) :: dt
real, dimension(imax,wlev), intent(out) :: km_o, ks_o, gammas_o
real, dimension(imax,wlev) :: gammas, km_hl, ks_hl
real, dimension(imax) :: d_zcr, t0
type(icedata), intent(in) :: ice
type(dgwaterdata), intent(inout) :: dgwater
type(waterdata), intent(inout) :: water
type(depthdata), intent(in) :: depth
type(turbdata), intent(inout) :: turb

if (diag>=2) write(6,*) "THREAD: Update diffuion coeff"
km_o = 0.
ks_o = 0.
gammas_o = 0.
if ( .not.mlo_active .or. .not.depth%data_allocated ) return

km_hl = 0.
ks_hl = 0.
gammas = 0.
where ( depth%dz(:,1)>1.e-4 )
  d_zcr = max(1.+max(water%eta,-delwater)/depth%depth_hl(:,wlev+1),minwater/depth%depth_hl(:,wlev+1))
elsewhere
  d_zcr = 1.
end where

call mlo_calc_k(km_hl,ks_hl,gammas,dt,d_zcr,depth,dgwater,water,turb)

do ii = 1,wlev
  where ( depth%dz(:,ii)>=1.e-4 )
    km_o(:,ii) = km_hl(:,ii) ! half levels
    ks_o(:,ii) = ks_hl(:,ii)
  end where  
end do
if ( oclosure==0 ) then
  if ( incradgam>0 ) then
    ! include radiation in counter-gradient term
    do iqw = 1,imax
      t0(iqw) = dgwater%wt0(iqw) + sum(dgwater%rad(iqw,1:dgwater%mixind(iqw)))
    end do
  else
    t0 = dgwater%wt0
  end if    
  do ii = 1,wlev
    where ( depth%dz(:,ii)>=1.e-4 )
      gammas(:,ii) = gammas(:,ii)*t0
      gammas_o(:,ii) = gammas(:,ii)
    end where
  end do  
else
  gammas_o = 0.
end if

return
end subroutine mlo_updatekm_imax

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Form seaice before flux calculations

subroutine mlonewice

implicit none

integer tile, is, ie

do tile = 1,ntiles
  is = (tile-1)*imax + 1
  ie = tile*imax
  call mlonewice_thread(0,depth_g(tile),ice_g(tile),water_g(tile))  
end do

end subroutine mlonewice

subroutine mlonewice_thread(diag,depth,ice,water)

implicit none

integer, intent(in) :: diag
type(icedata), intent(inout) :: ice
type(waterdata), intent(inout) :: water
type(depthdata), intent(in) :: depth

if (diag>=1) write(6,*) "Form new ice"
if ( .not.mlo_active .or. .not.depth%data_allocated ) return

call mlonewice_work(depth,ice,water)

end subroutine mlonewice_thread

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Update exponential weighted moving average

subroutine mlo_ema(dt,mode)

implicit none

integer tile, is, ie
real, intent(in) :: dt
character(len=*), intent(in) :: mode

do tile = 1,ntiles
  is = (tile-1)*imax + 1
  ie = tile*imax
  call mlo_ema_thread(dt,water_g(tile),depth_g(tile),mode)  
end do

end subroutine mlo_ema

subroutine mlo_ema_thread(dt,water,depth,mode)

implicit none

real, intent(in) :: dt
type(waterdata), intent(inout) :: water
type(depthdata), intent(in) :: depth
character(len=*), intent(in) :: mode

if ( .not.mlo_active .or. .not.depth%data_allocated ) return

select case(mode)
  case("uvw")
    call mlo_ema_uvw(dt,water,depth)
  case("ts")
    call mlo_ema_ts(dt,water,depth)  
  case("reset")
    call mlo_ema_reset(dt,water,depth)
  case default
    write(6,*) "ERROR: Unknown option for mlo_ema ",trim(mode)
    stop
end select
  
return
end subroutine mlo_ema_thread

end module mlo_ctrl
