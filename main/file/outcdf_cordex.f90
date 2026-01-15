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

! Output for CORDEX data.  Usually hourly.
    
module outcdf_cordex_m

use outcdf_common_m

private
public freqfile_cordex

contains
                    
!--------------------------------------------------------------
! HIGH FREQUENCY OUTPUT FILES
      
subroutine freqfile_cordex

use aerointerface                     ! Aerosol interface
use aerosol_arrays                    ! Aerosol arrays
use arrays_m                          ! Atmosphere dyamics prognostic arrays
use cc_mpi                            ! CC MPI routines
use const_phys                        ! Physical constants
use dates_m                           ! Date data
use dpsdt_m                           ! Vertical velocity
use extraout_m                        ! Additional diagnostics
use filnames_m                        ! Filenames
use histave_m                         ! Time average arrays
use infile                            ! Input file routines
use kuocom_m                          ! JLM convection
use liqwpar_m                         ! Cloud water mixing ratios
use morepbl_m                         ! Additional boundary layer diagnostics
use newmpar_m                         ! Grid parameters
use nharrs_m                          ! Non-hydrostatic atmosphere arrays
use nsibd_m                           ! Land-surface arrays
use parm_m                            ! Model configuration
use parmdyn_m                         ! Dynamics parameters
use parmgeom_m                        ! Coordinate data
use parmhdff_m                        ! Horizontal diffusion parameters
use parmhor_m                         ! Horizontal advection parameters
use pbl_m                             ! Boundary layer arrays
use prec_m                            ! Precipitation
use raddiag_m                         ! Radiation diagnostic
use screen_m                          ! Screen level diagnostics
use sigs_m                            ! Atmosphere sigma levels
use soil_m                            ! Soil and surface data
use soilsnow_m                        ! Soil, snow and surface data
use soilv_m                           ! Soil parameters
use tracers_m                         ! Tracer data
use uclem_ctrl, only :              & ! Urban
    uclem_avetemp, uclem_misc
use vvel_m                            ! Additional vertical velocity
use work2_m                           ! Diagnostic arrays
      
implicit none

include 'version.h'                   ! Model version data

integer, parameter :: freqvars = 45  ! number of variables to average
integer, dimension(:), allocatable :: vnode_dat
integer, dimension(:), allocatable :: procnode, procoffset
integer, dimension(5) :: adim
integer, dimension(4) :: sdim
integer, dimension(1) :: gpdim
integer, dimension(5) :: outdim
integer, dimension(1) :: msdim
integer idms,js,je,tile
integer ixp,iyp,izp,tlencd
integer icy,icm,icd,ich,icmi,ics
integer i,j,k,n,iq,fiarch
integer idnp, idgpn, idgpo
integer press_level, height_level
integer d4, ssize, fsize
integer, save :: fncid = -1
integer, save :: idnt = 0
integer, save :: idkdate = 0
integer, save :: idktime = 0
integer, save :: idmtimer = 0
real(kind=8), dimension(:,:), allocatable, save :: freqstore
real, dimension(ifull) :: umag, pmsl, outdata
real, dimension(ifull) :: ua_level, va_level, ta_level, hus_level, zg_level
real, dimension(ifull) :: wa_level
real, dimension(:,:), allocatable :: xpnt2
real, dimension(:,:), allocatable :: ypnt2
real, dimension(:), allocatable :: xpnt
real, dimension(:), allocatable :: ypnt
real, dimension(1) :: zpnt
real, dimension(kl) :: phi_local
real, dimension(ms) :: shallow_zse, zsoil
real xx, shallow_sum, new_sum
real press_level_pa
real(kind=8) tpnt
real, parameter :: shallow_max = 0.1 ! shallow soil depth (10cm)
logical, save :: first = .true.
logical local, lday, l6hr
logical cordex_core, cordex_tier1, cordex_tier2, cordex_urbrcc
logical cordex_tier2b
character(len=1024) ffile
character(len=80) lname
character(len=40) vname
character(len=33) grdtim
character(len=20) timorg

call START_LOG(outfile_begin)

! lprocformat mode is where one 'node' captian will
! write the output for that 'node' of processes.  Procformat supports virtual nodes, although
! they cannot be split across physical nodes.

! if myid==0 or local=.true., then this process needs to write to a file

local = localhist .and. vnode_myid==0
lday  = mod(ktau,nperday)==0.or.ktau==ntau
l6hr  = mod(ktau,nper6hr)==0.or.ktau==ntau

if ( localhist ) then
  d4    = 5
  ssize = 4
else
  d4    = 4
  ssize = 3
end if
fsize = ssize - 1 ! size of fixed variables

fiarch = ktau/tbave
cordex_core = .true.
cordex_tier1 = .true.
cordex_tier2 = .true.
cordex_tier2b = .true.
cordex_urbrcc = .true.

select case ( surf_cordex )
  case(0)
    cordex_core = .false.
    cordex_tier1 = .false.
    cordex_tier2 = .false.
    cordex_tier2b = .false.
    cordex_urbrcc = .false.
  case(1)
    cordex_core = .true.
    cordex_tier1 = .true.
    cordex_tier2 = .true.
    cordex_tier2b = .false.
    cordex_urbrcc = .false.
  case(2)
    cordex_core = .true.
    cordex_tier1 = .true.
    cordex_tier2 = .false.
    cordex_tier2b = .false.
    cordex_urbrcc = .false.
  case(3)
    cordex_core = .true.
    cordex_tier1 = .false.
    cordex_tier2 = .false.
    cordex_tier2b = .false.
    cordex_urbrcc = .false.
  case(11)
    cordex_core = .true.
    cordex_tier1 = .true.
    cordex_tier2 = .true.
    cordex_tier2b = .false.
    cordex_urbrcc = .true.
  case(12)
    cordex_core = .true.
    cordex_tier1 = .true.
    cordex_tier2 = .false.
    cordex_tier2b = .false.
    cordex_urbrcc = .true.  
  case(21)
    cordex_core = .true.
    cordex_tier1 = .true.
    cordex_tier2 = .true.
    cordex_tier2b = .true.
    cordex_urbrcc = .true.
  case default
    write(6,*) "ERROR: Invalid option for surf_cordex ",surf_cordex
    call ccmpi_abort(-1)
end select

! allocate arrays and open new file
if ( first ) then
  if ( myid==0 ) then
    write(6,*) "Initialise CORDEX output"
  end if
  allocate(freqstore(ifull,freqvars))
  freqstore(:,1:36) = 0._8
  freqstore(:,37) = -9.e9
  freqstore(:,38) =  9.e9
  if ( local ) then
    write(ffile,"(a,'.',i6.6)") trim(surfile), vnode_vleaderid
  else
    ffile = surfile
  end if
  if ( myid==0 .or. local ) then
    call ccnf_create(ffile,fncid)
    ! Create dimensions
    if ( local ) then
      call ccnf_def_dim(fncid,'longitude',il,adim(1))
      call ccnf_def_dim(fncid,'latitude',jl,adim(2))
    else
      call ccnf_def_dim(fncid,'longitude',il_g,adim(1))
      call ccnf_def_dim(fncid,'latitude',jl_g,adim(2))
    endif
    call ccnf_def_dim(fncid,'lev',1,adim(3))
    call ccnf_def_dim(fncid,'zsoil',ms,msdim(1))
    if ( local ) then
      call ccnf_def_dim(fncid,'processor',vnode_nproc,adim(4)) 
      if ( myid==0 ) then
        call ccnf_def_dim(fncid,'gprocessor',nproc,gpdim(1)) 
      else
        gpdim(1)=0
      end if
    end if
    tlencd = ntau/tbave
    call ccnf_def_dim(fncid,'time',tlencd,adim(d4))  
    ! Define coords.
    if ( local ) then
      outdim(1) = adim(1)
      outdim(2) = adim(4)
      call ccnf_def_var(fncid,'longitude','float',2,outdim(1:2),ixp)        
    else
      call ccnf_def_var(fncid,'longitude','float',1,adim(1:1),ixp)
    end if
    call ccnf_put_att(fncid,ixp,'point_spacing','even')
    call ccnf_put_att(fncid,ixp,'units','degrees_east')
    if ( local ) then
      outdim(1) = adim(2)
      outdim(2) = adim(4)
      call ccnf_def_var(fncid,'latitude','float',2,outdim(1:2),iyp)
    else
      call ccnf_def_var(fncid,'latitude','float',1,adim(2:2),iyp)
    end if
    call ccnf_put_att(fncid,iyp,'point_spacing','even')
    call ccnf_put_att(fncid,iyp,'units','degrees_north')
    call ccnf_def_var(fncid,'lev','float',1,adim(3:3),izp)
    call ccnf_put_att(fncid,izp,'positive','down')
    call ccnf_put_att(fncid,izp,'point_spacing','uneven')
    call ccnf_put_att(fncid,izp,'units','sigma_level')
    call ccnf_def_var(fncid,'zsoil','float',1,msdim(1:1),idms)
    call ccnf_put_att(fncid,idms,'point_spacing','uneven')
    call ccnf_put_att(fncid,idms,'units','m')    
    if ( local ) then
      call ccnf_def_var(fncid,'processor','float',1,adim(4:4),idnp)  
      if ( myid==0 ) then
        call ccnf_def_var(fncid,'gprocnode','int',1,gpdim(1:1),idgpn)
        call ccnf_def_var(fncid,'gprocoffset','int',1,gpdim(1:1),idgpo)
      end if
    end if
    call ccnf_def_var(fncid,'time','double',1,adim(d4:d4),idnt)
    call ccnf_put_att(fncid,idnt,'point_spacing','even')
    icy=kdate/10000
    icm=max(1,min(12,(kdate-icy*10000)/100))
    icd=max(1,min(31,(kdate-icy*10000-icm*100)))
    if ( icy<100 ) then
      icy=icy+1900
    end if
    ich=ktime/100
    icmi=(ktime-ich*100)
    ics=0
    write(timorg,'(i2.2,"-",a3,"-",i4.4,3(":",i2.2))') icd,month(icm),icy,ich,icmi,ics
    call ccnf_put_att(fncid,idnt,'time_origin',timorg)
    write(grdtim,'("minutes since ",i4.4,"-",i2.2,"-",i2.2," ",2(i2.2,":"),i2.2)') icy,icm,icd,ich,icmi,ics
    call ccnf_put_att(fncid,idnt,'units',grdtim)
    if ( leap==0 ) then
      call ccnf_put_att(fncid,idnt,'calendar','noleap')
    else if ( leap==2 ) then
      call ccnf_put_att(fncid,idnt,'calendar','360_day')  
    end if
    call ccnf_def_var(fncid,'kdate','int',1,adim(d4:d4),idkdate)
    call ccnf_def_var(fncid,'ktime','int',1,adim(d4:d4),idktime)
    call ccnf_def_var(fncid,'mtimer','int',1,adim(d4:d4),idmtimer)
    call ccnf_put_attg(fncid,'version',trim(version))        !   Model version

    ! Define global grid
    call ccnf_put_attg(fncid,'dt',dt)
    call ccnf_put_attg(fncid,'il_g',il_g)
    call ccnf_put_attg(fncid,'jl_g',jl_g)
    call ccnf_put_attg(fncid,'rlat0',rlat0)
    call ccnf_put_attg(fncid,'rlong0',rlong0)
    call ccnf_put_attg(fncid,'schmidt',schmidt)
    call ccnf_put_attg(fncid,'ms',ms)
    call ccnf_put_attg(fncid,'ntrac',ntrac)
    
    ! grid decomposition data
    if ( local ) then
      call ccnf_put_attg(fncid,'nproc',nproc)
      call ccnf_put_attg(fncid,'procmode',vnode_nproc)
      call ccnf_put_attg(fncid,'decomp','face')
    end if 
    
    ! ensemble data
    if ( driving_model_id /= ' ' ) then
      call ccnf_put_attg(fncid,'driving_model_id',trim(driving_model_id))
    end if
    if ( driving_model_ensemble_number /= ' ' ) then
      call ccnf_put_attg(fncid,'driving_model_ensemble_number',trim(driving_model_ensemble_number))
    end if
    if ( driving_experiment_name /= ' ' ) then
      call ccnf_put_attg(fncid,'driving_experiment_name',trim(driving_experiment_name))
    end if 
    
    ! solar data
    call ccnf_put_attg(fncid,'bpyear',bpyear)

    ! define variables
    if ( local ) then
      sdim(1:2) = adim(1:2) 
      sdim(3:4) = adim(4:5)
    else
      sdim(1:2) = adim(1:2)
      sdim(3)   = adim(4)
    end if
    if ( cordex_core ) then
      lname = 'Surface geopotential'
      call attrib(fncid,sdim(1:fsize),fsize,'zht',lname,'m-2 s-2',-1000.,90.e3,any_m,fixed_m,amean_m,float_m)
    end if
    if ( cordex_tier2 ) then
      lname = 'Soil type'        
      call attrib(fncid,sdim(1:fsize),fsize,'soilt',lname,'none',-650.,650.,any_m,fixed_m,anotdef_m,short_m)
      lname = 'Capacity of Soil to Store Water'
      call attrib(fncid,sdim(1:fsize),fsize,'mrsofc',lname,'kg m-2',0.,6500.,any_m,fixed_m,amean_m,short_m)
      lname = 'Urban fraction'
      call attrib(fncid,sdim(1:fsize),fsize,'sigmu',lname,'none',0.,3.25,any_m,fixed_m,land_m,short_m)
    end if  
    lname='x-component 10m wind'
    call attrib(fncid,sdim,ssize,'uas',lname,'m s-1',-130.,130.,any_m,point_m,amean_m,short_m)
    lname='y-component 10m wind'     
    call attrib(fncid,sdim,ssize,'vas',lname,'m s-1',-130.,130.,any_m,point_m,amean_m,short_m)
    lname='Near-Surface Air Temperature'     
    call attrib(fncid,sdim,ssize,'tscrn',lname,'K',100.,425.,any_m,point_m,amean_m,short_m)
    lname='Near-Surface Relative Humidity'     
    call attrib(fncid,sdim,ssize,'rhscrn',lname,'%',0.,200.,any_m,point_m,amean_m,short_m)
    lname='Precipitation'
    call attrib(fncid,sdim,ssize,'rnd',lname,'mm day-1',0.,1300.,any_m,tmean_m,amean_m,float_m)
    lname='Convective Precipitation'
    call attrib(fncid,sdim,ssize,'rnc',lname,'mm day-1',0.,1300.,any_m,tmean_m,amean_m,float_m)
    lname='Snowfall Flux'
    call attrib(fncid,sdim,ssize,'sno',lname,'mm day-1',0.,1300.,any_m,tmean_m,amean_m,float_m)
    lname='Graupelfall'
    call attrib(fncid,sdim,ssize,'grpl',lname,'mm day-1',0.,1300.,any_m,tmean_m,amean_m,float_m)
    lname ='Sea Level Pressure'
    call attrib(fncid,sdim,ssize,'pmsl',lname,'hPa',800.,1200.,any_m,point_m,amean_m,short_m)
    lname ='Surface Downwelling Shortwave Radiation'
    call attrib(fncid,sdim,ssize,'sgdn_ave',lname,'W m-2',-500.,2.e3,any_m,tmean_m,amean_m,float_m)
    if ( cordex_tier1 ) then
      lname ='Surface Direct Downwelling Shortwave Radiation'
      call attrib(fncid,sdim,ssize,'sgdndir_ave',lname,'W m-2',-500.,2.e3,any_m,tmean_m,amean_m,float_m)
    end if
    lname = 'Scaled Log Surface pressure'
    call attrib(fncid,sdim,ssize,'psf',lname,'none',-1.4,0.5,any_m,point_m,amean_m,short_m)
    lname = 'Screen mixing ratio'
    call attrib(fncid,sdim,ssize,'qgscrn',lname,'kg kg-1',0.,0.06,any_m,point_m,amean_m,short_m)
    lname = 'Total Cloud Fraction'
    call attrib(fncid,sdim,ssize,'cld',lname,'frac',0.,1.,any_m,tmean_m,amean_m,short_m)
    lname = 'Direct normal irradiance'
    call attrib(fncid,sdim,ssize,'dni',lname,'W m-2',-500.,2.e3,any_m,tmean_m,amean_m,float_m)
    if ( cordex_tier1 ) then
      do j = 1,3 ! 50m, 100m, 150m  
        height_level = height_level_data(j)
        call cordex_name(lname,"x-component ",height_level,"m wind")
        call cordex_name(vname,"ua",height_level,"m")
        call attrib(fncid,sdim,ssize,trim(vname),lname,'m s-1',-130.,130.,any_m,point_m,amean_m,short_m)
        call cordex_name(lname,"y-component ",height_level,"m wind")
        call cordex_name(vname,"va",height_level,"m")
        call attrib(fncid,sdim,ssize,trim(vname),lname,'m s-1',-130.,130.,any_m,point_m,amean_m,short_m)
      end do
    end if
    if ( cordex_tier2 ) then
      do j = 4,height_levels ! 200m, 250m, 300m  
        height_level = height_level_data(j)
        call cordex_name(lname,"x-component ",height_level,"m wind")
        call cordex_name(vname,"ua",height_level,"m")
        call attrib(fncid,sdim,ssize,trim(vname),lname,'m s-1',-130.,130.,any_m,point_m,amean_m,short_m)
        call cordex_name(lname,"y-component ",height_level,"m wind")
        call cordex_name(vname,"va",height_level,"m")
        call attrib(fncid,sdim,ssize,trim(vname),lname,'m s-1',-130.,130.,any_m,point_m,amean_m,short_m)
      end do
    end if
    if ( cordex_tier1 ) then
      height_level = height_level_data(1) ! 50m
      call cordex_name(lname,"Air temperature at ",height_level,"m")
      call cordex_name(vname,"ta",height_level,"m")
      call attrib(fncid,sdim,ssize,trim(vname),lname,'K',100.,400.,any_m,point_m,amean_m,short_m)
      call cordex_name(lname,"Specific Humidity at ",height_level,"m")
      call cordex_name(vname,"hus",height_level,"m")
      call attrib(fncid,sdim,ssize,trim(vname),lname,'1',0.,0.06,any_m,point_m,amean_m,short_m)
    end if  
    if ( cordex_core ) then
      lname = 'Daily Maximum Near-Surface Air Temperature'  
      call attrib(fncid,sdim,ssize,'tmaxscr',lname,'K',100.,425.,daily_m,max_m,amean_m,short_m)
      lname = 'Daily Minimum Near-Surface Air Temperature'
      call attrib(fncid,sdim,ssize,'tminscr',lname,'K',100.,425.,daily_m,min_m,amean_m,short_m)
    end if
    if ( cordex_tier1 ) then
      lname = 'Daily Maximum Hourly Precipitation Rate'
      call attrib(fncid,sdim,ssize,'prhmax',lname,'kg m-2 s-1',0.,2600.,daily_m,max_m,amean_m,float_m)
      lname = 'x-component max 10m wind (daily)'
      call attrib(fncid,sdim,ssize,'u10max',lname,'m s-1',-99.,99.,daily_m,max_m,amean_m,short_m)
      lname = 'y-component max 10m wind (daily)'
      call attrib(fncid,sdim,ssize,'v10max',lname,'m s-1',-99.,99.,daily_m,max_m,amean_m,short_m)
    end if  
    if ( cordex_core ) then
      lname = 'Surface Downwelling Longwave Radiation'
      call attrib(fncid,sdim,ssize,'rgdn_ave',lname,'W m-2',-500.,1.e3,any_m,tmean_m,amean_m,float_m)
    end if
    if ( cordex_tier1 ) then
      lname = 'Surface Upward Latent Heat Flux'
      call attrib(fncid,sdim,ssize,'eg_ave',lname,'W m-2',-3000.,3000.,any_m,tmean_m,amean_m,float_m)
      lname = 'Surface Upward Sensible Heat Flux'
      call attrib(fncid,sdim,ssize,'fg_ave',lname,'W m-2',-3000.,3000.,any_m,tmean_m,amean_m,float_m)
      lname = 'Solar net at ground (+ve down)'
      call attrib(fncid,sdim,ssize,'sgn_ave',lname,'W m-2',-500.,2000.,any_m,tmean_m,amean_m,float_m)
      lname = 'LW net at ground (+ve up)'
      call attrib(fncid,sdim,ssize,'rgn_ave',lname,'W m-2',-500.,1000.,any_m,tmean_m,amean_m,float_m)
    end if  
    if ( cordex_tier2 ) then
      lname = 'Avg potential evaporation'
      call attrib(fncid,sdim,ssize,'epot_ave',lname,'W m-2',-1000.,10.e3,any_m,tmean_m,amean_m,float_m)
    end if
    if ( cordex_tier1 ) then
      lname = 'Soil Frozen Water Content'
      call attrib(fncid,sdim,ssize,'mrfso',lname,'kg m-2',0.,6500.,sixhr_m,point_m,land_m,float_m)
      lname = 'Frozen Water Content in Upper Portion of Soil Column'
      call attrib(fncid,sdim,ssize,'mrfsos',lname,'kg m-2',0.,6500.,any_m,point_m,land_m,float_m)
      lname = 'Evaporation'
      call attrib(fncid,sdim,ssize,'evspsbl',lname,'mm day-1',-1300.,1300.,any_m,tmean_m,amean_m,float_m)
      lname = 'Surface runoff'
      call attrib(fncid,sdim,ssize,'mrros',lname,'mm day-1',0.,1300.,sixhr_m,tmean_m,land_m,float_m)
      lname = 'Runoff' ! mrro after pcc2hist
      call attrib(fncid,sdim,ssize,'runoff',lname,'mm day-1',0.,1300.,sixhr_m,tmean_m,land_m,float_m)
      lname = 'Total Soil Moisture Content'
      call attrib(fncid,sdim,ssize,'mrso',lname,'kg m-2',0.,6500.,sixhr_m,point_m,land_m,float_m)
      lname = 'Moisture in Upper Portion of Soil Column'
      call attrib(fncid,sdim,ssize,'mrsos',lname,'kg m-2',0.,6500.,any_m,point_m,land_m,float_m)
      lname = 'Snow melt' 
      call attrib(fncid,sdim,ssize,'snm',lname,'mm day-1',0.,1300.,sixhr_m,tmean_m,land_m,float_m)
      lname = 'TOA Outgoing Longwave Radiation'
      call attrib(fncid,sdim,ssize,'rtu_ave',lname,'W m-2',0.,800.,any_m,tmean_m,amean_m,float_m)
      lname = 'TOA Incident Shortwave Radiation'
      call attrib(fncid,sdim,ssize,'sint_ave',lname,'W m-2',0.,1600.,any_m,tmean_m,amean_m,float_m)
      lname = 'TOA Outgoing Shortwave Radiation'
      call attrib(fncid,sdim,ssize,'sot_ave',lname,'W m-2',0.,1000.,any_m,tmean_m,amean_m,float_m)
    end if 
    if ( cordex_tier2 .and. cordex_fix==0 ) then
      lname = 'High Level Cloud Fraction'
      call attrib(fncid,sdim,ssize,'clh',lname,'frac',0.,1.,any_m,tmean_m,amean_m,short_m)
      lname = 'Mid Level Cloud Fraction'
      call attrib(fncid,sdim,ssize,'clm',lname,'frac',0.,1.,any_m,tmean_m,amean_m,short_m)
      lname = 'Low Level Cloud Fraction'
      call attrib(fncid,sdim,ssize,'cll',lname,'frac',0.,1.,any_m,tmean_m,amean_m,short_m)
    else if ( cordex_tier2 ) then
      lname = 'High Level Cloud Fraction'
      call attrib(fncid,sdim,ssize,'clh',lname,'frac',0.,1.,sixhr_m,tmean_m,amean_m,short_m)
      lname = 'Mid Level Cloud Fraction'
      call attrib(fncid,sdim,ssize,'clm',lname,'frac',0.,1.,sixhr_m,tmean_m,amean_m,short_m)
      lname = 'Low Level Cloud Fraction'
      call attrib(fncid,sdim,ssize,'cll',lname,'frac',0.,1.,sixhr_m,tmean_m,amean_m,short_m)
    end if 
    if ( cordex_tier1 ) then
      lname = 'x-component wind stress'
      call attrib(fncid,sdim,ssize,'taux',lname,'N m-2',-50.,50.,any_m,tmean_m,amean_m,short_m)
      lname = 'y-component wind stress'
      call attrib(fncid,sdim,ssize,'tauy',lname,'N m-2',-50.,50.,any_m,tmean_m,amean_m,short_m)
    end if
    if ( cordex_tier2b ) then
      lname = 'Clear sky SW out at TOA'
      call attrib(fncid,sdim,ssize,'soc_ave',lname,'W m-2',0.,900.,sixhr_m,tmean_m,amean_m,float_m)
      lname = 'Clear sky SW at ground (+ve down)'
      call attrib(fncid,sdim,ssize,'sgc_ave',lname,'W m-2',-500.,2000.,sixhr_m,tmean_m,amean_m,float_m)
      lname = 'Surface Downwelling Clear-Sky Shortwave Radiation'
      call attrib(fncid,sdim,ssize,'sgdc_ave',lname,'W m-2',-500.,2000.,sixhr_m,tmean_m,amean_m,float_m)
      lname = 'Clear sky LW at TOA'
      call attrib(fncid,sdim,ssize,'rtc_ave',lname,'W m-2',0.,800.,sixhr_m,tmean_m,amean_m,float_m)
      lname = 'Clear sky LW at ground'
      call attrib(fncid,sdim,ssize,'rgc_ave',lname,'W m-2',-500.,1000.,sixhr_m,tmean_m,amean_m,float_m)
      lname = 'Surface Downwelling Clear-Sky Longwave Radiation'
      call attrib(fncid,sdim,ssize,'rgdc_ave',lname,'W m-2',-500.,2000.,sixhr_m,tmean_m,amean_m,float_m)
    end if  
    if ( cordex_tier2 ) then        
      if ( rescrn>0 ) then
        lname = 'Near-Surface Wind Speed of Gust'
        call attrib(fncid,sdim,ssize,'wsgs',lname,'m s-1',0.,350.,any_m,point_m,amean_m,short_m)
        lname = 'Daily Maximum Near-Surface Wind Speed of Gust'
        call attrib(fncid,sdim,ssize,'wsgsmax',lname,'m s-1',0.,350.,daily_m,max_m,amean_m,short_m)
        lname = 'Convective Available Potential Energy'
        call attrib(fncid,sdim,ssize,'CAPE',lname,'J kg-1',0.,20000.,any_m,point_m,amean_m,short_m)
        lname = 'Convective Inhibition'
        call attrib(fncid,sdim,ssize,'CIN',lname,'J kg-1',-20000.,0.,any_m,point_m,amean_m,short_m)
        lname = 'Lifted Index'
        call attrib(fncid,sdim,ssize,'LI',lname,'K',-100.,100.,any_m,point_m,amean_m,short_m)
      end if
    end if  
    if ( cordex_tier1 ) then
      lname = 'Surface Temperature'
      call attrib(fncid,sdim,ssize,'tsu',lname,'K',100.,425.,any_m,point_m,amean_m,short_m)
      lname = 'Height of Boundary Layer'
      call attrib(fncid,sdim,ssize,'pblh',lname,'m',0.,13000.,any_m,point_m,amean_m,short_m)
      lname = 'Water Vapor Path'
      call attrib(fncid,sdim,ssize,'prw',lname,'kg m-2',0.,130.,any_m,point_m,amean_m,short_m)
      lname = 'Condensed Water Path'
      call attrib(fncid,sdim,ssize,'clwvi',lname,'kg m-2',0.,130.,any_m,point_m,amean_m,short_m)
      lname = 'Ice Water Path'
      call attrib(fncid,sdim,ssize,'clivi',lname,'kg m-2',0.,130.,any_m,point_m,amean_m,short_m)
      lname = 'Snow Depth' ! liquid water
      call attrib(fncid,sdim,ssize,'snd',lname,'mm',0.,6500.,sixhr_m,point_m,land_m,short_m)
    end if
    if ( cordex_tier1 ) then
      ! fracice / siconca is supposed to be daily.  But we use 6hourly to make a sensible output for AXIOM
      lname = 'Sea ice fraction'
      call attrib(fncid,sdim,ssize,'fracice',lname,'none',0.,1.,sixhr_m,point_m,sea_m,short_m)
      lname = 'Sunshine hours per day'
      call attrib(fncid,sdim,ssize,'sunhours',lname,'hrs',0.,24.,daily_m,sum_m,amean_m,short_m)
    end if
    if ( cordex_tier2 ) then
      lname = 'Surface roughness'
      call attrib(fncid,sdim,ssize,'zolnd',lname,'m',0.,65.,daily_m,point_m,amean_m,float_m)
      if ( abs(iaero)>=2 ) then
        lname = 'atmosphere_optical_thickness_due_to_ambient_aerosol_particles'
        call attrib(fncid,sdim,ssize,'od550aer',lname,'1',0.,13.,daily_m,point_m,amean_m,short_m)
      end if  
    end if
    if ( cordex_tier1 ) then
      do k = 1,10 ! 1000, 925, 850, 700, 600, 500, 400, 300, 250, 200
        press_level = cordex_level_data(k)
        call cordex_name(lname,"x-component ",press_level,"hPa wind")
        call cordex_name(vname,"ua",press_level)
        call attrib(fncid,sdim,ssize,trim(vname),lname,'m s-1',-130.,130.,sixhr_m,point_m,amean_m,short_m)
        call cordex_name(lname,"y-component ",press_level,"hPa wind")
        call cordex_name(vname,"va",press_level)
        call attrib(fncid,sdim,ssize,trim(vname),lname,'m s-1',-130.,130.,sixhr_m,point_m,amean_m,short_m)
        lname = 'Air Temperature'     
        call cordex_name(vname,"ta",press_level)
        call attrib(fncid,sdim,ssize,trim(vname),lname,'K',100.,425.,sixhr_m,point_m,amean_m,short_m)
        lname = 'Specific Humidity'
        call cordex_name(vname,"hus",press_level)
        call attrib(fncid,sdim,ssize,trim(vname),lname,'1',0.,0.06,sixhr_m,point_m,amean_m,short_m)
        lname = 'Geopotential Height'
        call cordex_name(vname,"zg",press_level)
        call attrib(fncid,sdim,ssize,trim(vname),lname,'m',0.,130000.,sixhr_m,point_m,amean_m,short_m)
        lname = 'Upward Air Velocity'
        call cordex_name(vname,"wa",press_level)
        call attrib(fncid,sdim,ssize,trim(vname),lname,'m s-1',-130.,130.,sixhr_m,point_m,amean_m,short_m)
      end do 
    end if
    ! avaliable in std output
    if ( cordex_tier2b ) then
      do k = 11,cordex_levels ! 150, 100, 75, 50, 30, 20, 10
        press_level = cordex_level_data(k)
        call cordex_name(lname,"x-component ",press_level,"hPa wind")
        call cordex_name(vname,"ua",press_level)
        call attrib(fncid,sdim,ssize,trim(vname),lname,'m s-1',-130.,130.,sixhr_m,point_m,amean_m,short_m)
        call cordex_name(lname,"y-component ",press_level,"hPa wind")
        call cordex_name(vname,"va",press_level)
        call attrib(fncid,sdim,ssize,trim(vname),lname,'m s-1',-130.,130.,sixhr_m,point_m,amean_m,short_m)
        lname = 'Air Temperature'     
        call cordex_name(vname,"ta",press_level)
        call attrib(fncid,sdim,ssize,trim(vname),lname,'K',100.,425.,sixhr_m,point_m,amean_m,short_m)
        lname = 'Specific Humidity'
        call cordex_name(vname,"hus",press_level)
        call attrib(fncid,sdim,ssize,trim(vname),lname,'1',0.,0.06,sixhr_m,point_m,amean_m,short_m)
        lname = 'Geopotential Height'
        call cordex_name(vname,"zg",press_level)
        call attrib(fncid,sdim,ssize,trim(vname),lname,'m',0.,130000.,sixhr_m,point_m,amean_m,short_m)
        lname = 'Upward Air Velocity'
        call cordex_name(vname,"wa",press_level)
        call attrib(fncid,sdim,ssize,trim(vname),lname,'m s-1',-130.,130.,sixhr_m,point_m,amean_m,short_m)
      end do 
    end if
    
    if ( cordex_tier1 ) then   
      do k = 1,ms
        call cordex_name(vname,"tgg",k)  
        call cordex_name(lname,"Soil temperature lev ",k)
        call attrib(fncid,sdim,ssize,trim(vname),lname,'K',100.,425.,sixhr_m,point_m,land_m,short_m)
        call cordex_name(vname,"mrsol",k)  
        call cordex_name(lname,"Total Water Content of Soil Layer ",k)
        call attrib(fncid,sdim,ssize,trim(vname),lname,'kg m-2',0.,6500.,sixhr_m,point_m,land_m,short_m)
        call cordex_name(vname,"mrfsol",k)  
        call cordex_name(lname,"Frozen Water Content of Soil Layer ",k)
        call attrib(fncid,sdim,ssize,trim(vname),lname,'kg m-2',0.,6500.,sixhr_m,point_m,land_m,short_m)
      end do    
    end if   
    
    if ( cordex_urbrcc ) then
      lname = 'Urban anthropogenic flux'
      call attrib(fncid,sdim,ssize,'anth_ave',lname,'W m-2',0.,650.,any_m,tmean_m,land_m,short_m)
      lname = 'Skin temperature'
      call attrib(fncid,sdim,ssize,'tsskin',lname,'K',100.,425.,any_m,point_m,land_m,short_m)
      lname = 'Surface temperature pavements'
      call attrib(fncid,sdim,ssize,'tspav',lname,'K',100.,425.,any_m,point_m,land_m,short_m)
      lname = 'Surface temperature roof'
      call attrib(fncid,sdim,ssize,'tsroof',lname,'K',100.,425.,any_m,point_m,land_m,short_m)
      lname = 'Surface temperature green spaces'
      call attrib(fncid,sdim,ssize,'tsgree',lname,'K',100.,425.,any_m,point_m,land_m,short_m)
    end if    

    if ( output_windmax/=0 ) then
      lname = 'x-component max 10m wind (sub-daily)'
      call attrib(fncid,sdim,ssize,'u10m_max',lname,'m s-1',-99.,99.,any_m,max_m,amean_m,short_m) ! sub-daily
      lname = 'y-component max 10m wind (sub-daily)'
      call attrib(fncid,sdim,ssize,'v10m_max',lname,'m s-1',-99.,99.,any_m,max_m,amean_m,short_m) ! sub-daily
    end if
    
    lname = 'Updraft helicity (2-5km)'
    call attrib(fncid,sdim,ssize,'uh',lname,'m2 s-2',-520.,520.,any_m,point_m,amean_m,short_m)
    lname = 'Maximum updraft helicity (2-5km)'
    call attrib(fncid,sdim,ssize,'uhmax',lname,'m2 s-2',-520.,520.,any_m,max_m,amean_m,short_m)
    lname = 'Minimum updraft helicity (2-5km)'
    call attrib(fncid,sdim,ssize,'uhmin',lname,'m2 s-2',-520.,520.,any_m,min_m,amean_m,short_m)

    lname = 'Hail average diameter'
    call attrib(fncid,sdim,ssize,'hailradave',lname,'m',0.,1.3e-2,any_m,tmean_m,amean_m,float_m)
    lname = 'Hail maximum diameter'
    call attrib(fncid,sdim,ssize,'hailradmax',lname,'m',0.,1.3e-2,any_m,max_m,amean_m,float_m)
    
    ! end definition mode
    call ccnf_enddef(fncid)
    if ( local ) then
      ! procformat
      allocate(xpnt(il),xpnt2(il,vnode_nproc))
      do i = 1,ipan
        xpnt(i) = real(i + ioff)
      end do
      call ccmpi_gatherx(xpnt2,xpnt,0,comm_vnode)
      call ccnf_put_vara(fncid,ixp,(/1,1/),(/il,vnode_nproc/),xpnt2)
      deallocate(xpnt,xpnt2)
      allocate(ypnt(jl),ypnt2(jl,vnode_nproc))
      do n = 1,npan
        do j = 1,jpan
          i = j + (n-1)*jpan  
          ypnt(i) = real(j + joff + (n-noff)*il_g)
        end do
      end do
      call ccmpi_gatherx(ypnt2,ypnt,0,comm_vnode)
      call ccnf_put_vara(fncid,iyp,(/1,1/),(/jl,vnode_nproc/),ypnt2)
      deallocate(ypnt,ypnt2)
    else
      allocate(xpnt(il_g))
      do i=1,il_g
        xpnt(i) = real(i)
      end do
      call ccnf_put_vara(fncid,ixp,1,il_g,xpnt(1:il_g))
      deallocate(xpnt)
      allocate(ypnt(jl_g))
      do j=1,jl_g
        ypnt(j) = real(j)
      end do
      call ccnf_put_vara(fncid,iyp,1,jl_g,ypnt(1:jl_g))
      deallocate(ypnt)
    end if
    zpnt(1)=1.
    call ccnf_put_vara(fncid,izp,1,1,zpnt(1:1))
    zsoil(1)=0.5*zse(1)
    zsoil(2)=zse(1)+zse(2)*0.5
    do k = 3,ms
      zsoil(k)=sum(zse(1:k-1))+zse(k)*0.5
    end do
    call ccnf_put_vara(fncid,idms,1,ms,zsoil(1:ms))    
    
    if ( local ) then
      ! store local processor order in output file  
      allocate( vnode_dat(vnode_nproc) )  
      call ccmpi_gatherx(vnode_dat,myid,0,comm_vnode)
      call ccnf_put_vara(fncid,idnp,(/1/),(/vnode_nproc/),vnode_dat)
      deallocate( vnode_dat )
      ! store file id for a given processor number in output file number 000000
      if ( myid==0 ) then
        allocate( procnode(nproc) )
      else
        allocate( procnode(1) ) ! not used
      end if  
      call ccmpi_gatherx(procnode,vnode_vleaderid,0,comm_world) ! this is procnode_inv
      if ( myid==0 ) then
        call ccnf_put_vara(fncid,idgpn,(/1/),(/nproc/),procnode)  
      end if
      deallocate(procnode)
      ! store offset within a file for a given processor number in output file number 000000
      if ( myid==0 ) then
        allocate( procoffset(nproc) )
      else
        allocate( procoffset(1) ) ! not used
      end if
      call ccmpi_gatherx(procoffset,vnode_myid,0,comm_world) ! this is procoffset_inv
      if ( myid==0 ) then
        call ccnf_put_vara(fncid,idgpo,(/1/),(/nproc/),procoffset)  
      end if
      deallocate(procoffset)
    end if
    
  else if ( localhist ) then
    
    allocate(xpnt(il),xpnt2(il,vnode_nproc))
    do i = 1,ipan
      xpnt(i) = real(i + ioff)
    end do
    call ccmpi_gatherx(xpnt2,xpnt,0,comm_vnode)
    deallocate(xpnt,xpnt2)
    allocate(ypnt(jl),ypnt2(jl,vnode_nproc))
    do n = 1,npan
      do j = 1,jpan
        i = j + (n-1)*jpan  
        ypnt(i) = real(j + joff + (n-noff)*il_g)
      end do
    end do
    call ccmpi_gatherx(ypnt2,ypnt,0,comm_vnode)
    deallocate(ypnt,ypnt2)
    
    allocate( vnode_dat(vnode_nproc) )
    call ccmpi_gatherx(vnode_dat,myid,0,comm_vnode)
    deallocate( vnode_dat )
    allocate(procnode(1)) ! not used
    call ccmpi_gatherx(procnode,vnode_vleaderid,0,comm_world) ! this is procnode_inv
    deallocate(procnode)
    allocate(procoffset(1)) ! not used
    call ccmpi_gatherx(procoffset,vnode_myid,0,comm_world) ! this is procoffset_inv
    deallocate(procoffset)
    
  end if ! myid==0 .or. local ..else.. localhist
  
  if ( cordex_core ) then
    call histwrt(zs,'zht',fncid,fiarch,local,.true.)
    outdata(:) = real(isoilm_in(:))  ! use the raw soil data here to classify inland water bodies
    call histwrt(outdata,'soilt',fncid,fiarch,local,.true.) ! also defines land-sea mask
  end if
  if ( cordex_tier2 ) then
    outdata(:) = sfc(isoilm)*sum(zse)*1000.
    call histwrt(outdata,'mrsofc',fncid,fiarch,local,.true.)
    call histwrt(sigmu,'sigmu',fncid,fiarch,local,.true.)
  end if  
  
  first=.false.
  if ( myid==0 ) write(6,*) "Finished initialising CORDEX output"
 
end if

! store output
freqstore(1:ifull,1) = freqstore(1:ifull,1) + real(condx*(86400./dt/real(tbave)),8)
freqstore(1:ifull,2) = freqstore(1:ifull,2) + real(condc*(86400./dt/real(tbave)),8)
freqstore(1:ifull,3) = freqstore(1:ifull,3) + real(conds*(86400./dt/real(tbave)),8)
freqstore(1:ifull,4) = freqstore(1:ifull,4) + real(condg*(86400./dt/real(tbave)),8)
freqstore(1:ifull,5) = freqstore(1:ifull,5) + real(sgdn/real(tbave),8)
freqstore(1:ifull,6) = freqstore(1:ifull,6) + real(sgdndir/real(tbave),8)
freqstore(1:ifull,7) = freqstore(1:ifull,7) + real(cloudtot/real(tbave),8)
freqstore(1:ifull,8) = freqstore(1:ifull,8) + real(dni/real(tbave),8)
freqstore(1:ifull,9) = freqstore(1:ifull,9) + real(rgdn/real(tbave),8)
freqstore(1:ifull,10) = freqstore(1:ifull,10) + real(eg/real(tbave),8)
freqstore(1:ifull,11) = freqstore(1:ifull,11) + real(fg/real(tbave),8)
freqstore(1:ifull,12) = freqstore(1:ifull,12) + real(sgsave/real(tbave),8)
freqstore(1:ifull,13) = freqstore(1:ifull,13) + real(rgn/real(tbave),8)
freqstore(1:ifull,14) = freqstore(1:ifull,14) + real(epot/real(tbave),8)
freqstore(1:ifull,15) = freqstore(1:ifull,15) + real(rt/real(tbave),8)
freqstore(1:ifull,16) = freqstore(1:ifull,16) + real(sint/real(tbave),8)
freqstore(1:ifull,17) = freqstore(1:ifull,17) + real(sout/real(tbave),8)
freqstore(1:ifull,18) = freqstore(1:ifull,18) + real(cloudhi/real(tbave),8)
freqstore(1:ifull,19) = freqstore(1:ifull,19) + real(cloudmi/real(tbave),8)
freqstore(1:ifull,20) = freqstore(1:ifull,20) + real(cloudlo/real(tbave),8)
freqstore(1:ifull,21) = freqstore(1:ifull,21) + real(taux/real(tbave),8)
freqstore(1:ifull,22) = freqstore(1:ifull,22) + real(tauy/real(tbave),8)

freqstore(1:ifull,23) = freqstore(1:ifull,23) + real(soutclr/real(nperday),8)
freqstore(1:ifull,24) = freqstore(1:ifull,24) + real(sgclr/real(nperday),8)
freqstore(1:ifull,25) = freqstore(1:ifull,25) + real(sgdclr/real(nperday),8)
freqstore(1:ifull,26) = freqstore(1:ifull,26) + real(rtclr/real(nperday),8)
freqstore(1:ifull,27) = freqstore(1:ifull,27) + real(rgclr/real(nperday),8)
freqstore(1:ifull,28) = freqstore(1:ifull,28) + real(rgdclr/real(nperday),8)
if ( abs(iaero)>=2 ) then
  freqstore(1:ifull,29) = freqstore(1:ifull,29) + real(opticaldepth(:,4,1)/real(nperday),8)
end if
umag = sqrt(u(1:ifull,1)**2+v(1:ifull,1)**2)
where ( u10**2 > real(freqstore(1:ifull,30))**2 + real(freqstore(1:ifull,31))**2 )
  freqstore(1:ifull,30) = real(u10(:)*u(1:ifull,1)/max(0.001,umag),8)
  freqstore(1:ifull,31) = real(u10(:)*v(1:ifull,1)/max(0.001,umag),8)
end where
freqstore(1:ifull,32) = freqstore(1:ifull,32) + real(anthropogenic_flux/real(tbave),8)
freqstore(1:ifull,33) = freqstore(1:ifull,33) + real(runoff*(86400./dt/real(tbave)),8)
freqstore(1:ifull,34) = freqstore(1:ifull,34) + real(runoff_surface*(86400./dt/real(tbave)),8)
freqstore(1:ifull,35) = freqstore(1:ifull,35) + real(snowmelt/real(tbave),8)
freqstore(1:ifull,36) = freqstore(1:ifull,36) + real(evspsbl*(86400./dt/real(tbave)),8)
freqstore(1:ifull,37) = max( freqstore(1:ifull,37), real(updraft_helicity,8) )
freqstore(1:ifull,38) = max( freqstore(1:ifull,38), real(updraft_helicity,8) )
freqstore(1:ifull,39) = max( freqstore(1:ifull,39), real(dhail1,8) )
freqstore(1:ifull,40) = max( freqstore(1:ifull,40), real(dhail2,8) )
freqstore(1:ifull,41) = max( freqstore(1:ifull,41), real(dhail3,8) )
freqstore(1:ifull,42) = max( freqstore(1:ifull,42), real(dhail4,8) )
freqstore(1:ifull,43) = max( freqstore(1:ifull,43), real(dhail5,8) )
freqstore(1:ifull,44) = freqstore(1:ifull,44) + 0.2*( freqstore(:,39) + freqstore(:,40) + &
    freqstore(:,41) + freqstore(:,42) + freqstore(:,43) )
freqstore(1:ifull,45) = max( freqstore(:,39), freqstore(:,40), freqstore(:,41), freqstore(:,42), &
    freqstore(:,43), freqstore(:,45) )

shallow_zse(:) = 0.
shallow_sum = 0.
do j = 1,ms
  new_sum = min( shallow_sum + zse(j), shallow_max )
  shallow_zse(j) = new_sum - shallow_sum
  shallow_sum = new_sum
end do

! write data to file
if ( mod(ktau,tbave)==0 ) then
    
  if ( myid==0 .or. local ) then
    if ( myid==0 ) then
      write(6,*) "write CORDEX output"
    end if
    tpnt = real(ktau,8)*(real(dt,8)/60._8)
    call ccnf_put_vara(fncid,'time',fiarch,tpnt)
    call ccnf_put_vara(fncid,'kdate',fiarch,kdate)
    call ccnf_put_vara(fncid,'ktime',fiarch,ktime)
    call ccnf_put_vara(fncid,'mtimer',fiarch,mtimer)
  end if
  
  ! record output
  umag = sqrt(u(1:ifull,1)**2+v(1:ifull,1)**2)
  call mslp(pmsl,psl,zs,t)
  outdata = u10*u(1:ifull,1)/max(umag,1.E-6)
  call histwrt(outdata,"uas",fncid,fiarch,local,.true.)
  outdata = u10*v(1:ifull,1)/max(umag,1.E-6)
  call histwrt(outdata,"vas",fncid,fiarch,local,.true.)
  call histwrt(tscrn,"tscrn",fncid,fiarch,local,.true.)
  call histwrt(rhscrn,"rhscrn",fncid,fiarch,local,.true.)
  outdata = real(freqstore(:,1))
  call histwrt(outdata,"rnd",fncid,fiarch,local,.true.)
  outdata = real(freqstore(:,2))
  call histwrt(outdata,"rnc",fncid,fiarch,local,.true.)
  outdata = real(freqstore(:,3))
  call histwrt(outdata,"sno",fncid,fiarch,local,.true.)
  outdata = real(freqstore(:,4))
  call histwrt(outdata,"grpl",fncid,fiarch,local,.true.)
  outdata = pmsl/100.
  call histwrt(outdata,"pmsl",fncid,fiarch,local,.true.)
  outdata = real(freqstore(:,5))
  call histwrt(outdata,"sgdn_ave",fncid,fiarch,local,.true.)
  if ( cordex_tier1 ) then
    outdata = real(freqstore(:,6))  
    call histwrt(outdata,"sgdndir_ave",fncid,fiarch,local,.true.)
  end if  
  call histwrt(psl,"psf",fncid,fiarch,local,.true.)
  call histwrt(qgscrn,"qgscrn",fncid,fiarch,local,.true.)
  outdata = real(freqstore(:,7))
  call histwrt(outdata,"cld",fncid,fiarch,local,.true.)
  outdata = real(freqstore(:,8))
  call histwrt(outdata,"dni",fncid,fiarch,local,.true.)
  if ( cordex_tier1 ) then
    do j = 1,3  ! 50m, 100m, 150m
      height_level = height_level_data(j)  
      do iq = 1,ifull
        phi_local(1) = bet(1)*t(iq,1)
        do k = 2,kl
          phi_local(k) = phi_local(k-1) + bet(k)*t(iq,k) + betm(k)*t(iq,k-1)
        end do
        do k = 1,kl-1
          if ( phi_local(k)/grav<real(height_level) ) then
            n = k
          else
            exit
          end if
        end do
        xx = (real(height_level)*grav-phi_local(n))/(phi_local(n+1)-phi_local(n))
        ua_level(iq) = u(iq,n)*(1.-xx) + u(iq,n+1)*xx
        va_level(iq) = v(iq,n)*(1.-xx) + v(iq,n+1)*xx
      end do
      call cordex_name(vname,"ua",height_level,"m")
      call histwrt(ua_level,vname,fncid,fiarch,local,.true.)
      call cordex_name(vname,"va",height_level,"m")
      call histwrt(va_level,vname,fncid,fiarch,local,.true.)  
    end do 
  end if
  if ( cordex_tier2 ) then
    do j = 4,height_levels ! 200m, 250m, 300m 
      height_level = height_level_data(j)  
      do iq = 1,ifull
        phi_local(1) = bet(1)*t(iq,1)
        do k = 2,kl
          phi_local(k) = phi_local(k-1) + bet(k)*t(iq,k) + betm(k)*t(iq,k-1)
        end do
        do k = 1,kl-1
          if ( phi_local(k)/grav<real(height_level) ) then
            n = k
          else
            exit
          end if
        end do
        xx = (real(height_level)*grav-phi_local(n))/(phi_local(n+1)-phi_local(n))
        ua_level(iq) = u(iq,n)*(1.-xx) + u(iq,n+1)*xx
        va_level(iq) = v(iq,n)*(1.-xx) + v(iq,n+1)*xx
      end do
      call cordex_name(vname,"ua",height_level,"m")
      call histwrt(ua_level,vname,fncid,fiarch,local,.true.)
      call cordex_name(vname,"va",height_level,"m")
      call histwrt(va_level,vname,fncid,fiarch,local,.true.)  
    end do  
  end if
  if ( cordex_tier1 ) then
    height_level = height_level_data(1) ! 50m 
    do iq = 1,ifull
      phi_local(1) = bet(1)*t(iq,1)
      do k = 2,kl
        phi_local(k) = phi_local(k-1) + bet(k)*t(iq,k) + betm(k)*t(iq,k-1)
      end do
      do k = 1,kl-1
        if ( phi_local(k)/grav<real(height_level) ) then
          n = k
        else
          exit
        end if
      end do
      xx = (real(height_level)*grav-phi_local(n))/(phi_local(n+1)-phi_local(n))
      ta_level(iq) = t(iq,n)*(1.-xx) + t(iq,n+1)*xx
      hus_level(iq) = qg(iq,n)*(1.-xx) + qg(iq,n+1)*xx
      end do
    call cordex_name(vname,"ta",height_level,"m")
    call histwrt(ta_level,vname,fncid,fiarch,local,.true.)
    call cordex_name(vname,"hus",height_level,"m")
    call histwrt(hus_level,vname,fncid,fiarch,local,.true.)
  end if  
  if ( cordex_core ) then
    call histwrt(tmaxscr,'tmaxscr',fncid,fiarch,local,lday)
    call histwrt(tminscr,'tminscr',fncid,fiarch,local,lday)
  end if
  if ( cordex_tier1 ) then
    call histwrt(prhmax,'prhmax',fncid,fiarch,local,lday)    
    call histwrt(u10max,'u10max',fncid,fiarch,local,lday)
    call histwrt(v10max,'v10max',fncid,fiarch,local,lday)
  end if  
  if ( cordex_core ) then
    outdata = real(freqstore(:,9))  
    call histwrt(outdata,"rgdn_ave",fncid,fiarch,local,.true.)
  end if
  if ( cordex_tier1 ) then
    outdata = real(freqstore(:,10))  
    call histwrt(outdata,"eg_ave",fncid,fiarch,local,.true.)
    outdata = real(freqstore(:,11))
    call histwrt(outdata,"fg_ave",fncid,fiarch,local,.true.)
    outdata = real(freqstore(:,12))
    call histwrt(outdata,"sgn_ave",fncid,fiarch,local,.true.)
    outdata = real(freqstore(:,13))
    call histwrt(outdata,"rgn_ave",fncid,fiarch,local,.true.)
  end if
  if ( cordex_tier2 ) then
    outdata = real(freqstore(:,14))  
    call histwrt(outdata,"epot_ave",fncid,fiarch,local,.true.)
  end if  
  if ( cordex_tier1 ) then
    outdata = 0.
    do k = 1,ms
      outdata = outdata + wbice(:,k)*zse(k)*330.  
    end do
    call histwrt(outdata,"mrfso",fncid,fiarch,local,l6hr)
    !--
    outdata = 0.
    do k = 1,ms
      outdata = outdata + wbice(:,k)*shallow_zse(k)*330.  
    end do    
    call histwrt(outdata,"mrfsos",fncid,fiarch,local,.true.)
    !--
    outdata = real(freqstore(:,36))
    call histwrt(outdata,"evspsbl",fncid,fiarch,local,.true.)    
    !--
    outdata = real(freqstore(:,34))
    call histwrt(outdata,"mrros",fncid,fiarch,local,l6hr)
    !--
    outdata = real(freqstore(:,33))
    call histwrt(outdata,"runoff",fncid,fiarch,local,l6hr)
    !--
    outdata = 0.
    do k = 1,ms
      outdata = outdata + wb(:,k)*zse(k)*1000.  
    end do    
    call histwrt(outdata,"mrso",fncid,fiarch,local,l6hr)
    !--
    outdata = 0.
    do k = 1,ms
      outdata = outdata + wb(:,k)*shallow_zse(k)*1000.  
    end do    
    call histwrt(outdata,"mrsos",fncid,fiarch,local,.true.)   
    !--
    outdata = real(freqstore(:,35))
    call histwrt(outdata,"snm",fncid,fiarch,local,l6hr)
    !--
    outdata = real(freqstore(:,15))
    call histwrt(outdata,"rtu_ave",fncid,fiarch,local,.true.)
    outdata = real(freqstore(:,16))
    call histwrt(outdata,"sint_ave",fncid,fiarch,local,.true.)
    outdata = real(freqstore(:,17))
    call histwrt(outdata,"sot_ave",fncid,fiarch,local,.true.)
  end if
  if ( cordex_tier2 .and. cordex_fix==0 ) then
    outdata = real(freqstore(:,18))  
    call histwrt(outdata,"clh",fncid,fiarch,local,.true.)
    outdata = real(freqstore(:,19))
    call histwrt(outdata,"clm",fncid,fiarch,local,.true.)
    outdata = real(freqstore(:,20))
    call histwrt(outdata,"cll",fncid,fiarch,local,.true.)
  else if ( cordex_tier2 ) then
    outdata = real(freqstore(:,18))  
    call histwrt(outdata,"clh",fncid,fiarch,l6hr,.true.)
    outdata = real(freqstore(:,19))
    call histwrt(outdata,"clm",fncid,fiarch,l6hr,.true.)
    outdata = real(freqstore(:,20))
    call histwrt(outdata,"cll",fncid,fiarch,l6hr,.true.)
  end if
  if ( cordex_tier1 ) then
    outdata = real(freqstore(:,21))  
    call histwrt(outdata,"taux",fncid,fiarch,local,.true.)
    outdata = real(freqstore(:,22))
    call histwrt(outdata,"tauy",fncid,fiarch,local,.true.)
  end if
  if ( cordex_tier2b ) then
    outdata = real(freqstore(:,23))  
    call histwrt(outdata,"soc_ave",fncid,fiarch,local,l6hr)
    outdata = real(freqstore(:,24))
    call histwrt(outdata,"sgc_ave",fncid,fiarch,local,l6hr)
    outdata = real(freqstore(:,25))
    call histwrt(outdata,"sgdc_ave",fncid,fiarch,local,l6hr)    
    outdata = real(freqstore(:,26))
    call histwrt(outdata,"rtc_ave",fncid,fiarch,local,l6hr)
    outdata = real(freqstore(:,27))
    call histwrt(outdata,"rgc_ave",fncid,fiarch,local,l6hr)
    outdata = real(freqstore(:,28))
    call histwrt(outdata,"rgdc_ave",fncid,fiarch,local,l6hr)
  end if
  if ( cordex_tier2 ) then
    if ( rescrn>0 ) then
      call histwrt(wsgs,'wsgs',fncid,fiarch,local,.true.)  
      call histwrt(wsgsmax,'wsgsmax',fncid,fiarch,local,lday)
      call histwrt(cape_d,'CAPE',fncid,fiarch,local,.true.)
      call histwrt(cin_d,'CIN',fncid,fiarch,local,.true.)
      call histwrt(li_d,'LI',fncid,fiarch,local,.true.)
    end if
  end if
  if ( cordex_tier1 ) then
    call histwrt(tss,"tsu",fncid,fiarch,local,.true.)
    call histwrt(pblh,"pblh",fncid,fiarch,local,.true.)
    outdata = 0.
    do k = 1,kl
      outdata = outdata - dsig(k)*qg(1:ifull,k) ! sign of outdata defined so always positive 
    end do    
    outdata = outdata*ps(1:ifull)/grav
    call histwrt(outdata,"prw",fncid,fiarch,local,.true.)
    outdata = 0.
    do k = 1,kl
      outdata = outdata - dsig(k)*qlg(1:ifull,k) ! sign of outdata defined so always positive  
    end do    
    outdata = outdata*ps(1:ifull)/grav
    call histwrt(outdata,"clwvi",fncid,fiarch,local,.true.)
    outdata = 0.
    do k = 1,kl
      outdata = outdata - dsig(k)*qfg(1:ifull,k) ! sign of outdata defined so always positive 
    end do    
    outdata = outdata*ps(1:ifull)/grav
    call histwrt(outdata,"clivi",fncid,fiarch,local,.true.)
    call histwrt(snowd,"snd",fncid,fiarch,local,l6hr)
  end if
  if ( cordex_tier1 ) then
    call histwrt(fracice,"fracice",fncid,fiarch,local,lday)
    call histwrt(sunhours,'sunhours',fncid,fiarch,local,lday) 
  end if
  if ( cordex_tier2 ) then
    call histwrt(zo,'zolnd',fncid,fiarch,local,lday)  
    if ( abs(iaero)>=2 ) then
      outdata = real(freqstore(:,29))
      call histwrt(outdata,'od550aer',fncid,fiarch,local,lday)
    end if  
  end if
  if ( cordex_tier1 ) then
    do k = 1,10 ! 1000, 925, 850, 700, 600, 500, 400, 300, 250, 200
      press_level = cordex_level_data(k)
      press_level_pa = real(press_level)*100.  
      do iq = 1,ifull
        n = bisect(press_level_pa,ps(iq),sig(:)) 
        xx = (press_level_pa - ps(iq)*sig(n)) &
            /(ps(iq)*sig(n+1)-ps(iq)*sig(n))
        xx = min( max( xx, 0. ), 1. )
        ! special treatment for t
        if ( press_level_pa>ps(iq)*sig(1) ) then
          ta_level(iq) = t(iq,1)*(press_level_pa/(ps(iq)*sig(1)))**(6.5e-3*rdry/grav)
        else
          ta_level(iq) = t(iq,n)*(1.-xx) + t(iq,n+1)*xx
        end if
        hus_level(iq) = qg(iq,n)*(1.-xx) + qg(iq,n+1)*xx        
        hus_level(iq) = hus_level(iq)/(hus_level(iq)+1.) ! convert to specific humidity
        ua_level(iq) = u(iq,n)*(1.-xx) + u(iq,n+1)*xx
        va_level(iq) = v(iq,n)*(1.-xx) + v(iq,n+1)*xx
        zg_level(iq) = (phi(iq,n)*(1.-xx) + phi(iq,n+1)*xx)/grav
        wa_level(iq) = wvel(iq,n)*(1.-xx) + wvel(iq,n+1)*xx
      end do
      call cordex_name(vname,"ua",press_level)
      call histwrt(ua_level,trim(vname),fncid,fiarch,local,l6hr) 
      call cordex_name(vname,"va",press_level)
      call histwrt(va_level,trim(vname),fncid,fiarch,local,l6hr) 
      call cordex_name(vname,"ta",press_level)
      call histwrt(ta_level,trim(vname),fncid,fiarch,local,l6hr) 
      call cordex_name(vname,"hus",press_level)
      call histwrt(hus_level,trim(vname),fncid,fiarch,local,l6hr) 
      call cordex_name(vname,"zg",press_level)
      call histwrt(zg_level,trim(vname),fncid,fiarch,local,l6hr) 
      call cordex_name(vname,"wa",press_level)
      call histwrt(wa_level,trim(vname),fncid,fiarch,local,l6hr) 
    end do  
  end if
  ! avaliable in std output
  !if ( cordex_tier1 ) then
  !  do j = 11,cordex_levels ! 150, 100, 75, 50, 30, 20, 10
  !    press_level = cordex_level_data(j)
  !    press_level_pa = real(press_level)*100.
  !    do iq = 1,ifull
  !      n = bisect(press_level_pa,ps(iq),sig(:)) 
  !      xx = (press_level_pa - ps(iq)*sig(n)) &
  !          /(ps(iq)*sig(n+1)-ps(iq)*sig(n))
  !      xx = min( max( xx, 0. ), 1. )
  !      ! special treatment for t
  !      if ( press_level_pa>ps(iq)*sig(1) ) then
  !        ta_level(iq) = t(iq,1)*(press_level_pa/(ps(iq)*sig(1)))**(6.5e-3*rdry/grav)
  !      else
  !        ta_level(iq) = t(iq,n)*(1.-xx) + t(iq,n+1)*xx
  !      end if
  !      hus_level(iq) = qg(iq,n)*(1.-xx) + qg(iq,n+1)*xx
  !      hus_level(iq) = hus_level(iq)/(hus_level(iq)+1.)
  !      ua_level(iq) = u(iq,n)*(1.-xx) + u(iq,n+1)*xx
  !      va_level(iq) = v(iq,n)*(1.-xx) + v(iq,n+1)*xx
  !      zg_level(iq) = phi(iq,n)*(1.-xx) + phi(iq,n+1)*xx
  !      zg_level(iq) = zg_level(iq)/grav
  !      wa_level(iq) = wvel(iq,n)*(1.-xx) + wvel(iq,n+1)*xx
  !    end do
  !    call cordex_name(vname,"ua",press_level)
  !    call histwrt(ua_level,trim(vname),fncid,fiarch,local,l6hr) 
  !    call cordex_name(vname,"va",press_level)
  !    call histwrt(va_level,trim(vname),fncid,fiarch,local,l6hr) 
  !    call cordex_name(vname,"ta",press_level)
  !    call histwrt(ta_level,trim(vname),fncid,fiarch,local,l6hr) 
  !    call cordex_name(vname,"hus",press_level)
  !    call histwrt(hus_level,trim(vname),fncid,fiarch,local,l6hr) 
  !    call cordex_name(vname,"zg",press_level)
  !    call histwrt(zg_level,trim(vname),fncid,fiarch,local,l6hr) 
  !    call cordex_name(vname,"wa",press_level)
  !    call histwrt(wa_level,trim(vname),fncid,fiarch,local,l6hr) 
  !  end do  
  !end if  
    
  if ( cordex_tier1 ) then
    do k = 1,ms
      call cordex_name(vname,"tgg",k)  
      call histwrt(tgg(:,k),vname,fncid,fiarch,local,l6hr)
      outdata = wb(:,k)*zse(k)*1000.  
      call cordex_name(vname,"mrsol",k)  
      call histwrt(outdata,vname,fncid,fiarch,local,l6hr)
      outdata = wbice(:,k)*zse(k)*900.  
      call cordex_name(vname,"mrfsol",k)  
      call histwrt(outdata,vname,fncid,fiarch,local,l6hr)
    end do
  end if
  
  if ( cordex_urbrcc ) then
    outdata = real(freqstore(1:ifull,32))
    call histwrt(outdata,'anth_ave',fncid,fiarch,local,.true.) 
    call histwrt(urban_ts,'tsskin',fncid,fiarch,local,.true.)
    outdata = 999.
    call uclem_avetemp(outdata,"roadtemp1",0)  
    call histwrt(outdata,'tspav',fncid,fiarch,local,.true.)
    outdata = 999.
    call uclem_avetemp(outdata,"rooftemp1",0)  
    call histwrt(outdata,'tsroof',fncid,fiarch,local,.true.)
    outdata = 999.
    call uclem_misc(outdata,"vegt",0)  
    call histwrt(outdata,'tsgree',fncid,fiarch,local,.true.)    
  end if    

  if ( output_windmax/=0 ) then
    outdata = real(freqstore(1:ifull,30))  
    call histwrt(outdata,'u10m_max',fncid,fiarch,local,.true.)
    outdata = real(freqstore(1:ifull,31))
    call histwrt(outdata,'v10m_max',fncid,fiarch,local,.true.)
  end if
  
  call histwrt(updraft_helicity,'uh',fncid,fiarch,local,.true.)
  call histwrt(freqstore(:,37),'uhmax',fncid,fiarch,local,.true.)
  call histwrt(freqstore(:,38),'uhmin',fncid,fiarch,local,.true.)

  outdata = real(freqstore(:,44))
  call histwrt(outdata,'hailradave',fncid,fiarch,local,.true.)  
  outdata = real(freqstore(:,45))
  call histwrt(outdata,'hailradmax',fncid,fiarch,local,.true.)  
  
  freqstore(:,1:17) = 0._8
  if ( cordex_fix==0 ) then
    freqstore(:,18:20) = 0._8
  else if ( l6hr ) then
    freqstore(:,18:20) = 0._8
  end if
  freqstore(:,21:22) = 0._8
  if ( lday ) freqstore(:,23:29) = 0._8
  freqstore(:,30:32) = 0._8
  if ( l6hr ) freqstore(:,33:35) = 0._8
  freqstore(:,36) = 0._8
  freqstore(:,37) = -9.e9_8
  freqstore(:,38) = 9.e9_8
  freqstore(:,39:45) = 0._8
  
end if

if ( myid==0 .or. local ) then
  ! close file at end of run
  if ( ktau==ntau ) then
    call ccnf_close(fncid)
  end if
end if
      
call END_LOG(outfile_end)
      
return
end subroutine freqfile_cordex

end module outcdf_cordex_m
