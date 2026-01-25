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

! Output for high-frequency data.  Usually 10min
    
module outcdf_freq
    
use outcdf_common

implicit none

private
public freqfile_10

contains

subroutine freqfile_10

use arrays_m                          ! Atmosphere dyamics prognostic arrays
use cc_mpi                            ! CC MPI routines
use const_phys                        ! Physical constants
use dates_m                           ! Date data
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
use soilsnow_m                        ! Soil, snow and surface data
use soilv_m                           ! Soil parameters
use tracers_m                         ! Tracer data
use vvel_m                            ! Additional vertical velocity
      
implicit none

integer, parameter :: freqvars = 13  ! number of variables to average
integer, dimension(:), allocatable :: vnode_dat
integer, dimension(:), allocatable :: procnode, procoffset
integer, dimension(5) :: adim
integer, dimension(4) :: sdim
integer, dimension(1) :: gpdim
integer, dimension(5) :: outdim
integer ixp,iyp,izp
integer icy,icm,icd,ich,icmi,ics
integer i,j,n,fiarch,k,iq
integer idnp, idgpn, idgpo
integer press_level, tlenhf
integer d4, ssize, fsize
integer, save :: fncid = -1
integer, save :: idnt = 0
integer, save :: idkdate = 0
integer, save :: idktime = 0
integer, save :: idmtimer = 0
real(kind=8), dimension(:,:), allocatable, save :: freqstore
real, dimension(ifull) :: umag, outdata, pmsl
real, dimension(ifull) :: ua_level, va_level, ta_level, hus_level, zg_level
real, dimension(ifull) :: wa_level
real, dimension(:,:), allocatable :: xpnt2
real, dimension(:,:), allocatable :: ypnt2
real, dimension(:), allocatable :: xpnt
real, dimension(:), allocatable :: ypnt
real, dimension(1) :: zpnt
real(kind=8) tpnt
real press_level_pa, xx, sig_level
logical, save :: first = .true.
logical local
logical freq_core, freq_standard, freq_shep
character(len=1024) ffile
character(len=80) lname
character(len=40) vname
character(len=33) grdtim
character(len=20) timorg

call START_LOG(outfile_begin)

! procformat mode is where one 'node' captian will write the output for that
! 'node' of processes.  Procformat supports virtual nodes, although they
! cannot be split across physical nodes.

! if myid==0 or local=.true., then this process needs to write to a file

local = localhist .and. vnode_myid==0

if ( localhist ) then
  d4    = 5
  ssize = 4
else
  d4    = 4
  ssize = 3
end if
fsize = ssize - 1 ! size of fixed variables

freq_core = .false.
freq_standard = .false.
freq_shep = .false.

select case ( shep_cordex )
  case(0)
    freq_core = .true.
    freq_standard = .true.
    freq_shep = .false.
  case(1)
    freq_core = .true.
    freq_standard = .false.
    freq_shep = .true.
  case(2)
    freq_core = .true.
    freq_standard = .true.
    freq_shep = .true.
  case default
    write(6,*) "ERROR: Invalid option for shep_cordex ",shep_cordex
    call ccmpi_abort(-1)
end select

! allocate arrays and open new file
if ( first ) then
  if ( myid==0 ) then
    write(6,*) "Initialise sub hourly output"
  end if
  allocate(freqstore(ifull,freqvars))
  freqstore(:,:) = 0._8
  if ( local ) then
    write(ffile,"(a,'.',i6.6)") trim(freqfile), vnode_vleaderid
  else
    ffile = freqfile
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
    if ( local ) then
      call ccnf_def_dim(fncid,'processor',vnode_nproc,adim(4)) 
      if ( myid==0 ) then
        call ccnf_def_dim(fncid,'gprocessor',nproc,gpdim(1)) 
      else
        gpdim(1)=0
      end if
    end if
    tlenhf = ntau/tbave10
    call ccnf_def_dim(fncid,'time',tlenhf,adim(d4))  
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
    if ( local ) then
      call ccnf_def_var(fncid,'processor','float',1,adim(4:4),idnp)  
      if ( myid==0 ) then
        call ccnf_def_var(fncid,'gprocnode','int',1,gpdim(1:1),idgpn)
        call ccnf_def_var(fncid,'gprocoffset','int',1,gpdim(1:1),idgpo)
      end if
    end if
    call ccnf_def_var(fncid,'time','double',1,adim(d4:d4),idnt)
    call ccnf_put_att(fncid,idnt,'point_spacing','even')
    icy = kdate/10000
    icm = max(1,min(12,(kdate-icy*10000)/100))
    icd = max(1,min(31,(kdate-icy*10000-icm*100)))
    if ( icy<100 ) then
      icy = icy + 1900
    end if
    ich = ktime/100
    icmi = ktime - ich*100
    ics = 0
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
    
    if ( freq_shep ) then
      lname = 'Surface geopotential'
      call attrib(fncid,sdim(1:fsize),fsize,'zht',lname,'m-2 s-2',-1000.,90.e3,any_m,fixed_m,amean_m,float_m)
      lname = 'Soil type'        
      call attrib(fncid,sdim(1:fsize),fsize,'soilt',lname,'none',-650.,650.,any_m,fixed_m,anotdef_m,short_m)
    end if
    if ( freq_core ) then
      lname='x-component 10m wind'
      call attrib(fncid,sdim,ssize,'uas',lname,'m s-1',-130.,130.,any_m,point_m,amean_m,short_m)
      lname='y-component 10m wind'     
      call attrib(fncid,sdim,ssize,'vas',lname,'m s-1',-130.,130.,any_m,point_m,amean_m,short_m)
      lname='Near-Surface Air Temperature'     
      call attrib(fncid,sdim,ssize,'tscrn',lname,'K',100.,425.,any_m,point_m,amean_m,short_m)
      lname='Precipitation'
      call attrib(fncid,sdim,ssize,'rnd',lname,'mm day-1',0.,1300.,any_m,tmean_m,amean_m,float_m)
      lname = 'Scaled Log Surface pressure'
      call attrib(fncid,sdim,ssize,'psf',lname,'none',-1.4,0.5,any_m,point_m,amean_m,short_m)
    end if
    if ( freq_standard ) then
      lname='Near-Surface Relative Humidity'     
      call attrib(fncid,sdim,ssize,'rhscrn',lname,'%',0.,200.,any_m,point_m,amean_m,short_m)
      lname='Convective Precipitation'
      call attrib(fncid,sdim,ssize,'rnc',lname,'mm day-1',0.,1300.,any_m,tmean_m,amean_m,float_m)
    end if
    if ( freq_shep ) then
      lname = 'Evaporation'
      call attrib(fncid,sdim,ssize,'evspsbl',lname,'mm day-1',-1300.,1300.,any_m,tmean_m,amean_m,float_m)        
      lname = 'Screen mixing ratio'
      call attrib(fncid,sdim,ssize,'qgscrn',lname,'kg kg-1',0.,0.06,any_m,point_m,amean_m,short_m)
      lname ='Sea Level Pressure'
      call attrib(fncid,sdim,ssize,'pmsl',lname,'hPa',800.,1200.,any_m,point_m,amean_m,short_m)
      lname ='Surface Downwelling Shortwave Radiation'
      call attrib(fncid,sdim,ssize,'sgdn_ave',lname,'W m-2',-500.,2.e3,any_m,tmean_m,amean_m,float_m)
      lname = 'Surface Downwelling Longwave Radiation'
      call attrib(fncid,sdim,ssize,'rgdn_ave',lname,'W m-2',-500.,1.e3,any_m,tmean_m,amean_m,float_m)
      lname = 'Surface Temperature'
      call attrib(fncid,sdim,ssize,'tsu',lname,'K',100.,425.,any_m,point_m,amean_m,short_m)
      lname='Snowfall Flux'
      call attrib(fncid,sdim,ssize,'sno',lname,'mm day-1',0.,1300.,any_m,tmean_m,amean_m,float_m)
      lname = 'Surface runoff'
      call attrib(fncid,sdim,ssize,'mrros',lname,'mm day-1',0.,1300.,any_m,tmean_m,land_m,float_m)
      lname = 'Runoff' ! mrro after pcc2hist
      call attrib(fncid,sdim,ssize,'runoff',lname,'mm day-1',0.,1300.,any_m,tmean_m,land_m,float_m)
      lname = 'Snow melt' 
      call attrib(fncid,sdim,ssize,'snm',lname,'mm day-1',0.,1300.,any_m,tmean_m,land_m,float_m)
      lname = 'Solar net at ground (+ve down)'
      call attrib(fncid,sdim,ssize,'sgn_ave',lname,'W m-2',-500.,2000.,any_m,tmean_m,amean_m,float_m)
      lname = 'LW net at ground (+ve up)'
      call attrib(fncid,sdim,ssize,'rgn_ave',lname,'W m-2',-500.,1000.,any_m,tmean_m,amean_m,float_m)
      lname = 'Surface Upward Latent Heat Flux'
      call attrib(fncid,sdim,ssize,'eg_ave',lname,'W m-2',-3000.,3000.,any_m,tmean_m,amean_m,float_m)
      lname = 'Surface Upward Sensible Heat Flux'
      call attrib(fncid,sdim,ssize,'fg_ave',lname,'W m-2',-3000.,3000.,any_m,tmean_m,amean_m,float_m)
      lname = 'Height of Boundary Layer'
      call attrib(fncid,sdim,ssize,'pblh',lname,'m',0.,13000.,any_m,point_m,amean_m,short_m)
      lname = 'Convective Available Potential Energy'
      call attrib(fncid,sdim,ssize,'CAPE',lname,'J kg-1',0.,20000.,any_m,point_m,amean_m,short_m)
      lname = 'Convective Inhibition'
      call attrib(fncid,sdim,ssize,'CIN',lname,'J kg-1',-20000.,0.,any_m,point_m,amean_m,short_m)
      lname = 'Lifted Index'
      call attrib(fncid,sdim,ssize,'LI',lname,'K',-100.,100.,any_m,point_m,amean_m,short_m)      
      lname = 'Near-Surface Wind Speed of Gust'
      call attrib(fncid,sdim,ssize,'wsgs',lname,'m s-1',0.,350.,any_m,point_m,amean_m,short_m)
      do k = 1,10 ! 1000, 925, 850, 700, 600, 500, 400, 300, 250, 200
        press_level = cordex_level_data(k)
        call cordex_name(lname,"x-component ",press_level,"hPa wind")
        call cordex_name(vname,"ua",press_level)
        call attrib(fncid,sdim,ssize,trim(vname),lname,'m s-1',-130.,130.,any_m,point_m,amean_m,short_m)
        call cordex_name(lname,"y-component ",press_level,"hPa wind")
        call cordex_name(vname,"va",press_level)
        call attrib(fncid,sdim,ssize,trim(vname),lname,'m s-1',-130.,130.,any_m,point_m,amean_m,short_m)
        lname = 'Air Temperature'     
        call cordex_name(vname,"ta",press_level)
        call attrib(fncid,sdim,ssize,trim(vname),lname,'K',100.,425.,any_m,point_m,amean_m,short_m)
        lname = 'Specific Humidity'
        call cordex_name(vname,"hus",press_level)
        call attrib(fncid,sdim,ssize,trim(vname),lname,'1',0.,0.06,any_m,point_m,amean_m,short_m)
        lname = 'Geopotential Height'
        call cordex_name(vname,"zg",press_level)
        call attrib(fncid,sdim,ssize,trim(vname),lname,'m',0.,130000.,any_m,point_m,amean_m,short_m)
        lname = 'Upward Air Velocity'
        call cordex_name(vname,"wa",press_level)
        call attrib(fncid,sdim,ssize,trim(vname),lname,'m s-1',-130.,130.,any_m,point_m,amean_m,short_m)
      end do 
    end if

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
  
  if ( freq_shep ) then  
    call histwrt(zs,'zht',fncid,fiarch,local,.true.)
    outdata(:) = real(isoilm_in(:))  ! use the raw soil data here to classify inland water bodies
    call histwrt(outdata,'soilt',fncid,fiarch,local,.true.) ! also defines land-sea mask
  end if    
  
  first=.false.
  if ( myid==0 ) write(6,*) "Finished initialising sub hourly output"
 
end if

! store output
freqstore(1:ifull,1) = freqstore(1:ifull,1) + real(condx*(86400./dt/real(tbave10)),8)
freqstore(1:ifull,2) = freqstore(1:ifull,2) + real(condc*(86400./dt/real(tbave10)),8)
freqstore(1:ifull,3) = freqstore(1:ifull,3) + real(evspsbl*(86400./dt/real(tbave10)),8)
freqstore(1:ifull,4) = freqstore(1:ifull,4) + real(sgdn/real(tbave10),8)
freqstore(1:ifull,5) = freqstore(1:ifull,5) + real(rgdn/real(tbave10),8)
freqstore(1:ifull,6) = freqstore(1:ifull,6) + real(conds*(86400./dt/real(tbave10)),8)
freqstore(1:ifull,7) = freqstore(1:ifull,7) + real(runoff_surface*(86400./dt/real(tbave10)),8)
freqstore(1:ifull,8) = freqstore(1:ifull,8) + real(runoff*(86400./dt/real(tbave10)),8)
freqstore(1:ifull,9) = freqstore(1:ifull,9) + real(snowmelt/real(tbave10),8)
freqstore(1:ifull,10) = freqstore(1:ifull,10) + real(sgsave/real(tbave10),8)
freqstore(1:ifull,11) = freqstore(1:ifull,11) + real(rgn/real(tbave10),8)
freqstore(1:ifull,12) = freqstore(1:ifull,12) + real(eg/real(tbave10),8)
freqstore(1:ifull,13) = freqstore(1:ifull,13) + real(fg/real(tbave10),8)

fiarch = ktau/tbave10

! write data to file
if ( mod(ktau,tbave10)==0 ) then
    
  if ( myid==0 .or. local ) then
    if ( myid==0 ) then
      write(6,*) "write sub-hourly output"
    end if
    tpnt = real(ktau,8)*(real(dt,8)/60._8)
    call ccnf_put_vara(fncid,'time',fiarch,tpnt)
    call ccnf_put_vara(fncid,'kdate',fiarch,kdate)
    call ccnf_put_vara(fncid,'ktime',fiarch,ktime)
    call ccnf_put_vara(fncid,'mtimer',fiarch,mtimer)
  end if
  
  ! record output
  if ( freq_core ) then
    umag = sqrt(u(1:ifull,1)**2+v(1:ifull,1)**2)
    outdata = u10*u(1:ifull,1)/max(umag,1.E-6)
    call histwrt(outdata,"uas",fncid,fiarch,local,.true.)
    outdata = u10*v(1:ifull,1)/max(umag,1.E-6)
    call histwrt(outdata,"vas",fncid,fiarch,local,.true.)
    call histwrt(tscrn,"tscrn",fncid,fiarch,local,.true.)
    outdata = real(freqstore(:,1))
    call histwrt(outdata,"rnd",fncid,fiarch,local,.true.)
    call histwrt(psl,"psf",fncid,fiarch,local,.true.)
  end if
  if ( freq_standard ) then
    call histwrt(rhscrn,"rhscrn",fncid,fiarch,local,.true.)
    outdata = real(freqstore(:,2))
    call histwrt(outdata,"rnc",fncid,fiarch,local,.true.)
  end if
  if ( freq_shep ) then
    outdata = real(freqstore(:,3))  
    call histwrt(outdata,"evspsbl",fncid,fiarch,local,.true.)
    call histwrt(qgscrn,"qgscrn",fncid,fiarch,local,.true.)
    call mslp(pmsl,psl,zs,t)
    call histwrt(pmsl,"pmsl",fncid,fiarch,local,.true.)
    outdata = real(freqstore(:,4))
    call histwrt(outdata,"sgdn_ave",fncid,fiarch,local,.true.)
    outdata = real(freqstore(:,5))
    call histwrt(outdata,"rgdn_ave",fncid,fiarch,local,.true.)
    call histwrt(tss,"tsu",fncid,fiarch,local,.true.)
    outdata = real(freqstore(:,6))
    call histwrt(outdata,"sno",fncid,fiarch,local,.true.)
    outdata = real(freqstore(:,7))
    call histwrt(outdata,"mrros",fncid,fiarch,local,.true.)
    outdata = real(freqstore(:,8))
    call histwrt(outdata,"runoff",fncid,fiarch,local,.true.)
    outdata = real(freqstore(:,9))
    call histwrt(outdata,"snm",fncid,fiarch,local,.true.)
    outdata = real(freqstore(:,10))
    call histwrt(outdata,"sgn_ave",fncid,fiarch,local,.true.)
    outdata = real(freqstore(:,11))
    call histwrt(outdata,"rgn_ave",fncid,fiarch,local,.true.)
    outdata = real(freqstore(:,12))
    call histwrt(outdata,"eg_ave",fncid,fiarch,local,.true.)
    outdata = real(freqstore(:,13))
    call histwrt(outdata,"fg_ave",fncid,fiarch,local,.true.)
    call histwrt(pblh,"pblh",fncid,fiarch,local,.true.)
    call histwrt(cape_d,"CAPE",fncid,fiarch,local,.true.)
    call histwrt(cin_d,"CIN",fncid,fiarch,local,.true.)
    call histwrt(li_d,"LI",fncid,fiarch,local,.true.)
    call histwrt(wsgs,'wsgs',fncid,fiarch,local,.true.)  
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
      call histwrt(ua_level,trim(vname),fncid,fiarch,local,.true.) 
      call cordex_name(vname,"va",press_level)
      call histwrt(va_level,trim(vname),fncid,fiarch,local,.true.) 
      call cordex_name(vname,"ta",press_level)
      call histwrt(ta_level,trim(vname),fncid,fiarch,local,.true.) 
      call cordex_name(vname,"hus",press_level)
      call histwrt(hus_level,trim(vname),fncid,fiarch,local,.true.) 
      call cordex_name(vname,"zg",press_level)
      call histwrt(zg_level,trim(vname),fncid,fiarch,local,.true.) 
      call cordex_name(vname,"wa",press_level)
      call histwrt(wa_level,trim(vname),fncid,fiarch,local,.true.) 
    end do  
  end if
  
  freqstore(:,:) = 0._8

end if

if ( myid==0 .or. local ) then
  ! close file at end of run
  if ( ktau==ntau ) then
    call ccnf_close(fncid)
  end if
end if
      
call END_LOG(outfile_end)
      
return
end subroutine freqfile_10

end module outcdf_freq
