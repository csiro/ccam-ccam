! Conformal Cubic Atmospheric Model
    
! Copyright 2015 Commonwealth Scientific Industrial Research Organisation (CSIRO)
    
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

!      PE model on conformal-cubic grid
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!      input files are :namelist (via file called "input")
!                       "nrun.dat"
!      data input and output file names are specified in namelist 'datafile'
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!     sign convention:
!                      u+ve eastwards  (on the panel)
!                      v+ve northwards (on the panel)

program globpe

use aerointerface                          ! Aerosol interface
use aerosolldr, only : xtosav,xtg,naero  & ! LDR prognostic aerosols
    ,duste,dustwd,dustdd,dust_burden     &
    ,bce,bcwd,bcdd,bc_burden             &
    ,oce,ocwd,ocdd,oc_burden             &
    ,dmse,dmsso2o,dms_burden             &
    ,so2e,so2so4o,so2wd,so2dd,so2_burden &
    ,so4e,so4wd,so4dd,so4_burden         &
    ,Ch_dust,zvolcemi,aeroindir
use arrays_m                               ! Atmosphere dyamics prognostic arrays
use bigxy4_m                               ! Grid interpolation
use cable_ccam, only : proglai             ! CABLE
use carbpools_m, only : carbpools_init   & ! Carbon pools
    ,fpn,frs,frp
use cc_mpi                                 ! CC MPI routines
use cfrac_m                                ! Cloud fraction
use cloudmod                               ! Prognostic cloud fraction
use daviesnudge                            ! Far-field nudging
use diag_m                                 ! Diagnostic routines
use dpsdt_m                                ! Vertical velocity
use epst_m                                 ! Off-centre terms
use estab                                  ! Liquid saturation function
use extraout_m                             ! Additional diagnostics
use gdrag_m, only : gdrag_init             ! Gravity wave drag
use histave_m                              ! Time average arrays
use indata                                 ! Data initialisation
use indices_m                              ! Grid index arrays
use infile                                 ! Input file routines
use kuocomb_m                              ! JLM convection
use latlong_m                              ! Lat/lon coordinates
use leoncld_mod                            ! Prognostic cloud condensate
use liqwpar_m                              ! Cloud water mixing ratios
use map_m                                  ! Grid map arrays
use mlo, only : mlodiag,wlev,mxd,mindep  & ! Ocean physics and prognostic arrays
   ,minwater,zomode,zoseaice,factchseaice
use mlodynamics                            ! Ocean dynamics
use morepbl_m                              ! Additional boundary layer diagnostics
use nesting                                ! Nesting and assimilation
use nharrs_m, only : nharrs_init         & ! Non-hydrostatic atmosphere arrays
   ,lrestart
use nlin_m                                 ! Atmosphere non-linear dynamics
use nsibd_m                                ! Land-surface arrays
use outcdf                                 ! Output file routines
use parmhdff_m                             ! Horizontal diffusion parameters
use pbl_m                                  ! Boundary layer arrays
use permsurf_m, only : permsurf_init       ! Fixed surface arrays
use prec_m                                 ! Precipitation
use raddiag_m                              ! Radiation diagnostic
use river                                  ! River routing
use savuvt_m                               ! Saved dynamic arrays
use savuv1_m                               ! Saved dynamic arrays
use sbar_m                                 ! Saved dynamic arrays
use screen_m                               ! Screen level diagnostics
use seaesfrad_m                            ! SEA-ESF radiation
use sigs_m                                 ! Atmosphere sigma levels
use soil_m                                 ! Soil and surface data
use soilsnow_m                             ! Soil, snow and surface data
use tbar2d_m, only : tbar2d_init           ! Atmosphere dynamics reference temperature
use timeseries, only : write_ts            ! Tracer time series
use tkeeps                                 ! TKE-EPS boundary layer
use tracermodule, only : init_tracer     & ! Tracer routines
   ,trfiles,tracer_mass                  &
   ,interp_tracerflux,tracerlist
use tracers_m                              ! Tracer data
use unn_m                                  ! Saved dynamic arrays
use uvbar_m                                ! Saved dynamic arrays
use vecs_m, only : vecs_init               ! Eigenvectors for atmosphere dynamics
use vecsuv_m                               ! Map to cartesian coordinates
use vegpar_m                               ! Vegetation arrays
use vvel_m                                 ! Additional vertical velocity
use work2_m                                ! Diagnostic arrays
use work3_m                                ! Mk3 land-surface diagnostic arrays
use work3f_m                               ! Grid work arrays
use work3sav_m                             ! Water and tracer saved arrays
use workglob_m                             ! Additional grid interpolation
use xarrs_m                                ! Saved dynamic arrays
use xyzinfo_m                              ! Grid coordinate arrays

implicit none

include 'newmpar.h'                        ! Grid parameters
include 'const_phys.h'                     ! Physical constants
include 'darcdf.h'                         ! Netcdf data
include 'dates.h'                          ! Date data
include 'filnames.h'                       ! Filenames
include 'kuocom.h'                         ! Convection parameters
include 'parm.h'                           ! Model configuration
include 'parmdyn.h'                        ! Dynamics parameters
include 'parmgeom.h'                       ! Coordinate data
include 'parmhor.h'                        ! Horizontal advection parameters
include 'parmsurf.h'                       ! Surface parameters
include 'soilv.h'                          ! Soil parameters
include 'stime.h'                          ! File date data
include 'trcom2.h'                         ! Station data
include 'version.h'                        ! Model version data

#ifdef vampir
#include 'vt_user.inc'
#endif
      
integer leap
common/leap_yr/leap                        ! Leap year (1 to allow leap years)
integer nbarewet,nsigmf
common/nsib/nbarewet,nsigmf                ! Land-surface options

integer, dimension(8) :: tvals1, tvals2, nper3hr
integer ilx, io_nest, iq, irest, isoil
integer jalbfix, jlx, k, kktau
integer mins_dt, mins_gmt, mspeca, mtimer_in, nalpha
integer nlx, nmaxprsav, npa, npb, n3hr
integer nstagin, nstaguin, nwrite, nwtsav, mins_rad, secs_rad, mtimer_sav
integer nn, i, j, mstn, ierr, nperhr, nversion
integer ierr2, kmax, isoth, nsig, lapsbot
real, dimension(:,:), allocatable, save :: dums
real, dimension(:), allocatable, save :: spare1, spare2
real, dimension(:), allocatable, save :: spmean
real, dimension(9) :: temparray, gtemparray
real clhav, cllav, clmav, cltav, dsx, dtds, es
real gke, hourst, hrs_dt, evapavge, precavge, preccavge, psavge
real pslavge, pwater, rel_lat, rel_long, spavge, pwatr
real qtot, aa, bb, cc, bb_2, cc_2, rat
real targetlev
real, parameter :: con = 180./pi
character(len=60) comm, comment
character(len=47) header
character(len=10) timeval
character(len=8) rundate
logical odcalc

! version namelist
namelist/defaults/nversion
! main namelist
namelist/cardin/comment,dt,ntau,nwt,npa,npb,nhorps,nperavg,ia,ib, &
    ja,jb,id,jd,iaero,khdif,khor,nhorjlm,mex,mbd,nbd,ndi,ndi2,    &
    nhor,nlv,nmaxpr,nrad,ntaft,ntsea,ntsur,nvmix,restol,          &
    precon,kdate_s,ktime_s,leap,newtop,mup,lgwd,ngwd,rhsat,       &
    nextout,jalbfix,nalpha,nstag,nstagu,ntbar,nwrite,irest,nrun,  &
    nstn,rel_lat,rel_long,nrungcm,nsib,istn,jstn,iunp,slat,slon,  &
    zstn,name_stn,mh_bs,nritch_t,nt_adv,mfix,mfix_qg,namip,       &
    amipo3,nh,nhstest,nsemble,nspecial,panfg,panzo,nplens,rlatdn, &
    rlatdx,rlongdn,rlongdx,newrough,nglacier,newztsea,epsp,epsu,  &
    epsf,av_vmod,charnock,chn10,snmin,tss_sh,vmodmin,zobgin,      &
    rlong0,rlat0,schmidt,kbotdav,kbotu,nbox,nud_p,nud_q,nud_t,    &
    nud_uv,nud_hrs,nudu_hrs,nlocal,nbarewet,nsigmf,qgmin,io_in,   &
    io_nest,io_out,io_rest,tblock,tbave,localhist,unlimitedhist,  &
    m_fly,mstn,nqg,nurban,nmr,ktopdav,nud_sst,nud_sss,mfix_tr,    &
    mfix_aero,kbotmlo,ktopmlo,mloalpha,nud_ouv,nud_sfh,bpyear,    &
    rescrn,helmmeth,nmlo,ol,mxd,mindep,minwater,ocnsmag,ocneps,   &
    mlodiff,zomode,zoseaice,factchseaice,knh,ccycle,kblock,       &
    nud_aero,ch_dust,zvolcemi,aeroindir,helim,fc2,sigbot_gwd,     &
    alphaj,proglai,cgmap_offset,cgmap_scale,compression,filemode, &
    procformat,procmode,chunkoverride,pio
! radiation namelist
namelist/skyin/mins_rad,sw_resolution,sw_diff_streams
! file namelist
namelist/datafile/ifile,ofile,albfile,co2emfile,eigenv,hfile,     &
    icefile,mesonest,nmifile,o3file,radfile,restfile,rsmfile,     &
    scamfile,scrnfile,snowfile,so4tfile,soilfile,sstfile,surfile, &
    tmaxfile,tminfile,topofile,trcfil,vegfile,zofile,smoistfile,  &
    soil2file,radonemfile,co2_00,radon_00,surf_00,co2_12,         &
    radon_12,surf_12,laifile,albnirfile,urbanfile,bathfile,       &
    vegprev,vegnext,cnsdir,salfile,oxidantfile,casafile,phenfile
! convection and cloud microphysics namelist
namelist/kuonml/alflnd,alfsea,cldh_lnd,cldm_lnd,cldl_lnd,         &
    cldh_sea,cldm_sea,cldl_sea,convfact,convtime,shaltime,        &
    detrain,detrainx,dsig2,dsig4,entrain,fldown,iterconv,ksc,     &
    kscmom,kscsea,ldr,mbase,mdelay,methdetr,methprec,nbase,       &
    nclddia,ncvcloud,ncvmix,nevapcc,nevapls,nkuo,nrhcrit,         &
    nstab_cld,nuvconv,rhcv,rhmois,rhsat,sigcb,sigcll,sig_ct,      &
    sigkscb,sigksct,tied_con,tied_over,tied_rh,comm,acon,bcon,    &
    rcm,rcrit_l,rcrit_s,ncloud
! boundary layer turbulence namelist
namelist/turbnml/be,cm0,ce0,ce1,ce2,ce3,cq,ent0,dtrn0,dtrc0,m0,   &
    b1,b2,buoymeth,icm1,maxdts,mintke,mineps,minl,maxl

data nversion/0/
data comment/' '/,comm/' '/,irest/1/,jalbfix/1/,nalpha/1/
data mins_rad/-1/,nwrite/0/
data lapsbot/0/,io_nest/1/
      
#ifndef stacklimit
! For linux only - removes stacklimit on all processors
call setstacklimit(-1)
#endif

#ifdef i8r8
if ( kind(iq)/=8 .or. kind(es)/=8 ) then
  write(6,*) "ERROR: CCAM configured for double precision"
  stop
end if
#else
if ( kind(iq)/=4 .or. kind(es)/=4 ) then
  write(6,*) "ERROR: CCAM configured for single precision"
  stop
end if
#endif


!--------------------------------------------------------------
! INITALISE MPI ROUTINES
call ccmpi_init


!--------------------------------------------------------------
! INITALISE TIMING LOGS
call log_off()
call log_setup()
call START_LOG(model_begin)


!--------------------------------------------------------------
! READ NAMELISTS AND SET PARAMETER DEFAULTS
ia             = -1   ! diagnostic index
ib             = -1   ! diagnostic index
ntbar          = -1
rel_lat        = 0.
rel_long       = 0.
ktau           = 0
ol             = 20   ! default ocean levels
nhor           = -157
nhorps         = -1
khor           = -8
khdif          = 2
nhorjlm        = 1

! All processors read the namelist, so no MPI comms are needed
open(99,file="input",form="formatted",status="old")
read(99, defaults)
if ( myid == 0 ) then
  write(6,'(a20," running for nproc =",i7)') version,nproc
  write(6,*) 'Using defaults for nversion = ',nversion
end if
if ( nversion /= 0 ) then
  call change_defaults(nversion,mins_rad)
end if
read(99, cardin)
call ccmpi_shared_split
call ccmpi_node_leader
nperday = nint(24.*3600./dt)
nperhr  = nint(3600./dt)
do n3hr = 1,8
  nper3hr(n3hr) = nint(n3hr*3*3600/dt)
end do
if ( nwt == -99 )     nwt = nperday      ! set default nwt to 24 hours
if ( nperavg == -99 ) nperavg = nwt      ! set default nperavg to nwt
if ( nwrite == 0 )    nwrite = nperday   ! only used for outfile IEEE
if ( nmlo/=0 .and. abs(nmlo)<=9 ) then
  ol = max( ol, 1 )
else
  ol = 0
end if
wlev     = ol
mindep   = max( 0., mindep )
minwater = max( 0., minwater )
read(99, skyin)
read(99, datafile)
read(99, kuonml)
! try reading boundary layer turbulence namelist
read(99, turbnml, iostat=ierr)
if ( ierr /= 0 ) rewind(99)       ! rewind namelist if turbnml is not found
! try reading tracer namelist
ngas = 0
read(99, trfiles, iostat=ierr)        
if ( ierr /= 0 ) rewind(99)       ! rewind namelist if trfiles is not found
nagg = max( 10, naero )           ! maximum size of aggregation
nlx        = 0
mtimer_sav = 0


!--------------------------------------------------------------
! READ TOPOGRAPHY FILE TO DEFINE CONFORMAL CUBIC GRID
il_g    = 48 ! default global grid size
rlong0  = 0. ! default longitude
rlat0   = 0. ! default latitude
schmidt = 1. ! default schmidt factor for grid stretching
kl      = 18 ! default number of vertical levels
if ( myid==0 .and. io_in<=4 ) then
  ! open topo file and check its dimensions
  ! here used to supply rlong0,rlat0,schmidt
  ! Remander of topo file is read in indata.f90
  write(6,*) 'reading topofile header'
  call ccnf_open(topofile,ncidtopo,ierr)
  if ( ierr == 0 ) then
    ! Netcdf format
    lnctopo = 1 ! flag indicating netcdf file
    call ccnf_inq_dimlen(ncidtopo,'longitude',ilx)
    call ccnf_inq_dimlen(ncidtopo,'latitude',jlx)
    call ccnf_get_attg(ncidtopo,'lon0',rlong0)
    call ccnf_get_attg(ncidtopo,'lat0',rlat0)
    call ccnf_get_attg(ncidtopo,'schmidt',schmidt) 
  else
    ! ASCII format      
    lnctopo = 0 ! flag indicating ASCII file
    open(66,file=topofile,recl=2000,status='old',iostat=ierr)
    if ( ierr /= 0 ) then
      write(6,*) "Error opening topofile ",trim(topofile)
      call ccmpi_abort(-1)
    end if
    read(66,*) ilx,jlx,rlong0,rlat0,schmidt,dsx,header
  end if ! (ierr==0) ..else..
  il_g = ilx        
  write(6,*) 'ilx,jlx              ',ilx,jlx
  write(6,*) 'rlong0,rlat0,schmidt ',rlong0,rlat0,schmidt
end if      ! (myid==0.and.io_in<=4)
! store grid dimensions for broadcast below
temparray(1) = rlong0
temparray(2) = rlat0
temparray(3) = schmidt
temparray(4) = real(il_g)


!--------------------------------------------------------------
! READ EIGENV FILE TO DEFINE VERTICAL LEVELS
if ( myid == 0 ) then
  ! Remanded of file is read in indata.f
  open(28,file=eigenv,status='old',form='formatted',iostat=ierr)
  if ( ierr /= 0 ) then
    write(6,*) "Error opening eigenv file ",trim(eigenv)
    call ccmpi_abort(-1)
  end if
  read(28,*)kmax,lapsbot,isoth,nsig
  kl = kmax
  write(6,*)'kl,ol:              ',kl,ol
  write(6,*)'lapsbot,isoth,nsig: ',lapsbot,isoth,nsig
  temparray(5) = real(kl)
  temparray(6) = real(lapsbot)
  temparray(7) = real(isoth)
  temparray(8) = real(nsig)
end if
      
! Broadcast grid data to all processors
! (Since integers are smaller than 1e7, then they can be exactly
!  represented using real*4)
call ccmpi_bcast(temparray(1:8),0,comm_world)
rlong0  = temparray(1)
rlat0   = temparray(2)
schmidt = temparray(3)
il_g    = nint(temparray(4))
kl      = nint(temparray(5))
lapsbot = nint(temparray(6))
isoth   = nint(temparray(7))
nsig    = nint(temparray(8))

      
!--------------------------------------------------------------
! DEFINE newmpar VARIABLES AND DEFAULTS
! CCAM supports face and uniform grid decomposition over processes using preprocessor directives
! Face decomposition reduces MPI message passing, but only works for factors or multiples of six
! processes.  Uniform decomposition is less restrictive on the number of processes, but requires
! more MPI message passing.
#ifdef uniform_decomp
if ( myid == 0 ) then
  write(6,*) "Using uniform grid decomposition"
end if
#else
if ( myid == 0 ) then
  write(6,*) "Using face grid decomposition"
end if
if ( mod(nproc,6)/=0 .and. mod(6,nproc)/=0 ) then
  write(6,*) "ERROR: nproc must be a multiple of 6 or a factor of 6"
  call ccmpi_abort(-1)
end if
#endif
call proctest(npanels,il_g,nproc,nxp,nyp)
if ( nxp <= 0 ) then
  call badnproc(npanels,il_g,nproc)
end if
jl_g    = il_g + npanels*il_g                 ! size of grid along all panels (usually 6*il_g)
ifull_g = il_g*jl_g                           ! total number of global horizontal grid points
iquad   = 1 + il_g*((8*npanels)/(npanels+4))  ! grid size for interpolation calculations
il      = il_g/nxp                            ! local grid size on process in X direction
jl      = jl_g/nyp                            ! local grid size on process in Y direction
ifull   = il*jl                               ! total number of local horizontal grid points
! The perimeter of the processor region has length 2*(il+jl).
! The first row has 8 possible corner points per panel and the 
! second has 16. In practice these are not all distinct so there could
! be some optimisation.
#ifdef uniform_decomp
npan   = npanels + 1              ! number of panels on this process
iextra = (4*(il+jl)+24)*npan      ! size of halo for MPI message passing
#else      
npan   = max(1,(npanels+1)/nproc) ! number of panels on this process
iextra = 4*(il+jl) + 24*npan      ! size of halo for MPI message passing
#endif
! nrows_rad is a subgrid decomposition for radiation routines
nrows_rad = jl/6
do while( mod(jl,nrows_rad) /= 0 )
  nrows_rad = nrows_rad - 1
end do
if ( myid == 0 ) then
  write(6,*) "il_g,jl_g,il,jl   ",il_g,jl_g,il,jl
  write(6,*) "nxp,nyp,nrows_rad ",nxp,nyp,nrows_rad
end if

! some default values for unspecified parameters
if ( ia < 0 ) ia = il/2
if ( ib < 0 ) ib = ia + 3
if ( ldr == 0 ) mbase = 0
dsig4 = max(dsig2+.01,dsig4)
if( mbd/=0 .and. nbd/=0 ) then
  write(6,*) 'setting nbd=0 because mbd/=0'
  nbd = 0
endif
nud_hrs = abs(nud_hrs)  ! just for people with old -ves in namelist
if ( nudu_hrs == 0 ) nudu_hrs=nud_hrs


! **** do namelist fixes above this ***
      

!--------------------------------------------------------------
! DISPLAY NAMELIST
if ( myid == 0 ) then   
  write(6,*)'Dynamics options A:'
  write(6,*)'   mex   mfix  mfix_qg   mup    nh    precon' 
  write(6,'(i4,i6,i10,3i7)')mex,mfix,mfix_qg,mup,nh,precon
  write(6,*)'Dynamics options B:'
  write(6,*)'nritch_t ntbar  epsp    epsu   epsf   restol'
  write(6,'(i5,i7,1x,3f8.3,g9.2)')nritch_t,ntbar,epsp,epsu,epsf,restol
  write(6,*)'Dynamics options C:'
  write(6,*)'helmmeth mfix_aero mfix_tr'
  write(6,'(i8,i10,i8)') helmmeth,mfix_aero,mfix_tr
  write(6,*)'Horizontal advection/interpolation options:'
  write(6,*)' nt_adv mh_bs'
  write(6,'(i5,i7)') nt_adv,mh_bs
  write(6,*)'Horizontal wind staggering options:'
  write(6,*)'nstag nstagu'
  write(6,'(2i7)') nstag,nstagu
  write(6,*)'Horizontal mixing options:'
  write(6,*)' khdif  khor   nhor   nhorps nhorjlm'
  write(6,'(i5,11i7)') khdif,khor,nhor,nhorps,nhorjlm
  write(6,*)'Vertical mixing/physics options A:'
  write(6,*)' nvmix nlocal ncvmix  lgwd' 
  write(6,'(i5,6i7)') nvmix,nlocal,ncvmix,lgwd
  write(6,*)'Vertical mixing/physics options B:'
  write(6,*)' be   cm0  ce0  ce1  ce2  ce3  cq'
  write(6,'(7f5.2)') be,cm0,ce0,ce1,ce2,ce3,cq
  write(6,*)'Vertical mixing/physics options C:'
  write(6,*)' ent0  dtrn0 dtrc0   m0    b1    b2'
  write(6,'(6f6.2)') ent0,dtrn0,dtrc0,m0,b1,b2
  write(6,*)'Vertical mixing/physics options D:'
  write(6,*)' buoymeth icm1 maxdts'
  write(6,'(i9,i5,f7.1)') buoymeth,icm1,maxdts
  write(6,*)'Vertical mixing/physics options E:'
  write(6,*)'  mintke   mineps     minl     maxl'
  write(6,'(4g9.2)') mintke,mineps,minl,maxl
  write(6,*)'Vertical mixing/physics options F:'
  write(6,*)'  cgmap_offset   cgmap_scale'
  write(6,'(2f14.2)') cgmap_offset,cgmap_scale  
  write(6,*)'Gravity wave drag options:'
  write(6,*)' ngwd   helim     fc2  sigbot_gwd  alphaj'
  write(6,'(i5,2x,3f8.2,f12.6)') ngwd,helim,fc2,sigbot_gwd,alphaj
  write(6,*)'Cumulus convection options A:'
  write(6,*)' nkuo  sigcb sig_ct  rhcv  rhmois rhsat convfact convtime shaltime'
  write(6,'(i5,6f7.2,3x,9f8.2)') nkuo,sigcb,sig_ct,rhcv,rhmois,rhsat,convfact,convtime,shaltime
  write(6,*)'Cumulus convection options B:'
  write(6,*)' alflnd alfsea fldown iterconv ncvcloud nevapcc nevapls nuvconv'
  write(6,'(3f7.2,i6,i10,4i8)') alflnd,alfsea,fldown,iterconv,ncvcloud,nevapcc,nevapls,nuvconv
  write(6,*)'Cumulus convection options C:'
  write(6,*)' mbase mdelay methprec nbase detrain entrain methdetr detrainx dsig2  dsig4'
  write(6,'(3i6,i9,f8.2,f9.2,i8,4f8.2)') mbase,mdelay,methprec,nbase,detrain,entrain,methdetr,detrainx,dsig2,dsig4
  write(6,*)'Shallow convection options:'
  write(6,*)'  ksc  kscsea kscmom sigkscb sigksct tied_con tied_over tied_rh '
  write(6,'(i5,2i7,1x,3f8.3,2f10.3)') ksc,kscsea,kscmom,sigkscb,sigksct,tied_con,tied_over,tied_rh
  write(6,*)'Other moist physics options:'
  write(6,*)'  acon   bcon   qgmin      rcm    rcrit_l rcrit_s'
  write(6,'(2f7.2,2e10.2,2f7.2)') acon,bcon,qgmin,rcm,rcrit_l,rcrit_s
  write(6,*)'Radiation options A:'
  write(6,*)' nrad  mins_rad iaero  dt'
  write(6,'(i5,2i7,f10.2)') nrad,mins_rad,iaero,dt
  write(6,*)'Radiation options B:'
  write(6,*)' nmr bpyear sw_diff_streams sw_resolution'
  write(6,'(i4,f9.2,i4,a5)') nmr,bpyear,sw_diff_streams,sw_resolution
  write(6,*)'Aerosol options:'
  write(6,*)'  iaero ch_dust'
  write(6,'(i7,g9.2,f7.2)') iaero,ch_dust
  write(6,*)'  zvolcemi aeroindir'
  write(6,'(f7.2,i5)') zvolcemi,aeroindir
  write(6,*)'Cloud options A:'
  write(6,*)'  ldr nclddia nstab_cld nrhcrit sigcll '
  write(6,'(i5,i6,2i9,1x,f8.2)') ldr,nclddia,nstab_cld,nrhcrit,sigcll
  write(6,*)'Cloud options B:'
  write(6,*)'  ncloud'
  write(6,'(1i5)') ncloud
  write(6,*)'Soil, canopy and PBL options A:'
  write(6,*)' jalbfix nalpha nbarewet newrough nglacier nrungcm nsib  nsigmf'
  write(6,'(i5,9i8)') jalbfix,nalpha,nbarewet,newrough,nglacier,nrungcm,nsib,nsigmf
  write(6,*)'Soil, canopy and PBL options B:'
  write(6,*)' ntaft ntsea ntsur av_vmod tss_sh vmodmin  zobgin charnock chn10'
  write(6,'(i5,2i6,4f8.2,f8.3,f9.5)') ntaft,ntsea,ntsur,av_vmod,tss_sh,vmodmin,zobgin,charnock,chn10
  write(6,*)'Soil, canopy and PBL options C:'
  write(6,*)' nurban ccycle'
  write(6,'(2i7)') nurban,ccycle
  write(6,*)'Ocean/lake options:'
  write(6,*)' nmlo  ol      mxd   mindep minwater  ocnsmag   ocneps'
  write(6,'(i5,i4,5f9.2)') nmlo,ol,mxd,mindep,minwater,ocnsmag,ocneps
  write(6,*)' mlodiff  zomode zoseaice factchseaice'
  write(6,'(2i8,f9.6,f13.6)') mlodiff,zomode,zoseaice,factchseaice
  if ( mbd/=0 .or. nbd/=0 ) then
    write(6,*)'Nudging options A:'
    write(6,*)' nbd    nud_p  nud_q  nud_t  nud_uv nud_hrs nudu_hrs kbotdav  kbotu'
    write(6,'(i5,3i7,7i8)') nbd,nud_p,nud_q,nud_t,nud_uv,nud_hrs,nudu_hrs,kbotdav,kbotu
    write(6,*)'Nudging options B:'
    write(6,*)' mbd    ktopdav kblock'
    write(6,'(i5,2i8)') mbd,ktopdav,kblock
    write(6,*)'Nudging options C:'
    write(6,*)' nud_sst nud_sss nud_ouv nud_sfh ktopmlo kbotmlo mloalpha'
    write(6,'(6i8,i9)') nud_sst,nud_sss,nud_ouv,nud_sfh,ktopmlo,kbotmlo,mloalpha
  end if
  write(6,*)'Special and test options A:'
  write(6,*)' namip amipo3 newtop nhstest nplens nsemble nspecial panfg panzo'
  write(6,'(1i5,L7,4i7,i8,f9.1,f8.4)') namip,amipo3,newtop,nhstest,nplens,nsemble,nspecial,panfg,panzo
  write(6,*)'Special and test options B:'
  write(6,*)' knh rescrn'
  write(6,'(i4,i7)') knh,rescrn
  write(6,*)'I/O options:'
  write(6,*)' m_fly  io_in io_nest io_out io_rest  nwt  nperavg'
  write(6,'(i5,4i7,3i8)') m_fly,io_in,io_nest,io_out,io_rest,nwt,nperavg

  write(6, cardin)
  if ( nllp==0 .and. nextout>=4 ) then
    write(6,*) 'need nllp=3 for nextout>=4'
    call ccmpi_abort(-1)
  end if
  write(6, skyin)
  write(6, datafile)
  write(6, kuonml)
  write(6, turbnml)
end if ! myid=0
if ( newtop > 2 ) then
  write(6,*) 'newtop>2 no longer allowed'
  call ccmpi_abort(-1)
end if
if ( mfix_qg>0 .and. nkuo==4 ) then
  write(6,*) 'nkuo=4: mfix_qg>0 not allowed'
  call ccmpi_abort(-1)
end if
if ( mfix > 3 ) then
  write(6,*) 'mfix >3 not allowed now'
  call ccmpi_abort(-1)
end if
nstagin  = nstag    ! -ve nstagin gives swapping & its frequency
nstaguin = nstagu   ! only the sign of nstaguin matters (chooses scheme)
if ( nstagin==5 .or. nstagin<0 ) then
  nstag  = 4
  nstagu = 4
  if ( nstagin == 5 ) then  ! for backward compatability
    nstagin  = -1 
    nstaguin = 5  
  endif
endif
if ( kblock < 0 ) kblock = max( kl, ol ) ! must occur before indata
if ( mod(ntau,tblock*tbave)/=0 ) then
  write(6,*) "ERROR: tblock*tave must be a factor of ntau"
  write(6,*) "ntau,tblock,tbave ",ntau,tblock,tbave
  call ccmpi_abort(-1)
end if
if ( filemode.ge.2 .and. compression.gt.0 ) then
  write(6,*) "ERROR: NetCDF-3 file format cannot be used with compression"
  write(6,*) "filemode > 1 compression must equal 0"
  write(6,*) "filemode,compression ",filemode,compression
  call ccmpi_abort(-1)
end if
if ( procformat .and. procmode.gt.0 ) then
  if ( mod(nproc,procmode).ne.0 ) then
    write(6,*) "ERROR: procmode must be a multiple of the number of ranks"
    write(6,*) "nproc,procmode ",nproc,procmode
    call ccmpi_abort(-1)
  end if
end if
if ( pio .and. .not.procformat ) then
   write(6,*) "ERROR: Parallel I/O is not supported without procformat"
   write(6,*) "pio,procformat ",pio,procformat
   call ccmpi_abort(-1)
else
   if ( compression.gt.0 ) then
      write(6,*) "ERROR: Compression not supported with Parallel I/O"
      write(6,*) "pio,compression ",pio,compression
      call ccmpi_abort(-1)
   end if
end if


!--------------------------------------------------------------
! INITIALISE ifull_g ALLOCATABLE ARRAYS
call bigxy4_init(iquad)
call xyzinfo_init(ifull_g,ifull,iextra,myid,mbd,nud_uv)
call indices_init(ifull_g,ifull,iextra,npanels,npan)
call map_init(ifull_g,ifull,iextra,myid,mbd)
call latlong_init(ifull_g,ifull,iextra,myid)      
call vecsuv_init(ifull_g,ifull,iextra,myid)


!--------------------------------------------------------------
! SET UP CC GEOMETRY
! Only one process calls setxyz to save memory with large grids
if ( myid == 0 ) then
  write(6,*) "Calling setxyz"
  call workglob_init(ifull_g)
  call setxyz(il_g,rlong0,rlat0,schmidt,x_g,y_g,z_g,wts_g,ax_g,ay_g,az_g,bx_g,by_g,bz_g,xx4,yy4)
end if
! Broadcast the following global data.  xx4 and yy4 are used for calculating depature points.
call ccmpi_bcast(ds,0,comm_world)
call ccmpi_bcastr8(xx4,0,comm_world)
call ccmpi_bcastr8(yy4,0,comm_world)
! The following are only needed for the scale-selective filter
if ( mbd /= 0 ) then
  ! only need x_g, y_g and z_g for 2D filter.  1D filter recalculates
  ! these arrays from xx4 and yy4
  if ( nud_uv == 9 ) then
    call ccmpi_bcastr8(x_g,0,comm_world)
    call ccmpi_bcastr8(y_g,0,comm_world)
    call ccmpi_bcastr8(z_g,0,comm_world)
  end if
  ! both 1D and 2D filter need em_g
  call ccmpi_bcast(em_g,0,comm_world)
end if

if ( myid == 0 ) then
  write(6,*) "Calling ccmpi_setup"
end if
call ccmpi_setup(kblock)

      
!--------------------------------------------------------------
! DEALLOCATE ifull_g ARRAYS WHERE POSSIBLE
call worklocl_init(ifull)      
if ( myid == 0 ) then
  call ccmpi_distribute(rlong4_l,rlong4)
  call ccmpi_distribute(rlat4_l,rlat4)
  call workglob_end
  deallocate( wts_g, emu_g, emv_g )
  deallocate( ax_g, ay_g, az_g )
  deallocate( bx_g, by_g, bz_g )
  deallocate( f_g, fu_g, fv_g )
  deallocate( dmdx_g, dmdy_g )
  if ( mbd == 0 ) then
    deallocate( x_g, y_g, z_g, em_g )
  else if ( nud_uv /= 9 ) then
    deallocate( x_g, y_g, z_g )
  end if
  deallocate( rlatt_g, rlongg_g )
else
  call ccmpi_distribute(rlong4_l)
  call ccmpi_distribute(rlat4_l)
end if


!--------------------------------------------------------------
! INITIALISE LOCAL ARRAYS
allocate( spare1(ifull), spare2(ifull) )
allocate( dums(ifull,kl), spmean(kl) )
call arrays_init(ifull,iextra,kl)
call carbpools_init(ifull,iextra,kl,nsib,ccycle)
call cfrac_init(ifull,iextra,kl)
call cloudmod_init(ifull,iextra,kl,ncloud)
call dpsdt_init(ifull,iextra,kl)
call epst_init(ifull,iextra,kl)
call estab_init
call extraout_init(ifull,iextra,kl,nextout)
call gdrag_init(ifull,iextra,kl)
call histave_init(ifull,iextra,kl,ms)
call kuocomb_init(ifull,iextra,kl)
call liqwpar_init(ifull,iextra,kl)
call morepbl_init(ifull,iextra,kl)
call nharrs_init(ifull,iextra,kl)
call nlin_init(ifull,iextra,kl)
call nsibd_init(ifull,iextra,kl,nsib)
call parmhdff_init(ifull,iextra,kl)
call pbl_init(ifull,iextra,kl)
call permsurf_init(ifull,iextra,kl)
call prec_init(ifull,iextra,kl)
call raddiag_init(ifull,iextra,kl)
call savuvt_init(ifull,iextra,kl)
call savuv1_init(ifull,iextra,kl)
call sbar_init(ifull,iextra,kl)
call screen_init(ifull,iextra,kl)
call sigs_init(ifull,iextra,kl)
call soil_init(ifull,iextra,kl,iaero,nsib)
call soilsnow_init(ifull,iextra,kl,ms,nsib)
call tbar2d_init(ifull,iextra,kl)
call unn_init(ifull,iextra,kl)
call uvbar_init(ifull,iextra,kl)
call vecs_init(ifull,iextra,kl)
call vegpar_init(ifull,iextra,kl)
call vvel_init(ifull,iextra,kl)
call work2_init(ifull,iextra,kl,nsib)
call work3_init(ifull,iextra,kl,nsib)
call work3f_init(ifull,iextra,kl)
call xarrs_init(ifull,iextra,kl)
if ( nvmix == 6 ) then
  call tkeinit(ifull,iextra,kl,0)
end if
if ( tracerlist /= ' ' ) call init_tracer
call work3sav_init(ifull,iextra,kl,ilt,jlt,klt,ngasmax) ! must occur after tracers_init
if ( nbd/=0 .and. nud_hrs/=0 ) then
  if ( abs(iaero)>=2 .and. nud_aero/=0 ) then
    call dav_init(ifull,iextra,kl,naero)
  else
    call dav_init(ifull,iextra,kl,0)
  end if
end if
! Remaining arrays are allocated in indata.f90, since their
! definition requires additional input data (e.g, land-surface)

      
!--------------------------------------------------------------
! DISPLAY DIAGNOSTIC INDEX AND TIMER DATA
if ( mydiag ) then
  write(6,"('id,jd,rlongg,rlatt in degrees: ',2i4,2f8.2)") id,jd,con*rlongg(idjd),con*rlatt(idjd)
end if
call date_and_time(rundate)
call date_and_time(time=timeval)
if ( myid == 0 ) then
  write(6,*)'RUNDATE IS ',rundate
  write(6,*)'Starting time ',timeval
end if


!--------------------------------------------------------------
! READ INITIAL CONDITIONS
if ( myid == 0 ) then
  write(6,*) "Calling indata"
end if
call indataf(hourst,jalbfix,lapsbot,isoth,nsig,io_nest)

! fix nudging levels from pressure to level index
! this is done after indata has loaded sig
if ( kbotdav < 0 ) then
  targetlev = real(-kbotdav)/1000.
  do k = 1,kl
    if ( sig(k) <= targetlev ) then
      kbotdav = k
      if ( myid==0 ) then
        write(6,*) "kbotdav adjusted to ",kbotdav,"for sig ",sig(kbotdav)
      end if
      exit
    end if
  end do
  if ( kbotdav < 0 ) then
    write(6,*) "ERROR: Cannot locate nudging level for kbotdav ",kbotdav
    call ccmpi_abort(-1)
  end if
end if
if ( ktopdav == 0 ) then
  ktopdav = kl
else if ( ktopdav < 0 ) then
  targetlev = real(-ktopdav)/1000.
  do k = kl,1,-1
    if ( sig(k) >= targetlev ) then
      ktopdav = k
      if ( myid == 0 ) then
        write(6,*) "ktopdav adjusted to ",ktopdav,"for sig ",sig(ktopdav)
      end if
      exit
    end if
  end do
  if ( ktopdav < 0 ) then
    write(6,*) "ERROR: Cannot locate nudging level for ktopdav ",ktopdav
    call ccmpi_abort(-1)
  end if
end if
! fix ocean nuding levels
if ( kbotmlo == -1000 ) then
  kbotmlo = ol
else if ( kbotmlo < 0 )  then
  targetlev = real(-kbotmlo)/1000.
  do k = ol,1,-1
    if ( gosig(k) <= targetlev ) then
      kbotmlo = k
      if ( myid == 0 ) then
        write(6,*) "kbotmlo adjusted to ",kbotmlo,"for sig ",gosig(kbotmlo)
      end if
      exit
    end if
  end do
  if ( kbotmlo < 0 ) then
    write(6,*) "ERROR: Cannot locate nudging level for kbotmlo ",kbotmlo
    call ccmpi_abort(-1)
  end if   
end if
if ( ktopmlo < 0 ) then
  targetlev = real(-ktopmlo)/1000.
  do k = 1,ol
    if ( gosig(k) >= targetlev ) then
      ktopmlo = k
      if ( myid == 0 ) then
        write(6,*) "ktopmlo adjusted to ",ktopmlo,"for sig ",gosig(ktopmlo)
      end if
      exit
    end if
  end do
  if ( ktopmlo < 0 ) then
    write(6,*) "ERROR: Cannot locate nudging level for ktopmlo ",ktopmlo
    call ccmpi_abort(-1)
  end if
end if
if ( (ktopmlo<1.or.kbotmlo>ol.or.ktopmlo>kbotmlo) .and. nmlo/=0 ) then
  write(6,*) "ERROR: Invalid kbotmlo"
  write(6,*) "kbotmlo,ktopmlo ",kbotmlo,ktopmlo
  call ccmpi_abort(-1)
end if
if ( kbotdav<1 .or. ktopdav>kl .or. kbotdav>ktopdav ) then
  write(6,*) "ERROR: Invalid kbotdav and ktopdav"
  write(6,*) "kbotdav,ktopdav ",kbotdav,ktopdav
  call ccmpi_abort(-1)
end if
if ( kbotu == 0 ) kbotu = kbotdav
! identify reference level ntbar for temperature
if ( ntbar == -1 ) then
  ntbar = 1
  do while( sig(ntbar)>0.8 .and. ntbar<kl )
    ntbar = ntbar + 1
  end do
end if
! estimate radiation calling frequency
if ( mins_rad < 0 ) then
  ! automatic estimate for mins_rad
  secs_rad = min( nint((schmidt*112.*90./real(il_g))*8.*60.), 3600 )
  secs_rad = min( secs_rad, nint(real(nwt)*dt) )
  secs_rad = max( secs_rad, 1 )
  kountr   = nint(real(secs_rad)/dt)
  secs_rad = nint(real(kountr)*dt)
  do while ( mod( 3600, secs_rad )/=0 .and. kountr>1 )
    kountr = kountr - 1
    secs_rad = nint(real(kountr)*dt)
  end do
else
  ! user specified mins_rad
  kountr   = nint(real(mins_rad)*60./dt)  ! set default radiation to ~mins_rad m
  secs_rad = nint(real(kountr)*dt)        ! redefine to actual value
end if
if ( myid == 0 ) then
  write(6,*) "Radiation will use kountr ",kountr," for secs_rad ",secs_rad
end if
! for 6-hourly output of sint_ave etc, want 6*60*60 = N*secs_rad      
if ( (nrad==4.or.nrad==5) .and. mod(21600,secs_rad)/=0 ) then
  write(6,*) 'ERROR: CCAM would prefer 21600 = N*secs_rad ',secs_rad
  call ccmpi_abort(-1)
end if

! max/min diagnostics      
if ( nextout >= 4 ) call setllp
#ifdef debug
call maxmin(u,' u',ktau,1.,kl)
call maxmin(v,' v',ktau,1.,kl)
dums(:,:) = sqrt(u(1:ifull,:)**2+v(1:ifull,:)**2)  ! 3D 
call maxmin(dums,'sp',ktau,1.,kl)
call maxmin(t,' t',ktau,1.,kl)
call maxmin(qg,'qg',ktau,1.e3,kl)
call maxmin(qfg,'qf',ktau,1.e3,kl)
call maxmin(qlg,'ql',ktau,1.e3,kl)
call maxmin(wb,'wb',ktau,1.,ms)
call maxmin(tggsn,'tS',ktau,1.,3)
call maxmin(tgg,'tgg',ktau,1.,ms)
pwatr_l = 0.   ! in mm
do k = 1,kl
  pwatr_l = pwatr_l-sum(dsig(k)*wts(1:ifull)*(qg(1:ifull,k)+qlg(1:ifull,k)+qfg(1:ifull,k))*ps(1:ifull))
enddo
pwatr_l = pwatr_l/grav
temparray(1) = pwatr_l
call ccmpi_reduce( temparray(1:1), gtemparray(1:1), "sum", 0, comm_world )
pwatr = gtemparray(1)
if ( myid == 0 ) write (6,"('pwatr0 ',12f7.3)") pwatr
if ( ntrac > 0 ) then
  do ng = 1,ntrac
    write (text,'("g",i1)')ng
    call maxmin(tr(:,:,ng),text,ktau,1.,kl)
  end do
end if   ! (ntrac>0)
#endif


!--------------------------------------------------------------
! SETUP REMAINING PARAMETERS

! convection
! sig(kuocb) occurs for level just BELOW sigcb
kuocb = 1
do while( sig(kuocb+1)>=sigcb )
  kuocb = kuocb+1
end do
if ( myid == 0 ) write(6,*) 'convective cumulus scheme: kuocb,sigcb = ',kuocb,sigcb

! horizontal diffusion 
if ( khdif == -99 ) then   ! set default khdif appropriate to resolution
  khdif = 5
  if ( myid == 0 ) write(6,*) 'Model has chosen khdif =',khdif
endif
do k = 1,kl
  hdiff(k) = khdif*0.1
end do
if ( khor > 0 ) then
  do k = kl+1-khor,kl
    hdiff(k) = 2.*hdiff(k-1)
  end do
elseif ( khor < 0 ) then ! following needed +hdiff() (JLM 29/6/15)
  do k = 1,kl                    ! N.B. usually hdiff(k)=khdif*.1 
    ! increase hdiff between sigma=.15  and sigma=0., 0 to khor
    if ( sig(k) < 0.15 ) then
      hdiff(k) = .1*max(1.,(1.-sig(k)/.15)*abs(khor)) + hdiff(k)
    end if
  end do
  if ( myid == 0 ) write(6,*)'khor,hdiff: ',khor,hdiff
end if
if ( nud_p==0 .and. mfix==0 ) then
  write(6,*) "ERROR: Both nud_p=0 and mfix=0"
  write(6,*) "Model will not conserve mass"
  call ccmpi_abort(-1)
end if
if ( nud_q==0 .and. mfix_qg==0 ) then
  write(6,*) "ERROR: Both nud_q=0 and mfix_qg=0"
  write(6,*) "Model will not conserve moisture"
  call ccmpi_abort(-1)
end if
if ( nud_aero==0 .and. mfix_aero==0 .and. iaero/=0 ) then
  write(6,*) "ERROR: Both nud_aero=0 and mfix_aero=0"
  write(6,*) "Model will not conserve aerosols"
  call ccmpi_abort(-1)
end if
if ( mfix_tr==0 .and. ngas>0 ) then
  write(6,*) "ERROR: mfix_tr=0 and ngas>0"
  write(6,*) "Model will not conserve tracers"
  call ccmpi_abort(-1)
end if
      

call printa('zs  ',zs,0,0,ia,ib,ja,jb,0.,.01)
call printa('tss ',tss,0,0,ia,ib,ja,jb,200.,1.)
if ( mydiag ) write(6,*)'wb(idjd) ',(wb(idjd,k),k=1,6)
call printa('wb1   ',wb ,0,1,ia,ib,ja,jb,0.,100.)
call printa('wb6  ',wb,0,ms,ia,ib,ja,jb,0.,100.)

      
!--------------------------------------------------------------
! NRUN COUNTER
if ( myid == 0 ) then
  open(11, file='nrun.dat',status='unknown')
  if ( nrun == 0 ) then
    read(11,*,iostat=ierr2) nrun
    nrun = nrun + 1
  endif                  ! nrun==0
  write(6,*)'this is run ',nrun
  rewind 11
  write(11,*) nrun
  write(11,cardin)
  write(11,datafile)
  close(11)
end if


!-------------------------------------------------------------
! SETUP DIAGNOSTIC ARRAYS
rndmax(:)      = 0.
tmaxscr(:)     = 0.
tminscr(:)     = 400.
rhmaxscr(:)    = 0.
rhminscr(:)    = 400.
u10max(:)      = 0.
v10max(:)      = 0.
u1max(:)       = 0.
v1max(:)       = 0.
u2max(:)       = 0.
v2max(:)       = 0.
cape_max(:)    = 0.
cape_ave(:)    = 0.
u10mx(:)       = 0.
tscr_ave(:)    = 0.
qscrn_ave(:)   = 0.
dew_ave(:)     = 0.
epan_ave(:)    = 0.
epot_ave(:)    = 0.
eg_ave(:)      = 0.
fg_ave(:)      = 0.
ga_ave(:)      = 0.
rnet_ave(:)    = 0.
sunhours(:)    = 0.
riwp_ave(:)    = 0.
rlwp_ave(:)    = 0.
evap(:)        = 0.
precc(:)       = 0.
precip(:)      = 0.
convh_ave(:,:) = 0.
rnd_3hr(:,8)   = 0. ! i.e. rnd24(:)=0.
cbas_ave(:)    = 0.
ctop_ave(:)    = 0.
sno(:)         = 0.
hail(:)        = 0.
runoff(:)      = 0.
wb_ave(:,:)    = 0.
tsu_ave(:)     = 0.
alb_ave(:)     = 0.
fbeam_ave(:)   = 0.
psl_ave(:)     = 0.
mixdep_ave(:)  = 0.
koundiag       = 0
sint_ave(:)    = 0.  ! solar_in_top
sot_ave(:)     = 0.  ! solar_out_top
soc_ave(:)     = 0.  ! solar_out_top (clear sky)
sgdn_ave(:)    = 0.  ! solar_ground (down-welling) +ve down
sgn_ave(:)     = 0.  ! solar_ground (net) +ve down
rtu_ave(:)     = 0.  ! LW_out_top 
rtc_ave(:)     = 0.  ! LW_out_top (clear sky)
rgdn_ave(:)    = 0.  ! LW_ground (down-welling)  +ve down
rgn_ave(:)     = 0.  ! LW_ground (net)  +ve up
rgc_ave(:)     = 0.  ! LW_ground (clear sky)
sgc_ave(:)     = 0.  ! SW_ground (clear sky)
cld_ave(:)     = 0.
cll_ave(:)     = 0.
clm_ave(:)     = 0.
clh_ave(:)     = 0.
if ( ngas > 0 ) then
  traver       = 0.
end if
fpn_ave        = 0.
frs_ave        = 0.
frp_ave        = 0.
if ( abs(iaero) == 2 ) then
  duste        = 0.  ! Dust emissions
  dustdd       = 0.  ! Dust dry deposition
  dustwd       = 0.  ! Dust wet deposition
  dust_burden  = 0.  ! Dust burden
  bce          = 0.  ! Black carbon emissions
  bcdd         = 0.  ! Black carbon dry deposition
  bcwd         = 0.  ! Black carbon wet deposition
  bc_burden    = 0.  ! Black carbon burden
  oce          = 0.  ! Organic carbon emissions
  ocdd         = 0.  ! Organic carbon dry deposition
  ocwd         = 0.  ! Organic carbon wet deposition
  oc_burden    = 0.  ! Organic carbon burden
  dmse         = 0.  ! DMS emissions
  dmsso2o      = 0.  ! DMS -> SO2 oxidation
  so2e         = 0.  ! SO2 emissions
  so2so4o      = 0.  ! SO2 -> SO4 oxidation
  so2dd        = 0.  ! SO2 dry deposition
  so2wd        = 0.  ! SO2 wet deposiion
  so4e         = 0.  ! SO4 emissions
  so4dd        = 0.  ! SO4 dry deposition
  so4wd        = 0.  ! SO4 wet deposition
  dms_burden   = 0.  ! DMS burden
  so2_burden   = 0.  ! SO2 burden
  so4_burden   = 0.  ! SO4 burden
end if


!--------------------------------------------------------------
! OPEN OUTPUT FILES AND SAVE INITAL CONDITIONS
if ( nwt > 0 ) then
  ! write out the first ofile data set
  if ( myid == 0 ) then
    write(6,*)'calling outfile'
  end if
  call outfile(20,rundate,nwrite,nstagin,jalbfix,nalpha,mins_rad)  ! which calls outcdf
  if ( newtop < 0 ) then
    if ( myid == 0 ) then
      write(6,*) "newtop<0 requires a stop here"
    end if
    ! just for outcdf to plot zs  & write fort.22
    call ccmpi_abort(-1)
  end if
end if    ! (nwt>0)


!--------------------------------------------------------------
! INITIALISE DYNAMICS
dtin = dt
n3hr = 1   ! initial value at start of run
if ( myid == 0 ) then
  write(6,*) "number of time steps per day = ",nperday
  write(6,*) "nper3hr,nper6hr .. ",nper3hr(:)
end if
mspeca = 1
if ( mex/=1 .and. .not.lrestart ) then
  mspeca = 2
  dt     = dtin*.5
endif
call gettin(0)              ! preserve initial mass & T fields

nmaxprsav = nmaxpr
nwtsav    = nwt
hrs_dt    = dtin/3600.      ! time step in hours
mins_dt   = nint(dtin/60.)  ! time step in minutes
mtimer_in = mtimer
 
 
!--------------------------------------------------------------
! BEGIN MAIN TIME LOOP
if ( myid == 0 ) then
  call date_and_time(time=timeval,values=tvals1)
  write(6,*) "Start of loop time ", timeval
end if
call log_on()
call START_LOG(maincalc_begin)

do kktau = 1,ntau   ! ****** start of main time loop

  ktau     = kktau
  timer    = timer + hrs_dt                      ! timer now only used to give timeg
  timeg    = mod(timer+hourst,24.)
  mtimer   = mtimer_in + nint(ktau*dtin/60.)     ! 15/6/01 to allow dt < 1 minute
  mins_gmt = mod(mtimer+60*ktime/100,24*60)

  ! ***********************************************************************
  ! START ATMOSPHERE DYNAMICS
  ! ***********************************************************************

  
  ! NESTING ---------------------------------------------------------------
  if ( nbd /= 0 ) then
    ! Newtonian relaxiation
    call START_LOG(nestin_begin)
    call nestin
    call END_LOG(nestin_end)
  end if
      
  
  ! TRACERS ---------------------------------------------------------------
  ! interpolate tracer fluxes to current timestep
  if ( ngas > 0 ) then
    call interp_tracerflux(kdate,hrs_dt)
  end if

  
  ! DYNAMICS --------------------------------------------------------------
  if ( nstaguin>0 .and. ktau>=1 ) then   ! swapping here for nstaguin>0
    if ( nstagin<0 .and. mod(ktau-nstagoff,abs(nstagin))==0 ) then
      nstag  = 7-nstag  ! swap between 3 & 4
      nstagu = nstag
    endif
  endif

  do mspec = mspeca,1,-1    ! start of introductory time loop
    dtds = dt/ds
    un(1:ifull,:) = 0. 
    vn(1:ifull,:) = 0.
    tn(1:ifull,:) = 0.

    if ( mup/=1 .or. (ktau==1.and.mspec==mspeca.and..not.lrestart) ) then
      call bounds(psl)
      ! updps called first step or to permit clean restart option      
      call updps(0) 
    endif

    ! set up tau +.5 velocities in ubar, vbar
    if ( ktau<10 .and. nmaxpr==1 ) then
      if ( myid == 0 ) then
        write(6,*) 'ktau,mex,mspec,mspeca:',ktau,mex,mspec,mspeca
      end if
      call ccmpi_barrier(comm_world)
    endif
    sbar(:,2:kl) = sdot(:,2:kl)
    if ( ktau==1 .and. .not.lrestart ) then
      ! this sets (ubar,vbar) to ktau=1.5 values on 2nd time through
      ubar(:,:) = u(1:ifull,:)
      vbar(:,:) = v(1:ifull,:)
    elseif ( mex == 1) then
      ubar(:,:) = u(1:ifull,:)
      vbar(:,:) = v(1:ifull,:)
    elseif ( (ktau==2.and..not.lrestart) .or. mex==2 ) then        
      ! (tau+.5) from tau, tau-1
      ubar(:,:) = u(1:ifull,:)*1.5 - savu(:,:)*.5
      vbar(:,:) = v(1:ifull,:)*1.5 - savv(:,:)*.5
    elseif ( mex == 3 )then
      ! (tau+.5) from tau, tau-1, tau-2   ! ubar is savu1 here
      ubar(:,:) = u(1:ifull,:)+.5*(savu(:,:)-savu1(:,:))
      vbar(:,:) = v(1:ifull,:)+.5*(savv(:,:)-savv1(:,:))
    elseif ( mex==30 .and. (ktau>3.or.lrestart) ) then  ! using tau, tau-1, tau-2, tau-3
      do k = 1,kl
        do iq = 1,ifull
          bb = 1.5*u(iq,k) - 2.*savu(iq,k) + .5*savu1(iq,k)                             ! simple b
          bb_2 = (40.*u(iq,k) - 35.*savu(iq,k) - 16.*savu1(iq,k) + 11.*savu2(iq,k))/34. ! cwqls b
          cc = .5*u(iq,k) - savu(iq,k) + .5*savu1(iq,k)                                 ! simple c
          cc_2 = (10.*u(iq,k) - 13.*savu(iq,k) - 4.*savu1(iq,k) + 7.*savu2(iq,k))/34.   ! cwqls c
          aa = cc_2 - cc
          rat = max( 0., min( 1., cc_2/(aa+sign(1.e-9,aa)) ) )
          cc = rat*cc + (1.-rat)*cc_2 
          bb = rat*bb + (1.-rat)*bb_2 
          ubar(iq,k) = u(iq,k) + .5*bb + .25*cc
          bb = 1.5*v(iq,k) - 2.*savv(iq,k) + .5*savv1(iq,k)                           ! simple b
          bb_2 = (40.*v(iq,k)-35.*savv(iq,k)-16.*savv1(iq,k)+11.*savv2(iq,k))/34.     ! cwqls b
          cc = .5*v(iq,k) - savv(iq,k) + .5*savv1(iq,k)                               ! simple c
          cc_2 = (10.*v(iq,k)-13.*savv(iq,k)-4.*savv1(iq,k)+7.*savv2(iq,k))/34.       ! cwqls c
          aa = cc_2 - cc
          rat = max( 0., min( 1., cc_2/(aa+sign(1.e-9,aa)) ) )
          cc = rat*cc + (1.-rat)*cc_2 
          bb = rat*bb + (1.-rat)*bb_2 
          vbar(iq,k) = v(iq,k)+.5*bb+.25*cc
        enddo  ! iq loop
      enddo   ! k loop 
    else      ! i.e. mex >=4 and ktau>=3
      ! (tau+.5) from tau, tau-1, tau-2   ! ubar is savu1 here
      ubar(:,:) = (u(1:ifull,:)*15.-savu(:,:)*10.+savu1(:,:)*3.)/8.
      vbar(:,:) = (v(1:ifull,:)*15.-savv(:,:)*10.+savv1(:,:)*3.)/8.
    endif    ! (ktau==1) .. else ..
    if ( mod(ktau,nmaxpr)==0 .and. mydiag ) then
      nlx = max( 2, nlv )  ! as savs not defined for k=1
      write (6,"(i4,' savu2,savu1,savu,u,ubar',5f8.2)") ktau,savu2(idjd,nlv),savu1(idjd,nlv),savu(idjd,nlv), &
                                                        u(idjd,nlv),ubar(idjd,nlv)
      write (6,"(i4,' savv2,savv1,savv,v,vbar',5f8.2)") ktau,savv2(idjd,nlv),savv1(idjd,nlv),savv(idjd,nlv), &
                                                        v(idjd,nlv),vbar(idjd,nlv)
    end if
    if ( ktau>2 .and. epsp>1. .and. epsp<2. ) then
      if ( ktau==3 .and. nmaxpr==1 ) then
        if ( myid==0 ) then
          write(6,*)'using epsp= ',epsp
        end if
        call ccmpi_barrier(comm_world)
      end if
      where ( (sign(1.,dpsdt(1:ifull))/=sign(1.,dpsdtb(1:ifull))) .and. (sign(1.,dpsdtbb(1:ifull))/=sign(1.,dpsdtb(1:ifull))) )
        epst(1:ifull) = epsp - 1.
      elsewhere
        epst(1:ifull) = 0.
      end where
    endif ! (ktau>2.and.epsp>1..and.epsp<2.)

    if ( ktau<10 .and. mydiag ) then
      write(6,*)'savu,u,ubar ',ktau,savu(idjd,1),u(idjd,1),ubar(idjd,1)
    end if
    if ( ktau==1 .and. .not.lrestart .and. mspec==1 .and. mex/=1 ) then
      u(1:ifull,:) = savu(1:ifull,:)  ! reset u,v to original values
      v(1:ifull,:) = savv(1:ifull,:)
    end if
    savu2(1:ifull,:) = savu1(1:ifull,:)  
    savv2(1:ifull,:) = savv1(1:ifull,:)
    savs1(1:ifull,:) = savs(1:ifull,:)  
    savu1(1:ifull,:) = savu(1:ifull,:)  
    savv1(1:ifull,:) = savv(1:ifull,:)
    savs(1:ifull,:)  = sdot(1:ifull,2:kl)  
    savu(1:ifull,:)  = u(1:ifull,:)  ! before any time-splitting occurs
    savv(1:ifull,:)  = v(1:ifull,:)

    ! set diagnostic printout flag
    diag = ( ktau>=abs(ndi) .and. ktau<=ndi2 )
    if ( ndi < 0 ) then
      if ( ktau == (ktau/ndi)*ndi ) then
        diag = .true.
      end if
    endif

    ! update non-linear dynamic terms
    if ( nmaxpr == 1 ) then
      if ( myid == 0 ) then
        write(6,*) "Before nonlin"
      end if
      call ccmpi_barrier(comm_world)
    end if
    call nonlin
    if ( nmaxpr == 1 ) then
      if ( myid == 0 ) then
        write(6,*) "After nonlin"
      end if
      call ccmpi_barrier(comm_world)
    end if
    if ( diag ) then
      if ( mydiag ) write(6,*) 'before hadv'
      call printa('tx  ',tx,ktau,nlv,ia,ib,ja,jb,0.,1.)
      if ( mydiag ) then
        nlx = min( nlv, kl-8 )
        write(6,"('tx  ',9f8.2)") (tx(idjd,k),k=nlx,nlx+8)
        write(6,"('txe ',9f8.2)") (tx(ie(idjd),k),k=nlx,nlx+8)
        write(6,"('txw ',9f8.2)") (tx(iw(idjd),k),k=nlx,nlx+8)
        write(6,"('txn ',9f8.2)") (tx(in(idjd),k),k=nlx,nlx+8)
        write(6,"('txs ',9f8.2)") (tx(is(idjd),k),k=nlx,nlx+8)
        write(6,'(i2," qgv ",18f7.4)')ktau,(1000.*qg(idjd,k),k=1,kl)
      end if
      call printa('qgv ',qg,ktau,nlv,ia,ib,ja,jb,0.,1.e3)
    endif

    ! evaluate horizontal advection for combined quantities
    if ( nmaxpr == 1 ) then
      if ( myid == 0 ) then
        write(6,*) "Before upglobal"
      end if
      call ccmpi_barrier(comm_world)
    end if
    call upglobal
    if ( nmaxpr == 1 ) then
      if ( myid == 0 ) then
        write(6,*) "After upglobal"
      end if
      call ccmpi_barrier(comm_world)
    end if
    if ( diag ) then
      if ( mydiag ) then
        write(6,*) 'after hadv'
        write (6,"('tx  ',9f8.2)") (tx(idjd,k),k=nlx,nlx+8)
      end if
      call printa('tx  ',tx,ktau,nlv,ia,ib,ja,jb,200.,1.)
      if ( mydiag ) write(6,'(i2," qgh ",18f7.4)')ktau,1000.*qg(idjd,:)
    end if

    if ( nstaguin<0 .and. ktau>=1 ) then  ! swapping here (lower down) for nstaguin<0
      if ( nstagin<0 .and. mod(ktau-nstagoff,abs(nstagin))==0 ) then
        nstag  = 7 - nstag  ! swap between 3 & 4
        nstagu = nstag
      end if
    end if
    if ( nmaxpr == 1 ) then
      if ( myid == 0 ) then
        write(6,*) "Before adjust5"
      end if
      call ccmpi_barrier(comm_world)
    end if
    call adjust5
    if ( nmaxpr == 1 ) then
      if ( myid == 0 ) then
        write(6,*) "After adjust5"
      end if
      call ccmpi_barrier(comm_world)
    end if
 
    
    ! NESTING ---------------------------------------------------------------
    ! nesting now after mass fixers
    call START_LOG(nestin_begin)
    if ( nmaxpr == 1 ) then
      if ( myid == 0 ) then
        write(6,*) "Before nesting"
      end if
      call ccmpi_barrier(comm_world)
    end if
    if ( mspec == 1 ) then
      if ( mbd /= 0 ) then
        ! scale-selective filter
        call nestinb
      else if ( nbd /= 0 ) then
        ! Newtonian relaxiation
        call davies
      end if
    end if
    if ( nmaxpr == 1 ) then
      if ( myid == 0 ) then
        write(6,*) "After nesting"
      end if
      call ccmpi_barrier(comm_world)
    end if
    call END_LOG(nestin_end)

    
    ! DYNAMICS --------------------------------------------------------------
    if ( mspec == 2 ) then     ! for very first step restore mass & T fields
      call gettin(1)
    endif    !  (mspec==2) 
    if ( mfix_qg==0 .or. mspec==2 ) then
      qfg(1:ifull,:) = max(qfg(1:ifull,:),0.) 
      qlg(1:ifull,:) = max(qlg(1:ifull,:),0.)
      qg(1:ifull,:)  = max(qg(1:ifull,:),qgmin-qlg(1:ifull,:)-qfg(1:ifull,:),0.) 
    endif  ! (mfix_qg==0.or.mspec==2)

    dt = dtin
  end do ! ****** end of introductory time loop
  mspeca = 1

  
  ! HORIZONTAL DIFFUSION ----------------------------------------------------
  call START_LOG(hordifg_begin)
  if ( nmaxpr == 1 ) then
    if ( myid == 0 ) then
      write(6,*) "Before atm horizontal diffusion"
    end if
    call ccmpi_barrier(comm_world)
  end if
  if ( nhor < 0 ) then
    call hordifgt  ! now not tendencies
  end if
  if ( diag .and. mydiag ) then
    write(6,*) 'after hordifgt t ',t(idjd,:)
  end if
  if ( nmaxpr == 1 ) then
    if ( myid == 0 ) then
      write(6,*) "After atm horizontal diffusion"
    end if
    call ccmpi_barrier(comm_world)
  end if
  call END_LOG(hordifg_end)

  
  ! ***********************************************************************
  ! START OCEAN DYNAMICS
  ! ***********************************************************************

  ! nmlo=0   Prescriped SSTs and sea-ice with JLM skin enhancement
  ! nmlo=1   1D mixed-layer-ocean model
  ! nmlo=2   nmlo=1 plus river-routing and horiontal diffusion
  ! nmlo=3   nmlo=2 plus 3D dynamics
  ! nmlo>9   Use external PCOM ocean model

  if ( abs(nmlo) >= 2 ) then
    ! RIVER ROUTING ------------------------------------------------------
    ! This option can also be used with PCOM
    call START_LOG(river_begin)
    if ( nmaxpr == 1 ) then
      if ( myid == 0 ) then
        write(6,*) "Before river"
      end if
      call ccmpi_barrier(comm_world)
    end if
    call rvrrouter
    if ( nmaxpr == 1 ) then
      if ( myid == 0 ) then
        write(6,*) "After river"
      end if
      call ccmpi_barrier(comm_world)
    end if
    call END_LOG(river_end)
  end if

  
  call START_LOG(waterdynamics_begin)
  if ( abs(nmlo)>=3 .and. abs(nmlo)<=9 ) then
    ! DYNAMICS & DIFFUSION ------------------------------------------------
    if ( nmaxpr == 1 ) then
      if ( myid == 0 ) then
        write(6,*) "Before MLO dynamics"
      end if
      call ccmpi_barrier(comm_world)
    end if
    call mlohadv
    if ( nmaxpr == 1 ) then
      if ( myid == 0 ) then
        write(6,*) "After MLO dynamics"
      end if
      call ccmpi_barrier(comm_world)
    end if
  else if ( abs(nmlo)>=2 .and. abs(nmlo)<=9 ) then
    ! DIFFUSION ONLY ------------------------------------------------------
    if ( nmaxpr == 1 ) then
      if ( myid == 0 ) then
        write(6,*) "Before MLO diffusion"
      end if
      call ccmpi_barrier(comm_world)          
    end if
    call mlodiffusion
    if ( nmaxpr == 1 ) then
      if ( myid == 0 ) then
        write(6,*) "After MLO diffusion"
      end if
      call ccmpi_barrier(comm_world)          
    end if
  end if
  call END_LOG(waterdynamics_end)
      

  ! ***********************************************************************
  ! START PHYSICS 
  ! ***********************************************************************
  call START_LOG(phys_begin)
      
  
  ! GWDRAG ----------------------------------------------------------------
  call START_LOG(gwdrag_begin)
  if ( nmaxpr == 1 ) then
    if ( myid == 0 ) then
      write(6,*) "Before gwdrag"
    end if
    call ccmpi_barrier(comm_world)
  end if
  if ( ngwd < 0 ) call gwdrag  ! <0 for split - only one now allowed
  if ( nmaxpr == 1 ) then
    if ( myid == 0 ) then
      write(6,*) "After gwdrag"
    end if
    call ccmpi_barrier(comm_world)
  end if
  call END_LOG(gwdrag_end)
      
  
  ! CONVECTION ------------------------------------------------------------
  call START_LOG(convection_begin)
  if ( nmaxpr == 1 ) then
    if ( myid == 0 ) then
      write(6,*) "Before convection"
    end if
    call ccmpi_barrier(comm_world)
  end if
  convh_ave = convh_ave - t(1:ifull,:)*real(nperday)/real(nperavg)
  condc     = 0. ! default convective rainfall (assumed to be rain)
  condx     = 0. ! default total precip = rain + ice + snow + graupel (convection and large scale)
  conds     = 0. ! default total ice + snow (convection and large scale)
  condg     = 0. ! default total graupel (convection and large scale)
  ! Save aerosol concentrations for outside convective fraction of grid box
  if ( abs(iaero) >= 2 ) then
    xtosav(:,:,:) = xtg(1:ifull,:,:) ! Aerosol mixing ratio outside convective cloud
  end if
  ! Select convection scheme
  select case ( nkuo )
    case(5)
      call betts(t,qg,tn,land,ps) ! not called these days
    case(23,24)
      call convjlm                ! split convjlm 
    case(46)
      !call conjob                ! split Arakawa-Gordon scheme
      write(6,*) "ERROR: Conjob no longer supported with nkuo=46"
      call ccmpi_abort(-1)
  end select
  cbas_ave(:) = cbas_ave(:) + condc(:)*(1.1-sig(kbsav(:)))      ! diagnostic
  ctop_ave(:) = ctop_ave(:) + condc(:)*(1.1-sig(abs(ktsav(:)))) ! diagnostic
  if ( nmaxpr == 1 ) then
    if ( myid == 0 ) then
      write(6,*) "After convection"
    end if
    call ccmpi_barrier(comm_world)
  end if
  call END_LOG(convection_end)
     
  
  ! CLOUD MICROPHYSICS ----------------------------------------------------
  call START_LOG(cloud_begin)
  if ( nmaxpr == 1 ) then
    if ( myid == 0 ) then
      write(6,*) "Before cloud microphysics"
    end if
    call ccmpi_barrier(comm_world)
  end if
  if ( ldr /= 0 ) then
    ! LDR microphysics scheme
    call leoncld
  end if
  do k = 1,kl
    riwp_ave(1:ifull) = riwp_ave(1:ifull) - qfrad(:,k)*dsig(k)*ps(1:ifull)/grav ! ice water path
    rlwp_ave(1:ifull) = rlwp_ave(1:ifull) - qlrad(:,k)*dsig(k)*ps(1:ifull)/grav ! liq water path
  enddo
  convh_ave = convh_ave + t(1:ifull,:)*real(nperday)/real(nperavg)
  rnd_3hr(1:ifull,8) = rnd_3hr(1:ifull,8) + condx(:)  ! i.e. rnd24(:)=rnd24(:)+condx(:)
#ifdef debug
  if ( nmaxpr==1 .and. mydiag ) then
    write (6,"('qfrad',3p9f8.3/5x,9f8.3)") qfrad(idjd,:)
    write (6,"('qlrad',3p9f8.3/5x,9f8.3)") qlrad(idjd,:)
    write (6,"('qf   ',3p9f8.3/5x,9f8.3)") qfg(idjd,:)
  endif
#endif
  if ( nmaxpr == 1 ) then
    if ( myid == 0 ) then
      write(6,*) "After cloud microphysics"
    end if
    call ccmpi_barrier(comm_world)
  end if
  call END_LOG(cloud_end)

  
  ! RADIATION -------------------------------------------------------------
      
  ! nrad=4 Fels-Schwarzkopf radiation
  ! nrad=5 SEA-ESF radiation

  call START_LOG(radnet_begin)
  if ( nmaxpr == 1 ) then
    if ( myid == 0 ) then
      write(6,*) "Before radiation"
    end if
    call ccmpi_barrier(comm_world)
  end if
  if ( ncloud >= 4 ) then
    nettend = nettend + t(1:ifull,:)/dt
  end if
  odcalc = mod(ktau,kountr)==0 .or. ktau==1 ! ktau-1 better
  if ( nhstest < 0 ) then ! aquaplanet test -1 to -8  
    mtimer_sav = mtimer
    mtimer     = mins_gmt    ! so radn scheme repeatedly works thru same day
  end if    ! (nhstest<0)
  select case ( nrad )
    case(4)
      ! Fels-Schwarzkopf radiation
      call radrive(il*nrows_rad,odcalc)
    case(5)
      ! GFDL SEA-EFS radiation
      call seaesfrad(il*nrows_rad,odcalc)
    case DEFAULT
      ! use preset slwa array (use +ve nrad)
      slwa(:) = -10*nrad  
  end select
  if ( nhstest<0 ) then ! aquaplanet test -1 to -8  
    mtimer = mtimer_sav
  end if    ! (nhstest<0)
  if ( nmaxpr == 1 ) then
    call maxmin(slwa,'sl',ktau,.1,1)
    if ( myid == 0 ) then
      write(6,*) "After radiation"
    end if
    call ccmpi_barrier(comm_world)
  end if
  call END_LOG(radnet_end)


  ! HELD & SUAREZ ---------------------------------------------------------
  if ( ntsur<=1 .or. nhstest==2 ) then ! Held & Suarez or no surf fluxes
    eg(:)   = 0.
    fg(:)   = 0.
    cdtq(:) = 0.
    cduv(:) = 0.
  end if     ! (ntsur<=1.or.nhstest==2) 
  if ( nhstest == 2 ) then
    call hs_phys
  end if

  
  ! SURFACE FLUXES ---------------------------------------------
  ! (Includes ocean dynamics and mixing, as well as ice dynamics and thermodynamics)
  call START_LOG(sfluxnet_begin)
  if ( nmaxpr == 1 ) then
    if ( myid == 0 ) then
      write(6,*) "Before surface fluxes"
    end if
    call ccmpi_barrier(comm_world)
  end if
  if ( diag ) then
    call maxmin(u,'#u',ktau,1.,kl)
    call maxmin(v,'#v',ktau,1.,kl)
    call maxmin(t,'#t',ktau,1.,kl)
    call maxmin(qg,'qg',ktau,1.e3,kl)     
    call ccmpi_barrier(comm_world) ! stop others going past
  end if
  if ( ntsur > 1 ) then  ! should be better after convjlm
    call sflux(nalpha)
  endif   ! (ntsur>1)    
  if ( nmaxpr == 1 ) then
    if ( myid == 0 ) then
      write(6,*) "After surface fluxes"
    end if
    call ccmpi_barrier(comm_world)
  end if
  call END_LOG(sfluxnet_end)
  

  ! AEROSOLS --------------------------------------------------------------
  ! MJT notes - aerosols called before vertical mixing so that convective
  ! and strat cloud can be separated consistently with cloud microphysics
  if ( abs(iaero) >= 2 ) then
    call START_LOG(aerosol_begin)
    if ( nmaxpr == 1 ) then
      if ( myid == 0 ) then
        write(6,*) "Before aerosols"
      end if
      call ccmpi_barrier(comm_world)
    end if
    call aerocalc
    if ( nmaxpr == 1 ) then
      if ( myid == 0 ) then
        write(6,*) "After aerosols"
      end if
      call ccmpi_barrier(comm_world)
    end if
    call END_LOG(aerosol_end)
  end if

  
  ! VERTICAL MIXING ------------------------------------------------------
  if ( nmaxpr == 1 ) then
    if ( myid == 0 ) then
      write(6,*) "Before PBL mixing"
    end if
    call ccmpi_barrier(comm_world)
  end if
  if ( ntsur >= 1 ) then ! calls vertmix but not sflux for ntsur=1
    call START_LOG(vertmix_begin)
    if ( nmaxpr==1 .and. mydiag ) then
      write (6,"('pre-vertmix t',9f8.3/13x,9f8.3)") t(idjd,:)
    end if
    call vertmix
    if ( nmaxpr==1 .and. mydiag ) then
      write (6,"('aft-vertmix t',9f8.3/13x,9f8.3)") t(idjd,:)
    end if
    call END_LOG(vertmix_end)
  endif  ! (ntsur>=1)
  if ( ncloud >= 4 ) then
    nettend = (nettend-t(1:ifull,:)/dt)
  end if
  if ( nmaxpr == 1 ) then
    if ( myid == 0 ) then
      write(6,*) "After PBL mixing"
    end if
    call ccmpi_barrier(comm_world)
  end if
      
  
  ! Update diagnostics for consistancy in history file
  if ( rescrn > 0 ) then
    call autoscrn
  end if

  
  ! PHYSICS LOAD BALANCING ------------------------------------------------
  ! This is the end of the physics. The next routine makes the load imbalance
  ! overhead explicit rather than having it hidden in one of the diagnostic
  ! calls.
#ifdef loadbal
  call phys_loadbal
#endif
  call END_LOG(phys_end)

  
  ! ***********************************************************************
  ! TRACER OUTPUT
  ! ***********************************************************************
  ! rml 16/02/06 call tracer_mass, write_ts
  if ( ngas > 0 ) then
    call tracer_mass !also updates average tracer array
    call write_ts(ktau,ntau,dt)
  endif

  
  ! ***********************************************************************
  ! DIAGNOSTICS AND OUTPUT
  ! ***********************************************************************

  ! STATION OUTPUT ---------------------------------------------
  if ( nstn > 0 ) then
    call stationa ! write every time step
  end if
       
  ! DIAGNOSTICS ------------------------------------------------
  if ( mod(ktau,nmaxpr)==0 .and. mydiag ) then
    write(6,*)
    write (6,"('ktau =',i5,' gmt(h,m):',f6.2,i5,' runtime(h,m):',f7.2,i6)") ktau,timeg,mins_gmt,timer,mtimer
    ! some surface (or point) diagnostics
    isoil = isoilm(idjd)
    write(6,*) 'land,isoil,ivegt,isflag ',land(idjd),isoil,ivegt(idjd),isflag(idjd)
    write (6,"('snage,snowd,alb   ',f8.4,2f8.2)") snage(idjd),snowd(idjd),albvisnir(idjd,1)
    write (6,"('sicedep,fracice,runoff ',3f8.2)") sicedep(idjd),fracice(idjd),runoff(idjd)
    write (6,"('tgg(1-6)   ',9f8.2)") (tgg(idjd,k),k=1,6)
    write (6,"('tggsn(1-3) ',9f8.2)") (tggsn(idjd,k),k=1,3)
    write (6,"('wb(1-6)    ',9f8.3)") (wb(idjd,k),k=1,6)
    write (6,"('wbice(1-6) ',9f8.3)") (wbice(idjd,k),k=1,6)
    write (6,"('smass(1-3) ',9f8.2)") (smass(idjd,k),k=1,3) ! as mm of water
    write (6,"('ssdn(1-3)  ',9f8.2)") (ssdn(idjd,k),k=1,3)
    iq = idjd
    pwater = 0.   ! in mm
    do k = 1,kl
      qtot   = qg(iq,k)+qlg(iq,k)+qfg(iq,k)
      pwater = pwater-dsig(k)*qtot*ps(iq)/grav
    enddo
    write (6,"('pwater,condc,condx,rndmax,rmc',9f8.3)") pwater,condc(idjd),condx(idjd),rndmax(idjd),cansto(idjd)
    write (6,"('wetfac,sno,evap,precc,precip',6f8.2)") wetfac(idjd),sno(idjd),evap(idjd),precc(idjd),precip(idjd)
    write (6,"('tmin,tmax,tscr,tss,tpan',9f8.2)") tminscr(idjd),tmaxscr(idjd),tscrn(idjd),tss(idjd),tpan(idjd)
    write (6,"('u10,ustar,pblh',9f8.2)") u10(idjd),ustar(idjd),pblh(idjd)
    write (6,"('ps,qgscrn',5f8.2,f8.3)") .01*ps(idjd),1000.*qgscrn(idjd)
    write (6,"('dew_,eg_,epot,epan,eg,fg,ga',9f8.2)") dew_ave(idjd),eg_ave(idjd),epot(idjd),epan(idjd),eg(idjd),fg(idjd),ga(idjd)
    write (6,"('zo,cduv',2f8.5)") zo(idjd),cduv(idjd)/vmod(idjd)
    write (6,"('slwa,sint,sg,rt,rg    ',9f8.2)") slwa(idjd),sintsave(idjd),sgsave(idjd),rtsave(idjd),rgsave(idjd)
    write (6,"('cll,clm,clh,clt ',9f8.2)") cloudlo(idjd),cloudmi(idjd),cloudhi(idjd),cloudtot(idjd)
    write (6,"('u10max,v10max,rhmin,rhmax   ',9f8.2)") u10max(iq),v10max(iq),rhminscr(iq),rhmaxscr(iq)
    write (6,"('kbsav,ktsav,convpsav ',2i3,f8.4,9f8.2)") kbsav(idjd),ktsav(idjd),convpsav(idjd)
    write (6,"('t   ',9f8.3/4x,9f8.3)") t(idjd,:)
    write (6,"('u   ',9f8.3/4x,9f8.3)") u(idjd,:)
    write (6,"('v   ',9f8.3/4x,9f8.3)") v(idjd,:)
    write (6,"('qg  ',9f8.3/4x,9f8.3)") qg(idjd,:)
    write (6,"('qf  ',9f8.3/4x,9f8.3)") qfg(idjd,:)
    write (6,"('ql  ',9f8.3/4x,9f8.3)") qlg(idjd,:)
    write (6,"('cfrac',9f8.3/5x,9f8.3)") cfrac(idjd,:)
    do k = 1,kl
      es        = establ(t(idjd,k))
      spmean(k) = 100.*qg(idjd,k)*max(ps(idjd)*sig(k)-es,1.)/(.622*es) ! max as for convjlm
    enddo
    write (6,"('rh  ',9f8.3/4x,9f8.3)") spmean(:)
    write (6,"('omgf ',9f8.3/5x,9f8.3)") ps(idjd)*dpsldt(idjd,:) ! in Pa/s
    write (6,"('sdot ',9f8.3/5x,9f8.3)") sdot(idjd,1:kl)
    if ( nextout >= 4 ) then
      write (6,"('xlat,long,pres ',3f8.2)") tr(idjd,nlv,ngas+1),tr(idjd,nlv,ngas+2),tr(idjd,nlv,ngas+3)
    end if
  endif  ! (mod(ktau,nmaxpr)==0.and.mydiag)
  
  if ( ndi == -ktau ) then
    nmaxpr = 1         ! diagnostic prints; reset 6 lines on
    if ( ndi2 == 0 ) ndi2 = ktau + 40
  endif
  if ( ktau == ndi2 ) then
    if ( myid == 0 ) write(6,*)'reset nmaxpr'
    nmaxpr = nmaxprsav
  endif
  if ( mod(ktau,nmaxpr)==0 .or. ktau==ntau ) then
    call maxmin(u,' u',ktau,1.,kl)
    call maxmin(v,' v',ktau,1.,kl)
    dums(:,:) = u(1:ifull,:)**2 + v(1:ifull,:)**2 ! 3D
    call average(dums,spmean,spavge)
    do k = 1,kl
      spmean(k) = sqrt(spmean(k))
    enddo
    dums(1:ifull,1:kl) = sqrt(dums(1:ifull,1:kl)) ! 3D
    spavge = sqrt(spavge)
    call maxmin(dums,'sp',ktau,1.,kl)
    call maxmin(t,' t',ktau,1.,kl)
    call maxmin(qg,'qg',ktau,1.e3,kl)
    call maxmin(qfg,'qf',ktau,1.e3,kl)
    call maxmin(qlg,'ql',ktau,1.e3,kl)
    call maxmin(sdot,'sd',ktau,1.,kl)  ! grid length units 
    if ( myid == 0 ) then
      write(6,'("spmean ",9f8.3)') spmean
      write(6,'("spavge ",f8.3)') spavge
    end if
    dums = qg(1:ifull,:)
    call average(dums,spmean,spavge)
    if ( myid == 0 ) then
      write(6,'("qgmean ",9f8.5)') spmean
      write(6,'("qgavge ",f8.5)') spavge
    end if
    call maxmin(wb,'wb',ktau,1.,ms)
    call maxmin(tggsn,'tggsn',ktau,1.,3)
    call maxmin(tgg,'tg',ktau,1.,ms)
    call maxmin(tss,'ts',ktau,1.,1)
    call maxmin(pblh,'pb',ktau,1.,1)
    call maxmin(precip,'pr',ktau,1.,1)
    call maxmin(precc,'pc',ktau,1.,1)
    call maxmin(convpsav,'co',ktau,1.,1)
    call maxmin(sno,'sn',ktau,1.,1)        ! as mm during timestep
    call maxmin(rhscrn,'rh',ktau,1.,1)
    call maxmin(ps,'ps',ktau,.01,1)
    psavge    = sum(ps(1:ifull)*wts(1:ifull))
    pslavge   = sum(psl(1:ifull)*wts(1:ifull))
    preccavge = sum(precc(1:ifull)*wts(1:ifull))
    precavge  = sum(precip(1:ifull)*wts(1:ifull))
    ! KE calculation, not taking into account pressure weighting
    gke = 0.
    do k = 1,kl
      gke = gke - sum( 0.5 * wts(1:ifull) * dsig(k) * ( u(1:ifull,k)**2 + v(1:ifull,k)**2 ) )
    end do
    cllav = sum(wts(1:ifull)*cloudlo(1:ifull))
    clmav = sum(wts(1:ifull)*cloudmi(1:ifull))
    clhav = sum(wts(1:ifull)*cloudhi(1:ifull))
    cltav = sum(wts(1:ifull)*cloudtot(1:ifull))

    ! All this combined into a single reduction
    temparray = (/ psavge, pslavge, preccavge, precavge, gke, cllav, clmav,clhav, cltav /)
    call ccmpi_reduce(temparray,gtemparray,"sum",0,comm_world)
    if ( myid == 0 ) then
      write(6,97) gtemparray(1:5) ! psavge,pslavge,preccavge,precavge,gke
97    format(' average ps, psl, precc, prec, gke: ',f10.2,f10.6,2f6.2,f7.2)
      write(6,971) gtemparray(6:9) ! cllav,clmav,clhav,cltav
971   format(' global_average cll, clm, clh, clt: ',4f6.2)
    end if
    if ( mydiag ) then
      write(6,98) ktau,diagvals(ps)
98    format(i7,' ps diag:',9f9.1)
      if ( t(idjd,kl) > 258. ) then
        write(6,*) 't(idjd,kl) > 258. for idjd = ',idjd
        write(6,91) ktau,(t(idjd,k),k=kl-8,kl)
91      format(i7,'    t',9f7.2)
        write(6,92) ktau,(sdot(idjd,k),k=kl-8,kl)
92      format(i7,' sdot',9f7.3)
      end if             ! (t(idjd,kl)>258.)
    end if               ! myid==0
  endif                  ! (mod(ktau,nmaxpr)==0)

  ! update diag_averages and daily max and min screen temps 
  ! N.B. runoff is accumulated in sflux
  tmaxscr(1:ifull)     = max( tmaxscr(1:ifull), tscrn )
  tminscr(1:ifull)     = min( tminscr(1:ifull), tscrn )
  rhmaxscr(1:ifull)    = max( rhmaxscr(1:ifull), rhscrn )
  rhminscr(1:ifull)    = min( rhminscr(1:ifull), rhscrn )
  rndmax(1:ifull)      = max( rndmax(1:ifull), condx )
  cape_max(1:ifull)    = max( cape_max(1:ifull), cape )
  cape_ave(1:ifull)    = cape_ave(1:ifull) + cape
  u10mx(1:ifull)       = max( u10mx(1:ifull), u10 )  ! for hourly scrnfile
  dew_ave(1:ifull)     = dew_ave(1:ifull) - min( 0., eg )    
  epan_ave(1:ifull)    = epan_ave(1:ifull) + epan
  epot_ave(1:ifull)    = epot_ave(1:ifull) + epot 
  eg_ave(1:ifull)      = eg_ave(1:ifull) + eg    
  fg_ave(1:ifull)      = fg_ave(1:ifull) + fg
  ga_ave(1:ifull)      = ga_ave(1:ifull) + ga
  rnet_ave(1:ifull)    = rnet_ave(1:ifull) + rnet
  tscr_ave(1:ifull)    = tscr_ave(1:ifull) + tscrn 
  qscrn_ave(1:ifull)   = qscrn_ave(1:ifull) + qgscrn 
  wb_ave(1:ifull,1:ms) = wb_ave(1:ifull,1:ms) + wb
  tsu_ave(1:ifull)     = tsu_ave(1:ifull) + tss
  call mslp(spare2,psl,zs,t) ! calculate MSLP from psl
  spare2 = spare2/100.       ! convert MSLP to hPa
  psl_ave(1:ifull)     = psl_ave(1:ifull) + spare2(1:ifull)
  spare1(1:ifull)      = 0.
  call mlodiag(spare1,0)     ! obtain ocean mixed level depth
  mixdep_ave(1:ifull)  = mixdep_ave(1:ifull) + spare1(1:ifull)
  spare1(:) = u(1:ifull,1)**2 + v(1:ifull,1)**2
  spare2(:) = u(1:ifull,2)**2 + v(1:ifull,2)**2
  do iq = 1,ifull
    if ( u10(iq)**2 > u10max(iq)**2 +v10max(iq)**2 ) then
      u10max(iq) = u10(iq)*u(iq,1)/max(.001,sqrt(spare1(iq)))
      v10max(iq) = u10(iq)*v(iq,1)/max(.001,sqrt(spare1(iq)))
    end if
    if ( spare1(iq) > u1max(iq)**2+v1max(iq)**2 ) then
      u1max(iq) = u(iq,1)
      v1max(iq) = v(iq,1)
    end if
    if ( spare2(iq) > u2max(iq)**2+v2max(iq)**2 ) then
      u2max(iq) = u(iq,2)
      v2max(iq) = v(iq,2)
    end if
  end do
  if ( ngas > 0 ) then
    traver(:,:,1:ngas) = traver(:,:,1:ngas) + tr(1:ilt*jlt,:,1:ngas)
  end if
  fpn_ave(1:ifull) = fpn_ave(1:ifull) + fpn
  frs_ave(1:ifull) = frs_ave(1:ifull) + frs
  frp_ave(1:ifull) = frp_ave(1:ifull) + frp

  ! rnd03 to rnd21 are accumulated in mm     
  if ( myid == 0 ) then
    write(6,*) 'ktau,mod,nper3hr ',ktau,mod(ktau-1,nperday)+1,nper3hr(n3hr)
  end if
  if ( mod(ktau-1,nperday)+1 == nper3hr(n3hr) ) then
    rnd_3hr(1:ifull,n3hr) = rnd_3hr(1:ifull,8)
    if ( nextout >= 2 ) then
      spare1(:) = max( .001, sqrt(u(1:ifull,1)**2+v(1:ifull,1)**2) )
      u10_3hr(:,n3hr) = u10(:)*u(1:ifull,1)/spare1(:)
      v10_3hr(:,n3hr) = u10(:)*v(1:ifull,1)/spare1(:)
      tscr_3hr(:,n3hr) = tscrn(:)
      spare1(:) = establ(t(1:ifull,1)) ! spare1 = es
      rh1_3hr(1:ifull,n3hr) = 100.*qg(1:ifull,1)*(ps(1:ifull)*sig(1)-spare1(:))/(.622*spare1(:))
    end if    ! (nextout==2)
    n3hr = n3hr+1
    if ( n3hr > 8 ) n3hr = 1
  endif    ! (mod(ktau,nperday)==nper3hr(n3hr))

  if ( ktau==ntau .or. mod(ktau,nperavg)==0 ) then
    cape_ave(1:ifull)   = cape_ave(1:ifull)/min(ntau,nperavg)
    dew_ave(1:ifull)    = dew_ave(1:ifull)/min(ntau,nperavg)
    epan_ave(1:ifull)   = epan_ave(1:ifull)/min(ntau,nperavg)
    epot_ave(1:ifull)   = epot_ave(1:ifull)/min(ntau,nperavg)
    eg_ave(1:ifull)     = eg_ave(1:ifull)/min(ntau,nperavg)
    fg_ave(1:ifull)     = fg_ave(1:ifull)/min(ntau,nperavg)
    ga_ave(1:ifull)     = ga_ave(1:ifull)/min(ntau,nperavg)    
    rnet_ave(1:ifull)   = rnet_ave(1:ifull)/min(ntau,nperavg)
    sunhours(1:ifull)   = sunhours(1:ifull)/min(ntau,nperavg)
    riwp_ave(1:ifull)   = riwp_ave(1:ifull)/min(ntau,nperavg)
    rlwp_ave(1:ifull)   = rlwp_ave(1:ifull)/min(ntau,nperavg)
    tscr_ave(1:ifull)   = tscr_ave(1:ifull)/min(ntau,nperavg)
    qscrn_ave(1:ifull)  = qscrn_ave(1:ifull)/min(ntau,nperavg)
    do k=1,ms
      wb_ave(1:ifull,k) = wb_ave(1:ifull,k)/min(ntau,nperavg)
    end do
    tsu_ave(1:ifull)    = tsu_ave(1:ifull)/min(ntau,nperavg)
    psl_ave(1:ifull)    = psl_ave(1:ifull)/min(ntau,nperavg)
    mixdep_ave(1:ifull) = mixdep_ave(1:ifull)/min(ntau,nperavg)
    sgn_ave(1:ifull)    = sgn_ave(1:ifull)/min(ntau,nperavg)  ! Dec07 because of solar fit
    sint_ave(1:ifull)   = sint_ave(1:ifull)/max(koundiag,1)
    sot_ave(1:ifull)    = sot_ave(1:ifull)/max(koundiag,1)
    soc_ave(1:ifull)    = soc_ave(1:ifull)/max(koundiag,1)
    sgdn_ave(1:ifull)   = sgdn_ave(1:ifull)/max(koundiag,1)
    rtu_ave(1:ifull)    = rtu_ave(1:ifull)/max(koundiag,1)
    rtc_ave(1:ifull)    = rtc_ave(1:ifull)/max(koundiag,1)
    rgdn_ave(1:ifull)   = rgdn_ave(1:ifull)/max(koundiag,1)
    rgn_ave(1:ifull)    = rgn_ave(1:ifull)/max(koundiag,1)
    rgc_ave(1:ifull)    = rgc_ave(1:ifull)/max(koundiag,1)
    sgc_ave(1:ifull)    = sgc_ave(1:ifull)/max(koundiag,1)
    cld_ave(1:ifull)    = cld_ave(1:ifull)/max(koundiag,1)
    cll_ave(1:ifull)    = cll_ave(1:ifull)/max(koundiag,1)
    clm_ave(1:ifull)    = clm_ave(1:ifull)/max(koundiag,1)
    clh_ave(1:ifull)    = clh_ave(1:ifull)/max(koundiag,1)
    alb_ave(1:ifull)    = alb_ave(1:ifull)/max(koundiag,1)
    fbeam_ave(1:ifull)  = fbeam_ave(1:ifull)/max(koundiag,1)
    cbas_ave(1:ifull)   = 1.1 - cbas_ave(1:ifull)/max(1.e-4,precc(:))  ! 1.1 for no precc
    ctop_ave(1:ifull)   = 1.1 - ctop_ave(1:ifull)/max(1.e-4,precc(:))  ! 1.1 for no precc
    if ( ngas > 0 ) then
      traver(1:ifull,1:kl,1:ngas) = traver(1:ifull,1:kl,1:ngas)/min(ntau,nperavg)
    end if
    fpn_ave(1:ifull)    = fpn_ave(1:ifull)/min(ntau,nperavg)
    frs_ave(1:ifull)    = frs_ave(1:ifull)/min(ntau,nperavg)
    frp_ave(1:ifull)    = frp_ave(1:ifull)/min(ntau,nperavg)
    if ( abs(iaero) == 2 ) then
      duste        = duste/min(ntau,nperavg)       ! Dust emissions
      dustdd       = dustdd/min(ntau,nperavg)      ! Dust dry deposition
      dustwd       = dustwd/min(ntau,nperavg)      ! Dust wet deposition
      dust_burden  = dust_burden/min(ntau,nperavg) ! Dust burden
      bce          = bce/min(ntau,nperavg)         ! Black carbon emissions
      bcdd         = bcdd/min(ntau,nperavg)        ! Black carbon dry deposition
      bcwd         = bcwd/min(ntau,nperavg)        ! Black carbon wet deposition
      bc_burden    = bc_burden/min(ntau,nperavg)   ! Black carbon burden
      oce          = oce/min(ntau,nperavg)         ! Organic carbon emissions
      ocdd         = ocdd/min(ntau,nperavg)        ! Organic carbon dry deposition
      ocwd         = ocwd/min(ntau,nperavg)        ! Organic carbon wet deposition
      oc_burden    = oc_burden/min(ntau,nperavg)   ! Organic carbon burden
      dmse         = dmse/min(ntau,nperavg)        ! DMS emissions
      dmsso2o      = dmsso2o/min(ntau,nperavg)     ! DMS -> SO2 oxidation
      so2e         = so2e/min(ntau,nperavg)        ! SO2 emissions
      so2so4o      = so2so4o/min(ntau,nperavg)     ! SO2 -> SO4 oxidation
      so2dd        = so2dd/min(ntau,nperavg)       ! SO2 dry deposition
      so2wd        = so2wd/min(ntau,nperavg)       ! SO2 wet deposiion
      so4e         = so4e/min(ntau,nperavg)        ! SO4 emissions
      so4dd        = so4dd/min(ntau,nperavg)       ! SO4 dry deposition
      so4wd        = so4wd/min(ntau,nperavg)       ! SO4 wet deposition
      dms_burden   = dms_burden/min(ntau,nperavg)  ! DMS burden
      so2_burden   = so2_burden/min(ntau,nperavg)  ! SO2 burden
      so4_burden   = so4_burden/min(ntau,nperavg)  ! SO4 burden
    end if
  end if    ! (ktau==ntau.or.mod(ktau,nperavg)==0)

  call log_off()
  if ( ktau==ntau .or. mod(ktau,nwt)==0 ) then
    call outfile(20,rundate,nwrite,nstagin,jalbfix,nalpha,mins_rad)  ! which calls outcdf

    if ( ktau==ntau .and. irest==1 ) then
      ! Don't include the time for writing the restart file
      call END_LOG(maincalc_end)
      ! write restart file
      call outfile(19,rundate,nwrite,nstagin,jalbfix,nalpha,mins_rad)
      if ( myid == 0 ) then
        write(6,*)'finished writing restart file in outfile'
      end if
      call START_LOG(maincalc_begin)
    endif  ! (ktau==ntau.and.irest==1)
  endif    ! (ktau==ntau.or.mod(ktau,nwt)==0)
      
  ! write high temporal frequency fields
  if ( surfile /= ' ' ) then
    call freqfile
  end if
  call log_on()
 
  if ( mod(ktau,nperavg) == 0 ) then   
    ! produce some diags & reset most averages once every nperavg
    if ( nmaxpr == 1 ) then
      precavge = sum(precip(1:ifull)*wts(1:ifull))
      evapavge = sum(evap(1:ifull)*wts(1:ifull))   ! in mm/day
      pwatr    = 0.   ! in mm
      do k = 1,kl
        pwatr = pwatr - sum(dsig(k)*wts(1:ifull)*(qg(1:ifull,k)+qlg(1:ifull,k)+qfg(1:ifull,k))*ps(1:ifull))/grav
      enddo
      temparray(1:3) = (/ precavge, evapavge, pwatr /)
      call ccmpi_reduce(temparray,gtemparray,"max",0,comm_world)
      if ( myid == 0 ) then
        precavge = gtemparray(1)
        evapavge = gtemparray(2)
        pwatr    = gtemparray(3)
        write(6,985) pwatr,precavge,evapavge ! MJT bug fix
985     format(' average pwatr,precc,prec,evap: ',4f7.3)
      end if
    end if
    ! also zero most averaged fields every nperavg
    convh_ave(:,:) = 0.
    cbas_ave(:)    = 0.
    ctop_ave(:)    = 0.
    dew_ave(:)     = 0.
    epan_ave(:)    = 0.
    epot_ave(:)    = 0.
    eg_ave(:)      = 0.
    fg_ave(:)      = 0.
    rnet_ave(:)    = 0.
    sunhours(:)    = 0.
    riwp_ave(:)    = 0.
    rlwp_ave(:)    = 0.
    qscrn_ave(:)   = 0.
    tscr_ave(:)    = 0.
    wb_ave(:,:)    = 0.
    tsu_ave(:)     = 0.
    alb_ave(:)     = 0.
    fbeam_ave(:)   = 0.
    psl_ave(:)     = 0.
    mixdep_ave(:)  = 0.
    koundiag       = 0
    sint_ave(:)    = 0.
    sot_ave(:)     = 0.
    soc_ave(:)     = 0.
    sgdn_ave(:)    = 0.
    sgn_ave(:)     = 0.
    rtu_ave(:)     = 0.
    rtc_ave(:)     = 0.
    rgdn_ave(:)    = 0.
    rgn_ave(:)     = 0.
    rgc_ave(:)     = 0.
    sgc_ave(:)     = 0.
    cld_ave(:)     = 0.
    cll_ave(:)     = 0.
    clm_ave(:)     = 0.
    clh_ave(:)     = 0.
    ! zero evap, precip, precc, sno, runoff fields each nperavg (3/12/04) 
    evap(:)        = 0.  
    precip(:)      = 0.  ! converted to mm/day in outcdf
    precc(:)       = 0.  ! converted to mm/day in outcdf
    sno(:)         = 0.  ! converted to mm/day in outcdf
    hail(:)        = 0.  ! converted to mm/day in outcdf
    runoff(:)      = 0.  ! converted to mm/day in outcdf
    u10mx(:)       = 0.
    cape_max(:)    = 0.
    cape_ave(:)    = 0.
    if ( ngas > 0 ) then
      traver = 0.
    end if
    fpn_ave = 0.
    frs_ave = 0.
    frp_ave = 0.
    if ( abs(iaero) == 2 ) then
      duste        = 0.  ! Dust emissions
      dustdd       = 0.  ! Dust dry deposition
      dustwd       = 0.  ! Dust wet deposition
      dust_burden  = 0.  ! Dust burden
      bce          = 0.  ! Black carbon emissions
      bcdd         = 0.  ! Black carbon dry deposition
      bcwd         = 0.  ! Black carbon wet deposition
      bc_burden    = 0.  ! Black carbon burden
      oce          = 0.  ! Organic carbon emissions
      ocdd         = 0.  ! Organic carbon dry deposition
      ocwd         = 0.  ! Organic carbon wet deposition
      oc_burden    = 0.  ! Organic carbon burden
      dmse         = 0.  ! DMS emissions
      dmsso2o      = 0.  ! DMS -> SO2 oxidation
      so2e         = 0.  ! SO2 emissions
      so2so4o      = 0.  ! SO2 -> SO4 oxidation
      so2dd        = 0.  ! SO2 dry deposition
      so2wd        = 0.  ! SO2 wet deposiion
      so4e         = 0.  ! SO4 emissions
      so4dd        = 0.  ! SO4 dry deposition
      so4wd        = 0.  ! SO4 wet deposition
      dms_burden   = 0.  ! DMS burden
      so2_burden   = 0.  ! SO2 burden
      so4_burden   = 0.  ! SO4 burden
    end if
  endif  ! (mod(ktau,nperavg)==0)

  if ( mod(ktau,nperday) == 0 ) then   ! re-set at the end of each 24 hours
    if ( ntau<10*nperday .and. nstn>0 ) then     ! print stn info
      do nn = 1,nstn
        if ( .not.mystn(nn) ) cycle
        i = istn(nn)
        j = jstn(nn)
        iq = i+(j-1)*il
        write(6,956) ktau,iunp(nn),name_stn(nn),rnd_3hr(iq,4),rnd_3hr(iq,8), &
                     tmaxscr(iq)-273.16+(zs(iq)/grav-zstn(nn))*stdlapse,     &
                     tminscr(iq)-273.16+(zs(iq)/grav-zstn(nn))*stdlapse,     &
                     tmaxscr(iq)-273.16,tminscr(iq)-273.16
956     format(i5,i3,a5,6f7.1)
      end do
    end if  ! (ntau<10*nperday)
    rndmax (:)  = 0.
    tmaxscr(:)  = tscrn(:) 
    tminscr(:)  = tscrn(:) 
    rhmaxscr(:) = rhscrn(:) 
    rhminscr(:) = rhscrn(:) 
    u10max(:)   = 0.
    v10max(:)   = 0.
    u1max(:)    = 0.
    v1max(:)    = 0.
    u2max(:)    = 0.
    v2max(:)    = 0.
    rnd_3hr(:,8)= 0.       ! i.e. rnd24(:)=0.
    if ( nextout >= 4 ) then
      call setllp ! from Nov 11, reset once per day
    end if
  endif   ! (mod(ktau,nperday)==0)
  
  if ( namip /= 0 ) then
    if ( nmlo == 0 ) then
      if ( mod(ktau,nperday) == 0 ) then
        if ( myid == 0 ) then
          write(6,*) 'amipsst called at end of day for ktau,mtimer,namip ',ktau,mtimer,namip  
        end if
        call amipsst
      end if
    else
      ! call evey time-step for nudging
      call amipsst
    end if
  end if

#ifdef vampir
  ! Flush vampir trace information to disk to save memory.
  VT_BUFFER_FLUSH()
#endif

end do                  ! *** end of main time loop
call END_LOG(maincalc_end)
call log_off()

! Report timings of run
if ( myid == 0 ) then
  call date_and_time(time=timeval,values=tvals2)
  write(6,*) "End of time loop ", timeval
  write(6,*) "normal termination of run"
  call date_and_time(time=timeval)
  write(6,*) "End time ", timeval
  aa = 3600.*(tvals2(5)-tvals1(5)) + 60.*(tvals2(6)-tvals1(6)) + (tvals2(7)-tvals1(7)) + 0.001*(tvals2(8)-tvals1(8))
  if ( aa <= 0. ) aa = aa + 86400.
  write(6,*) "Model time in main loop",aa
end if
call END_LOG(model_end)

! close mesonest files
if ( mbd/=0 .or. nbd/=0 ) then
  call histclose
end if

#ifdef simple_timer
! report subroutine timings
call simple_timer_finalize
#endif

! finalize MPI comms
call ccmpi_finalize

end


!--------------------------------------------------------------
! PREPARE SPECIAL TRACER ARRAYS
! sets tr arrays for lat, long, pressure if nextout>=4 &(nllp>=3)
subroutine setllp
      
use arrays_m           ! Atmosphere dyamics prognostic arrays
use cc_mpi             ! CC MPI routines
use latlong_m          ! Lat/lon coordinates
use sigs_m             ! Atmosphere sigma levels
use tracers_m          ! Tracer data
      
implicit none
      
include 'newmpar.h'    ! Grid parameters
include 'const_phys.h' ! Physical constants
      
integer k
      
if ( nllp < 3 ) then
  write(6,*) "ERROR: Incorrect setting of nllp",nllp
  call ccmpi_abort(-1)
end if
      
do k = 1,klt
  tr(1:ifull,k,ngas+1) = rlatt(1:ifull)*180./pi
  tr(1:ifull,k,ngas+2) = rlongg(1:ifull)*180./pi
  tr(1:ifull,k,ngas+3) = .01*ps(1:ifull)*sig(k)  ! in HPa
enddo
if ( nllp >= 4 ) then   ! theta
  do k = 1,klt
    tr(1:ifull,k,ngas+4) = t(1:ifull,k)*(1.e-5*ps(1:ifull)*sig(k))**(-rdry/cp)
  enddo
endif   ! (nllp>=4)
if ( nllp >= 5 ) then   ! mixing_ratio (g/kg)
  do k = 1,klt
    tr(1:ifull,k,ngas+5) = 1000.*qg(1:ifull,k)
  enddo
endif   ! (nllp>=5)
      
return
end subroutine setllp


!--------------------------------------------------------------
! INTIAL PARAMETERS
blockdata main_blockdata

implicit none

include 'newmpar.h'          ! Grid parameters
include 'dates.h'            ! Date data
include 'filnames.h'         ! Filenames
include 'kuocom.h'           ! Convection parameters
include 'parm.h'             ! Model configuration
include 'parmdyn.h'          ! Dynamics parmaters
include 'parmgeom.h'         ! Coordinate data
include 'parmhor.h'          ! Horizontal advection parameters
include 'parmsurf.h'         ! Surface parameters
include 'soilv.h'            ! Soil parameters
include 'stime.h'            ! File date data
include 'trcom2.h'           ! Station data

integer leap
common/leap_yr/leap          ! Leap year (1 to allow leap years)
integer nbarewet,nsigmf
common/nsib/nbarewet,nsigmf  ! Land-surface options

! for cardin
data ia/1/,ib/3/,id/2/,ja/1/,jb/10/,jd/5/,nlv/1/
data ndi/1/,ndi2/0/,nmaxpr/99/     
data kdate_s/-1/,ktime_s/-1/,leap/0/
data mbd/0/,nbd/0/,nbox/1/,kbotdav/4/,kbotu/0/
data nud_p/0/,nud_q/0/,nud_t/0/,nud_uv/1/,nud_hrs/24/,nudu_hrs/0/
data ktopdav/0/,kblock/-1/
data nud_aero/0/
data nud_sst/0/,nud_sss/0/,nud_ouv/0/,nud_sfh/0/
data mloalpha/10/,kbotmlo/-1000/,ktopmlo/1/
      
! Dynamics options A & B      
data mex/30/,mfix/3/,mfix_qg/1/,mup/1/,nh/0/
data nritch_t/300/,epsp/-15./,epsu/0./,epsf/0./
data precon/-2900/,restol/4.e-7/
data schmidt/1./,rlong0/0./,rlat0/90./,nrun/0/
data helmmeth/0/,mfix_tr/0/,mfix_aero/0/
! Horiz advection options
data nt_adv/7/,mh_bs/4/
! Horiz wind staggering options
data nstag/-10/,nstagu/-1/,nstagoff/0/
! Horizontal mixing options (now in globpe)
! data khdif/2/,khor/-8/,nhor/-157/,nhorps/-1/,nhorjlm/1/
! Vertical mixing options
data nvmix/3/,nlocal/6/,ncvmix/0/,lgwd/0/,ngwd/-5/
data helim/800./,fc2/1./,sigbot_gwd/0./,alphaj/1.e-6/
data cgmap_offset/0./,cgmap_scale/1./
! Cumulus convection options
data nkuo/23/,sigcb/1./,sig_ct/1./,rhcv/0./,rhmois/.1/,rhsat/1./
data convfact/1.02/,convtime/.33/,shaltime/0./
data alflnd/1.1/,alfsea/1.1/,fldown/.6/,iterconv/3/,ncvcloud/0/
data nevapcc/0/,nevapls/-4/,nuvconv/0/
data mbase/101/,mdelay/-1/,methprec/8/,nbase/-4/,detrain/.15/
data entrain/.05/,methdetr/2/,detrainx/0./,dsig2/.15/,dsig4/.4/
! Shallow convection options
data ksc/-95/,kscsea/0/,kscmom/1/,sigkscb/.95/,sigksct/.8/
data tied_con/2./,tied_over/0./,tied_rh/.75/
! Other moist physics options
data acon/.2/,bcon/.07/,qgmin/1.e-6/,rcm/.92e-5/
data rcrit_l/.75/,rcrit_s/.85/ 
! Radiation options
data nrad/4/
data nmr/0/,bpyear/0./
! Cloud options
data ldr/1/,nclddia/1/,nstab_cld/0/,nrhcrit/10/,sigcll/.95/ 
data cldh_lnd/95./,cldm_lnd/85./,cldl_lnd/75./
data cldh_sea/95./,cldm_sea/90./,cldl_sea/80./
data ncloud/0/
! Soil, canopy, PBL options
data nbarewet/0/,newrough/0/,nglacier/1/
data nrungcm/-1/,nsib/3/,nsigmf/1/
data ntaft/2/,ntsea/6/,ntsur/6/,av_vmod/.7/,tss_sh/1./
data vmodmin/.2/,zobgin/.02/,charnock/.018/,chn10/.00125/
data newztsea/1/,newtop/1/                
data snmin/.11/  ! 1000. for 1-layer; ~.11 to turn on 3-layer snow
data nurban/0/,ccycle/0/
! Special and test options
data namip/0/,amipo3/.false./,nhstest/0/,nsemble/0/,nspecial/0/
data panfg/4./,panzo/.001/,nplens/0/,rlatdx/0./,rlatdn/0./
data rlongdn/0./,rlongdx/0./
data rescrn/0/,knh/-1/
! I/O options
data m_fly/4/,io_in/1/,io_out/1/,io_rest/1/
data nperavg/-99/,nwt/-99/,tblock/1/,tbave/1/
data nextout/3/,localhist/.false./,unlimitedhist/.true./
data nstn/0/  
data slat/nstnmax*-89./,slon/nstnmax*0./,iunp/nstnmax*6/
data zstn/nstnmax*0./,name_stn/nstnmax*'   '/ 
data compression/1/,filemode/0/,procformat/.false./,procmode/0/
data chunkoverride/0/
data pio/.false./
! Ocean options
data nmlo/0/
! Aerosol options
data iaero/0/      

! initialize file names to something
data albfile/' '/,icefile/' '/,maskfile/' '/
data snowfile/' '/,sstfile/' '/,topofile/' '/,zofile/' '/
data rsmfile/' '/,scamfile/' '/,soilfile/' '/,vegfile/' '/
data co2emfile/' '/,so4tfile/' '/
data smoistfile/' '/,soil2file/' '/,restfile/' '/
data radonemfile/' '/,surfile/' '/,surf_00/'s_00a '/
data surf_12/'s_12a '/,co2_00/' '/,co2_12/' '/,radon_00/' '/
data radon_12/' '/,ifile/' '/,ofile/' '/,nmifile/' '/
data eigenv/' '/,radfile/' '/,o3file/' '/,hfile/' '/,mesonest/' '/
data scrnfile/' '/,tmaxfile/' '/,tminfile/' '/,trcfil/' '/
data laifile/' '/,albnirfile/' '/,urbanfile/' '/,bathfile/' '/
data vegprev/' '/,vegnext/' '/,cnsdir/' '/,salfile/' '/
data oxidantfile/' '/,casafile/' '/,phenfile/' '/
! floating point:
data timer/0./,mtimer/0/

! stuff from insoil  for soilv.h
data rlaim44/4.8, 6.3, 5., 3.75, 2.78, 2.5, 3.9, 2.77, 2.04, 2.6,         & ! 1-10
             1.69, 1.9, 1.37, 1.5, 1.21, 1.58, 1.41, 2.3, 1.2, 1.71,      & ! 11-20
             1.21, 2.3, 2.3, 1.2, 1.2, 1.87, 1., 3., .01, .01, 1.2,       & ! 21-31
             6., 5.5, 5., 4.5, 5., 4., 3., 3.5, 1., 4., .5, 4., 0./         ! 32-44
data rlais44/1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,                      & ! 1-10
             1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,                      & ! 11-20
             1., 1., 1., 1., .6, .6, .5, 1., 0., 0., 1.,                  & ! 21-31
             2., 2., 2., 2., 2., 1.5, 1.5, 1.5, 1., .5, .5, .5, 0./         ! 32-44
data rsunc44/370., 330., 260., 200., 150., 130., 200., 150., 110., 160.,  & ! 1-10
             100., 120.,  90.,  90.,  80.,  90.,  90., 150.,  80., 100.,  & ! 11-20
             80.,  80.,  80.,  60.,  60., 120.,  80., 180., 2*995., 80.,  & ! 21-31
             350., 4*300., 3*230., 150., 230., 995., 150., 9900./           ! 32-44
data scveg44/0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,                      & ! 1-10
             0., 0., 0., 0., 0., .1, .1, .1, .1, .1,                      & ! 11-20
             .1, .2, .4, .2, .1, .1, .1, 0., 0., 0., 0.,                  & ! 21-31
             .05, 0., 0., 0., 0., .05, .05, .05, .1, 0., 0., .4, 0./        ! 32-44
data slveg44/0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,                      & ! 1-10
             0., 0., 0., 0., 0., .1, .1, .1, .1, .1,                      & ! 11-20
             .1, .2, .4, .2, .1, .1, .1, 0., 0., 0., 0.,                  & ! 21-31
             1., 5.5, 3., 1., 3., 3., 3.5, 3., .5, 3.5, .1, 3.5, 0./        ! 32-44
data froot/.05, .10, .35, .40, .10/       ! 10/02/99 veg. root distr.

data silt/.08, .33, .17, .2, .06, .25, .15, .70, .33, .2, .33, .33, .17/    ! with mxst=13
data clay/.09, .3, .67, .2, .42, .48, .27, .17, .30, .2, .3, .3, .67/       ! with mxst=13
data sand/.83, .37, .16, .6, .52, .27, .58, .13, .37, .6, .37, .37, .17/    ! with mxst=13
data swilt/0., .072, .216, .286, .135, .219, .283, .175, .395, .216, .1142, .1547, .2864, .2498/
data sfc/1.,  .143, .301, .367, .218, .31 , .37 , .255, .45, .301, .22 , .25 , .367, .294/
data ssat/2., .398, .479, .482, .443, .426, .482, .420, .451, .479, .435, .451, .482, .476/
data bch/4.2, 7.1, 11.4, 5.15, 10.4, 10.4, 7.12, 5.83, 7.1, 4.9, 5.39, 11.4, 8.52/ ! bch for gravity term
data hyds/166.e-6, 4.e-6, 1.e-6, 21.e-6, 2.e-6, 1.e-6, 6.e-6,800.e-6, 1.e-6, 34.e-6, 7.e-6, 1.3e-6, 2.5e-6/
data sucs/-.106, -.591, -.405, -.348, -.153, -.49, -.299,-.356, -.153, -.218, -.478, -.405, -.63/ ! phisat (m)
data rhos/7*2600., 1300.,  910., 4*2600. /     ! soil density
data  css/7* 850., 1920., 2100., 4*850./       ! heat capacity

data zse/.022, .058, .154, .409, 1.085, 2.872/ ! layer thickness
! so depths of centre of layers: .011, .051, .157, .4385, 1.1855, 3.164
! with base at 4.6     

end
      
!--------------------------------------------------------------
! WRITE STATION DATA
subroutine stationa

use arrays_m           ! Atmosphere dyamics prognostic arrays
use cc_mpi             ! CC MPI routines
use diag_m             ! Diagnostic routines
use estab              ! Liquid saturation function
use extraout_m         ! Additional diagnostics
use map_m              ! Grid map arrays
use morepbl_m          ! Additional boundary layer diagnostics
use nsibd_m            ! Land-surface arrays
use pbl_m              ! Boundary layer arrays
use prec_m             ! Precipitation
use screen_m           ! Screen level diagnostics
use sigs_m             ! Atmosphere sigma levels
use soil_m             ! Soil and surface data
use soilsnow_m         ! Soil, snow and surface data
use tracers_m          ! Tracer data
use vecsuv_m           ! Map to cartesian coordinates
use vegpar_m           ! Vegetation arrays
use work2_m            ! Diagnostic arrays
use work3_m            ! Mk3 land-surface diagnostic arrays
use xyzinfo_m          ! Grid coordinate arrays

implicit none

include 'newmpar.h'    ! Grid parameters
include 'const_phys.h' ! Physical constants
include 'dates.h'      ! Date data
include 'parm.h'       ! Model configuration
include 'parmgeom.h'   ! Coordinate data
include 'soilv.h'      ! Soil parameters
include 'trcom2.h'     ! Station data

integer leap
common/leap_yr/leap    ! Leap year (1 to allow leap years)

integer i, j, iq, iqt, isoil, k2, nn
real coslong, sinlong, coslat, sinlat, polenx, poleny, polenz
real zonx, zony, zonz, den, costh, sinth, uzon, vmer, rh1, rh2
real es, wbav, rh_s

coslong = cos(rlong0*pi/180.)   ! done here, where work2 has arrays
sinlong = sin(rlong0*pi/180.)
coslat  = cos(rlat0*pi/180.)
sinlat  = sin(rlat0*pi/180.)
polenx  = -coslat
poleny  = 0.
polenz  = sinlat
do nn = 1,nstn
  ! Check if this station is in this processors region
  if ( .not. mystn(nn) ) cycle 
  if ( ktau == 1 ) then
    write (iunp(nn),950) kdate,ktime,leap
  end if
950 format("#",i9,2i5)
  i = istn(nn)
  j = jstn(nn)
  iq = i + (j-1)*il
  zonx  = real(            -polenz*y(iq))
  zony  = real(polenz*x(iq)-polenx*z(iq))
  zonz  = real(polenx*y(iq)             )
  den   = sqrt( max(zonx**2+zony**2+zonz**2,1.e-7) ) 
  costh =  (zonx*ax(iq)+zony*ay(iq)+zonz*az(iq))/den
  sinth = -(zonx*bx(iq)+zony*by(iq)+zonz*bz(iq))/den
  uzon  = costh*u(iq,1)-sinth*v(iq,1)
  vmer  = sinth*u(iq,1)+costh*v(iq,1)
  es   = establ(t(iq,1))
  rh1  = 100.*qg(iq,1)*(ps(iq)*sig(1)-es)/(.622*es)
  es   = establ(t(iq,2))
  rh2  = 100.*qg(iq,2)*(ps(iq)*sig(2)-es)/(.622*es)
  es   = establ(tscrn(iq))
  rh_s = 100.*qgscrn(iq)*(ps(iq)-es)/(.622*es)
  wbav = (zse(1)*wb(iq,1)+zse(2)*wb(iq,2)+zse(3)*wb(iq,3)+zse(4)*wb(iq,4))/(zse(1)+zse(2)+zse(3)+zse(4))
  iqt = min( iq, ilt*jlt ) ! Avoid bounds problems if there are no tracers
  k2  = min( 2, klt )
  write (iunp(nn),951) ktau,tscrn(iq)-273.16,rnd_3hr(iq,8),      &
        tss(iq)-273.16,tgg(iq,1)-273.16,tgg(iq,2)-273.16,        &
        tgg(iq,3)-273.16,t(iq,1)-273.16,0.,wb(iq,1),wb(iq,2),    &
        cloudlo(iq),cloudmi(iq)+1.,cloudhi(iq)+2.,               &
        cloudtot(iq)+3.,fg(iq),eg(iq),0.,0.,rnet(iq),sgsave(iq), &
        qg(iq,1)*1.e3,uzon,vmer,precc(iq),qg(iq,2)*1.e3,rh1,rh2, &
        0.,0.,0.,0.,.01*ps(iq),wbav,epot(iq),qgscrn(iq)*1.e3,    &
        rh_s,u10(iq),uscrn(iq),condx(iq)
  ! N.B. qgscrn formula needs to be greatly improved
951 format(i4,6f7.2,                                             &
           2f7.2, 2f6.3, 4f5.2,                                  & ! t1 ... cld
           5f7.1,f6.1,f5.1,                                      & ! fg ... qg1
           2f6.1,f7.2, f5.1,2f6.1, 2(1x,f5.1),                   & ! uu ... co2_2
           2(1x,f5.1) ,f7.1,f6.3,f7.1,5f6.1,                     & ! rad_1 ... rh_s
           f7.2)                                                   ! condx
  if ( ktau == ntau ) then
    write (iunp(nn),952)
952 format("#   2tscrn 3precip 4tss  5tgg1  6tgg2  7tgg3",               &
           "   8t1    9tgf  10wb1 11wb2 cldl cldm cldh  cld",            &
           "   16fg   17eg  18fgg  19egg  20rnet 21sg 22qg1",            &
           " 23uu   24vv 25precc qg2  rh1 28rh2 29co2_1 co2_2",          &
           " rad_1 rad_2  ps 34wbav 35epot qgscrn 37rh_s 38u10 uscrn",   &
           " 40condx")
    write (iunp(nn),953) land(iq),isoilm(iq),ivegt(iq),zo(iq),zs(iq)/grav
953 format("# land,isoilm,ivegt,zo,zs/g: ",l2,2i3,2f9.3)
    isoil = max(1,isoilm(iq))
    write (iunp(nn),954) sigmf(iq),swilt(isoil),sfc(isoil),ssat(isoil),0.5*sum(albvisnir(iq,:))
954 format("#sigmf,swilt,sfc,ssat,alb: ",5f7.3)
  end if
end do
return
end subroutine stationa
      
!--------------------------------------------------------------
! PREVIOUS VERSION DEFAULT PARAMETERS
subroutine change_defaults(nversion,mins_rad)

use parmhdff_m              ! Horizontal diffusion parameters

implicit none

include 'newmpar.h'         ! Grid parameters
include 'kuocom.h'          ! Convection parameters
include 'parm.h'            ! Model configuration
include 'parmdyn.h'         ! Dynamics parmaters
include 'parmhor.h'         ! Horizontal advection parameters
include 'parmsurf.h'        ! Surface parameters

integer nbarewet,nsigmf
common/nsib/nbarewet,nsigmf ! Land-surface options

integer, intent(in) :: nversion
integer, intent(inout) :: mins_rad

if ( nversion < 1510 ) then
  mins_rad = 60
end if
if ( nversion < 907 ) then
  mfix = 1         ! new is 3
  newrough = 2     ! new is 0
  newtop = 0       ! new is 1
  nvmix = 5        ! new is 3
  ksc = 0          ! new is -95
  sig_ct = .8      ! new is 1.
end if
if ( nversion < 904 ) then
  newtop = 1       ! new is 0
  nvmix = 3        ! new is 5
  ksc = -95        ! new is 0
  sig_ct = -.8     ! new is .8
end if
if( nversion < 809 ) then
  nvmix = 5        ! new is 3
  ksc = 0          ! new is -95
  sig_ct = .8      ! new is -.8
end if
if ( nversion < 806 ) then
  nvmix = 3        ! new is 5
  ksc = -95        ! new is 0
  nclddia = 5      ! new is 1
end if
if ( nversion == 803 ) then
  restol = 2.e-7   ! new is 4.e-7
end if
if ( nversion < 803 ) then
  restol = 5.e-7   ! new is 2.e-7
  alflnd = 1.15    ! new is 1.1
  alfsea = 1.05    ! new is 1.1
  entrain = 0.     ! new is .05
  ksc = 0          ! new is -95
endif
if ( nversion == 709 ) then
  ksc = 99
end if
if ( nversion < 709 ) then
  precon = 0       ! new is -2900
  restol = 2.e-7   ! new is 5.e-7
  mbase = 2000     ! new is 101
  mdelay = 0       ! new is -1
  nbase = -2       ! new is -4
  sigkscb = -.2    ! new is .95
  sigksct = .75    ! new is .8
  tied_con = 6.    ! new is 2.
  tied_over = 2.   ! new is 0.
  tied_rh = .99    ! new is .75
end if
if ( nversion < 705 ) then
  nstag = 5        ! new is -10
  nstagu = 5       ! new is -1.
  detrain = .3     ! new is .15
end if
if ( nversion < 704 ) then
  mex = 4          ! new is 30.
  ntsur = 2        ! new is 6
end if
if ( nversion < 703 ) then
  ntbar = 4        ! new is 6
  ntsur = 7        ! new is 2
  vmodmin = 2.     ! new is .2
  nbase = 1        ! new is -2
end if
if ( nversion < 701 ) then
  nbase = 0        ! new is 1  new variable
end if
if ( nversion < 608 ) then
  epsp = -20.      ! new is -15.
end if
if ( nversion < 606 ) then
  epsp = 1.1       ! new is -20.
  newrough = 0     ! new is 2
  nstag = -10      ! new is 5
  nstagu = 3       ! new is 5
  ntsur = 6        ! new is 7
  mbase = 10       ! new is 2000
end if
if ( nversion < 604 ) then
  mh_bs = 3        ! new is 4
end if
if ( nversion < 602 ) then
  ntbar = 9        ! new is 4
end if
if ( nversion < 601 ) then
  epsp = 1.2       ! new is 1.1
  newrough = 2     ! new is 0
  restol = 1.e-6   ! new is 2.e-7
end if
if ( nversion < 511 ) then
  nstag = 3        ! new is -10
  ! mins_rad = 120   ! new is 72   not in common block, so can't assign here
  detrain = .1     ! new is .3
  mbase = 1        ! new is 10
  nuvconv = 5      ! new is 0
  sigcb = .97      ! new is 1.
end if
if ( nversion < 510 ) then
  epsp = .1        ! new is 1.2
  epsu = .1        ! new is 0.
  khdif = 5        ! new is 2
  khor = 0         ! new is -8
  nbarewet = 7     ! new is 0
  newrough = 0     ! new is 2
  nhor = 0         ! new is -157
  nhorps = 1       ! new is -1
  nlocal = 5       ! new is 6
  ntsur = 7        ! new is 6
  ! jalbfix = 0      ! new is 1  not in common block, so can't assign here
  tss_sh = 0.      ! new is 1.
  zobgin = .05     ! new is .02
  detrain = .4     ! new is .1
  convtime = .3    ! new is .33
  iterconv = 2     ! new is 3
  mbase = 0        ! new is 1
  sigcb = 1.       ! new is .97
  sigkscb = .98    ! new is -2.
  tied_rh = .75    ! new is .99
  ldr = 2          ! new is 1
  rcm = 1.e-5      ! new is .92e-5
end if
if ( nversion < 509 ) then
  ntsur = 6        ! new is 7
end if
if ( nversion < 508 ) then
  mh_bs = 1        ! new is 3
  nvmix = 4        ! new is 3
  entrain = .3     ! new is 0.
endif
if ( nversion < 506 ) then
  mh_bs = 4        ! new is 1
end if
if ( nversion < 503 ) then
  ntsur = 5        ! new is 6
end if
if ( nversion < 411 ) then
  nstag = -3       ! new is 3
  nstagu = -3      ! new is 3
  nhor = 155       ! new is 0
  nlocal = 1       ! new is 5
  ngwd = 0         ! new is -5
  nevapls = 5      ! new is -4
  nuvconv = 0      ! new is 5
  detrain = .05    ! new is .4
  entrain = 0.     ! new is .3
  detrainx = 1.    ! new is 0.
  dsig2 = .1       ! new is .15
  dsig4 = .55      ! new is .4
  kscmom = 0       ! new is 1
  ldr = 1          ! new is 2
  nbarewet = 2     ! new is 7
  av_vmod = 1.     ! new is .7
  chn10 = .00137   ! new is .00125
end if

return
end subroutine change_defaults

!--------------------------------------------------------------
! ERROR MESSAGE FOR INVALID nproc
subroutine badnproc(npanels,il_g,nproc)

use cc_mpi                                 ! CC MPI routines

implicit none

integer, intent(in) :: il_g, nproc, npanels
integer nxpa, nxpb, nyp, nproc_low, nproc_high
integer ilg_low, ilg_high

if ( myid == 0 ) then
  write(6,*)
  write(6,*) "ERROR: Invalid number of processors for this grid"
  do nproc_low = nproc,1,-1
    call proctest(npanels,il_g,nproc_low,nxpa,nyp)
    if ( nxpa > 0 ) exit
  end do
  do nproc_high = nproc,10*nproc
    call proctest(npanels,il_g,nproc_high,nxpb,nyp)
    if ( nxpb > 0 ) exit
  end do
  if ( nxpb > 0 ) then
    write(6,*) "Consider using processor numbers ",nproc_low," or ",nproc_high
  else
    write(6,*) "Consider using processor number  ",nproc_low  
  end if
  do ilg_low = il_g,1,-1
    call proctest(npanels,ilg_low,nproc,nxpa,nyp)
    if ( nxpa > 0 ) exit
  end do
  do ilg_high = il_g,10*il_g
    call proctest(npanels,ilg_high,nproc,nxpb,nyp)
    if ( nxpb > 0 ) exit
  end do    
  if ( nxpa>0 .and. nxpb>0 ) then
    write(6,*) "Alternatively, try grid sizes    ",ilg_low," or ",ilg_high
  else if ( nxpa > 0 ) then
    write(6,*) "Alternatively, try grid size     ",ilg_low
  else if ( nxpb > 0 ) then
    write(6,*) "Alternatively, try grid size     ",ilg_high  
  end if
  write(6,*)
end if
call ccmpi_barrier(comm_world)
call ccmpi_abort(-1)

return
end subroutine badnproc

!--------------------------------------------------------------
! TEST GRID DECOMPOSITION    
subroutine proctest(npanels,il_g,nproc,nxp,nyp)

implicit none

integer, intent(in) :: il_g, nproc, npanels
integer, intent(out) :: nxp, nyp
integer jl_g

#ifdef uniform_decomp
jl_g = il_g + npanels*il_g     ! size of grid along all panels (usually 6*il_g)
nxp = nint(sqrt(real(nproc)))  ! number of processes in X direction
nyp = nproc/nxp                ! number of processes in Y direction
! search for vaild process decomposition.  CCAM enforces the same grid size on each process
do while ( (mod(il_g,max(nxp,1))/=0.or.mod(nproc,max(nxp,1))/=0.or.mod(il_g,nyp)/=0) .and. nxp>0 )
  nxp = nxp - 1
  nyp = nproc/max(nxp,1)
end do
#else
if ( mod(nproc,6)/=0 .and. mod(6,nproc)/=0 ) then
  nxp = -1
else
  jl_g = il_g + npanels*il_g                 ! size of grid along all panels (usually 6*il_g)
  nxp = max( 1, nint(sqrt(real(nproc)/6.)) ) ! number of processes in X direction
  nyp = nproc/nxp                            ! number of processes in Y direction
  ! search for valid process decomposition.  CCAM enforces the same grid size on each process
  do while ( (mod(il_g,max(nxp,1))/=0.or.mod(nproc/6,max(nxp,1))/=0.or.mod(jl_g,max(nyp,1))/=0) .and. nxp>0 )
    nxp = nxp - 1
    nyp = nproc/max(nxp,1)
  end do
end if
#endif

return
end subroutine proctest
    
