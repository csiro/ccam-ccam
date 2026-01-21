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

module config

private
public globpe_init, globpe_finalize

contains

!--------------------------------------------------------------
! INITIALISE CCAM
subroutine globpe_init

use aerointerface, only : aeroindir      & ! Aerosol interface
    ,aero_split,aerosol_u10              &
    ,ch_dust,zvolcemi,so4mtn,carbmtn     &
    ,saltsmallmtn,saltlargemtn           &
    ,enhanceu10,naero
use arrays_m                               ! Atmosphere dyamics prognostic arrays
use bigxy4_m                               ! Grid interpolation
use cc_acc                                 ! CC ACC routines
use cc_mpi                                 ! CC MPI routines
use cfrac_m                                ! Cloud fraction
use const_phys                             ! Physical constants
use darcdf_m                               ! Netcdf data
use daviesnudge                            ! Far-field nudging
use diag_m                                 ! Diagnostic routines
use dpsdt_m                                ! Vertical velocity
use epst_m                                 ! Off-centre terms
use estab                                  ! Liquid saturation function
use extraout_m                             ! Additional diagnostics
use filnames_m                             ! Filenames
use gdrag_m, only : gdrag_init             ! Gravity wave drag
use getopt_m                               ! Command option parsing
use histave_m                              ! Time average arrays
use indata                                 ! Data initialisation
use indices_m                              ! Grid index arrays
use infile                                 ! Input file routines
use kuocom_m                               ! JLM convection
use latlong_m                              ! Lat/lon coordinates
use liqwpar_m                              ! Cloud water mixing ratios
use map_m                                  ! Grid map arrays
use mlodynamics                            ! Ocean dynamics
use module_aux_rad                         ! Additional cloud and radiation routines
use module_ctrl_microphysics               ! Interface for cloud microphysics
use morepbl_m                              ! Additional boundary layer diagnostics
use newmpar_m                              ! Grid parameters
use nharrs_m                               ! Non-hydrostatic atmosphere arrays
use nlin_m                                 ! Atmosphere non-linear dynamics
use nsibd_m                                ! Land-surface arrays
use ozoneread                              ! Ozone input routines
use parm_m                                 ! Model configuration
use parmdyn_m                              ! Dynamics parameters
use parmgeom_m                             ! Coordinate data
use parmhdff_m                             ! Horizontal diffusion parameters
use parmhor_m                              ! Horizontal advection parameters
use parmvert_m                             ! Vertical advection parameters
use pbl_m                                  ! Boundary layer arrays
use permsurf_m, only : permsurf_init       ! Fixed surface arrays
use prec_m                                 ! Precipitation
use raddiag_m                              ! Radiation diagnostic
use river                                  ! River routing
use riverarrays_m                          ! River data
use savuvt_m                               ! Saved dynamic arrays
use savuv1_m                               ! Saved dynamic arrays
use sbar_m                                 ! Saved dynamic arrays
use screen_m                               ! Screen level diagnostics
use seaesfrad_m                            ! SEA-ESF radiation
use setxyz_m                               ! Define CCAM grid
use sflux_m                                ! Surface flux routines
use sigs_m                                 ! Atmosphere sigma levels
use soil_m                                 ! Soil and surface data
use soilsnow_m                             ! Soil, snow and surface data
use staguvmod                              ! Reversible grid staggering 
use stime_m                                ! File date data
use tbar2d_m, only : tbar2d_init           ! Atmosphere dynamics reference temperature
use tkeeps                                 ! TKE-EPS boundary layer
use tracermodule, only : tracerlist      & ! Tracer routines
    ,sitefile,shipfile,writetrpm         &
    ,init_tracer
use tracers_m                              ! Tracer data
use unn_m                                  ! Saved dynamic arrays
use usage_m                                ! Usage message
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

include 'version.h'                        ! Model version data

integer, dimension(:), allocatable, save :: dumi
integer, dimension(3) :: shsize ! for share_ifullg
integer ierr, k, new_nproc, ilx, jlx, i, ng
integer isoth, nsig, lapsbot
integer secs_rad, nversion
integer mstn, mbd_min
integer opt, nopt, secs_nud
integer ateb_intairtmeth, ateb_intmassmeth
integer npa, npb, tkecduv, tblock  ! depreciated namelist options
integer o3_time_interpolate        ! depreciated namelist options
integer kmlo, calcinloop           ! depreciated namelist options
integer fnproc_bcast_max, nriver   ! depreciated namelist options
integer ateb_conductmeth           ! depreciated namelist options
integer ateb_useonewall            ! depreciated namelist options
integer cable_climate              ! depreciated namelist options
integer surf_windfarm, adv_precip  ! depreciated namelist options
real, dimension(:,:), allocatable, save :: dums
real, dimension(:), allocatable, save :: dumr, gosig_in
real, dimension(8) :: temparray
real, dimension(1) :: gtemparray
real targetlev, dsx, pwatr_l, pwatr, tscale
real ateb_zocanyon, ateb_zoroof, ateb_energytol
real cgmap_offset, cgmap_scale      ! depreciated namelist options
real ateb_ac_smooth, ateb_ac_copmax ! depreciated namelist options
real ateb_alpha, cable_version      ! depreciated namelist options 
real zimax,mlomaxuv                 ! depreciated namelist options
real plume_alpha                    ! depreciated namelist options
real ocnlap                         ! depreciated namelist options
logical procformat                  ! depreciated namelist options
logical unlimitedhist               ! depreciated namelist options
character(len=1024) nmlfile
character(len=MAX_ARGLEN) optarg
character(len=60) comm, comment
character(len=47) header
character(len=10) timeval
character(len=8) text, rundate
character(len=1024) vegprev, vegnext, vegnext2 ! depreciated namelist options

! version namelist
namelist/defaults/nversion
! main namelist
namelist/cardin/comment,dt,ntau,nwt,nhorps,nperavg,ia,ib,         &
    ja,jb,id,jd,iaero,khdif,khor,nhorjlm,mex,mbd,nbd,             &
    mbd_maxscale,mbd_maxgrid,ndi,ndi2,nhor,nlv,nmaxpr,nrad,ntaft, &
    ntsea,ntsur,nvmix,restol,precon,kdate_s,ktime_s,leap,newtop,  &
    mup,lgwd,ngwd,rhsat,nextout,jalbfix,nalpha,nstag,nstagu,      &
    ntbar,nwrite,irest,nrun,nstn,nrungcm,nsib,istn,jstn,iunp,     &
    slat,slon,zstn,name_stn,mh_bs,nritch_t,nt_adv,mfix,mfix_qg,   &
    namip,amipo3,nh,nhstest,nsemble,nspecial,panfg,panzo,         &
    rlatdn,rlatdx,rlongdn,rlongdx,newrough,nglacier,newztsea,     &
    epsp,epsu,epsf,epsh,av_vmod,charnock,chn10,snmin,tss_sh,      &
    vmodmin,zobgin,rlong0,rlat0,schmidt,kbotdav,kbotu,nud_p,      &
    nud_q,nud_t,nud_uv,nud_hrs,nudu_hrs,sigramplow,sigramphigh,   &
    nlocal,nbarewet,nsigmf,io_in,io_nest,io_out,io_rest,          &
    tblock,tbave,localhist,synchist,m_fly,                        &
    nurban,ktopdav,mbd_mlo,mbd_maxscale_mlo,nud_sst,nud_sss,      &
    mfix_tr,mfix_aero,kbotmlo,ktopmlo,mloalpha,nud_ouv,nud_sfh,   &
    rescrn,helmmeth,nmlo,ol,knh,kblock,nud_aero,                  &
    nud_period,mfix_t,zo_clearing,intsch_mode,qg_fix,             &
    always_mspeca,ntvd,tbave10,maxuv,maxcolour,estab_bug_fix,     &
    procmode,compression,hp_output,pil_single,process_rate_mode,  & ! file io
    chunk_time,                                                   &
    maxtilesize,async_length,nagg,                                & ! MPI, OMP & ACC
    ensemble_mode,ensemble_period,ensemble_rsfactor,              & ! ensemble
    ch_dust,helim,fc2,sigbot_gwd,alphaj,nmr,qgmin,mstn,           & ! backwards compatible
    npa,npb,cgmap_offset,cgmap_scale,procformat,fnproc_bcast_max, & ! depreciated
    nriver,unlimitedhist,adv_precip                                 ! depreciated
! radiation and aerosol namelist
namelist/skyin/mins_rad,sw_resolution,sw_diff_streams,            & ! radiation
    liqradmethod,iceradmethod,so4radmethod,carbonradmethod,       &
    dustradmethod,seasaltradmethod,bpyear,qgmin,lwem_form,        & 
    siglow,sigmid,linecatalog_form,continuum_form,do_co2_10um,    &
    do_quench,remain_rayleigh_bug,use_rad_year,rad_year,          &
    ch_dust,zvolcemi,aeroindir,so4mtn,carbmtn,saltsmallmtn,       & ! aerosols
    saltlargemtn,enhanceu10,aerosol_u10,aero_split,               &
    o3_vert_interpolate,                                          & ! ozone
    o3_time_interpolate                                             ! depreciated
! file namelist
namelist/datafile/ifile,ofile,albfile,eigenv,icefile,mesonest,    &
    o3file,radfile,restfile,rsmfile,so4tfile,soilfile,sstfile,    &
    surfile,topofile,vegfile,zofile,surf_00,surf_12,laifile,      &
    albnirfile,urbanfile,bathfile,freqfile,                       &
    cnsdir,salfile,oxidantfile,casafile,phenfile,casapftfile,     &
    ensembleoutfile,solarfile,ch4file,n2ofile,cfc11file,          &
    cfc12file,cfc113file,hcfc22file,                              &
    save_aerosols,save_pbl,save_cloud,save_land,save_maxmin,      &
    save_ocean,save_radiation,save_urban,save_carbon,save_river,  &
    diaglevel_aerosols,diaglevel_pbl,diaglevel_cloud,             &
    diaglevel_land,diaglevel_maxmin,diaglevel_ocean,              &
    diaglevel_radiation,diaglevel_urban,diaglevel_carbon,         &
    diaglevel_river,diaglevel_pop,                                &
    surf_cordex,surf_windfarm,output_windmax,cordex_fix,          &
    wbclimfile,shep_cordex,                                       &
    vegprev,vegnext,vegnext2                                        ! depreciated
! convection and cloud microphysics namelist
namelist/kuonml/alflnd,alfsea,cldh_lnd,cldm_lnd,cldl_lnd,         & ! convection
    cldh_sea,cldm_sea,cldl_sea,convfact,convtime,shaltime,        &
    detrain,detrainx,dsig2,dsig4,entrain,fldown,iterconv,ksc,     &
    kscmom,kscsea,ldr,mbase,mdelay,methdetr,methprec,nbase,       &
    ncvcloud,ncvmix,nevapcc,nkuo,nrhcrit,                         &
    nstab_cld,nuvconv,rhcv,rhmois,rhsat,sigcb,sigcll,sig_ct,      &
    sigkscb,sigksct,tied_con,tied_over,tied_rh,comm,acon,bcon,    &
    rcm,nscheme,                                                  &
    rcrit_l,rcrit_s,ncloud,nclddia,nmr,nevapls,cld_decay,         & ! cloud
    vdeposition_mode,tiedtke_form,cloud_aerosol_mode,             &
    cloud_ice_method,leon_snowmeth,lin_aerosolmode,maxlintime,    &
    lin_adv,qlg_max,qfg_max                                                       
! boundary layer turbulence and gravity wave namelist
namelist/turbnml/be,cm0,ce0,ce1,ce2,ce3,cqmix,ent0,ent1,entc0,    & ! EDMF PBL scheme
    dtrc0,m0,b1,b2,buoymeth,maxdts,mintke,mineps,minl,maxl,       &
    stabmeth,tkemeth,qcmf,ezmin,ent_min,mfbeta,                   &
    tke_timeave_length,plume_alpha,tcalmeth,                      &
    wg_tau,wg_prob,ugs_meth,                                      & ! wind gusts
    amxlsq,dvmodmin,                                              & ! JH PBL scheme
    ngwd,helim,fc2,sigbot_gwd,alphaj,                             & ! GWdrag
    tkecduv,zimax                                                   ! depreciated
! land, urban and carbon namelist
namelist/landnml/proglai,ccycle,soil_struc,cable_pop,             & ! CABLE
    progvcmax,fwsoil_switch,cable_litter,                         &
    gs_switch,cable_climate,smrf_switch,strf_switch,              &
    cable_gw_model,cable_roughness,cable_potev,                   &
    wt_transport, cable_enablefao,                                &
    ateb_energytol,ateb_resmeth,ateb_zohmeth,                     & ! urban
    ateb_acmeth,ateb_nrefl,                                       &
    ateb_scrnmeth,ateb_wbrelaxc,ateb_wbrelaxr,                    &
    ateb_nfgits,ateb_tol,                                         &
    ateb_zosnow,ateb_snowemiss,ateb_maxsnowalpha,                 &
    ateb_minsnowalpha,ateb_maxsnowden,ateb_minsnowden,            &
    ateb_refheight,ateb_zomratio,ateb_zocanyon,ateb_zoroof,       &
    ateb_maxrfwater,ateb_maxrdwater,ateb_maxrfsn,ateb_maxrdsn,    &
    ateb_maxvwatf,ateb_intairtmeth,ateb_intmassmeth,              &
    ateb_cvcoeffmeth,ateb_statsmeth,ateb_lwintmeth,               &
    ateb_infilmeth,ateb_ac_heatcap,ateb_ac_coolcap,               &
    ateb_ac_deltat,ateb_acfactor,ateb_soilunder,                  &
    siburbanfrac,freshwaterlake_fix,                              & ! special options
    wbclim_lonn,wbclim_lonx,wbclim_latn,wbclim_latx,              & 
    ateb_ac_smooth,ateb_ac_copmax,ateb_conductmeth,ateb_ncyits,   & ! depreciated
    ateb_useonewall,ateb_alpha,cable_version
! ocean namelist
namelist/mlonml/mlodiff,ocnsmag,ocneps,usetide,zomode,zoseaice,   & ! MLO
    factchseaice,minwater,mxd,mindep,otaumode,alphavis_seaice,    &
    alphanir_seaice,mlojacobi,usepice,mlosigma,nodrift,           &
    kmlo,mlontvd,alphavis_seasnw,alphanir_seasnw,mlodiff_numits,  &
    ocnlap,mlo_adjeta,mstagf,mlodps,mlo_limitsal,nxtrrho,mlo_bs,  &
    mlo_step,mlo_uvcoupl,fluxwgt,mlointschf,ocnepr,delwater,      &
    mloiceadv,minsal,maxsal,                                      &
    pdl,pdu,k_mode,eps_mode,limitL,fixedce3,nops,nopb,            & ! k-e
    fixedstabfunc,omink,omineps,oclosure,ominl,omaxl,             &
    mlo_timeave_length,kemaxdt,                                   &
    rivermd,basinmd,rivercoeff,                                   & ! River
    mlomfix,calcinloop,mlomaxuv                                     ! Depreciated
! tracer namelist
namelist/trfiles/tracerlist,sitefile,shipfile,writetrpm

! some defaults to avoid confusion
tblock = 0
kmlo = 0
calcinloop = 0
ateb_ac_copmax = 0.
ateb_ac_smooth = 0.
zimax = 0.
tkecduv = 0.
procformat = .true.
unlimitedhist = .false.
adv_precip = 0

!--------------------------------------------------------------
! READ COMMAND LINE OPTIONS
nmlfile = "input"
do
  call getopt("hc:",nopt,opt,optarg)
  if ( opt==-1 ) exit  ! End of options
  select case ( char(opt) )
    case ( "h" )
      call help
    case ( "c" )
      nmlfile = optarg
    case default
      if ( myid==0 ) write(6,*) "ERROR: Unknown command line option ",char(opt)
      call usage
  end select
end do


!--------------------------------------------------------------
! READ NAMELISTS AND SET PARAMETER DEFAULTS
nversion            = 0
comm                = ' '
comment             = ' '
ia                  = -1   ! diagnostic index
ib                  = -1   ! diagnostic index
ntbar               = -1
ktau                = 0
ol                  = 20   ! default ocean levels
nhor                = -157
nhorps              = -1
khor                = -8
khdif               = 2
nhorjlm             = 1
ngas                = 0
ateb_energytol      = 0.1
ateb_intairtmeth    = 0
ateb_intmassmeth    = 0
ateb_zocanyon       = zocanyon
ateb_zoroof         = zoroof
lapsbot             = 0
npa                 = 0   ! depreciated
npb                 = 0   ! depreciated
cgmap_offset        = 0.  ! depreciated
cgmap_scale         = 0.  ! depreciated
o3_time_interpolate = 0   ! depreciated
vegprev             = ' ' ! depreciated
vegnext             = ' ' ! depreciated
vegnext2            = ' ' ! depreciated


! All processors read the namelist, so no MPI comms are needed
if ( myid==0 ) then
  open(99,file=trim(nmlfile),form="formatted",status="old",iostat=ierr)
  if ( ierr/=0 ) then
    write(6,*) "ERROR: Cannot open namelist ",trim(nmlfile)  
    call ccmpi_abort(-1)
  end if
  read(99, defaults)
end if
call ccmpi_bcast(nversion,0,comm_world)
if ( nversion/=0 ) then
  call change_defaults(nversion)
end if
if ( myid==0 ) then
  read(99, cardin)
end if
call broadcast_cardin
if ( myid==0 ) then
  read(99, skyin)
end if
call broadcast_skyin
if ( myid==0 ) then
  read(99, datafile)
end if
call broadcast_datafile
if ( myid==0 ) then
  read(99, kuonml)
end if
call broadcast_kuonml
if ( myid==0 ) then
  read(99, turbnml, iostat=ierr)  ! try reading PBL and GWdrag namelist
  if ( ierr/=0 ) then
    rewind(99)
    ! if namelist is not missing, then trigger an error message
    if ( .not.is_iostat_end(ierr) ) read(99, turbnml)
  end if
end if
call broadcast_turbnml
if ( myid==0 ) then
  read(99, landnml, iostat=ierr)  ! try reading land/carbon namelist
  if ( ierr/=0 ) then
    rewind(99)
    ! if namelist is not missing, then trigger an error message
    if ( .not.is_iostat_end(ierr) ) read(99, landnml)
  end if
  energytol = real(ateb_energytol,8)
  zocanyon = ateb_zocanyon
  zoroof = ateb_zoroof
  intairtmeth = ateb_intairtmeth
  intmassmeth = ateb_intmassmeth
end if
call broadcast_landnml
if ( myid==0 ) then
  read(99, mlonml, iostat=ierr)   ! try reading ocean namelist
  if ( ierr/=0 ) then
    rewind(99)
    ! if namelist is not missing, then trigger an error message
    if ( .not.is_iostat_end(ierr) ) read(99, mlonml)
  end if
end if
call broadcast_mlonml
if ( myid==0 ) then
  rewind(99)  
  read(99, trfiles, iostat=ierr)  ! try reading tracer namelist
  if ( ierr/=0 ) then
    rewind(99)
    ! if namelist is not missing, then trigger an error message
    if ( .not.is_iostat_end(ierr) ) read(99, trfiles)
  end if
end if
call broadcast_trfiles
if ( myid==0 ) then
  close(99)
end if

if ( dt<=0. ) then
  write(6,*) "ERROR: dt must be greather than zero"
  call ccmpi_abort(-1)
end if
if ( dt>3600. ) then
  write(6,*) "ERROR: dt must be less or equal to 3600."
  call ccmpi_abort(-1)
end if
if ( nvmix==9 .and. (nmlo==0.or.nhstest<0) ) then
  write(6,*) "ERROR: nvmix=9 requires nmlo/=0 and nhstest>=0"
  call ccmpi_abort(-1)
end if
nagg = max( nagg, 4 ) ! use 4 for two staguv u & v arrays
nperday = nint(24.*3600./dt)           ! time-steps in one day
nperhr  = nint(3600./dt)               ! time-steps in one hour
nper6hr = nint(6.*3600./dt)            ! time-steps in six hours
if ( nwt==-99 )     nwt = nperday      ! set default nwt to 24 hours
if ( nperavg==-99 ) nperavg = nwt      ! set default nperavg to nwt
if ( nwrite==0 )    nwrite = nperday   ! only used for outfile IEEE
if ( nwt<=0 ) then
  write(6,*) "ERROR: nwt must be greater than zero or nwt=-99"
  call ccmpi_abort(-1)
end if
if ( nhstest<0 .and. nmlo/=0 ) then
  write(6,*) "ERROR: nhstest<0 requires nmlo=0"
  call ccmpi_abort(-1)
end if
if ( nmlo/=0 .and. abs(nmlo)<=9 ) then ! set ocean levels if required
  ol = max( ol, 1 )
else
  ol = 0
end if
wlev     = ol                   ! set nmlo and nmlodynamics ocean levels
mindep   = max( 0., mindep )    ! limit ocean minimum depth below sea-level
minwater = max( 0., minwater )  ! limit ocean minimum water level
! Update radiation sea-ice albedo to match MLO albedo values
seaice_albvis = alphavis_seaice
seaice_albnir = alphanir_seaice


!--------------------------------------------------------------
! READ TOPOGRAPHY FILE TO DEFINE CONFORMAL CUBIC GRID

il_g    = 48 ! default global grid size (replaced with size in topography file)
rlong0  = 0. ! default longitude (replaced with longitude in topography file)
rlat0   = 0. ! default latitude (replaced with latitiude in topography file)
schmidt = 1. ! default schmidt factor for grid stretching (replaced with schmidt in topography file)
kl      = 18 ! default number of vertical levels (replaced with levels in eigen file)

if ( myid==0 ) then
  ! open topo file and check its dimensions
  ! here used to supply rlong0,rlat0,schmidt
  ! Remander of topo file is read in indata.f90
  lnctopo = -1 ! flag indicating file not yet identified
  ! NetCDF format
  if ( lnctopo==-1 ) then
    call ccnf_open(topofile,ncidtopo,ierr)
    if ( ierr==0 ) then
      lnctopo = 1 ! flag indicating netcdf file
      call ccnf_inq_dimlen(ncidtopo,'longitude',ilx)
      call ccnf_inq_dimlen(ncidtopo,'latitude',jlx)
      call ccnf_get_attg(ncidtopo,'lon0',rlong0)
      call ccnf_get_attg(ncidtopo,'lat0',rlat0)
      call ccnf_get_attg(ncidtopo,'schmidt',schmidt) 
    end if
  end if
  ! ASCII format      
  if ( lnctopo==-1 ) then
    open(66,file=topofile,recl=2000,status='old',iostat=ierr)
    if ( ierr==0 ) then
      lnctopo = 0 ! flag indicating ASCII file  
      read(66,*) ilx,jlx,rlong0,rlat0,schmidt,dsx,header
    end if  
  end if
  ! Failed to read topo file
  if ( lnctopo==-1 ) then
    write(6,*) "Error opening topofile ",trim(topofile)
    call ccmpi_abort(-1)
  end if
  ! specify grid size based on topography file dimensions
  il_g = ilx        
  ! store grid dimensions for broadcast below
  temparray(1) = rlong0
  temparray(2) = rlat0
  temparray(3) = schmidt
  temparray(4) = real(il_g)
end if      ! (myid==0)


!--------------------------------------------------------------
! READ EIGENV FILE TO DEFINE VERTICAL LEVELS

if ( myid==0 ) then
  ! Remanded of file is read in indata.f90
  open(28,file=eigenv,status='old',form='formatted',iostat=ierr)
  if ( ierr/=0 ) then
    write(6,*) "Error opening eigenv file ",trim(eigenv)
    call ccmpi_abort(-1)
  end if
  read(28,*)kl,lapsbot,isoth,nsig
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
! Face decomposition reduces the number of MPI messages, but only works for factors or multiples
! of six processes.
! Uniform decomposition has been depreciated due to performance issues (more MPI messages)
! that provides misleading scaling results with increasing nproc.
call reducenproc(npanels,il_g,nproc,new_nproc,nxp,nyp)
call ccmpi_reinit(new_nproc) 

jl_g    = il_g + npanels*il_g                 ! size of grid along all panels (usually 6*il_g)
ifull_g = il_g*jl_g                           ! total number of global horizontal grid points
iquad   = 1 + il_g*((8*npanels)/(npanels+4))  ! grid size for interpolation
il      = il_g/nxp                            ! local grid size on process in X direction
jl      = jl_g/nyp                            ! local grid size on process in Y direction
ifull   = il*jl                               ! total number of local horizontal grid points
! The perimeter of the processor region has length 2*(il+jl).
! The first row has 8 possible corner points per panel and the 
! second has 16. In practice these are not all distinct so there could
! be some optimisation.
npan = max(1, (npanels+1)/nproc)   ! number of panels on this process
iextra = (4*(il+jl)+24*npan) + 4   ! size of halo for MPI message passing (jl includes npan)
! nrows_rad is a subgrid decomposition for older radiation routines
nrows_rad = max( min( maxtilesize/il, jl ), 1 ) 
do while( mod(jl, nrows_rad) /= 0 )
  nrows_rad = nrows_rad - 1
end do
! tiles for newer physics routines
call calc_phys_tiles(ntiles,maxtilesize,ifull)
imax = ifull/ntiles

! Initalise OpenACC for GPUs and OpenMP
#ifdef usempi3
! since processes might have been remapped, then use node_myid
! to determine GPU assigned to each process
call ccacc_init(node_myid,ngpus)
#else
call ccacc_init(myid,ngpus)
#endif

! Display model configuration information in log file
if ( myid==0 ) then
  write(6,'(" ",A)') trim(version)
  write(6,*) 'Running for nproc                        = ',nproc
  write(6,*) 'Using defaults for nversion              = ',nversion
#ifdef usempi3
#ifdef share_ifullg
  write(6,*) 'Using shared memory with number of nodes = ',nodecaptain_nproc
#else
  write(6,*) 'Node aware with number of nodes          = ',nodecaptain_nproc
#endif
#endif
#ifdef i8r8
  write(6,*) 'Using double precision mode'
#endif
#ifdef _OPENACC
  write(6,*) 'Using OpenACC with GPUs per node         = ',ngpus  
#endif
  write(6,*) 'Reading namelist from ',trim(nmlfile)
  write(6,*) 'rlong0,rlat0,schmidt ',rlong0,rlat0,schmidt
  write(6,*) 'kl,ol                ',kl,ol
  write(6,*) 'lapsbot,isoth,nsig   ',lapsbot,isoth,nsig
  write(6,*) 'ntiles,imax          ',ntiles,ifull/ntiles
  write(6,*) 'il_g,jl_g,il,jl      ',il_g,jl_g,il,jl
  write(6,*) 'nxp,nyp              ',nxp,nyp
end if

! some default values for unspecified parameters
if ( ia<0 ) ia = il/2          ! diagnostic point
if ( ib<0 ) ib = ia + 3        ! diagnostic point
dsig4 = max(dsig2+.01, dsig4)  ! convection

! check nudging settings - adjust mbd scale parameter to satisfy mbd_maxscale and mbd_maxgrid settings
if ( ensemble_mode>0 .and. (mbd/=0.or.nbd/=0.or.mbd_mlo/=0) ) then
  write(6,*) "ERROR: mbd=0, nbd=0 and mbd_mlo=0 are required for ensemble_mode>0"
  call ccmpi_abort(-1)
end if
if ( mbd/=0 .and. nbd/=0 ) then
  if ( myid==0 ) then  
    write(6,*) 'WARN: setting nbd=0 because mbd/=0'
  end if  
  nbd = 0
end if
if ( mbd<0 ) then
  write(6,*) "ERROR: mbd<0 is invalid"
  call ccmpi_abort(-1)
end if
if ( mbd/=0 ) then
  if ( mbd_maxscale==0 ) then
    write(6,*) "ERROR: mbd_maxscale must be >0 when mbd/=0"
    call ccmpi_abort(-1)
  end if
  mbd_min = int(20.*112.*90.*schmidt/real(mbd_maxscale))
  if ( mbd<mbd_min .and. mbd/=0 ) then
    if ( myid==0 ) then
      write(6,*) "Satisfy mbd_maxscale by increasing mbd = ",mbd_min
    end if
    mbd = mbd_min
  end if
  if ( mbd_maxgrid==0 ) then
    write(6,*) "ERROR: mbd_maxgrid must be >0 when mbd/=0"
    call ccmpi_abort(-1)
  end if
  mbd_min = int(20.*real(il_g)/real(mbd_maxgrid))
  if ( mbd<mbd_min .and. mbd/=0 ) then
    if ( myid==0 ) then
      write(6,*) "Satisfy mbd_maxgrid by adjusting mbd = ",mbd_min
    end if
    mbd = mbd_min
  end if
  nud_hrs = abs(nud_hrs)  ! just for people with old -ves in namelist
  if ( nudu_hrs==0 ) then
    nudu_hrs = nud_hrs
  end if
end if
if ( mbd_mlo<0 ) then
  write(6,*) "ERROR: mbd_mlo<0 is invalid"
  call ccmpi_abort(-1)
end if
if ( mbd_mlo/=0 ) then
  if ( mbd_maxscale_mlo==0 ) then
    write(6,*) "ERROR: mbd_maxscale_mlo must be >0 when mbd_mlo/=0"
    call ccmpi_abort(-1)
  end if
  mbd_min = int(20.*112.*90.*schmidt/real(mbd_maxscale_mlo))
  if ( mbd_mlo<mbd_min ) then
    if ( myid==0 ) then
      write(6,*) "Satisfy mbd_maxscale_mlo by adjusting mbd_mlo = ",mbd_min
    end if
    mbd_mlo = mbd_min
  end if
end if
! number of vertical levels in spectral nudging for MPI.
if ( kblock<0 ) then
  kblock = max(kl, ol) ! must occur before indata
  if ( myid==0 ) then
    write(6,*) "Adjusting vertical kblock = ",kblock
  end if  
end if
if ( wgcoeff<0. ) then
  tscale = max( 3600., wg_tau )
  ! Schreur et al (2008) "Theory of a TKE based parameterisation of wind gusts" HIRLAM newsletter 54.
  wgcoeff = sqrt(max(0.,2.*log((tscale/wg_tau)*(1./sqrt(2.*pi))*(1./log(1./wg_prob)))))
  if ( myid==0 ) then
    write(6,*) "Adjusting wgcoeff = ",wgcoeff
  end if
end if


! **** do namelist fixes above this line ***

!--------------------------------------------------------------
! REMAP MPI PROCESSES

! Optimise the MPI process ranks to reduce inter-node message passing
call ccmpi_remap

! the grid size is defined by il_g grid-points
!   il_g is the number of grid-points for the grid along the X-axis
!   jl_g is the number of grid-points for the grid along the Y-axis,
!     where jl_g=6*il_g due to the six panels of the cube
!   total number of horizontal grid-points for the grid is ifull_g=il_g*jl_g
! the grid is divided into node_dx*node_dy nodes
!   node_dx is the number of nodes for the grid along the X-axis
!   node_dy is the number of nodes for the grid along the Y-axis
!   Usually node_dx*node_dy is equal to the total number of physical nodes equal to
!   nodecaptian_nproc.  However if processes on a node are not fully allocated
!   then a node can be decomposed into smaller 'virtual' nodes until all processes on
!   a virtual node are fully allocated.
! each node is divided into node_nx*node_ny processes
!   nxp is the number of processes for the grid along the X-axis
!   nyp is the number of processes for the grid along the Y-axis
!   node_nx is the number of processes for a node along the X-axis,
!     where node_nx=nxp/node_dx
!   node_ny is the number of processes for a node along the Y-axis,
!     where node_ny=nyp/node_dy
!   il is the number of horizontal grid-points for a process along the X-axis,
!     where il=il_g/nxp
!   jl is the number of horizontal grid-points for a process along the Y-axis,
!     where jl=jl_g/nyp
!   total number of processes for the grid is nproc=nxp*nyp
!   total number of horizontal grid-points for a process is ifull=il*jl
! each process is divided into ntiles (only for physics and chemistry)
!   the number of grid-points per tile is imax=ifull/ntiles
! MPI routines use ipan, jpan and npan to decompose the grid on a process
!   npan is the number of panels on a process ( 1>=npan>=6 )
!   ipan=il is the numnber of grid-points for a process along the X-axis
!   jpan=jl/npan is the number of grid-points per panel for a process along the Y-axis
!   nxproc is the number of processes per panel along the X-axis (nxproc=nxp)
!   nyproc is the number of processes per panel along the Y-axis (6*nyproc=nyp)
! CCAM will optimise nxp, node_nx and npan (constrained by the number of processes,
! number of nodes and number of cubic panels, respectively) to reduce MPI message
! size and number between processes and nodes.


!--------------------------------------------------------------
! PARAMETER TESTS

if ( nextout>=4 ) then
  if ( nllp<3 ) then
    if ( myid==0 ) then
      write(6,*) "WARN: Increase nllp=3 for nextout>=4"
    end if
    nllp = 3
  end if
end if
if ( newtop>2 ) then
  write(6,*) 'newtop>2 no longer allowed'
  call ccmpi_abort(-1)
end if
if ( mfix_qg>0 .and. nkuo==4 ) then
  write(6,*) 'nkuo=4: mfix_qg>0 not allowed'
  call ccmpi_abort(-1)
end if
nstagin  = nstag    ! -ve nstagin gives swapping & its frequency
nstaguin = nstagu   ! only the sign of nstaguin matters (chooses scheme)
if ( nstagin==5 .or. nstagin<0 ) then
  nstag  = 4
  nstagu = 4
  if ( nstagin==5 ) then  ! for backward compatability
    nstagin  = -1 
    nstaguin = 5  
  endif
endif
if ( surfile /= ' ' ) then
  if ( tbave<=0 ) then
    write(6,*) "ERROR: tbave must be greater than zero"
    write(6,*) "tbave ",tbave
    call ccmpi_abort(-1)  
  end if
  if ( mod(ntau, tbave)/=0 ) then
    write(6,*) "ERROR: tave must be a factor of ntau"
    write(6,*) "ntau,tbave ",ntau,tbave
    call ccmpi_abort(-1)
  end if
end if
if ( freqfile /= ' ' ) then
  if ( tbave10<=0 ) then
    write(6,*) "ERROR: tbave10 must be greater than zero"
    write(6,*) "tbave10 ",tbave10
    call ccmpi_abort(-1)  
  end if
  if ( mod(ntau, tbave10)/=0 ) then
    write(6,*) "ERROR: tave must be a factor of ntau"
    write(6,*) "ntau,tbave10 ",ntau,tbave10
    call ccmpi_abort(-1)
  end if
end if


!--------------------------------------------------------------
! SHARED MEMORY AND FILE IO CONFIGURATION

! This is the procformat IO system where a single output file is
! written per (virtual) node
call ccmpi_procformat_init(localhist,procmode) 


!--------------------------------------------------------------
! DISPLAY NAMELIST

if ( myid==0 ) then   
  write(6,*)'Dynamics options:'
  write(6,*)'   mex   mfix  mfix_qg   mup    nh    precon' 
  write(6,'(i4,i6,i10,3i7)')mex,mfix,mfix_qg,mup,nh,precon
  write(6,*)'nritch_t ntbar  epsp    epsu   epsf   restol'
  write(6,'(i5,i7,1x,3f8.3,g9.2)')nritch_t,ntbar,epsp,epsu,epsf,restol
  write(6,*)'helmmeth mfix_aero mfix_tr'
  write(6,'(i8,i10,i8)') helmmeth,mfix_aero,mfix_tr
  write(6,*)'epsh'
  write(6,'(f8.3)') epsh
  write(6,*)'Horizontal advection/interpolation options:'
  write(6,*)' nt_adv mh_bs'
  write(6,'(i5,i7)') nt_adv,mh_bs
  write(6,*)'Horizontal wind staggering options:'
  write(6,*)'nstag nstagu'
  write(6,'(2i7)') nstag,nstagu
  write(6,*)'Horizontal mixing options:'
  write(6,*)' khdif  khor   nhor   nhorps nhorjlm'
  write(6,'(i5,11i7)') khdif,khor,nhor,nhorps,nhorjlm
  write(6,*)'Vertical mixing/physics options:'
  write(6,*)' nvmix nlocal ncvmix  lgwd' 
  write(6,'(i5,6i7)') nvmix,nlocal,ncvmix,lgwd
  write(6,*)' be   cm0  ce0  ce1  ce2  ce3  cqmix'
  write(6,'(7f5.2)') be,cm0,ce0,ce1,ce2,ce3,cqmix
  write(6,*)' ent0  ent1  entc0  dtrc0   m0    b1    b2'
  write(6,'(7f6.2)') ent0,ent1,entc0,dtrc0,m0,b1,b2
  write(6,*)' buoymeth stabmeth maxdts qcmf'
  write(6,'(2i9,f8.2,g9.2)') buoymeth,stabmeth,maxdts,qcmf
  write(6,*)'  mintke   mineps     minl     maxl'
  write(6,'(4g9.2)') mintke,mineps,minl,maxl
  write(6,*) ' tkemeth ezmin ent_min'
  write(6,'(i5,2f8.2)') tkemeth,ezmin,ent_min
  write(6,*) ' amxlsq'
  write(6,'(f8.2)') amxlsq
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
  write(6,*)' nrad  mins_rad  dt'
  write(6,'(i5,i7,f10.2)') nrad,mins_rad,dt
  write(6,*)'Radiation options B:'
  write(6,*)' nmr bpyear sw_diff_streams sw_resolution'
  write(6,'(i4,f9.2,i4," ",a5,i4)') nmr,bpyear,sw_diff_streams,sw_resolution
  write(6,*)'Radiation options C:'
  write(6,*)' liqradmethod iceradmethod carbonradmethod'
  write(6,'(3i4)') liqradmethod,iceradmethod,carbonradmethod
  write(6,*)'Aerosol options:'
  write(6,*)'  iaero ch_dust zvolcemi aeroindir'
  write(6,'(i7,g9.2,f7.2,i5)') iaero,ch_dust,zvolcemi,aeroindir
  write(6,*)'Cloud options:'
  write(6,*)'  ldr nclddia nstab_cld nrhcrit sigcll '
  write(6,'(i5,i6,2i9,1x,f8.2)') ldr,nclddia,nstab_cld,nrhcrit,sigcll
  write(6,*)'  ncloud'
  write(6,'(i5)') ncloud
  write(6,*)'Soil and canopy options:'
  write(6,*)' jalbfix nalpha nbarewet newrough nglacier nrungcm nsib  nsigmf'
  write(6,'(i5,9i8)') jalbfix,nalpha,nbarewet,newrough,nglacier,nrungcm,nsib,nsigmf
  write(6,*)' ntaft ntsea ntsur av_vmod tss_sh vmodmin  zobgin charnock chn10'
  write(6,'(i5,2i6,4f8.2,f8.3,f9.5)') ntaft,ntsea,ntsur,av_vmod,tss_sh,vmodmin,zobgin,charnock,chn10
  write(6,*)' ccycle proglai soil_struc cable_pop progvcmax fwsoil_switch cable_litter'
  write(6,'(7i7)') ccycle,proglai,soil_struc,cable_pop,progvcmax,fwsoil_switch,cable_litter
  write(6,*)' gs_switch smrf_switch strf_switch'
  write(6,'(3i7)') gs_switch,smrf_switch,strf_switch
  write(6,*)' nurban siburbanfrac'
  write(6,'(i7,f8.4)') nurban,siburbanfrac
  write(6,*)'Ocean/lake options:'
  write(6,*)' nmlo  ol      mxd   mindep minwater  ocnsmag   ocneps'
  write(6,'(i5,i4,5f9.2)') nmlo,ol,mxd,mindep,minwater,ocnsmag,ocneps
  write(6,*)' mlodiff  zomode zoseaice factchseaice otaumode'
  write(6,'(2i8,f9.6,f13.6,i8)') mlodiff,zomode,zoseaice,factchseaice,otaumode
  write(6,*)' usetide mlojacobi alphavis_seaice alphanir_seaice'
  write(6,'(2i8,2f8.4)') usetide,mlojacobi,alphavis_seaice,alphanir_seaice
  write(6,*)'River options:'
  write(6,*)' rivermd basinmd rivercoeff'
  write(6,'(2i8,g9.2)') rivermd,basinmd,rivercoeff
  write(6,*)'Nudging options:'
  write(6,*)' nbd    nud_p  nud_q  nud_t  nud_uv nud_hrs nudu_hrs kbotdav  kbotu'
  write(6,'(i5,3i7,7i8)') nbd,nud_p,nud_q,nud_t,nud_uv,nud_hrs,nudu_hrs,kbotdav,kbotu
  write(6,*)' mbd    mbd_maxscale mbd_maxgrid mbd_maxscale_mlo ktopdav kblock'
  write(6,'(i5,2i12,i16,2i8)') mbd,mbd_maxscale,mbd_maxgrid,mbd_maxscale_mlo,ktopdav,kblock
  write(6,*)' nud_sst nud_sss nud_ouv nud_sfh ktopmlo kbotmlo mloalpha'
  write(6,'(6i8,i9)') nud_sst,nud_sss,nud_ouv,nud_sfh,ktopmlo,kbotmlo,mloalpha
  write(6,*)' sigramplow sigramphigh nud_period'
  write(6,*)'Ensemble options:'
  write(6,*)' ensemble_mode ensemble_period ensemble_rsfactor'
  write(6,'(2i5,f8.4)') ensemble_mode,ensemble_period,ensemble_rsfactor
  write(6,'(2f10.6,i9)') sigramplow,sigramphigh,nud_period
  write(6,*)'Special and test options A:'
  write(6,*)' namip amipo3 newtop nhstest nsemble nspecial panfg panzo'
  write(6,'(1i5,L7,3i7,i8,f9.1,f8.4)') namip,amipo3,newtop,nhstest,nsemble,nspecial,panfg,panzo
  write(6,*)'Special and test options B:'
  write(6,*)' knh rescrn maxtilesize'
  write(6,'(i4,2i7)') knh,rescrn,maxtilesize
  write(6,*)'I/O options:'
  write(6,*)' m_fly  io_in io_nest io_out io_rest  nwt  nperavg'
  write(6,'(i5,4i7,3i8)') m_fly,io_in,io_nest,io_out,io_rest,nwt,nperavg
  write(6,*)' hp_output procmode compression'
  write(6,'(i5,2i5)') hp_output,procmode,compression

  write(6, cardin)
  write(6, skyin)
  write(6, datafile)
  write(6, kuonml)
  write(6, turbnml)
  write(6, landnml)
  write(6, mlonml)
end if ! myid=0


!--------------------------------------------------------------
! INITIALISE ifull_g ALLOCATABLE ARRAYS

#ifdef share_ifullg
! Allocate xx4, yy4, em_g, x_g, y_g and z_g as shared
! memory within a node.  The node captain is responsible
! for updating these arrays.
shsize(1:2) = (/ iquad, iquad /)
call ccmpi_allocshdatar8(xx4,shsize(1:2),xx4_win)
call ccmpi_allocshdatar8(yy4,shsize(1:2),yy4_win)
shsize(1) = ifull_g
call ccmpi_allocshdata(em_g,shsize(1:1),em_g_win)
call ccmpi_allocshdatar8(x_g,shsize(1:1),x_g_win)
call ccmpi_allocshdatar8(y_g,shsize(1:1),y_g_win)
call ccmpi_allocshdatar8(z_g,shsize(1:1),z_g_win)
#else
! Allocate xx4, yy4, em_g, x_g, y_g and z_g for each process
allocate( xx4(iquad,iquad), yy4(iquad,iquad) )
allocate( em_g(ifull_g) )
allocate( x_g(ifull_g), y_g(ifull_g), z_g(ifull_g) )
#endif
call xyzinfo_init(ifull_g,ifull,iextra,myid)
call map_init(ifull_g,ifull,iextra,myid)
call latlong_init(ifull_g,ifull,myid)      
call vecsuv_init(ifull_g,ifull,iextra,myid)
call workglob_init(ifull_g,ifull,myid)
call indices_init(ifull,npan)


!--------------------------------------------------------------
! SET UP CC GEOMETRY

! Only one process calls setxyz to save memory with large grids
if ( myid==0 ) then
  write(6,*) "Calling setxyz"
  call setxyz(il_g,rlong0,rlat0,schmidt,x_g,y_g,z_g,wts_g,ax_g,ay_g,az_g,bx_g,by_g,bz_g,xx4,yy4, &
              id,jd,ktau,ds)
end if
! Broadcast the following global data
! xx4 and yy4 are used for calculating depature points
! em_g, x_g, y_g and z_g are for the scale-selective filter (1D and 2D versions)
#ifdef share_ifullg
if ( myid==0 ) then
  write(6,*) "Update global arrays with shared memory"
end if  
if ( node_myid==0 ) then
  call ccmpi_bcastr8(xx4,0,comm_nodecaptain)
  call ccmpi_bcastr8(yy4,0,comm_nodecaptain)
  call ccmpi_bcast(em_g,0,comm_nodecaptain)
  call ccmpi_bcastr8(x_g,0,comm_nodecaptain)
  call ccmpi_bcastr8(y_g,0,comm_nodecaptain)
  call ccmpi_bcastr8(z_g,0,comm_nodecaptain)
end if
call ccmpi_barrier(comm_node)
#else
if ( myid==0 ) then
  write(6,*) "Update global arrays"
end if
! make copies of global arrays on all processes
call ccmpi_bcastr8(xx4,0,comm_world)
call ccmpi_bcastr8(yy4,0,comm_world)
call ccmpi_bcast(em_g,0,comm_world)
call ccmpi_bcastr8(x_g,0,comm_world)
call ccmpi_bcastr8(y_g,0,comm_world)
call ccmpi_bcastr8(z_g,0,comm_world)
#endif
call ccmpi_bcast(ds,0,comm_world)

if ( myid==0 ) then
  write(6,*) "Calling ccmpi_setup"
end if
call ccmpi_setup(id,jd,idjd,dt)

!--------------------------------------------------------------
! DEALLOCATE ifull_g ARRAYS WHERE POSSIBLE
if ( myid==0 ) then
  deallocate( wts_g, emu_g, emv_g )
  deallocate( ax_g, ay_g, az_g )
  deallocate( bx_g, by_g, bz_g )
  deallocate( f_g, fu_g, fv_g )
  deallocate( rlatt_g, rlongg_g )
end if


!--------------------------------------------------------------
! INITIALISE LOCAL ARRAYS
allocate( dums(ifull,kl) )
call arrays_init(ifull,iextra,kl)
call carbpools_init(ifull,nsib,ccycle)
call cfrac_init(ifull,iextra,kl,ncloud)
call dpsdt_init(ifull,epsp)
call epst_init(ifull)
call extraout_init(ifull,nextout)
call gdrag_init(ifull)
call histave_init(ifull,kl,ms,ccycle,output_windmax)
call kuocom_init(ifull,kl)
call liqwpar_init(ifull,iextra,kl,process_rate_mode)
call morepbl_init(ifull,kl)
call nharrs_init(ifull,iextra,kl)
call nlin_init(ifull,kl)
call nsibd_init(ifull,nsib)
call parmhdff_init(kl)
call pbl_init(ifull)
call permsurf_init(ifull)
call prec_init(ifull)
call raddiag_init(ifull,kl)
call riverarrays_init(ifull,iextra)
call savuvt_init(ifull,kl)
call savuv1_init(ifull,kl)
call sbar_init(ifull,kl)
call screen_init(ifull)
call sigs_init(kl)
call soil_init(ifull,iaero,nsib)
call soilsnow_init(ifull,ms,nsib)
call tbar2d_init(ifull)
call unn_init(ifull,kl)
call uvbar_init(ifull,kl)
call vecs_init(kl)
call vegpar_init(ifull)
call vvel_init(ifull,iextra,kl)
call work2_init(ifull,nsib)
call work3_init(ifull,nsib)
call work3f_init(ifull,kl)
call xarrs_init(ifull,iextra,kl)
if ( nvmix==6 .or. nvmix==9 ) then
  call tkeinit(ifull,iextra,kl)
end if
call init_tracer
call work3sav_init(ifull,kl,ntrac) ! must occur after tracers_init
if ( nbd/=0 .or. mbd/=0 ) then
  if ( abs(iaero)>=2 .and. nud_aero/=0 ) then
    call dav_init(ifull,kl,naero,nbd,mbd)
  else
    call dav_init(ifull,kl,0,nbd,mbd)
  end if
end if
! Remaining arrays are allocated in indata.f90, since their
! dimension size requires additional input data (e.g, land-surface)
 
!--------------------------------------------------------------
! DISPLAY DIAGNOSTIC INDEX AND TIMER DATA
if ( mydiag ) then
  write(6,"(' id,jd,rlongg,rlatt in degrees: ',2i4,2f8.2)") id,jd,180./pi*rlongg(idjd),180./pi*rlatt(idjd)
end if
call date_and_time(rundate)
call date_and_time(time=timeval)
if ( myid==0 ) then
  write(6,*)'RUNDATE IS ',rundate
  write(6,*)'Starting time ',timeval
end if


!--------------------------------------------------------------
! READ INITIAL CONDITIONS
if ( myid==0 ) then
  write(6,*) "Calling indata"
end if
call indataf(lapsbot,isoth,nsig,nmlfile)


!--------------------------------------------------------------
! SETUP REMAINING PARAMETERS
if ( myid==0 ) then
  write(6,*) "Setup remaining parameters"
end if
  
! fix nudging levels from pressure to level index
! this is done after indata has loaded sig
if ( kbotdav<0 ) then
  targetlev = real(-kbotdav)/1000.
  do k = 1,kl
    if ( sig(k)<=targetlev ) then
      kbotdav = k
      if ( myid==0 ) then
        write(6,*) "Nesting kbotdav adjusted to ",kbotdav," for sig ",sig(kbotdav)
      end if
      exit
    end if
  end do
  if ( kbotdav<0 ) then
    write(6,*) "ERROR: Cannot locate nudging level for kbotdav ",kbotdav
    call ccmpi_abort(-1)
  end if
end if
if ( ktopdav==0 ) then
  ktopdav = kl
else if ( ktopdav<0 ) then
  targetlev = real(-ktopdav)/1000.
  do k = kl,1,-1
    if ( sig(k)>=targetlev ) then
      ktopdav = k
      if ( myid == 0 ) then
        write(6,*) "Nesting ktopdav adjusted to ",ktopdav," for sig ",sig(ktopdav)
      end if
      exit
    end if
  end do
  if ( ktopdav<0 ) then
    write(6,*) "ERROR: Cannot locate nudging level for ktopdav ",ktopdav
    call ccmpi_abort(-1)
  end if
end if
if ( kbotdav<1 .or. ktopdav>kl .or. kbotdav>ktopdav ) then
  write(6,*) "ERROR: Invalid kbotdav and ktopdav"
  write(6,*) "kbotdav,ktopdav ",kbotdav,ktopdav
  call ccmpi_abort(-1)
end if
if ( kbotu==0 ) kbotu = kbotdav

! fix ocean nuding levels
if ( nmlo/=0 .and. abs(nmlo)<9 ) then
  allocate( gosig_in(ol) )
  call mlovlevels(gosig_in,sigma=.true.)
  if ( kbotmlo<0 )  then
    targetlev = real(-kbotmlo)/1000.   
    do k = ol,1,-1
      if ( gosig_in(k)<=targetlev ) then
        kbotmlo = k
        if ( myid==0 ) then
          write(6,*) "Nesting kbotmlo adjusted to ",kbotmlo," for sig ",gosig_in(kbotmlo)
        end if
        exit
      end if
    end do
    if ( kbotmlo<0 ) then
      write(6,*) "ERROR: Cannot locate nudging level for kbotmlo ",kbotmlo
      call ccmpi_abort(-1)
    end if   
  end if
  if ( ktopmlo<0 ) then
    targetlev = real(-ktopmlo)/1000.
    do k = 1,ol
      if ( gosig_in(k)>=targetlev ) then
        ktopmlo = k
        if ( myid==0 ) then
          write(6,*) "Nesting ktopmlo adjusted to ",ktopmlo," for sig ",gosig_in(ktopmlo)
        end if
        exit
      end if
    end do
    if ( ktopmlo<0 ) then
      write(6,*) "ERROR: Cannot locate nudging level for ktopmlo ",ktopmlo
      call ccmpi_abort(-1)
    end if
  end if
  if ( ktopmlo<1 .or. kbotmlo>ol .or. ktopmlo>kbotmlo ) then
    write(6,*) "ERROR: Invalid kbotmlo"
    write(6,*) "kbotmlo,ktopmlo ",kbotmlo,ktopmlo
    call ccmpi_abort(-1)
  end if
  deallocate(gosig_in)
end if  

! identify reference level ntbar for temperature
if ( ntbar==-1 ) then
  ntbar = 1
  do while( sig(ntbar)>0.8 .and. ntbar<kl )
    ntbar = ntbar + 1
  end do
end if

! estimate radiation calling frequency
if ( mins_rad<0 ) then
  ! automatic estimate for mins_rad
  secs_rad = min(nint((schmidt*112.*90./real(il_g))*8.*60.), nint(real(nwt)*dt), 3600)
  kountr   = max(nint(real(secs_rad)/dt), 1)
  secs_rad = nint(real(kountr)*dt)
  do while ( (mod(3600, secs_rad)/=0 .or. mod(nint(real(nwt)*dt), secs_rad)/=0) .and. kountr>1 )
    kountr = kountr - 1
    secs_rad = nint(real(kountr)*dt)
  end do
else
  ! user specified mins_rad
  kountr   = nint(real(mins_rad)*60./dt)  ! set default radiation to ~mins_rad m
  secs_rad = nint(real(kountr)*dt)        ! redefine to actual value
end if
if ( myid==0 ) then
  write(6,*) "Radiation will use kountr ",kountr," for secs_rad ",secs_rad
end if
! for 6-hourly output of sint_ave etc, want 6*60*60 = N*secs_rad      
if ( (nrad==4.or.nrad==5) .and. mod(21600,secs_rad)/=0 ) then
  write(6,*) 'ERROR: CCAM would prefer 21600 = N*secs_rad ',secs_rad
  call ccmpi_abort(-1)
end if

! max/min diagnostics      
if ( nextout>=4 ) call setllp

if ( nmaxpr<=ntau ) then
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
    pwatr_l = pwatr_l - sum(dsig(k)*wts(1:ifull)*(qg(1:ifull,k)+qlg(1:ifull,k)+qfg(1:ifull,k))*ps(1:ifull))
  enddo
  pwatr_l = pwatr_l/grav
  temparray(1) = pwatr_l
  call ccmpi_reduce( temparray(1:1), gtemparray(1:1), "sum", 0, comm_world )
  pwatr = gtemparray(1)
  if ( myid==0 ) write (6,"('pwatr0 ',12f7.3)") pwatr
  if ( ntrac>0 ) then
    do ng = 1,ntrac
      write (text,'("g",i1)')ng
      call maxmin(tr(:,:,ng),text,ktau,1.,kl)
    end do
  end if   ! (ntrac>0)
end if  

! convection ( for vertmix (nvmix==3) and radriv90 )
! sig(kuocb) occurs for level just BELOW sigcb
kuocb = 1
do while( sig(kuocb+1)>=sigcb )
  kuocb = kuocb + 1
end do
if ( myid==0 ) write(6,*) 'convective cumulus scheme: kuocb,sigcb = ',kuocb,sigcb

! horizontal diffusion 
if ( khdif==-99 ) then   ! set default khdif appropriate to resolution
  khdif = 5
  if ( myid==0 ) write(6,*) 'Model has chosen khdif =',khdif
endif
do k = 1,kl
  hdiff(k) = khdif*0.1
end do
if ( khor>0 ) then
  do k = kl+1-khor,kl
    hdiff(k) = 2.*hdiff(k-1)
  end do
elseif ( khor<0 ) then ! following needed +hdiff() (JLM 29/6/15)
  do k = 1,kl          ! N.B. usually hdiff(k)=khdif*.1 
    ! increase hdiff between sigma=.15  and sigma=0., 0 to khor
    if ( sig(k)<0.15 ) then
      hdiff(k) = .1*max(1.,(1.-sig(k)/.15)*abs(khor)) + hdiff(k)
    end if
  end do
  if ( myid==0 ) write(6,*)'khor,hdiff: ',khor,hdiff
end if

! nudging
! nud_period=-1 uses input data time-step (e.g., 6-hourly)
! nud_period=0 uses resolution dependent time-step
if ( nud_period==0 ) then
  ! automatic estimate for nudging period
  secs_nud = min(nint((schmidt*112.*90./real(il_g))*8.*60.), nint(real(nwt)*dt), 3600)
  nud_period = max(nint(real(secs_nud)/60.),1)
  do while ( (mod(60, nud_period)/=0 .or. mod(nint(real(nwt)*dt/60.), nud_period)/=0 .or. &
              mod(nud_period*60, nint(dt))/=0 ) .and. nud_period>1 )
    nud_period = nud_period - 1
  end do
end if
if ( myid==0 ) then
  write(6,*) "Nudging will use nud_period ",nud_period
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
  write(6,*) "ERROR: ngas>0 and mfix_tr=0"
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
if ( myid==0 ) then
  open(11, file='nrun.dat', status='unknown')
  if ( nrun==0 ) then
    read(11,*,iostat=ierr) nrun
    nrun = nrun + 1
  end if   ! nrun==0
  write(6,*) 'this is run ',nrun
  rewind 11
  write(11,*) nrun
  write(11,cardin)
  write(11,skyin)
  write(11,datafile)
  write(11,kuonml)
  write(11,turbnml)
  write(11,landnml)
  write(11,mlonml)
  close(11)
end if

deallocate( dums )
  
return
end subroutine globpe_init
    
!--------------------------------------------------------------
! PREVIOUS VERSION DEFAULT PARAMETERS
subroutine change_defaults(nversion)

use kuocom_m                ! JLM convection
use newmpar_m               ! Grid parameters
use parm_m                  ! Model configuration
use parmdyn_m               ! Dynamics parmaters
use parmhor_m               ! Horizontal advection parameters
use parmhdff_m              ! Horizontal diffusion parameters

implicit none

integer, intent(in) :: nversion

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
  nbase = 0        ! new is 1
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
  mins_rad = 120   ! new is 72
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
  jalbfix = 0      ! new is 1
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
! Broadcast cardin namelist
subroutine broadcast_cardin

use aerointerface, only : ch_dust        & ! Aerosol arrays
    ,zvolcemi,so4mtn,carbmtn             &
    ,saltsmallmtn,saltlargemtn           &
    ,enhanceu10
use cc_acc                                 ! CC ACC routines
use cc_mpi                                 ! CC MPI routines
use estab                                  ! Liquid saturation function
use indata                                 ! Data initialisation
use infile                                 ! Input file routines
use kuocom_m                               ! JLM convection
use module_ctrl_microphysics
use newmpar_m                              ! Grid parameters
use nharrs_m                               ! Non-hydrostatic atmosphere arrays
use parm_m                                 ! Model configuration
use parmdyn_m                              ! Dynamics parameters
use parmgeom_m                             ! Coordinate data
use parmhdff_m                             ! Horizontal diffusion parameters
use parmhor_m                              ! Horizontal advection parameters
use parmvert_m                             ! Vertical advection parameters
use staguvmod                              ! Reversible grid staggering 
use stime_m                                ! File date data

implicit none

integer i
integer, dimension(123) :: dumi
real, dimension(34) :: dumr
    
dumr(:) = 0.
dumi(:) = 0
if ( myid==0 ) then
  dumr(1)   = dt
  dumr(2)   = restol
  dumr(3)   = panfg
  dumr(4)   = panzo
  dumr(5)   = rlatdn
  dumr(6)   = rlatdx
  dumr(7)   = rlongdn
  dumr(8)   = rlongdx
  dumr(9)   = epsp
  dumr(10)  = epsu
  dumr(11)  = epsf
  dumr(12)  = epsh
  dumr(13)  = av_vmod
  dumr(14)  = charnock
  dumr(15)  = chn10
  dumr(16)  = snmin
  dumr(17)  = tss_sh
  dumr(18)  = vmodmin
  dumr(19)  = zobgin
  dumr(20)  = rlong0
  dumr(21)  = rlat0
  dumr(22)  = schmidt
  dumr(23)  = sigramplow
  dumr(24)  = sigramphigh
  dumr(25)  = ch_dust
  dumr(26)  = helim
  dumr(27)  = fc2
  dumr(28)  = sigbot_gwd
  dumr(29)  = alphaj
  dumr(30)  = qgmin
  dumr(31)  = rhsat
  dumr(32)  = ensemble_rsfactor
  dumr(33)  = zo_clearing
  dumr(34)  = maxuv
  dumi(1)   = ntau
  dumi(2)   = nwt
  dumi(3)   = nhorps
  dumi(4)   = nperavg
  dumi(5)   = ia
  dumi(6)   = ib
  dumi(7)   = ja
  dumi(8)   = jb
  dumi(9)   = id
  dumi(10)  = jd
  dumi(11)  = iaero
  dumi(12)  = khdif
  dumi(13)  = khor
  dumi(14)  = nhorjlm
  dumi(15)  = mex
  dumi(16)  = mbd
  dumi(17)  = nbd
  dumi(18)  = mbd_maxscale
  dumi(19)  = mbd_maxgrid
  dumi(20)  = ndi
  dumi(21)  = ndi2
  dumi(22)  = nhor
  dumi(23)  = nlv
  dumi(24)  = nmaxpr
  dumi(25)  = nrad
  dumi(26)  = ntaft
  dumi(27)  = ntsea
  dumi(28)  = ntsur
  dumi(29)  = nvmix
  dumi(30)  = precon
  dumi(31)  = kdate_s
  dumi(32)  = ktime_s
  dumi(33)  = leap
  dumi(34)  = newtop
  dumi(35)  = mup
  dumi(36)  = lgwd
  dumi(37)  = ngwd
  dumi(38)  = nextout
  dumi(39)  = jalbfix
  dumi(40)  = nalpha
  dumi(41)  = nstag
  dumi(42)  = nstagu
  dumi(43)  = ntbar
  dumi(44)  = nwrite
  dumi(45)  = irest
  dumi(46)  = nrun
  dumi(47)  = nstn
  dumi(48)  = nrungcm
  dumi(49)  = nsib
  dumi(50)  = mh_bs
  dumi(51)  = nritch_t
  dumi(52)  = nt_adv
  dumi(53)  = mfix
  dumi(54)  = mfix_qg
  dumi(55)  = namip
  if ( amipo3 ) dumi(56) = 1
  dumi(57)  = nh
  dumi(58)  = nhstest
  dumi(59)  = nsemble
  dumi(60)  = nspecial
  dumi(61)  = newrough
  dumi(62)  = nglacier
  dumi(63)  = newztsea
  dumi(64)  = kbotdav
  dumi(65)  = kbotu
  dumi(66)  = nud_p
  dumi(67)  = nud_q
  dumi(68)  = nud_t
  dumi(69)  = nud_uv
  dumi(70)  = nud_hrs
  dumi(71)  = nudu_hrs
  dumi(72)  = nlocal
  dumi(73)  = nbarewet
  dumi(74)  = nsigmf
  dumi(75)  = io_in
  dumi(76)  = io_nest
  dumi(77)  = io_out
  dumi(78)  = io_rest
  dumi(79)  = tbave
  if ( synchist ) dumi(80) = 1
  dumi(81)  = m_fly
  dumi(82)  = nurban
  dumi(83)  = ktopdav
  dumi(84)  = mbd_mlo
  dumi(85)  = mbd_maxscale_mlo
  dumi(86)  = nud_sst
  dumi(87)  = nud_sss
  dumi(88)  = mfix_tr
  dumi(89)  = mfix_aero
  dumi(90)  = kbotmlo
  dumi(91)  = ktopmlo
  dumi(92)  = mloalpha
  dumi(93)  = nud_ouv
  dumi(94)  = nud_sfh
  dumi(95)  = rescrn
  dumi(96) = helmmeth
  dumi(97) = nmlo
  dumi(98) = ol
  dumi(99) = knh
  dumi(100) = kblock
  dumi(101) = nud_aero
  dumi(102) = nud_period
  dumi(103) = procmode
  dumi(104) = compression
  dumi(105) = nmr
  dumi(106) = maxtilesize
  dumi(107) = mfix_t
  dumi(108) = ensemble_mode
  dumi(109) = ensemble_period
  dumi(110) = hp_output
  dumi(111) = intsch_mode
  dumi(112) = qg_fix
  if ( always_mspeca ) dumi(113) = 1
  dumi(114) = ntvd
  dumi(115) = tbave10  
  dumi(116) = async_length
  dumi(117) = nagg
  dumi(118) = pil_single
  if ( localhist ) dumi(119) = 1
  dumi(120) = maxcolour
  dumi(121) = process_rate_mode
  dumi(122) = chunk_time
  dumi(123) = estab_bug_fix
end if
call ccmpi_bcast(dumr,0,comm_world)
call ccmpi_bcast(dumi,0,comm_world)
dt                = dumr(1)
restol            = dumr(2)
panfg             = dumr(3)
panzo             = dumr(4)
rlatdn            = dumr(5)
rlatdx            = dumr(6)
rlongdn           = dumr(7)
rlongdx           = dumr(8)
epsp              = dumr(9)
epsu              = dumr(10)
epsf              = dumr(11)
epsh              = dumr(12)
av_vmod           = dumr(13)
charnock          = dumr(14)
chn10             = dumr(15)
snmin             = dumr(16)
tss_sh            = dumr(17)
vmodmin           = dumr(18)
zobgin            = dumr(19)
rlong0            = dumr(20)
rlat0             = dumr(21)
schmidt           = dumr(22)
sigramplow        = dumr(23)
sigramphigh       = dumr(24)
ch_dust           = dumr(25)
helim             = dumr(26)
fc2               = dumr(27)
sigbot_gwd        = dumr(28)
alphaj            = dumr(29)
qgmin             = dumr(30)
rhsat             = dumr(31)
ensemble_rsfactor = dumr(32)
zo_clearing       = dumr(33)
maxuv             = dumr(34)
ntau              = dumi(1)
nwt               = dumi(2)
nhorps            = dumi(3)
nperavg           = dumi(4)
ia                = dumi(5)
ib                = dumi(6)
ja                = dumi(7)
jb                = dumi(8)
id                = dumi(9)
jd                = dumi(10)
iaero             = dumi(11)
khdif             = dumi(12)
khor              = dumi(13)
nhorjlm           = dumi(14)
mex               = dumi(15)
mbd               = dumi(16)
nbd               = dumi(17)
mbd_maxscale      = dumi(18)
mbd_maxgrid       = dumi(19)
ndi               = dumi(20)
ndi2              = dumi(21)
nhor              = dumi(22)
nlv               = dumi(23)
nmaxpr            = dumi(24)
nrad              = dumi(25)
ntaft             = dumi(26)
ntsea             = dumi(27)
ntsur             = dumi(28)
nvmix             = dumi(29)
precon            = dumi(30)
kdate_s           = dumi(31)
ktime_s           = dumi(32)
leap              = dumi(33)
newtop            = dumi(34)
mup               = dumi(35)
lgwd              = dumi(36)
ngwd              = dumi(37)
nextout           = dumi(38)
jalbfix           = dumi(39)
nalpha            = dumi(40)
nstag             = dumi(41)
nstagu            = dumi(42)
ntbar             = dumi(43)
nwrite            = dumi(44)
irest             = dumi(45)
nrun              = dumi(46)
nstn              = dumi(47)
nrungcm           = dumi(48)
nsib              = dumi(49)
mh_bs             = dumi(50)
nritch_t          = dumi(51)
nt_adv            = dumi(52)
mfix              = dumi(53)
mfix_qg           = dumi(54)
namip             = dumi(55)
amipo3            = dumi(56)==1
nh                = dumi(57)
nhstest           = dumi(58)
nsemble           = dumi(59)
nspecial          = dumi(60)
newrough          = dumi(61)
nglacier          = dumi(62)
newztsea          = dumi(63)
kbotdav           = dumi(64)
kbotu             = dumi(65)
nud_p             = dumi(66)
nud_q             = dumi(67)
nud_t             = dumi(68)
nud_uv            = dumi(69)
nud_hrs           = dumi(70)
nudu_hrs          = dumi(71)
nlocal            = dumi(72)
nbarewet          = dumi(73)
nsigmf            = dumi(74)
io_in             = dumi(75)
io_nest           = dumi(76)
io_out            = dumi(77)
io_rest           = dumi(78)
tbave             = dumi(79)
synchist          = dumi(80)==1
m_fly             = dumi(81)
nurban            = dumi(82)
ktopdav           = dumi(83)
mbd_mlo           = dumi(84)
mbd_maxscale_mlo  = dumi(85)
nud_sst           = dumi(86)
nud_sss           = dumi(87)
mfix_tr           = dumi(88)
mfix_aero         = dumi(89)
kbotmlo           = dumi(90)
ktopmlo           = dumi(91)
mloalpha          = dumi(92)
nud_ouv           = dumi(93)
nud_sfh           = dumi(94)
rescrn            = dumi(95)
helmmeth          = dumi(96)
nmlo              = dumi(97)
ol                = dumi(98)
knh               = dumi(99)
kblock            = dumi(100)
nud_aero          = dumi(101)
nud_period        = dumi(102)
procmode          = dumi(103)
compression       = dumi(104)
nmr               = dumi(105)
maxtilesize       = dumi(106)
mfix_t            = dumi(107)
ensemble_mode     = dumi(108)
ensemble_period   = dumi(109)
hp_output         = dumi(110)
intsch_mode       = dumi(111)
qg_fix            = dumi(112)
always_mspeca     = dumi(113)==1
ntvd              = dumi(114)
tbave10           = dumi(115)
async_length      = dumi(116)
nagg              = dumi(117)
pil_single        = dumi(118)
localhist         = dumi(119)==1
maxcolour         = dumi(120)
process_rate_mode = dumi(121)
chunk_time        = dumi(122)
estab_bug_fix     = dumi(123)
if ( nstn>0 ) then
  call ccmpi_bcast(istn(1:nstn),0,comm_world)
  call ccmpi_bcast(jstn(1:nstn),0,comm_world)
  call ccmpi_bcast(iunp(1:nstn),0,comm_world)
  call ccmpi_bcast(slat(1:nstn),0,comm_world)
  call ccmpi_bcast(slon(1:nstn),0,comm_world)
  call ccmpi_bcast(zstn(1:nstn),0,comm_world)
  call ccmpi_bcast(name_stn(1:nstn),0,comm_world)
end if

return
end subroutine broadcast_cardin

!--------------------------------------------------------------
! Broadcast skyin namelist
subroutine broadcast_skyin

use aerointerface, only : aeroindir      & ! Aerosol interface
    ,aero_split,aerosol_u10              &
    ,ch_dust,zvolcemi                    &
    ,so4mtn,carbmtn,saltsmallmtn         &
    ,saltlargemtn,enhanceu10
use cc_mpi                                 ! CC MPI routines
use module_aux_rad                         ! Additional cloud and radiation routines
use ozoneread                              ! Ozone input routines
use parm_m                                 ! Model configuration
use seaesfrad_m                            ! SEA-ESF radiation

implicit none

integer, dimension(17) :: dumi
real, dimension(10) :: dumr
    
dumr = 0.
dumi = 0
if ( myid==0 ) then
  dumr(1)  = bpyear
  dumr(2)  = qgmin
  dumr(3)  = ch_dust
  dumr(4)  = zvolcemi
  dumr(5)  = so4mtn
  dumr(6)  = carbmtn
  dumr(7)  = saltsmallmtn
  dumr(8)  = saltlargemtn
  dumr(9)  = siglow
  dumr(10) = sigmid
  dumi(1)  = mins_rad
  dumi(2)  = liqradmethod
  dumi(3)  = iceradmethod
  dumi(4)  = so4radmethod
  dumi(5)  = carbonradmethod
  dumi(6)  = dustradmethod
  dumi(7)  = seasaltradmethod
  dumi(8)  = aeroindir
  dumi(9)  = o3_vert_interpolate
  if ( do_co2_10um ) dumi(10) = 1
  dumi(11) = aerosol_u10  
  dumi(12) = aero_split
  dumi(13) = enhanceu10
  if ( do_quench ) dumi(14) = 1
  if ( remain_rayleigh_bug ) dumi(15) = 1
  if ( use_rad_year ) dumi(16) = 1
  dumi(17) = rad_year
end if
call ccmpi_bcast(dumr,0,comm_world)
call ccmpi_bcast(dumi,0,comm_world)
call ccmpi_bcast(sw_resolution,0,comm_world)
call ccmpi_bcast(lwem_form,0,comm_world)
call ccmpi_bcast(linecatalog_form,0,comm_world)
call ccmpi_bcast(continuum_form,0,comm_world)
bpyear              = dumr(1)
qgmin               = dumr(2)
ch_dust             = dumr(3)
zvolcemi            = dumr(4)
so4mtn              = dumi(5)
carbmtn             = dumr(6)
saltsmallmtn        = dumr(7)
saltlargemtn        = dumr(8)
siglow              = dumr(9)
sigmid              = dumr(10)
mins_rad            = dumi(1)
liqradmethod        = dumi(2)
iceradmethod        = dumi(3)
so4radmethod        = dumi(4)
carbonradmethod     = dumi(5)
dustradmethod       = dumi(6)
seasaltradmethod    = dumi(7)
aeroindir           = dumi(8)
o3_vert_interpolate = dumi(9)
do_co2_10um         = dumi(10)==1
aerosol_u10         = dumi(11)
aero_split          = dumi(12)
enhanceu10          = dumi(13)
do_quench           = dumi(14)==1
remain_rayleigh_bug = dumi(15)==1
use_rad_year        = dumi(16)==1
rad_year            = dumi(17)

return
end subroutine broadcast_skyin

!--------------------------------------------------------------
! Broadcast datafile namelist
subroutine broadcast_datafile

use cc_mpi                                 ! CC MPI routines
use filnames_m                             ! Filenames
use parm_m                                 ! Model configuration

implicit none

integer, dimension(25) :: dumi
    
dumi = 0
if ( myid==0 ) then
  if ( save_aerosols ) dumi(1)=1
  if ( save_pbl ) dumi(2)=1
  if ( save_cloud ) dumi(3)=1
  if ( save_land ) dumi(4)=1
  if ( save_maxmin ) dumi(5)=1
  if ( save_ocean ) dumi(6)=1
  if ( save_radiation ) dumi(7)=1
  if ( save_urban ) dumi(8)=1
  if ( save_carbon ) dumi(9)=1
  if ( save_river ) dumi(10)=1
  dumi(11) = diaglevel_aerosols
  dumi(12) = diaglevel_pbl
  dumi(13) = diaglevel_cloud
  dumi(14) = diaglevel_land
  dumi(15) = diaglevel_maxmin
  dumi(16) = diaglevel_ocean
  dumi(17) = diaglevel_radiation
  dumi(18) = diaglevel_urban
  dumi(19) = diaglevel_carbon
  dumi(20) = diaglevel_river
  dumi(21) = diaglevel_pop
  dumi(22) = surf_cordex
  dumi(23) = output_windmax
  dumi(24) = cordex_fix
  dumi(25) = shep_cordex
end if
call ccmpi_bcast(dumi,0,comm_world)
call ccmpi_bcast(ifile,0,comm_world)
call ccmpi_bcast(ofile,0,comm_world)
call ccmpi_bcast(mesonest,0,comm_world)
call ccmpi_bcast(restfile,0,comm_world)
call ccmpi_bcast(surfile,0,comm_world)
call ccmpi_bcast(surf_00,0,comm_world)
call ccmpi_bcast(surf_12,0,comm_world)
call ccmpi_bcast(cnsdir,0,comm_world)
call ccmpi_bcast(ensembleoutfile,0,comm_world)
!call ccmpi_bcast(albfile,0,comm_world)
!call ccmpi_bcast(eigenv,0,comm_world)
!call ccmpi_bcast(icefile,0,comm_world)
!call ccmpi_bcast(rsmfile,0,comm_world)
!call ccmpi_bcast(so4tfile,0,comm_world)
!call ccmpi_bcast(soilfile,0,comm_world)
!call ccmpi_bcast(sstfile,0,comm_world)
!call ccmpi_bcast(topofile,0,comm_world)
!call ccmpi_bcast(vegfile,0,comm_world)
!call ccmpi_bcast(zofile,0,comm_world)
!call ccmpi_bcast(laifile,0,comm_world)
!call ccmpi_bcast(albnirfile,0,comm_world)
!call ccmpi_bcast(urbanfile,0,comm_world)
!call ccmpi_bcast(bathfile,0,comm_world)
!call ccmpi_bcast(salfile,0,comm_world)
!call ccmpi_bcast(oxidantfile,0,comm_world)
!call ccmpi_bcast(casafile,0,comm_world)
!call ccmpi_bcast(phenfile,0,comm_world)
call ccmpi_bcast(solarfile,0,comm_world)
call ccmpi_bcast(radfile,0,comm_world)
call ccmpi_bcast(ch4file,0,comm_world)
call ccmpi_bcast(n2ofile,0,comm_world)
call ccmpi_bcast(cfc11file,0,comm_world)
call ccmpi_bcast(cfc12file,0,comm_world)
call ccmpi_bcast(cfc113file,0,comm_world)
call ccmpi_bcast(hcfc22file,0,comm_world)
!call ccmpi_bcast(o3file,0,comm_world)
call ccmpi_bcast(freqfile,0,comm_world)
call ccmpi_bcast(wbclimfile,0,comm_world)
save_aerosols  = dumi(1)==1
save_pbl       = dumi(2)==1
save_cloud     = dumi(3)==1
save_land      = dumi(4)==1
save_maxmin    = dumi(5)==1
save_ocean     = dumi(6)==1
save_radiation = dumi(7)==1
save_urban     = dumi(8)==1
save_carbon    = dumi(9)==1
save_river     = dumi(10)==1
diaglevel_aerosols  = dumi(11)
diaglevel_pbl       = dumi(12)
diaglevel_cloud     = dumi(13)
diaglevel_land      = dumi(14)
diaglevel_maxmin    = dumi(15)
diaglevel_ocean     = dumi(16)
diaglevel_radiation = dumi(17)
diaglevel_urban     = dumi(18)
diaglevel_carbon    = dumi(19)
diaglevel_river     = dumi(20)
diaglevel_pop       = dumi(21)
surf_cordex         = dumi(22)
output_windmax      = dumi(23)
cordex_fix          = dumi(24)
shep_cordex         = dumi(25)

return
end subroutine broadcast_datafile

!--------------------------------------------------------------
! Broadcast kuonml namelist
subroutine broadcast_kuonml

use cc_mpi                                 ! CC MPI routines
use kuocom_m                               ! JLM convection
use module_ctrl_microphysics               ! Interface for cloud microphysics
use parm_m                                 ! Model configuration

implicit none

integer, dimension(29) :: dumi
real, dimension(37) :: dumr
    
dumr = 0.
dumi = 0
if ( myid==0 ) then
  dumr(1)  = alflnd
  dumr(2)  = alfsea
  dumr(3)  = cldh_lnd
  dumr(4)  = cldm_lnd
  dumr(5)  = cldl_lnd
  dumr(6)  = cldh_sea
  dumr(7)  = cldm_sea
  dumr(8)  = cldl_sea
  dumr(9)  = convfact
  dumr(10) = convtime
  dumr(11) = shaltime
  dumr(12) = detrain
  dumr(13) = detrainx
  dumr(14) = dsig2
  dumr(15) = dsig4
  dumr(16) = entrain
  dumr(17) = fldown
  dumr(18) = rhcv
  dumr(19) = rhmois
  dumr(20) = rhsat
  dumr(21) = sigcb
  dumr(22) = sigcll
  dumr(23) = sig_ct
  dumr(24) = sigkscb
  dumr(25) = sigksct
  dumr(26) = tied_con
  dumr(27) = tied_over
  dumr(28) = tied_rh
  dumr(29) = acon
  dumr(30) = bcon
  dumr(31) = rcm
  dumr(32) = rcrit_l
  dumr(33) = rcrit_s
  dumr(34) = cld_decay
  dumr(35) = maxlintime
  dumr(36) = qlg_max
  dumr(37) = qfg_max
  dumi(1)  = iterconv
  dumi(2)  = ksc
  dumi(3)  = kscmom
  dumi(4)  = kscsea
  dumi(5)  = ldr
  dumi(6)  = mbase
  dumi(7)  = mdelay
  dumi(8)  = methdetr
  dumi(9)  = methprec
  dumi(10) = nbase
  dumi(11) = ncvcloud
  dumi(12) = ncvmix
  dumi(13) = nevapcc
  dumi(14) = nkuo
  dumi(15) = nrhcrit
  dumi(16) = nstab_cld
  dumi(17) = nuvconv
  dumi(18) = ncloud
  dumi(19) = nclddia
  dumi(20) = nmr
  dumi(21) = nevapls
  dumi(22) = vdeposition_mode
  dumi(23) = tiedtke_form
  dumi(24) = cloud_aerosol_mode
  dumi(25) = lin_aerosolmode  
  dumi(26) = cloud_ice_method
  dumi(27) = leon_snowmeth
  dumi(28) = lin_adv
  dumi(29) = nscheme
end if
call ccmpi_bcast(dumr,0,comm_world)
call ccmpi_bcast(dumi,0,comm_world)
alflnd             = dumr(1)
alfsea             = dumr(2)
cldh_lnd           = dumr(3)
cldm_lnd           = dumr(4) 
cldl_lnd           = dumr(5)
cldh_sea           = dumr(6) 
cldm_sea           = dumr(7)
cldl_sea           = dumr(8)
convfact           = dumr(9)
convtime           = dumr(10)
shaltime           = dumr(11) 
detrain            = dumr(12)
detrainx           = dumr(13)
dsig2              = dumr(14)
dsig4              = dumr(15)
entrain            = dumr(16)
fldown             = dumr(17)
rhcv               = dumr(18)
rhmois             = dumr(19)
rhsat              = dumr(20)
sigcb              = dumr(21)
sigcll             = dumr(22)
sig_ct             = dumr(23)
sigkscb            = dumr(24)
sigksct            = dumr(25)
tied_con           = dumr(26)
tied_over          = dumr(27)
tied_rh            = dumr(28)
acon               = dumr(29)
bcon               = dumr(30)
rcm                = dumr(31)
rcrit_l            = dumr(32)
rcrit_s            = dumr(33)
cld_decay          = dumr(34)
maxlintime         = dumr(35)
qlg_max            = dumr(36)
qfg_max            = dumr(37)
iterconv           = dumi(1) 
ksc                = dumi(2)
kscmom             = dumi(3)
kscsea             = dumi(4)
ldr                = dumi(5)
mbase              = dumi(6)
mdelay             = dumi(7)
methdetr           = dumi(8) 
methprec           = dumi(9)
nbase              = dumi(10)
ncvcloud           = dumi(11)
ncvmix             = dumi(12)
nevapcc            = dumi(13)
nkuo               = dumi(14)
nrhcrit            = dumi(15)
nstab_cld          = dumi(16)
nuvconv            = dumi(17)
ncloud             = dumi(18)
nclddia            = dumi(19) 
nmr                = dumi(20)
nevapls            = dumi(21)
vdeposition_mode   = dumi(22)
tiedtke_form       = dumi(23)
cloud_aerosol_mode = dumi(24)
lin_aerosolmode    = dumi(25)
cloud_ice_method   = dumi(26)
leon_snowmeth      = dumi(27)
lin_adv            = dumi(28)
nscheme            = dumi(29)

return
end subroutine broadcast_kuonml

!--------------------------------------------------------------
! Broadcast turbnml namelist
subroutine broadcast_turbnml

use cc_mpi                                 ! CC MPI routines
use parm_m                                 ! Model configuration
use tkeeps                                 ! TKE-EPS boundary layer

implicit none

integer, dimension(6) :: dumi
real, dimension(32) :: dumr
    
dumr = 0.
dumi = 0
if ( myid==0 ) then
  dumr(1)  = be
  dumr(2)  = cm0
  dumr(3)  = ce0
  dumr(4)  = ce1
  dumr(5)  = ce2
  dumr(6)  = ce3
  dumr(7)  = cqmix
  dumr(8)  = ent0
  dumr(9)  = ent1
  dumr(10) = entc0
  dumr(11) = dtrc0
  dumr(12) = m0
  dumr(13) = b1
  dumr(14) = b2
  dumr(15) = maxdts
  dumr(16) = mintke
  dumr(17) = mineps
  dumr(18) = minl
  dumr(19) = maxl
  dumr(20) = qcmf
  dumr(21) = ezmin
  dumr(22) = amxlsq
  dumr(23) = helim
  dumr(24) = fc2
  dumr(25) = sigbot_gwd
  dumr(26) = alphaj
  dumr(27) = ent_min
  dumr(28) = mfbeta
  dumr(29) = dvmodmin
  dumr(30) = tke_timeave_length
  dumr(31) = wg_tau
  dumr(32) = wg_prob
  dumi(1)  = buoymeth
  dumi(2)  = stabmeth
  dumi(3)  = tkemeth
  dumi(4)  = ngwd
  dumi(5)  = ugs_meth
  dumi(6)  = tcalmeth
end if
call ccmpi_bcast(dumr,0,comm_world)
call ccmpi_bcast(dumi,0,comm_world)
be                 = dumr(1)
cm0                = dumr(2)
ce0                = dumr(3)
ce1                = dumr(4)
ce2                = dumr(5)
ce3                = dumr(6)
cqmix              = dumr(7)
ent0               = dumr(8)
ent1               = dumr(9)
entc0              = dumr(10)
dtrc0              = dumr(11)
m0                 = dumr(12)
b1                 = dumr(13)
b2                 = dumr(14)
maxdts             = dumr(15)
mintke             = dumr(16)
mineps             = dumr(17) 
minl               = dumr(18)
maxl               = dumr(19)
qcmf               = dumr(20)
ezmin              = dumr(21)
amxlsq             = dumr(22)
helim              = dumr(23)
fc2                = dumr(24)
sigbot_gwd         = dumr(25)
alphaj             = dumr(26)
ent_min            = dumr(27)
mfbeta             = dumr(28)
dvmodmin           = dumr(29)
tke_timeave_length = dumr(30)
wg_tau             = dumr(31)
wg_prob            = dumr(32)
buoymeth           = dumi(1)
stabmeth           = dumi(2)
tkemeth            = dumi(3)
ngwd               = dumi(4)
ugs_meth           = dumi(5)
tcalmeth           = dumi(6)

return
end subroutine broadcast_turbnml

!--------------------------------------------------------------
! Broadcast landnml namelist
subroutine broadcast_landnml

use cc_mpi                                 ! CC MPI routines
use parm_m                                 ! Model configuration
use river                                  ! River routing
use sflux_m                                ! Surface flux routines

implicit none

integer, dimension(32) :: dumi
real, dimension(26) :: dumr
    
dumr = 0.
dumi = 0
if ( myid==0 ) then
  dumr(1)  = real(energytol) 
  dumr(2)  = ateb_tol
  dumr(3)  = ateb_zosnow
  dumr(4)  = ateb_snowemiss
  dumr(5)  = ateb_maxsnowalpha
  dumr(6)  = ateb_minsnowalpha
  dumr(7)  = ateb_maxsnowden
  dumr(8)  = ateb_minsnowden
  dumr(9)  = ateb_refheight
  dumr(10) = ateb_zomratio
  dumr(11) = zocanyon
  dumr(12) = zoroof
  dumr(13) = ateb_maxrfwater
  dumr(14) = ateb_maxrdwater
  dumr(15) = ateb_maxrfsn
  dumr(16) = ateb_maxrdsn
  dumr(17) = ateb_maxvwatf
  dumr(18) = ateb_ac_heatcap
  dumr(19) = ateb_ac_coolcap
  dumr(20) = ateb_ac_deltat
  dumr(21) = ateb_acfactor
  dumr(22) = siburbanfrac
  dumr(23) = wbclim_lonn
  dumr(24) = wbclim_lonx
  dumr(25) = wbclim_latn
  dumr(26) = wbclim_latx
  dumi(1)  = proglai
  dumi(2)  = ccycle
  dumi(3)  = soil_struc
  dumi(4)  = cable_pop
  dumi(5)  = progvcmax
  dumi(6)  = fwsoil_switch
  dumi(7)  = cable_litter
  dumi(8)  = gs_switch
  dumi(9)  = smrf_switch
  dumi(10) = strf_switch
  dumi(11) = ateb_resmeth
  dumi(12) = ateb_zohmeth
  dumi(13) = ateb_acmeth
  dumi(14) = ateb_nrefl
  dumi(15) = ateb_scrnmeth
  dumi(16) = ateb_wbrelaxc
  dumi(17) = ateb_wbrelaxr
  dumi(18) = ateb_nfgits
  dumi(19) = intairtmeth
  dumi(20) = intmassmeth
  dumi(21) = ateb_cvcoeffmeth
  dumi(22) = ateb_statsmeth
  dumi(23) = ateb_lwintmeth
  dumi(24) = ateb_infilmeth
  dumi(25) = cable_roughness
  dumi(26) = cable_potev
  dumi(27) = ateb_soilunder
  dumi(28) = wt_transport
  dumi(29) = cable_gw_model
  dumi(30) = freshwaterlake_fix
  dumi(31) = cable_enablefao
  dumi(32) = ateb_ncyits
end if
call ccmpi_bcast(dumr,0,comm_world)
call ccmpi_bcast(dumi,0,comm_world)
energytol          = real(dumr(1),8)
ateb_tol           = dumr(2)
ateb_zosnow        = dumr(3)
ateb_snowemiss     = dumr(4)
ateb_maxsnowalpha  = dumr(5)
ateb_minsnowalpha  = dumr(6)
ateb_maxsnowden    = dumr(7)
ateb_minsnowden    = dumr(8)
ateb_refheight     = dumr(9) 
ateb_zomratio      = dumr(10)
zocanyon           = dumr(11)
zoroof             = dumr(12)
ateb_maxrfwater    = dumr(13)
ateb_maxrdwater    = dumr(14)
ateb_maxrfsn       = dumr(15)
ateb_maxrdsn       = dumr(16)
ateb_maxvwatf      = dumr(17) 
ateb_ac_heatcap    = dumr(18)
ateb_ac_coolcap    = dumr(19)
ateb_ac_deltat     = dumr(20)
ateb_acfactor      = dumr(21)
siburbanfrac       = dumr(22) 
wbclim_lonn        = dumr(23)
wbclim_lonx        = dumr(24)
wbclim_latn        = dumr(25)
wbclim_latx        = dumr(26)
proglai            = dumi(1)
ccycle             = dumi(2)
soil_struc         = dumi(3)
cable_pop          = dumi(4)
progvcmax          = dumi(5)
fwsoil_switch      = dumi(6)
cable_litter       = dumi(7)
gs_switch          = dumi(8)
smrf_switch        = dumi(9)
strf_switch        = dumi(10)
ateb_resmeth       = dumi(11)
ateb_zohmeth       = dumi(12)
ateb_acmeth        = dumi(13)
ateb_nrefl         = dumi(14) 
ateb_scrnmeth      = dumi(15)
ateb_wbrelaxc      = dumi(16) 
ateb_wbrelaxr      = dumi(17) 
ateb_nfgits        = dumi(18) 
intairtmeth        = dumi(19) 
intmassmeth        = dumi(20)
ateb_cvcoeffmeth   = dumi(21) 
ateb_statsmeth     = dumi(22) 
ateb_lwintmeth     = dumi(23) 
ateb_infilmeth     = dumi(24)
cable_roughness    = dumi(25)
cable_potev        = dumi(26)
ateb_soilunder     = dumi(27)
wt_transport       = dumi(28)
cable_gw_model     = dumi(29)
freshwaterlake_fix = dumi(30)
cable_enablefao    = dumi(31)
ateb_ncyits        = dumi(32)

return
end subroutine broadcast_landnml

!--------------------------------------------------------------
! Broadcast mlonml namelist
subroutine broadcast_mlonml

use cc_mpi                                 ! CC MPI routines
use mlo_ctrl                               ! Ocean physics control layer
use mlodynamics                            ! Ocean dynamics
use river                                  ! River routing

implicit none

integer, dimension(31) :: dumi
real, dimension(25) :: dumr    

dumr = 0.
dumi = 0
if ( myid==0 ) then
  dumr(1)  = ocnsmag
  dumr(2)  = ocneps
  dumr(3)  = zoseaice
  dumr(4)  = factchseaice
  dumr(5)  = minwater
  dumr(6)  = mxd
  dumr(7)  = mindep
  dumr(8)  = alphavis_seaice
  dumr(9)  = alphanir_seaice
  dumr(10) = rivercoeff
  dumr(11) = pdl
  dumr(12) = pdu
  dumr(13) = omink
  dumr(14) = omineps
  dumr(15) = ominl
  dumr(16) = omaxl
  dumr(17) = mlo_timeave_length
  dumr(18) = kemaxdt
  dumr(19) = alphavis_seasnw
  dumr(20) = alphanir_seasnw
  dumr(21) = fluxwgt
  dumr(22) = ocnepr
  dumr(23) = delwater
  dumr(24) = minsal
  dumr(25) = maxsal
  dumi(1)  = mlodiff
  dumi(2)  = usetide
  dumi(3)  = zomode
  dumi(4)  = otaumode
  dumi(5)  = rivermd
  dumi(6)  = basinmd
  dumi(7)  = mlojacobi
  dumi(8)  = usepice
  dumi(9)  = mlosigma
  dumi(10) = oclosure
  dumi(11) = k_mode
  dumi(12) = eps_mode
  dumi(13) = limitL
  dumi(14) = fixedce3
  dumi(15) = nops
  dumi(16) = nopb
  dumi(17) = fixedstabfunc
  dumi(18) = mlomfix
  dumi(19) = nodrift
  dumi(20) = mlontvd
  dumi(21) = mlodiff_numits
  dumi(22) = mlo_adjeta
  dumi(23) = mstagf
  dumi(24) = mlodps
  dumi(25) = mlo_limitsal
  dumi(26) = mlo_bs
  dumi(27) = mlo_step
  dumi(28) = mlo_uvcoupl
  dumi(29) = mlointschf
  dumi(30) = nxtrrho
  dumi(31) = mloiceadv
end if
call ccmpi_bcast(dumr,0,comm_world)
call ccmpi_bcast(dumi,0,comm_world)
ocnsmag            = dumr(1) 
ocneps             = dumr(2) 
zoseaice           = dumr(3) 
factchseaice       = dumr(4)
minwater           = dumr(5) 
mxd                = dumr(6)
mindep             = dumr(7)
alphavis_seaice    = dumr(8)
alphanir_seaice    = dumr(9)
rivercoeff         = dumr(10)
pdl                = dumr(11)
pdu                = dumr(12)
omink              = dumr(13)
omineps            = dumr(14)
ominl              = dumr(15)
omaxl              = dumr(16)
mlo_timeave_length = dumr(17)
kemaxdt            = dumr(18)
alphavis_seasnw    = dumr(19)
alphanir_seasnw    = dumr(20)
fluxwgt            = dumr(21)
ocnepr             = dumr(22)
delwater           = dumr(23)
minsal             = dumr(24)
maxsal             = dumr(25)
mlodiff            = dumi(1)
usetide            = dumi(2) 
zomode             = dumi(3) 
otaumode           = dumi(4) 
rivermd            = dumi(5)
basinmd            = dumi(6)
mlojacobi          = dumi(7)
usepice            = dumi(8)
mlosigma           = dumi(9)
oclosure           = dumi(10)
k_mode             = dumi(11)
eps_mode           = dumi(12)
limitL             = dumi(13)
fixedce3           = dumi(14)
nops               = dumi(15)
nopb               = dumi(16)
fixedstabfunc      = dumi(17)
mlomfix            = dumi(18)
nodrift            = dumi(19)
mlontvd            = dumi(20)
mlodiff_numits     = dumi(21)
mlo_adjeta         = dumi(22)
mstagf             = dumi(23)
mlodps             = dumi(24)
mlo_limitsal       = dumi(25)
mlo_bs             = dumi(26)
mlo_step           = dumi(27)
mlo_uvcoupl        = dumi(28)
mlointschf         = dumi(29)
nxtrrho            = dumi(30)
mloiceadv          = dumi(31)
    
return
end subroutine broadcast_mlonml
    
!--------------------------------------------------------------
! Broadcast trfiles namelist
subroutine broadcast_trfiles

use cc_mpi                                 ! CC MPI routines
use kuocom_m                               ! JLM convection
use tracermodule, only : tracerlist      & ! Tracer routines
    ,sitefile,shipfile,writetrpm         &
    ,init_tracer

implicit none

integer, dimension(1) :: dumi
    
dumi = 0
if ( myid==0 ) then
  if ( writetrpm ) dumi(1) = 1
end if
call ccmpi_bcast(tracerlist,0,comm_world)
if ( tracerlist/=' ' ) then
  call ccmpi_bcast(dumi,0,comm_world)
  call ccmpi_bcast(sitefile,0,comm_world)
  call ccmpi_bcast(shipfile,0,comm_world)
  writetrpm = dumi(1)==1
end if  
    
return
end subroutine broadcast_trfiles   
    
!--------------------------------------------------------------
! Find valid nproc
subroutine reducenproc(npanels,il_g,nproc,newnproc,nxp,nyp)

implicit none

integer, intent(in) :: il_g, nproc, npanels
integer, intent(out) :: newnproc, nxp, nyp
integer nproc_low, nxp_test, nyp_test

nxp_test = 0
nyp_test = 0

! try face decompositoin
do nproc_low = nproc,1,-1
  call proctest_face(npanels,il_g,nproc_low,nxp_test,nyp_test)
  if ( nxp_test>0 ) exit
end do
newnproc = nproc_low
nxp = nxp_test
nyp = nyp_test

return
end subroutine reducenproc

!--------------------------------------------------------------
! Find valid ntiles for physics
subroutine calc_phys_tiles(ntiles,maxtilesize,ifull)    

implicit none

integer, intent(in) :: maxtilesize, ifull
integer, intent(out) :: ntiles
integer i, tmp, imax

!find imax if maxtilesize isn't already a factor of ifull
imax = min( max( maxtilesize, 1 ), ifull )
tmp = imax
imax = -1 ! missing flag
! first attempt to find multiple of 8
do i = tmp,8,-1
  if ( mod(ifull,i)==0 .and. mod(i,8)==0 ) then
    imax = i
    exit
  end if
end do
if ( imax<1 ) then
  ! second attempt if multiple of 8 is not possible
  do i = tmp,1,-1
    if ( mod(ifull,i)==0 ) then
      imax = i
      exit
    end if
  end do
end if

!find the number of tiles
ntiles = ifull/imax

return
end subroutine calc_phys_tiles
    
!--------------------------------------------------------------
! TEST GRID DECOMPOSITION - FACE   
subroutine proctest_face(npanels,il_g,nproc,nxp,nyp)

implicit none

integer, intent(in) :: il_g, nproc, npanels
integer, intent(out) :: nxp, nyp
integer jl_g

if ( mod(nproc,6)/=0 .and. mod(6,nproc)/=0 ) then
  nxp = -1
  nyp = -1
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

return
end subroutine proctest_face

subroutine globpe_finalize

use cc_mpi                                 ! CC MPI routines
use bigxy4_m                               ! Grid interpolation
use map_m                                  ! Grid map arrays
use parm_m                                 ! Model configuration
use xyzinfo_m                              ! Grid coordinate arrays

implicit none

#ifdef share_ifullg
call ccmpi_freeshdata(xx4_win)
call ccmpi_freeshdata(yy4_win)
call ccmpi_freeshdata(em_g_win)
call ccmpi_freeshdata(x_g_win)
call ccmpi_freeshdata(y_g_win)
call ccmpi_freeshdata(z_g_win)
#else
deallocate(xx4, yy4)
deallocate(em_g)
deallocate(x_g, y_g, z_g)
#endif
call ccmpi_filewinfinalize_exit
if ( mbd/=0 .and. nud_uv/=9 ) then
  call deallocateglobalpack
end if
nullify(xx4, yy4)
nullify(em_g)
nullify(x_g, y_g, z_g)

! finalize MPI comms
call ccmpi_finalize

return
end subroutine globpe_finalize
    
end module config
