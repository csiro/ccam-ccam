program ccamscm

use aerointerface                          ! Aerosol interface
use aerosolldr, only : xtosav,xtg,naero  & ! LDR prognostic aerosols
    ,duste,dustwd,dustdd,dust_burden     &
    ,bce,bcwd,bcdd,bc_burden             &
    ,oce,ocwd,ocdd,oc_burden             &
    ,dmse,dmsso2o,dms_burden             &
    ,so2e,so2so4o,so2wd,so2dd,so2_burden &
    ,so4e,so4wd,so4dd,so4_burden         &
    ,Ch_dust,zvolcemi,aeroindir          &
    ,so4mtn,carbmtn                      &
    ,saltsmallmtn,saltlargemtn
use arrays_m                               ! Atmosphere dyamics prognostic arrays
use ateb, only : atebnmlfile             & ! Urban
    ,ateb_energytol=>energytol           &
    ,ateb_resmeth=>resmeth               &
    ,ateb_useonewall=>useonewall         &
    ,ateb_zohmeth=>zohmeth               &
    ,ateb_acmeth=>acmeth                 &
    ,ateb_nrefl=>nrefl                   &
    ,ateb_vegmode=>vegmode               &
    ,ateb_soilunder=>soilunder           &
    ,ateb_conductmeth=>conductmeth       &
    ,ateb_scrnmeth=>scrnmeth             &
    ,ateb_wbrelaxc=>wbrelaxc             &
    ,ateb_wbrelaxr=>wbrelaxr             &
    ,ateb_lweff=>lweff                   &
    ,ateb_ncyits=>ncyits                 &
    ,ateb_nfgits=>nfgits                 &
    ,ateb_tol=>tol                       &
    ,ateb_alpha=>alpha                   &
    ,ateb_zosnow=>zosnow                 &
    ,ateb_snowemiss=>snowemiss           &
    ,ateb_maxsnowalpha=>maxsnowalpha     &
    ,ateb_minsnowalpha=>minsnowalpha     &
    ,ateb_maxsnowden=>maxsnowden         &
    ,ateb_minsnowden=>minsnowden         &
    ,ateb_refheight=>refheight           &
    ,ateb_zomratio=>zomratio             &
    ,ateb_zocanyon=>zocanyon             &
    ,ateb_zoroof=>zoroof                 &
    ,ateb_maxrfwater=>maxrfwater         &
    ,ateb_maxrdwater=>maxrdwater         &
    ,ateb_maxrfsn=>maxrfsn               &
    ,ateb_maxrdsn=>maxrdsn               &
    ,ateb_maxvwatf=>maxvwatf             &
    ,ateb_intairtmeth=>intairtmeth       &
    ,ateb_intmassmeth=>intmassmeth       &
    ,ateb_cvcoeffmeth=>cvcoeffmeth       &
    ,ateb_statsmeth=>statsmeth           &
    ,ateb_behavmeth=>behavmeth           &
    ,ateb_infilmeth=>infilmeth           &
    ,ateb_ac_heatcap=>ac_heatcap         &
    ,ateb_ac_coolcap=>ac_coolcap         &
    ,ateb_ac_heatprop=>ac_heatprop       &
    ,ateb_ac_coolprop=>ac_coolprop       &
    ,ateb_ac_smooth=>ac_smooth           &
    ,ateb_ac_deltat=>ac_deltat           &
    ,ateb_acfactor=>acfactor             &
    ,ateb_ac_copmax=>ac_copmax
use cable_ccam, only : proglai           & ! CABLE
    ,soil_struc,cable_pop,progvcmax      &
    ,fwsoil_switch,cable_litter          &
    ,gs_switch,cable_climate
use carbpools_m, only : carbpools_init   & ! Carbon pools
    ,fpn,frs,frp
use cc_mpi                                 ! CC MPI routines
use cc_omp                                 ! CC OpenMP routines
use cfrac_m                                ! Cloud fraction
use cloudmod                               ! Prognostic cloud fraction
use const_phys                             ! Physical constants
use convjlm_m                              ! Convection
use convjlm22_m                            ! Convection v2
use dates_m                                ! Date data
use estab                                  ! Liquid saturation function
use extraout_m                             ! Additional diagnostics
use filnames_m                             ! Filenames
use getopt_m                               ! Command option parsing
use gdrag_m, only : gdrag_init,gwdrag    & ! Gravity wave drag
    ,gdrag_sbl
use histave_m                              ! Time average arrays
use infile                                 ! Input file routines
use kuocomb_m                              ! JLM convection
use latlong_m                              ! Lat/lon coordinates
use leoncld_mod                            ! Prognostic cloud condensate
use liqwpar_m                              ! Cloud water mixing ratios
use map_m                                  ! Grid map arrays
use mlo, only : mlodiag,wlev,mxd,mindep  & ! Ocean physics and prognostic arrays
   ,minwater,zomode,zoseaice             &
   ,factchseaice,minwater,mxd,mindep     &
   ,alphavis_seaice,alphanir_seaice      &
   ,otaumode
use morepbl_m                              ! Additional boundary layer diagnostics
use newmpar_m                              ! Grid parameters
use nharrs_m, only : nharrs_init         & ! Non-hydrostatic atmosphere arrays
   ,lrestart
use nsibd_m                                ! Land-surface arrays
use parm_m                                 ! Model configuration
use parmdyn_m                              ! Dynamics parameters
use parmgeom_m                             ! Coordinate data
use parmhdff_m                             ! Horizontal diffusion parameters
use parmhor_m                              ! Horizontal advection parameters
use pbl_m                                  ! Boundary layer arrays
use prec_m                                 ! Precipitation
use raddiag_m                              ! Radiation diagnostic
use riverarrays_m                          ! River data
use radisw_m
use savuvt_m                               ! Saved dynamic arrays
use screen_m                               ! Screen level diagnostics
use seaesfrad_m                            ! SEA-ESF radiation
use sflux_m                                ! Surface flux routines
use sigs_m                                 ! Atmosphere sigma levels
use soil_m                                 ! Soil and surface data
use soilsnow_m                             ! Soil, snow and surface data
use stime_m                                ! File date data
use tkeeps                                 ! TKE-EPS boundary layer
use vegpar_m                               ! Vegetation arrays
use vertmix_m                              ! Boundary layer turbulent mixing
use vvel_m                                 ! Additional vertical velocity
use work2_m                                ! Diagnostic arrays
use work3_m                                ! Mk3 land-surface diagnostic arrays
use work3f_m                               ! Grid work arrays

implicit none
    
include 'kuocom.h'                         ! Convection parameters
include 'version.h'                        ! Model version data

integer io_nest, npa, npb, mstn
integer secs_rad
integer iq, k
integer ivegt_in, isoil_in, gablsflux
integer jyear, jmonth, jday, jhour, jmin, mins
integer spinup, spinup_start, ntau_end, ntau_spinup
integer opt, nopt
integer, save :: iarch_nudge = 0
real, dimension(1000) :: press_in
real press_surf, gridres
real es
real rlong_in, rlat_in, z_in
real ateb_bldheight, ateb_hwratio, ateb_sigvegc, ateb_sigmabld
real ateb_industryfg, ateb_trafficfg, ateb_vegalphac
real ateb_wallalpha, ateb_roadalpha, ateb_roofalpha
real ateb_wallemiss, ateb_roademiss, ateb_roofemiss
real ateb_infilach,  ateb_intgains
real, dimension(4) :: ateb_roof_thick, ateb_roof_cp, ateb_roof_cond
real, dimension(4) :: ateb_wall_thick, ateb_wall_cp, ateb_wall_cond
real, dimension(4) :: ateb_road_thick, ateb_road_cp, ateb_road_cond
real, dimension(4) :: ateb_slab_thick, ateb_slab_cp, ateb_slab_cond
real, dimension(:,:), allocatable :: t_save, qg_save, u_save, v_save
real, dimension(:), allocatable :: psl_save
character(len=60) comm, comment
character(len=80) metforcing, timeoutput, profileoutput
character(len=80) lsmforcing, lsmoutput
character(len=80) scm_mode
character(len=1024) nmlfile
character(len=MAX_ARGLEN) :: optarg
logical oxidant_update
logical fixtsurf, nolatent, noradiation
logical nogwdrag, noconvection, nocloud, noaerosol, novertmix
logical lsm_only

namelist/scmnml/rlong_in,rlat_in,kl,press_in,press_surf,gridres,  &
    z_in,ivegt_in,isoil_in,metforcing,lsmforcing,lsmoutput,       &
    fixtsurf,nolatent,timeoutput,profileoutput,noradiation,       &
    nogwdrag,noconvection,nocloud,noaerosol,novertmix,            &
    lsm_only,                                                     &
    gablsflux,scm_mode,spinup_start,ntau_spinup,                  &
    ateb_bldheight,ateb_hwratio,ateb_sigvegc,ateb_sigmabld,       &
    ateb_industryfg,ateb_trafficfg,ateb_vegalphac,                &
    ateb_wallalpha,ateb_roadalpha,ateb_roofalpha,                 &
    ateb_wallemiss,ateb_roademiss,ateb_roofemiss,                 &
    ateb_roof_thick,ateb_roof_cp,ateb_roof_cond,                  &
    ateb_wall_thick,ateb_wall_cp,ateb_wall_cond,                  &
    ateb_road_thick,ateb_road_cp,ateb_road_cond,                  &
    ateb_slab_thick,ateb_slab_cp,ateb_slab_cond,                  &
    ateb_infilach,ateb_intgains
! main namelist
namelist/cardin/comment,dt,ntau,nwt,npa,npb,nhorps,nperavg,ia,ib, &
    ja,jb,id,jd,iaero,khdif,khor,nhorjlm,mex,mbd,nbd,             &
    mbd_maxscale,ndi,ndi2,nhor,nlv,nmaxpr,nrad,ntaft,ntsea,ntsur, &
    nvmix,restol,precon,kdate_s,ktime_s,leap,newtop,mup,lgwd,     &
    ngwd,rhsat,nextout,jalbfix,nalpha,nstag,nstagu,ntbar,nwrite,  &
    irest,nrun,nrungcm,nsib,                                      &
    mh_bs,nritch_t,nt_adv,mfix,mfix_qg,                           &
    namip,amipo3,nh,nhstest,nsemble,nspecial,panfg,panzo,         &
    rlatdn,rlatdx,rlongdn,rlongdx,newrough,nglacier,newztsea,     &
    epsp,epsu,epsf,epsh,av_vmod,charnock,chn10,snmin,tss_sh,      &
    vmodmin,zobgin,rlong0,rlat0,schmidt,kbotdav,kbotu,nud_p,      &
    nud_q,nud_t,nud_uv,nud_hrs,nudu_hrs,sigramplow,sigramphigh,   &
    nlocal,nbarewet,nsigmf,qgmin,io_in,io_nest,io_out,io_rest,    &
    tblock,tbave,localhist,m_fly,nurban,nmr,ktopdav,              &
    nud_sst,nud_sss,mfix_tr,mfix_aero,kbotmlo,ktopmlo,mloalpha,   &
    nud_ouv,nud_sfh,bpyear,rescrn,helmmeth,nmlo,ol,mxd,mindep,    &
    minwater,zomode,zoseaice,factchseaice,                        &
    knh,ccycle,kblock,nud_aero,ch_dust,zvolcemi,aeroindir,helim,  &
    fc2,sigbot_gwd,alphaj,proglai,cgmap_offset,cgmap_scale,       &
    nriver,amxlsq
! radiation and aerosol namelist
namelist/skyin/mins_rad,sw_resolution,sw_diff_streams,            & ! radiation
    liqradmethod,iceradmethod,so4radmethod,carbonradmethod,       &
    dustradmethod,seasaltradmethod,bpyear,qgmin,lwem_form,        & 
    ch_dust,zvolcemi,aeroindir,so4mtn,carbmtn,saltsmallmtn,       & ! aerosols
    saltlargemtn
! file namelist
namelist/datafile/ifile,ofile,albfile,eigenv,icefile,mesonest,    &
    o3file,radfile,restfile,rsmfile,so4tfile,soilfile,sstfile,    &
    surfile,topofile,vegfile,zofile,surf_00,surf_12,laifile,      &
    albnirfile,urbanfile,bathfile,vegprev,vegnext,vegnext2,       &
    cnsdir,salfile,oxidantfile,casafile,phenfile,                 &
    save_aerosols,save_pbl,save_cloud,save_land,save_maxmin,      &
    save_ocean,save_radiation,save_urban,save_carbon,save_river
! convection and cloud microphysics namelist
namelist/kuonml/alflnd,alfsea,cldh_lnd,cldm_lnd,cldl_lnd,         &
    cldh_sea,cldm_sea,cldl_sea,convfact,convtime,shaltime,        &
    detrain,detrainx,dsig2,dsig4,entrain,fldown,iterconv,ksc,     &
    kscmom,kscsea,ldr,mbase,mdelay,methdetr,methprec,nbase,       &
    nclddia,ncvcloud,ncvmix,nevapcc,nevapls,nkuo,nrhcrit,         &
    nstab_cld,nuvconv,rhcv,rhmois,rhsat,sigcb,sigcll,sig_ct,      &
    sigkscb,sigksct,tied_con,tied_over,tied_rh,comm,acon,bcon,    &
    rcm,rcrit_l,rcrit_s,ncloud
! boundary layer turbulence and gravity wave namelist
namelist/turbnml/be,cm0,ce0,ce1,ce2,ce3,cq,ent0,ent1,entc0,dtrc0, & !EDMF PBL scheme
    m0,b1,b2,buoymeth,maxdts,mintke,mineps,minl,maxl,             &
    stabmeth,tke_umin,tkemeth,qcmf,ezmin,ent_min,                 &
    amxlsq,                                                       & !JH PBL scheme
    ngwd,helim,fc2,sigbot_gwd,alphaj                                !GWdrag
! land, urban and carbon namelist
namelist/landnml/proglai,ccycle,soil_struc,cable_pop,             & ! CABLE
    progvcmax,fwsoil_switch,cable_litter,                         &
    gs_switch,cable_climate,                                      &
    ateb_energytol,ateb_resmeth,ateb_useonewall,ateb_zohmeth,     & ! urban
    ateb_acmeth,ateb_nrefl,ateb_vegmode,ateb_soilunder,           &
    ateb_conductmeth,ateb_scrnmeth,ateb_wbrelaxc,ateb_wbrelaxr,   &
    ateb_lweff,ateb_ncyits,ateb_nfgits,ateb_tol,ateb_alpha,       &
    ateb_zosnow,ateb_snowemiss,ateb_maxsnowalpha,                 &
    ateb_minsnowalpha,ateb_maxsnowden,ateb_minsnowden,            &
    ateb_refheight,ateb_zomratio,ateb_zocanyon,ateb_zoroof,       &
    ateb_maxrfwater,ateb_maxrdwater,ateb_maxrfsn,ateb_maxrdsn,    &
    ateb_maxvwatf,ateb_intairtmeth,ateb_intmassmeth,              &
    ateb_cvcoeffmeth,ateb_statsmeth,ateb_behavmeth,               &
    ateb_infilmeth,ateb_ac_heatcap,ateb_ac_coolcap,               &
    ateb_ac_heatprop,ateb_ac_coolprop,ateb_ac_smooth,             &
    ateb_ac_deltat,ateb_acfactor,ateb_ac_copmax
! ocean namelist
namelist/mlonml/zomode,zoseaice,                                  &
    factchseaice,minwater,mxd,mindep,otaumode,                    &
    alphavis_seaice,alphanir_seaice

data comment/' '/,comm/' '/

write(6,*) "CCAM SCM"
write(6,*) version
write(6,*)

nh = 5
gablsflux = 0
nproc = 1
maxthreads = 1
maxtilesize = 1
ntiles = 1
imax = 1
scm_mode = "gabls4"
spinup_start = 1
ntau_spinup = 0
fixtsurf = .false.
nolatent = .false.
noradiation = .false.
nogwdrag = .false.
noconvection = .false.
nocloud = .false.
noaerosol = .false.
novertmix = .false.
lsm_only = .false.
ateb_bldheight = -999.
ateb_hwratio = -999.
ateb_sigvegc = -999.
ateb_sigmabld = -999.
ateb_industryfg = -999.
ateb_trafficfg = -999.
ateb_roofalpha = -999.
ateb_wallalpha = -999.
ateb_roadalpha = -999.
ateb_roofemiss = -999.
ateb_wallemiss = -999.
ateb_roademiss = -999.
ateb_vegalphac = -999.
ateb_roof_thick = -999.
ateb_roof_cp = -999.
ateb_roof_cond = -999.
ateb_wall_thick = -999. 
ateb_wall_cp = -999.
ateb_wall_cond = -999.
ateb_road_thick = -999.
ateb_road_cp = -999.
ateb_road_cond = -999.
ateb_slab_thick = -999.
ateb_slab_cp = -999.
ateb_slab_cond = -999.
ateb_infilach = -999.
ateb_intgains = -999.
ateb_energytol = 0.005_8

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
! READ COMMAND LINE OPTIONS
nmlfile = "input"
do
  call getopt("c:",nopt,opt,optarg)
  if ( opt==-1 ) exit  ! End of options
  select case ( char(opt) )
    case ( "c" )
      nmlfile = optarg
    case default
      if ( myid==0 ) write(6,*) "ERROR: Unknown command line option ",char(opt)
      stop
  end select
end do

write(6,*) "Reading namelist"
open(99,file=nmlfile,form="formatted",status="old")
read(99,scmnml)
read(99,cardin)
read(99,skyin)
read(99,datafile)
read(99,kuonml)
read(99,turbnml)
read(99,landnml)

write(6,*) "scm_mode ",trim(scm_mode)

myid = 0
mydiag = .false.

il_g = 1
jl_g = 1
ifull_g = 1
il = 1
jl = 1
ifull = 1
iextra = 0
nrows_rad = 1

procformat=.false.
procmode = nproc
comm_vnode = comm_world
vnode_nproc = 1
vnode_myid = 0
comm_vleader = comm_world
vleader_nproc = nproc
vleader_myid = myid
vnode_vleaderid = myid

maxtilesize = 1
ntiles = 1

nperday = nint(24.*3600./dt)           ! time-steps in one day
nperhr  = nint(3600./dt)               ! time-steps in one hour

schmidt = gridres*real(il_g)/(90.*112.)
write(6,*) "gridres,schmidt ",gridres,schmidt

nriver=0
if ( abs(nmlo)>=2 ) then
  write(6,*) "ERROR: Cannot use dynamical MLO with nmlo ",nmlo
  stop -1
end if

call map_init(ifull_g,ifull,iextra,myid)

call arrays_init(ifull,iextra,kl)
call carbpools_init(ifull,nsib,ccycle)
call cfrac_init(ifull,kl)
call cloudmod_init(ifull,iextra,kl,ncloud)
call estab_init
call extraout_init(ifull,nextout)
call gdrag_init(ifull)
call histave_init(ifull,kl,ms,ccycle)
call kuocomb_init(ifull,kl)
call latlong_init(ifull_g,ifull,myid)
call liqwpar_init(ifull,iextra,kl)
call morepbl_init(ifull,kl)
call nharrs_init(ifull,iextra,kl)
call nsibd_init(ifull,nsib,cable_climate)
call pbl_init(ifull)
call prec_init(ifull)
call raddiag_init(ifull,kl)
call riverarrays_init(ifull,iextra,nriver)
call savuvt_init(ifull,kl)
call screen_init(ifull)
call sigs_init(kl)
call soil_init(ifull,iaero,nsib)
call soilsnow_init(ifull,ms,nsib)
call vvel_init(ifull,kl)
call vegpar_init(ifull)
call work2_init(ifull,nsib)
call work3_init(ifull,nsib)
call work3f_init(ifull,kl)
if ( nvmix==6 ) then
  call tkeinit(ifull,iextra,kl,0)
end if

rlatt(1) = rlat_in*pi/180.
rlongg(1) = rlong_in*pi/180.
dtin = dt
rrvco2 = 330./1.e6

write(6,*) "Calling initialscm"
call initialscm(scm_mode,metforcing,lsmforcing,press_in(1:kl),press_surf,z_in,ivegt_in, &
                isoil_in,ateb_bldheight,ateb_hwratio,ateb_sigvegc,ateb_sigmabld,        &
                ateb_industryfg,ateb_trafficfg,ateb_vegalphac,                          &
                ateb_roofalpha,ateb_wallalpha,ateb_roadalpha,                           &
                ateb_roofemiss,ateb_wallemiss,ateb_roademiss,                           &
                ateb_roof_thick,ateb_roof_cp,ateb_roof_cond,                            &
                ateb_wall_thick,ateb_wall_cp,ateb_wall_cond,                            &
                ateb_road_thick,ateb_road_cp,ateb_road_cond,                            &
                ateb_slab_thick,ateb_slab_cp,ateb_slab_cond,                            &
                ateb_infilach,ateb_intgains)

allocate( t_save(ifull,kl), qg_save(ifull,kl), u_save(ifull,kl), v_save(ifull,kl) )
allocate( psl_save(ifull) )
t_save = t
u_save = u
v_save = v
qg_save = qg
psl_save = psl

call nantest("after initialisation",1,ifull)

kountr = 1
if ( myid==0 ) then
  write(6,*) "Radiation will use kountr ",kountr
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
grpl(:)        = 0.
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
if ( ccycle>0 ) then
  fnee_ave = 0.  
  fpn_ave  = 0.
  frd_ave  = 0.
  frp_ave  = 0.
  frpw_ave = 0.
  frpr_ave = 0.
  frs_ave  = 0.
end if
if ( abs(iaero)==2 ) then
  duste         = 0.  ! Dust emissions
  dustdd        = 0.  ! Dust dry deposition
  dustwd        = 0.  ! Dust wet deposition
  dust_burden   = 0.  ! Dust burden
  bce           = 0.  ! Black carbon emissions
  bcdd          = 0.  ! Black carbon dry deposition
  bcwd          = 0.  ! Black carbon wet deposition
  bc_burden     = 0.  ! Black carbon burden
  oce           = 0.  ! Organic carbon emissions
  ocdd          = 0.  ! Organic carbon dry deposition
  ocwd          = 0.  ! Organic carbon wet deposition
  oc_burden     = 0.  ! Organic carbon burden
  dmse          = 0.  ! DMS emissions
  dmsso2o       = 0.  ! DMS -> SO2 oxidation
  so2e          = 0.  ! SO2 emissions
  so2so4o       = 0.  ! SO2 -> SO4 oxidation
  so2dd         = 0.  ! SO2 dry deposition
  so2wd         = 0.  ! SO2 wet deposiion
  so4e          = 0.  ! SO4 emissions
  so4dd         = 0.  ! SO4 dry deposition
  so4wd         = 0.  ! SO4 wet deposition
  dms_burden    = 0.  ! DMS burden
  so2_burden    = 0.  ! SO2 burden
  so4_burden    = 0.  ! SO4 burden
end if

! NUDGING
ktau = 0
call nudgescm(scm_mode,metforcing,fixtsurf,iarch_nudge)
call nantest("after nudging",1,ifull)

savu(1:ifull,:) = u(1:ifull,:)
savv(1:ifull,:) = v(1:ifull,:)
savs(1:ifull,:) = sdot(1:ifull,2:kl)

! INITIAL OUTPUT
call outputscm(scm_mode,timeoutput,profileoutput,lsmoutput)

call gdrag_sbl
select case ( nkuo )
  case(21,22)
    call convjlm22_init
  case(23,24)
    call convjlm_init
  end select
select case(nrad)
  case(5)
    call seaesfrad_init
end select
call sflux_init
call vertmix_init

! spin-up?
do spinup = spinup_start,1,-1

  if ( spinup>1 ) then
    ntau_end = ntau_spinup
    write(6,*) "Spin-up Phase ",spinup
  else
    ntau_end = ntau
    write(6,*) "Main simulation"
  end if  
  
  iarch_nudge = 0 ! reset nudging
  t = t_save
  qg = qg_save
  u = u_save
  v = v_save
  psl = psl_save

  ! main simulation
  do ktau = 1,ntau_end
    write(6,*) "Time ",ktau,ntau
    if ( isoilm(1)==9 ) then
      write(6,*) "tggsn ",tggsn,isflag
      !write(6,*) "ssdn ",ssdn
      !write(6,*) "snowd ",snowd
    end if
    if ( abs(nurban)>0 ) then
      write(6,*) "rnet,fg,eg ",rnet(1),fg(1),eg(1)
      write(6,*) "storage,anthro ",urban_storage_flux(1),anthropogenic_flux(1)
    end if
    if ( nvmix==6 ) then
      write(6,*) "tke,eps,kh ",tke(1,1),eps(1,1),cm0*tke(1,1)*tke(1,1)/eps(1,1)
    end if
    if ( abs(nurban)==1 ) then
      write(6,*) "tss,t1,pblh ",tss(1),t(1,1),pblh(1)
    end if
  
    ! RESET AVERAGES
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
    grpl(:)        = 0.
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
    if ( ccycle>0 ) then
      fnee_ave = 0.  
      fpn_ave  = 0.
      frd_ave  = 0.
      frp_ave  = 0.
      frpw_ave = 0.
      frpr_ave = 0.
      frs_ave  = 0.
    end if  
    if ( abs(iaero)==2 ) then
      duste         = 0.  ! Dust emissions
      dustdd        = 0.  ! Dust dry deposition
      dustwd        = 0.  ! Dust wet deposition
      dust_burden   = 0.  ! Dust burden
      bce           = 0.  ! Black carbon emissions
      bcdd          = 0.  ! Black carbon dry deposition
      bcwd          = 0.  ! Black carbon wet deposition
      bc_burden     = 0.  ! Black carbon burden
      oce           = 0.  ! Organic carbon emissions
      ocdd          = 0.  ! Organic carbon dry deposition
      ocwd          = 0.  ! Organic carbon wet deposition
      oc_burden     = 0.  ! Organic carbon burden
      dmse          = 0.  ! DMS emissions
      dmsso2o       = 0.  ! DMS -> SO2 oxidation
      so2e          = 0.  ! SO2 emissions
      so2so4o       = 0.  ! SO2 -> SO4 oxidation
      so2dd         = 0.  ! SO2 dry deposition
      so2wd         = 0.  ! SO2 wet deposiion
      so4e          = 0.  ! SO4 emissions
      so4dd         = 0.  ! SO4 dry deposition
      so4wd         = 0.  ! SO4 wet deposition
      dms_burden    = 0.  ! DMS burden
      so2_burden    = 0.  ! SO2 burden
      so4_burden    = 0.  ! SO4 burden
    end if

  
    mtimer = nint(real(ktau)*dt/60.)

    ! calculate shear
    call calcshear
  
    ! Update saved arrays
    savu(1:ifull,:) = u(1:ifull,:)
    savv(1:ifull,:) = v(1:ifull,:)
    savs(1:ifull,:) = sdot(1:ifull,2:kl)

    ! MISC
    if ( nrad==5 ) then
      call seaesfrad_settime
    end if    
    ! aerosol timer calculations
    call getzinp(jyear,jmonth,jday,jhour,jmin,mins)
    oxidant_update = oxidant_timer<=mins-updateoxidant
    ! initialse surface rainfall to zero
    condc(:) = 0. ! default convective rainfall (assumed to be rain)
    condx(:) = 0. ! default total precip = rain + ice + snow + graupel (convection and large scale)
    conds(:) = 0. ! default total ice + snow (convection and large scale)
    condg(:) = 0. ! default total graupel (convection and large scale)
    ! Save aerosol concentrations for outside convective fraction of grid box
    if ( abs(iaero)>=2 ) then
      xtosav(:,:,:) = xtg(:,:,:) ! Aerosol mixing ratio outside convective cloud
    end if
    call nantest("start of physics",1,ifull)
  
    ! GWDRAG
    if ( ngwd<0 .and. .not.nogwdrag ) then
      call gwdrag  ! <0 for split - only one now allowed
    end if
    call nantest("after gravity wave drag",1,ifull)

    ! CONVECTION
    convh_ave(1:ifull,1:kl) = convh_ave(1:ifull,1:kl) - t(1:ifull,1:kl)*real(nperday)/real(nperavg)
    ! Select convection scheme
    if ( .not.noconvection ) then
      select case ( nkuo )
        case(21,22)
          call convjlm22              ! split convjlm 
        case(23,24)
          call convjlm                ! split convjlm 
      end select
    end if    
    call fixqg(1,ifull)
    call nantest("after convection",1,ifull)    

    ! CLOUD MICROPHYSICS
    if ( ldr/=0 .and. .not.nocloud ) then
      ! LDR microphysics scheme
      call leoncld
    end if
    convh_ave(1:ifull,1:kl) = convh_ave(1:ifull,1:kl) + t(1:ifull,1:kl)*real(nperday)/real(nperavg)    
    call nantest("after cloud microphysics",1,ifull) 

    ! RADIATON
    if ( ncloud>=4 ) then
      nettend(1:ifull,1:kl) = nettend(1:ifull,1:kl) + t(1:ifull,1:kl)/dt
    end if  
    if ( .not.noradiation ) then
      select case ( nrad )
        case(5)
          ! GFDL SEA-EFS radiation
          call seaesfrad
      end select
    end if
    call nantest("after radiation",1,ifull)   

    ! REPLACE SCM WITH INPUT DATA, PRIOR TO SFLUX
    if ( lsm_only ) then
      call replace_scm(lsmforcing)
      call nantest("after replace_scm",1,ifull)   
    end if
    
    ! SURFACE FLUXES
    if ( ntsur>1 ) then  ! should be better after convjlm
      call sflux
    endif   ! (ntsur>1) 
    call nantest("after surface fluxes",1,ifull)
  
    call nudgescm(scm_mode,metforcing,fixtsurf,iarch_nudge)
  
    !if ( gablsflux>0 ) then
    !  ppa(:) = sig(1)*ps(:) 
    !  prhoa(:) = ppa(:)/(rdry*t(1:ifull,1))
    !  zref(:) = (bet(1)*t(1:ifull,1)+phi_nh(:,1))/grav
    !  call gabls_flux(tss,ps,t(1:ifull,1),ppa,u(1:ifull,1),v(1:ifull,1),prhoa,zref,zolnd, &
    !                    psfth,ustar,vmod)
    !  rhos(:) = ps(1:ifull)/(rdry*tss(1:ifull))
    !  fg(:) = rhos(:)*cp*psfth(:)
    !  eg(:) = 0.
    !  cduv(:) = ustar(:)**2/vmod(:)  ! ustar = sqrt(cd)*vmod = sqrt(cduv*vmod)
    !end if
  
    if ( nolatent ) then
      eg(:) = 0.
      cdtq(:) = 0.
    end if

    ! AEROSOLS
    if ( abs(iaero)>=2 .and. .not.noaerosol ) then
      call aerocalc(oxidant_update,mins)
    end if
    call nantest("after aerosols",1,ifull)

    ! VERTICAL MIXING
    if ( ntsur>=1 .and. .not.novertmix ) then ! calls vertmix but not sflux for ntsur=1
      call vertmix
    endif  ! (ntsur>=1)
    call fixqg(1,ifull)
    call nantest("after PBL mixing",1,ifull)
  
    ! DIAGNOSTICS
    if ( rescrn > 0 ) then
      call autoscrn(1,ifull)
    end if
    ! Convection diagnostic output
    cbas_ave(1:ifull) = cbas_ave(1:ifull) + condc(1:ifull)*(1.1-sig(kbsav(1:ifull)))      ! diagnostic
    ctop_ave(1:ifull) = ctop_ave(1:ifull) + condc(1:ifull)*(1.1-sig(abs(ktsav(1:ifull)))) ! diagnostic
    ! Microphysics diagnostic output
    do k = 1,kl
      riwp_ave(1:ifull) = riwp_ave(1:ifull) - qfrad(1:ifull,k)*dsig(k)*ps(1:ifull)/grav ! ice water path
      rlwp_ave(1:ifull) = rlwp_ave(1:ifull) - qlrad(1:ifull,k)*dsig(k)*ps(1:ifull)/grav ! liq water path
    end do
    rnd_3hr(1:ifull,8) = rnd_3hr(1:ifull,8) + condx(1:ifull)  ! i.e. rnd24(:)=rnd24(:)+condx(:)

    ! MISC    
    ! Update aerosol timer
    if ( oxidant_update ) then
      oxidant_timer = mins
    end if
  
    ! OUTPUT
    if ( spinup==1 ) then
      call outputscm(scm_mode,timeoutput,profileoutput,lsmoutput)
    end if
  
  end do ! ktau
    
end do   ! spinup

write(6,*) "SCM complete"

end

subroutine calcshear

use arrays_m
use const_phys
use newmpar_m
use nharrs_m
use parm_m
use savuvt_m
use sigs_m
use tkeeps, only : shear

implicit none

integer k
real, dimension(ifull,kl) :: uav, vav
real, dimension(ifull,kl) :: dudz, dvdz
real, dimension(ifull,kl) :: zg, tnhs
real, dimension(ifull,2) :: zgh
real, dimension(ifull) :: r1, r2

! calculate height on full levels and non-hydrostatic temp correction
tnhs(:,1)=phi_nh(:,1)/bet(1)
zg(:,1) = (zs(1:ifull)+bet(1)*t(1:ifull,1))/grav
do k=2,kl
  tnhs(:,k)=(phi_nh(:,k)-phi_nh(:,k-1)-betm(k)*tnhs(:,k-1))/bet(k)
  zg(:,k) = zg(:,k-1) + (bet(k)*t(1:ifull,k)+betm(k)*t(1:ifull,k-1))/grav
end do ! k  loop
zg(:,:) = zg(:,:) + phi_nh(:,:)/grav

do k=1,kl        
  ! weighted horizontal velocities
  uav(1:ifull,k)=av_vmod*u(1:ifull,k)+(1.-av_vmod)*savu(1:ifull,k)
  vav(1:ifull,k)=av_vmod*v(1:ifull,k)+(1.-av_vmod)*savv(1:ifull,k)
end do

! calculate vertical gradients
zgh(:,2)=ratha(1)*zg(:,2)+rathb(1)*zg(:,1) ! upper half level
r1=uav(1:ifull,1)
r2=ratha(1)*uav(1:ifull,2)+rathb(1)*uav(1:ifull,1)          
dudz(1:ifull,1)=(r2-r1)/(zgh(1:ifull,2)-zg(1:ifull,1))
r1=vav(1:ifull,1)
r2=ratha(1)*vav(1:ifull,2)+rathb(1)*vav(1:ifull,1)          
dvdz(1:ifull,1)=(r2-r1)/(zgh(1:ifull,2)-zg(1:ifull,1))
do k=2,kl-1
  zgh(:,1)=zgh(:,2) ! lower half level
  zgh(:,2)=ratha(k)*zg(:,k+1)+rathb(k)*zg(:,k) ! upper half level
  r1=ratha(k-1)*uav(1:ifull,k)+rathb(k-1)*uav(1:ifull,k-1)
  r2=ratha(k)*uav(1:ifull,k+1)+rathb(k)*uav(1:ifull,k)          
  dudz(1:ifull,k)=(r2-r1)/(zgh(1:ifull,2)-zgh(1:ifull,1))
  r1=ratha(k-1)*vav(1:ifull,k)+rathb(k-1)*vav(1:ifull,k-1)
  r2=ratha(k)*vav(1:ifull,k+1)+rathb(k)*vav(1:ifull,k)          
  dvdz(1:ifull,k)=(r2-r1)/(zgh(1:ifull,2)-zgh(1:ifull,1))
end do
zgh(:,1)=zgh(:,2) ! lower half level
r1=ratha(kl-1)*uav(1:ifull,kl)+rathb(kl-1)*uav(1:ifull,kl-1)
r2=uav(1:ifull,kl)          
dudz(1:ifull,kl)=(r2-r1)/(zg(1:ifull,kl)-zgh(1:ifull,1))
r1=ratha(kl-1)*vav(1:ifull,kl)+rathb(kl-1)*vav(1:ifull,kl-1)
r2=vav(1:ifull,kl)          
dvdz(1:ifull,kl)=(r2-r1)/(zg(1:ifull,kl)-zgh(1:ifull,1))

if ( nvmix==6 ) then
  do k = 1,kl
    shear(:,k) = dudz(:,k)**2 + dvdz(:,k)**2
  end do
end if  

end subroutine calcshear
    
subroutine initialscm(scm_mode,metforcing,lsmforcing,press_in,press_surf,z_in,ivegt_in, &
                      isoil_in,ateb_bldheight,ateb_hwratio,ateb_sigvegc,ateb_sigmabld,  &
                      ateb_industryfg,ateb_trafficfg,ateb_vegalphac,                    &
                      ateb_roofalpha,ateb_wallalpha,ateb_roadalpha,                     &
                      ateb_roofemiss,ateb_wallemiss,ateb_roademiss,                     &
                      ateb_roof_thick,ateb_roof_cp,ateb_roof_cond,                      &
                      ateb_wall_thick,ateb_wall_cp,ateb_wall_cond,                      &
                      ateb_road_thick,ateb_road_cp,ateb_road_cond,                      &
                      ateb_slab_thick,ateb_slab_cp,ateb_slab_cond,                      &
                      ateb_infilach,ateb_intgains)

use aerointerface                          ! Aerosol interface
use aerosolldr                             ! LDR prognostic aerosols
use arrays_m                               ! Atmosphere dyamics prognostic arrays
use ateb                                   ! Urban
use cable_ccam                             ! CABLE
use carbpools_m                            ! Carbon pools
use cfrac_m                                ! Cloud fraction
use const_phys                             ! Physical constants
use dates_m                                ! Date data
use filnames_m                             ! Filenames
use gdrag_m                                ! Gravity wave drag
use infile                                 ! Input file routines
use latlong_m                              ! Lat/lon coordinates
use liqwpar_m                              ! Cloud water mixing ratios
use map_m                                  ! Grid map arrays
use mlo                                    ! Ocean physics and prognostic arrays
use newmpar_m                              ! Grid parameters
use nharrs_m                               ! Non-hydrostatic atmosphere arrays
use nsibd_m                                ! Land-surface arrays
use parm_m                                 ! Model configuration
use parmdyn_m                              ! Dynamics parameters
use parmgeom_m                             ! Coordinate data
use pbl_m                                  ! Boundary layer arrays
use sigs_m                                 ! Atmosphere sigma levels
use soil_m                                 ! Soil and surface data
use soilsnow_m                             ! Soil, snow and surface data
use stime_m                                ! File date data
use vegpar_m                               ! Vegetation arrays
use vvel_m                                 ! Additional vertical velocity

use scmarrays_m

implicit none

integer, intent(in) :: ivegt_in, isoil_in
integer iq, i, k, lapsbot
integer ncid, ncstatus, nlev, slev, tlev
integer, dimension(ifull,5) :: ivs
integer, dimension(271,mxvt) :: greenup, fall, phendoy1
integer, dimension(3) :: spos, npos
integer, dimension(ifull) :: urbantype
real, dimension(kl), intent(in) :: press_in
real, intent(in) :: press_surf, z_in
real, dimension(ifull,5) :: svs, vlin, vlinprev, vlinnext, vlinnext2
real, dimension(ifull,5) :: casapoint
real, dimension(ifull) :: aa
real, dimension(ifull,wlev,4) :: mlodwn
real, dimension(ifull,2) :: ocndwn
real, dimension(ifull) :: depth
real, dimension(kl) :: datout
real, dimension(:), allocatable, save :: dat_in, sig_in, dep_in
real, dimension(:), allocatable, save :: dephl_in, thick_in
real, dimension(:), allocatable, save :: new_in
real, dimension(:), allocatable, save :: height_in, height_model
real, dimension(1) :: psurf_in
real, intent(in) :: ateb_bldheight, ateb_hwratio, ateb_sigvegc, ateb_sigmabld
real, intent(in) :: ateb_industryfg, ateb_trafficfg, ateb_vegalphac
real, intent(in) :: ateb_roofalpha, ateb_wallalpha, ateb_roadalpha
real, intent(in) :: ateb_roofemiss, ateb_wallemiss, ateb_roademiss
real, intent(in) :: ateb_infilach,  ateb_intgains
real, dimension(4), intent(in) :: ateb_roof_thick, ateb_roof_cp, ateb_roof_cond
real, dimension(4), intent(in) :: ateb_wall_thick, ateb_wall_cp, ateb_wall_cond
real, dimension(4), intent(in) :: ateb_road_thick, ateb_road_cp, ateb_road_cond
real, dimension(4), intent(in) :: ateb_slab_thick, ateb_slab_cp, ateb_slab_cond
real, dimension(8) :: atebparm
real gwdfac, hefact, c
character(len=*), intent(in) :: metforcing, lsmforcing, scm_mode
character(len=20) vname
logical, save :: lsmonly = .false.

allocate( sdepth(ifull,3) )

kdate = kdate_s
ktime = ktime_s
lapsbot = 2
write(6,*) "kdate,ktime = ",kdate,ktime

! Model data
sig(1:kl) = press_in(1:kl)/press_surf
write(6,*) "sig ",sig(1:kl)
sigmh(1) = 1.
do k=2,kl
  sigmh(k) = 0.5*(sig(k-1)+sig(k))
end do
dsig(1:kl-1) = sigmh(2:kl)-sigmh(1:kl-1)
dsig(kl) = -sigmh(kl)
!     rata and ratb are used to interpolate half level values to full levels
!     ratha and rathb are used to interpolate full level values to half levels
ratha(kl) = 0. ! not used
rathb(kl) = 0. ! not used
rata(kl)=(sigmh(kl)-sig(kl))/sigmh(kl)
ratb(kl)=sig(kl)/sigmh(kl)
do k=1,kl-1
  bet(k+1)=rdry*log(sig(k)/sig(k+1))*.5
  rata(k)=(sigmh(k)-sig(k))/(sigmh(k)-sigmh(k+1))
  ratb(k)=(sig(k)-sigmh(k+1))/(sigmh(k)-sigmh(k+1))
  ratha(k)=(sigmh(k+1)-sig(k))/(sig(k+1)-sig(k))
  rathb(k)=(sig(k+1)-sigmh(k+1))/(sig(k+1)-sig(k))
enddo
c=grav/stdlapse
bet(1)=c *(sig(1)**(-rdry/c)-1.)
if(lapsbot==1)bet(1)=-rdry*log(sig(1))
betm(1:kl)=bet(1:kl)
if(lapsbot==2)then     ! may need refinement for non-equal spacing
  do k=2,kl
    bet(k)=.5*rdry*(sig(k-1)-sig(k))/sig(k)
    betm(k)=.5*rdry*(sig(k-1)-sig(k))/sig(k-1)
  enddo
  bet(1)=rdry*(1.-sig(1))/sig(1)
elseif(lapsbot==3)then ! possibly suits nh
  betm(:)=0.
  do k=2,kl
    bet(k)=rdry*log(sig(k-1)/sig(k))
  enddo
  bet(1)=-rdry*log(sig(1))
endif

zmin = 280.*bet(1)/grav
write(6,*) 'zmin = ',zmin

! define map factors
ds = 2.*pi*rearth/(4.*real(il_g))
em(:) = ds*real(il_g)/(schmidt*90.*112.*1000.)

f(:) = real(2._8*2._8*real(pi,8)*sin(real(rlatt(:),8))/86400._8)

! Main prognostic variables
write(6,*) "Initialise main prognostic variables"
t(:,:) = 280.
u(:,:) = 0.
v(:,:) = 0.
qg(:,:) = 0.
qlg(:,:) = 0.
qfg(:,:) = 0.
qsng(:,:) = 0.
qgrg(:,:) = 0.
phi_nh(:,:) = 0.
cfrac(:,:) = 0.
rfrac(:,:) = 0.
sfrac(:,:) = 0.
gfrac(:,:) = 0.
psl(:) = 0.
ps(:) = press_surf
psl(:) = log(ps(:)/1.e5)
land(:) = .true.
tss(:) = 300.
tgg(:,:) = 300.
wb(:,:) = 0.
wbice(:,:) = 0.
tggsn(:,:) = 0.
smass(:,:) = 0.
ssdn(:,:) = 0.
snowd(:) = 0.
ssdnn(:) = 0.
isflag(:) = 0
snage(:) = 0.
fracice(:) = 0.
sicedep(:) = 0.
if ( ccycle>0 ) then
  cplant(:,:) = 0.
  csoil(:,:) = 0.
end if    
isflag(:) = 0
snowd(:) = 0.
tggsn(:,:) = 240.
ssdn(:,:) = 300.
ssdnn(:) = 300.
smass(:,:) = 10.*ssdn(:,:)
!tss(:) = tggsn(:,1)
!tgg(:,:) = 240. 

! Assume constant surface pressure for now
dpsldt(:,:) = 0.

! LOAD HE
zs(:) = grav*z_in
he(:) = 0.

! LOAD INITIAL CONDITIONS
if ( scm_mode=="sublime" ) then
    
  write(6,*) "Loading MET initial conditions"
  call ccnf_open(metforcing,ncid,ncstatus)
    
  call ccnf_inq_dimlen(ncid,'atmospheric_layers',nlev)
 
  allocate( height_in(nlev), height_model(kl) )
  spos(1:3) = (/ 1, 1, 1 /)
  npos(1:3) = (/ 1, 1, nlev /)
  call ccnf_get_vara(ncid,'height',spos(1:3),npos(1:3),height_in)
  
  ! surface pressure
  spos(1:3) = (/ 1, 1, 1 /)
  npos(1:3) = (/ 1, 1, 1 /)
  call ccnf_get_vara(ncid,'PSFC',spos(1:3),npos(1:3),ps)
  psl(:) = log(ps(:)/1.e5)
  
  ! terrain height
  spos(1:3) = (/ 1, 1, 1 /)
  npos(1:3) = (/ 1, 1, 1 /)
  call ccnf_get_vara(ncid,'HGT',spos(1:3),npos(1:3),zs)
  zs(1) = zs(1)*grav ! convert to geopotential
  
  allocate( dat_in(nlev), new_in(nlev) )

  ! potential temperature
  spos(1:3) = (/ 1, 1, 1 /)
  npos(1:3) = (/ 1, 1, nlev /)
  call ccnf_get_vara(ncid,'Theta',spos(1:3),npos(1:3),dat_in)

  ! need to iterate to get consistent solution with temperature
  do i = 1,5
    height_model(1) = bet(1)*t(1,1)/grav
    do k = 2,kl
      height_model(k) = height_model(k-1) + (bet(k)*t(1,k) + betm(k)*t(1,k-1))/grav
    end do
    call vinterp2m(height_in,height_model,dat_in,datout,nlev,kl)
    do k=1,kl
      t(1,k) = datout(k)*(1.e5/(ps(1)*sig(k)))**(-rdry/cp)
    end do
  end do ! iterative i loop
  write(6,*) "height_model ",height_model(1:kl)
  
  ! U wind
  spos(1:3) = (/ 1, 1, 1 /)
  npos(1:3) = (/ 1, 1, nlev /)
  call ccnf_get_vara(ncid,'U',spos(1:3),npos(1:3),dat_in)
  call vinterp2m(height_in,height_model,dat_in,datout,nlev,kl)
  u(1,1:kl) = datout
    
  ! V wind
  spos(1:3) = (/ 1, 1, 1 /)
  npos(1:3) = (/ 1, 1, nlev /)
  call ccnf_get_vara(ncid,'V',spos(1:3),npos(1:3),dat_in)
  call vinterp2m(height_in,height_model,dat_in,datout,nlev,kl)
  v(1,1:kl) = datout
  
  ! qv mixing ratio
  spos(1:3) = (/ 1, 1, 1 /)
  npos(1:3) = (/ 1, 1, nlev /)
  call ccnf_get_vara(ncid,'Qvapor',spos(1:3),npos(1:3),dat_in)
  call vinterp2m(height_in,height_model,dat_in,datout,nlev,kl)
  qg(1,1:kl) = max(datout,2.e-7)
  
  ! ignore tke initial conditions for now
  
  call ccnf_close(ncid)
  
  deallocate(dat_in, new_in )
  deallocate(height_in, height_model)  
  
  tss = t(:,1)
    
else if ( scm_mode=="gabls4" ) then
  if ( .not.lsmonly ) then
    write(6,*) "Loading MET initial conditions"
    call ccnf_open(metforcing,ncid,ncstatus)
  
    call ccnf_inq_dimlen(ncid,'lev',nlev)
    spos(1:2) = (/ 1, 1 /)
    npos(1:2) = (/ nlev, 1 /)
    allocate( dat_in(nlev), sig_in(nlev), new_in(nlev) )
    call ccnf_get_vara(ncid,'psurf',psurf_in(1:1))
    call ccnf_get_vara(ncid,'pf',spos(1:1),npos(1:1),dat_in)
    do k=1,nlev
      sig_in(nlev-k+1) = dat_in(k)/psurf_in(1)
    end do
    write(6,*) "sig_in ",sig_in(:)
    call ccnf_get_vara(ncid,'t',spos(1:1),npos(1:1),new_in)
    do k=1,nlev
      dat_in(nlev-k+1) = new_in(k)
    end do
    call vinterpolate(dat_in,sig_in,datout,sig,nlev,kl)
    t(1,:) = datout(:)
    call ccnf_get_vara(ncid,'qv',spos(1:1),npos(1:1),new_in)
    do k=1,nlev
      dat_in(nlev-k+1) = new_in(k)
    end do
    call vinterpolate(dat_in,sig_in,datout,sig,nlev,kl)
    qg(1,:) = datout(:)/(1.-datout(:))
    call ccnf_get_vara(ncid,'u',spos(1:1),npos(1:1),new_in)
    do k=1,nlev
      dat_in(nlev-k+1) = new_in(k)
    end do
    call vinterpolate(dat_in,sig_in,datout,sig,nlev,kl)
    u(1,:) = datout(:)
    call ccnf_get_vara(ncid,'v',spos(1:1),npos(1:1),new_in)
    do k=1,nlev
      dat_in(nlev-k+1) = new_in(k)
    end do
    call vinterpolate(dat_in,sig_in,datout,sig,nlev,kl)
    v(1,:) = datout(:)
    deallocate( dat_in, sig_in, new_in )
  
    call ccnf_inq_dimlen(ncid,'levs',slev)
    spos(1:2) = (/ 1, 1 /)
    npos(1:2) = (/ slev, 1 /)  
    allocate( dat_in(slev), dep_in(slev), dephl_in(0:slev), thick_in(slev) )  
    call ccnf_get_vara(ncid,'depth_sn',spos(1:1),npos(1:1),dep_in)
    write(6,*) "dep_in ",dep_in(:)
    dephl_in(0) = 0.
    dephl_in(1:slev-1) = 0.5*(dep_in(1:slev-1)+dep_in(2:slev))
    dephl_in(slev) = 2.*dep_in(slev) - dephl_in(slev-1)
    thick_in(1:slev) = dephl_in(1:slev) - dephl_in(0:slev-1)
    sdepth(:,:) = 0.
    tlev=15
    do k=1,3
      sdepth(:,1) = sdepth(:,1) + thick_in(k)
    end do
    do k = 4,tlev
      sdepth(:,2) = sdepth(:,2) + thick_in(k)  
    end do
    do k = tlev+1,slev
      sdepth(:,3) = sdepth(:,3) + thick_in(k)  
    end do
    write(6,*) "sdepth ",sdepth
    call ccnf_get_vara(ncid,'Tsnow',spos(1:1),npos(1:1),dat_in)
    tggsn(:,:) = 0.
    do k=1,3
      tggsn(:,1) = tggsn(:,1) + dat_in(k)*thick_in(k)/sdepth(:,1)
    end do
    do k = 4,tlev
      tggsn(:,2) = tggsn(:,2) + dat_in(k)*thick_in(k)/sdepth(:,2)
    end do
    do k = tlev+1,slev
      tggsn(:,3) = tggsn(:,3) + dat_in(k)*thick_in(k)/sdepth(:,3)
    end do
    write(6,*) "tggsn ",tggsn
    call ccnf_get_vara(ncid,'snow_density',spos(1:1),npos(1:1),dat_in)
    ssdn(:,:) = 0.
    do k=1,3
      ssdn(:,1) = ssdn(:,1) + dat_in(k)*thick_in(k)/sdepth(:,1)
    end do
    do k = 4,tlev
      ssdn(:,2) = ssdn(:,2) + dat_in(k)*thick_in(k)/sdepth(:,2)
    end do
    do k = tlev+1,slev
      ssdn(:,3) = ssdn(:,3) + dat_in(k)*thick_in(k)/sdepth(:,3)
    end do
    write(6,*) "ssdn ",ssdn
    smass(:,:) = sdepth(:,:)*ssdn(:,:)
    write(6,*) "smass ",smass
    deallocate( dat_in, dep_in, dephl_in, thick_in )
    call ccnf_close(ncid)
    isflag(:) = 1
    snowd(:) = sum(sdepth*ssdn,2)
    write(6,*) "snowd ",snowd
    ssdnn(:) = sum(ssdn(:,:)*smass(:,:),2)/snowd
    write(6,*) "ssdnn ",ssdnn  
    tss(:) = tggsn(:,1)
    do k=1,ms
      tgg(:,k) = tggsn(:,3)
    end do
  
  else
    write(6,*) "Loading LSM initial conditions"
    call ccnf_open(lsmforcing,ncid,ncstatus)
    call ccnf_inq_dimlen(ncid,'nsnow_layer',slev)
    spos(1:2) = (/ 1, 1 /)
    npos(1:2) = (/ slev, 1 /) 
    allocate( dat_in(slev), dep_in(slev) )
    call ccnf_get_vara(ncid,'mid_layer_depth',spos(1:1),npos(1:1),dep_in)
    write(6,*) "dep_in ",dep_in(:)
    call ccnf_get_vara(ncid,'snow_temperature',spos(1:1),npos(1:1),dat_in)
    ! interpolate .... unfinished
    call ccnf_get_vara(ncid,'snow_density',spos(1:1),npos(1:1),dat_in)
    ! interpolate .... unfinished
    deallocate( dat_in, dep_in )  
    call ccnf_close(ncid)
    isflag(:) = 1
    snowd(:) = 10.
    tggsn(:,:) = 240. ! fix with initial conditions
    ssdn(:,:) = 300.
    ssdnn(:) = 300.
    smass(:,:) = 10.*ssdn(:,:)
    tss(:) = tggsn(:,1)
    tgg(:,:) = 240.
  end if

else
    
  write(6,*) "ERROR: Unknown option for scm_mode ",trim(scm_mode)
  stop
end if

! UPDATE GRAVITY WAVE DRAG DATA (lgwd)
write(6,*) "Initialise gravity wave drag"
gwdfac = 0.01*lgwd       ! most runs used .02 up to fri  10-10-1997
hefact = 0.1*abs(ngwd)   ! hal used hefact=1. (equiv to ngwd=10)
write(6,*)'hefact,helim,gwdfac: ',hefact,helim,gwdfac
helo(:) = 0.
if ( lgwd>0 ) then
  do iq=1,ifull
    if(land(iq))then
      if(he(iq)==0.)write(6,*)'zero he over land for iq = ',iq
      aa(iq)=min(gwdfac*max(he(iq),.01),.8*zmin)   ! already in meters
      ! replace he by square of corresponding af value
      helo(iq)=( .4/log(zmin/aa(iq)) )**2
    endif
  enddo   ! iq loop
end if ! lgwd>0
if ( ngwd/=0 ) then
!****    limit launching height : Palmer et al use limit on variance of
!****    (400 m)**2. we use launching height = std dev. we limit
!****    launching height to  2*400=800 m. this may be a bit severe.
!****    according to Palmer this prevents 2-grid noise at steep edge of
!****    himalayas etc.
  he(1:ifull) = min(hefact*he(1:ifull),helim)
endif     ! (ngwd/=0)

!--------------------------------------------------------------
! READ SURFACE DATA (nsib and nspecial)
! nsib=3 (original land surface scheme with original 1deg+Dean's datasets)
! nsib=5 (original land surface scheme with MODIS datasets)
! nsib=6 (CABLE land surface scheme with internal screen diagnostics)
! nsib=7 (CABLE land surface scheme with CCAM screen diagnostics)
if (nsib>=1) then
  write(6,*) "Initialise land-surface"
  call insoil
  albvisnir(:,1) = 0.81
  albvisnir(:,2) = 0.81
  albvisdir(:) = 0.81
  albvisdif(:) = 0.81
  albnirdir(:) = 0.81
  albnirdif(:) = 0.81
  rsmin(:) = 300.
  zolnd(:) = 0.001
  vlai(:) = 0.
  ivegt(:) = ivegt_in
  isoilm(:) = isoil_in
  isoilm_in(:) = isoilm(:)
  if (nsib==6.or.nsib==7) then
    ! albvisnir at this point holds soil albedo for cable initialisation
    ivs(:,1) = ivegt(:)
    ivs(:,2:) = 14
    svs(:,1) = 1.
    svs(:,2:) = 0.
    vlinprev(:,1) = vlai(:)
    vlinprev(:,2:) = 0.
    vlin(:,1) = vlai(:)
    vlin(:,2:) = 0.
    vlinnext(:,1) = vlai(:)
    vlinnext(:,2:) = 0.
    vlinnext2(:,1) = vlai(:)
    vlinnext2(:,2:) = 0.
    casapoint(:,:) = 0.
    greenup  = -50
    fall     = 367
    phendoy1 = 2    
    call cbmparm(ivs,svs,vlinprev,vlin,vlinnext,vlinnext2,casapoint,greenup,fall,phendoy1)
    ! albvisnir at this point holds net albedo
  end if
end if

!--------------------------------------------------------------
! UPDATE BIOSPHERE DATA (nsib)
if ( nsib==6 .or. nsib==7 ) then
  ! Load CABLE data
  write(6,*) 'Importing CABLE data'
  call defaulttile
end if


!-----------------------------------------------------------------
! INITIALISE URBAN SCHEME (nurban)
! nurban=0 no urban
! nurban=1 urban (save in restart file)
! nurban=-1 urban (save in history and restart files)
if (nurban/=0) then
  write(6,*) 'Initialising ateb urban scheme'
  sigmu(:) = 1.
  where ( .not.land(1:ifull) .or. sigmu<0.01 )
    sigmu(:) = 0.
  end where
  call atebinit(ifull,sigmu(:),0)
  urbantype = 1
  if ( ateb_bldheight>-900. ) then
    atebparm(1:8) = ateb_bldheight
    call atebdeftype(atebparm(1:8),urbantype,'bldheight',0)
  end if
  if ( ateb_hwratio>-900. ) then
    atebparm(1:8) = ateb_hwratio
    call atebdeftype(atebparm(1:8),urbantype,'hwratio',0)
  end if
  if ( ateb_sigvegc>-900. ) then
    atebparm(1:8) = ateb_sigvegc
    call atebdeftype(atebparm(1:8),urbantype,'sigvegc',0)
  end if
  if ( ateb_sigmabld>-900. ) then
    atebparm(1:8) = ateb_sigmabld
    call atebdeftype(atebparm(1:8),urbantype,'sigmabld',0)
  end if
  if ( ateb_industryfg>-900. ) then
    atebparm(1:8) = ateb_industryfg
    call atebdeftype(atebparm(1:8),urbantype,'industryfg',0)
  end if
  if ( ateb_trafficfg>-900. ) then
    atebparm(1:8) = ateb_trafficfg
    call atebdeftype(atebparm(1:8),urbantype,'trafficfg',0)
  end if
  if ( ateb_wallemiss>-900. ) then
    atebparm(1:8) = ateb_wallemiss
    call atebdeftype(atebparm(1:8),urbantype,'wallemiss',0)
  end if
  if ( ateb_roofemiss>-900. ) then
    atebparm(1:8) = ateb_roofemiss
    call atebdeftype(atebparm(1:8),urbantype,'roofemiss',0)
  end if
  if ( ateb_roademiss>-900. ) then
    atebparm(1:8) = ateb_roademiss
    call atebdeftype(atebparm(1:8),urbantype,'roademiss',0)
  end if
  if ( ateb_wallalpha>-900. ) then
    atebparm(1:8) = ateb_wallalpha
    call atebdeftype(atebparm(1:8),urbantype,'wallalpha',0)
  end if
  if ( ateb_roofalpha>-900. ) then
    atebparm(1:8) = ateb_roofalpha
    call atebdeftype(atebparm(1:8),urbantype,'roofalpha',0)
  end if
  if ( ateb_roadalpha>-900. ) then
    atebparm(1:8) = ateb_roadalpha
    call atebdeftype(atebparm(1:8),urbantype,'roadalpha',0)
  end if
  if ( ateb_vegalphac>-900. ) then
    atebparm(1:8) = ateb_vegalphac
    call atebdeftype(atebparm(1:8),urbantype,'vegalphac',0)
  end if
  if ( ateb_infilach>-900. ) then
    atebparm(1:8) = ateb_infilach
    call atebdeftype(atebparm(1:8),urbantype,'infilach',0)
  end if
  if ( ateb_intgains>-900. ) then
    atebparm(1:8) = ateb_intgains
    call atebdeftype(atebparm(1:8),urbantype,'intgains',0)
  end if
  if ( all( ateb_roof_thick>-900. ) ) then
    do k = 1,4  
      write(vname,'("roofthick",(I1.1))') k  
      atebparm(1:8) = ateb_roof_thick(k)
      call atebdeftype(atebparm(1:8),urbantype,vname,0)
    end do
  end if
  if ( all( ateb_roof_cp>-900. ) ) then
    do k = 1,4  
      write(vname,'("roofcp",(I1.1))') k  
      atebparm(1:8) = ateb_roof_cp(k)
      call atebdeftype(atebparm(1:8),urbantype,vname,0)
    end do
  end if
  if ( all( ateb_roof_cond>-900. ) ) then
    do k = 1,4  
      write(vname,'("roofcond",(I1.1))') k  
      atebparm(1:8) = ateb_roof_cond(k)
      call atebdeftype(atebparm(1:8),urbantype,vname,0)
    end do
  end if
  if ( all( ateb_wall_thick>-900. ) ) then
    do k = 1,4  
      write(vname,'("wallthick",(I1.1))') k  
      atebparm(1:8) = ateb_wall_thick(k)
      call atebdeftype(atebparm(1:8),urbantype,vname,0)
    end do
  end if
  if ( all( ateb_wall_cp>-900. ) ) then
    do k = 1,4  
      write(vname,'("wallcp",(I1.1))') k  
      atebparm(1:8) = ateb_wall_cp(k)
      call atebdeftype(atebparm(1:8),urbantype,vname,0)
    end do
  end if
  if ( all( ateb_wall_cond>-900. ) ) then
    do k = 1,4  
      write(vname,'("wallcond",(I1.1))') k  
      atebparm(1:8) = ateb_wall_cond(k)
      call atebdeftype(atebparm(1:8),urbantype,vname,0)
    end do
  end if
  if ( all( ateb_road_thick>-900. ) ) then
    do k = 1,4  
      write(vname,'("roadthick",(I1.1))') k  
      atebparm(1:8) = ateb_road_thick(k)
      call atebdeftype(atebparm(1:8),urbantype,vname,0)
    end do
  end if
  if ( all( ateb_road_cp>-900. ) ) then
    do k = 1,4  
      write(vname,'("roadcp",(I1.1))') k  
      atebparm(1:8) = ateb_road_cp(k)
      call atebdeftype(atebparm(1:8),urbantype,vname,0)
    end do
  end if
  if ( all( ateb_road_cond>-900. ) ) then
    do k = 1,4  
      write(vname,'("roadcond",(I1.1))') k  
      atebparm(1:8) = ateb_road_cond(k)
      call atebdeftype(atebparm(1:8),urbantype,vname,0)
    end do
  end if
  if ( all( ateb_slab_thick>-900. ) ) then
    do k = 1,4  
      write(vname,'("slabthick",(I1.1))') k  
      atebparm(1:8) = ateb_slab_thick(k)
      call atebdeftype(atebparm(1:8),urbantype,vname,0)
    end do
  end if
  if ( all( ateb_slab_cp>-900. ) ) then
    do k = 1,4  
      write(vname,'("slabcp",(I1.1))') k  
      atebparm(1:8) = ateb_slab_cp(k)
      call atebdeftype(atebparm(1:8),urbantype,vname,0)
    end do
  end if
  if ( all( ateb_slab_cond>-900. ) ) then
    do k = 1,4  
      write(vname,'("slabcond",(I1.1))') k  
      atebparm(1:8) = ateb_slab_cond(k)
      call atebdeftype(atebparm(1:8),urbantype,vname,0)
    end do
  end if
else
  sigmu(:) = 0.
  call atebdisable(0) ! disable urban
end if

!-----------------------------------------------------------------
! INTIALISE MIXED LAYER OCEAN (nmlo)
! nmlo<0 same as below, but save all data in history file
! nmlo=0 no mixed layer ocean
! nmlo=1 mixed layer ocean (KPP)
! nmlo=2 same as 1, but with Smag horz diffusion and river routing
! nmlo=3 same as 2, but with horizontal and vertical advection
if (nmlo/=0.and.abs(nmlo)<=9) then
  write(6,*) 'Initialising MLO'
  depth=200.
  where (land)
    depth=0.
  elsewhere
    depth=max(depth,2.*minwater)
  end where
  call mloinit(ifull,depth,0)
end if

!-----------------------------------------------------------------
! INITIALISE PROGNOSTIC AEROSOLS (iaero)
! iaero=1 prescribed aerosol direct effect
! iaero=2 LDR prognostic aerosol (direct only with nrad=4, direct+indirect with nrad=5)
select case(abs(iaero))
  case(2)
    write(6,*) "Initialise aerosols"
    write(6,*) "ERROR: Prognostic aerosols not implemented in SCM"
    stop -1
!    call load_aerosolldr(so4tfile,oxidantfile,kdate)
end select


!-----------------------------------------------------------------
! UPDATE MIXED LAYER OCEAN DATA (nmlo)
if ( nmlo/=0 .and. abs(nmlo)<=9 ) then
  write(6,*) 'Using MLO defaults'
  allocate( micdwn(ifull,11) )
  ocndwn(:,2) = 0.
  do k = 1,wlev
    call mloexpdep(0,depth,k,0)
    ! This polynomial fit is from MOM3, based on Levitus
    where (depth(:)<2000.)
      mlodwn(:,k,1) = 18.4231944       &
        - 0.43030662E-1*depth(:)       &
        + 0.607121504E-4*depth(:)**2   &
        - 0.523806281E-7*depth(:)**3   &
        + 0.272989082E-10*depth(:)**4  &
        - 0.833224666E-14*depth(:)**5  &
        + 0.136974583E-17*depth(:)**6  &
        - 0.935923382E-22*depth(:)**7
      mlodwn(:,k,1) = mlodwn(:,k,1)*(tss-273.16)/18.4231944 + 273.16 - wrtemp
    elsewhere
      mlodwn(:,k,1)= 275.16 - wrtemp
    end where
    mlodwn(:,k,2) = 34.72
    mlodwn(:,k,3:4) = 0.
  end do
  micdwn(:,1) = tss
  micdwn(:,2) = tss
  micdwn(:,3) = tss
  micdwn(:,4) = tss
  micdwn(:,5:6) = 0.
  where ( fracice(:)>0. )
    micdwn(:,1) = min(tss,271.)
    micdwn(:,2) = min(tss,271.)
    micdwn(:,3) = min(tss,271.)
    micdwn(:,4) = min(tss,271.)
    micdwn(:,5) = 1.
    micdwn(:,6) = 2.
  end where
  micdwn(:,7:10) = 0.
  micdwn(:,11) = 10.
  where ( .not.land )
    fracice = micdwn(:,5)
    sicedep = micdwn(:,6)
    snowd = micdwn(:,7)*1000.
  end where          
  !mlodwn(1:ifull,1:wlev,1)=max(mlodwn(1:ifull,1:wlev,1),271.-wrtemp)
  mlodwn(1:ifull,1:wlev,2) = max(mlodwn(1:ifull,1:wlev,2),0.)
  micdwn(1:ifull,1:4) = min(max(micdwn(1:ifull,1:4),0.),300.)
  micdwn(1:ifull,11) = max(micdwn(1:ifull,11),0.)
  call mloload(mlodwn,ocndwn(:,2),micdwn,0)
  deallocate( micdwn )
  do k = 1,ms
    call mloexport(0,tgg(:,k),k,0)
    where ( tgg(:,k)<100. )
      tgg(:,k) = tgg(:,k) + wrtemp
    end where    
  end do
  do k = 1,3
    call mloexpice(tggsn(:,k),k,0)
  end do 
end if

!-----------------------------------------------------------------
! UPDATE URBAN DATA (nurban)
if (nurban/=0) then
  write(6,*) 'Importing ateb urban data'
  allocate( atebdwn(ifull,5) )
  atebdwn(:,1)=289.                     ! roof temp 1
  atebdwn(:,2)=(0.75*289.)+(0.25*293.)  ! roof temp 2
  atebdwn(:,3)=(0.50*289.)+(0.50*293.)  ! roof temp 3
  atebdwn(:,4)=(0.25*289.)+(0.75*293.)  ! roof temp 4
  atebdwn(:,5)=293.                     ! roof temp 5
  do k = 1,5
    write(vname,'("rooftemp",I1.1)') k  
    atebdwn(:,k) = atebdwn(:,k) - urbtemp
    call atebloadd(atebdwn(:,k),vname,0)
  end do 

  atebdwn(:,1)=290.5                    ! walleast temp 1
  atebdwn(:,2)=(0.75*290.5)+(0.25*293.) ! walleast temp 2
  atebdwn(:,3)=(0.50*290.5)+(0.50*293.) ! walleast temp 3
  atebdwn(:,4)=(0.25*290.5)+(0.75*293.) ! walleast temp 4
  atebdwn(:,5)=293.                     ! walleast temp 5
  do k = 1,5
    write(vname,'("walletemp",I1.1)') k  
    atebdwn(:,k) = atebdwn(:,k) - urbtemp
    call atebloadd(atebdwn(:,k),vname,0)
  end do

  atebdwn(:,1)=290.5                    ! wallwest temp 1
  atebdwn(:,2)=(0.75*290.5)+(0.25*293.) ! wallwest temp 2
  atebdwn(:,3)=(0.50*290.5)+(0.50*293.) ! wallwest temp 3
  atebdwn(:,4)=(0.25*290.5)+(0.75*293.) ! wallwest temp 4
  atebdwn(:,5)=293.                     ! wallwest temp 5
  do k = 1,5
    write(vname,'("wallwtemp",I1.1)') k  
    atebdwn(:,k) = atebdwn(:,k) - urbtemp
    call atebloadd(atebdwn(:,k),vname,0)
  end do

  atebdwn(:,1)=289.                      ! road temp 1
  atebdwn(:,2)=(0.966*289.)+(0.034*287.) ! road temp 2
  atebdwn(:,3)=(0.833*289.)+(0.167*287.) ! road temp 3
  atebdwn(:,4)=(0.533*289.)+(0.467*287.) ! road temp 4
  atebdwn(:,5)=287.                      ! road temp 5
  do k = 1,5
    write(vname,'("roadtemp",I1.1)') k  
    atebdwn(:,k) = atebdwn(:,k) - urbtemp
    call atebloadd(atebdwn(:,k),vname,0)
  end do

  atebdwn(:,1)=293.          ! slab temp 1 (internal air temp)
  atebdwn(:,2)=293.          ! slab temp 2 (internal air temp)
  atebdwn(:,3)=293.          ! slab temp 3 (internal air temp) 
  atebdwn(:,4)=293.          ! slab temp 4 (internal air temp)
  atebdwn(:,5)=293.          ! slab temp 5 (internal air temp)
  do k = 1,5
    write(vname,'("slabtemp",I1.1)') k  
    atebdwn(:,k) = atebdwn(:,k) - urbtemp
    call atebloadd(atebdwn(:,k),vname,0)
  end do

  atebdwn(:,1)=293.          ! intm temp 1 (internal air temp)
  atebdwn(:,2)=293.          ! intm temp 2 (internal air temp)
  atebdwn(:,3)=293.          ! intm temp 3 (internal air temp)
  atebdwn(:,4)=293.          ! intm temp 4 (internal air temp)
  atebdwn(:,5)=293.          ! intm temp 5 (internal air temp)
  do k = 1,5
    write(vname,'("intmtemp",I1.1)') k  
    atebdwn(:,k) = atebdwn(:,k) - urbtemp
    call atebloadd(atebdwn(:,k),vname,0)
  end do
  atebdwn(:,1) = 293.             ! room air temp
  call atebloadd(atebdwn(:,1),"roomtemp",0)
  atebdwn(:,1)=0.5*0.26+0.5*0.34  ! Soil water road
  call atebloadd(atebdwn(:,1),"canyonsoilmoisture",0)
  atebdwn(:,1)=0.18               ! Green roof water
  call atebloadd(atebdwn(:,1),"roofsoilmoisture",0)
  atebdwn(:,1)=0.   ! road water
  call atebloadd(atebdwn(:,1),"roadsurfacewater",0)
  atebdwn(:,1)=0.   ! roof water
  call atebloadd(atebdwn(:,1),"roofsurfacewater",0)
  atebdwn(:,1)=0.   ! canyon leaf water
  call atebloadd(atebdwn(:,1),"canyonleafwater",0)
  atebdwn(:,1)=0.   ! roof leaf water
  call atebloadd(atebdwn(:,1),"roofleafwater",0)
  atebdwn(:,1)=0.   ! road snow
  call atebloadd(atebdwn(:,1),"roadsnowdepth",0)
  atebdwn(:,1)=0.   ! roof snow
  call atebloadd(atebdwn(:,1),"roofsnowdepth",0)
  atebdwn(:,1)=100. ! road snow density
  call atebloadd(atebdwn(:,1),"roadsnowdensity",0)
  atebdwn(:,1)=100. ! roof snow density
  call atebloadd(atebdwn(:,1),"roofsnowdensity",0)
  atebdwn(:,1)=0.85 ! road snow albedo
  call atebloadd(atebdwn(:,1),"roadsnowalbedo",0)
  atebdwn(:,1)=0.85 ! roof snow albedo
  call atebloadd(atebdwn(:,1),"roofsnowalbedo",0)
  deallocate( atebdwn )
end if


if ( abs(iaero)>=2 ) then
   write(6,*) "Aerosol initial conditions"
   xtg(:,:,:) = 0.
end if
    
write(6,*) "Finised initialisation"

return
end subroutine initialscm
    
subroutine nudgescm(scm_mode,metforcing,fixtsurf,iarch_nudge)

use arrays_m                               ! Atmosphere dyamics prognostic arrays
use cable_ccam
use const_phys                             ! Physical constants
use infile                                 ! Input file routines
use map_m                                  ! Grid map arrays
use newmpar_m                              ! Grid parameters
use parm_m                                 ! Model configuration
use pbl_m                                  ! Boundary layer arrays
use sigs_m                                 ! Atmosphere sigma levels
use soilsnow_m                             ! Soil, snow and surface data

use scmarrays_m

implicit none

integer, intent(inout) :: iarch_nudge
integer, save :: ncid
integer, save :: nlev, ntimes
integer :: ncstatus, k, l
integer, dimension(4) :: spos, npos
real, save :: time_a = -1.
real, save :: time_b = -1.
real :: time_ktau, x
real, dimension(1) :: tset
real, dimension(:), allocatable, save :: t_force_a, q_force_a, u_force_a, v_force_a
real, dimension(:), allocatable, save :: t_force_b, q_force_b, u_force_b, v_force_b
real, dimension(:), allocatable, save :: new_in, dat_in, sig_in
real, dimension(:), allocatable, save :: ug_file_a, vg_file_a, time_file, height_file, height_file_a, height_file_b
real, dimension(:), allocatable, save :: ug_file_b, vg_file_b, ug_file, vg_file
real, dimension(:), allocatable, save :: theta_adv, qv_adv, u_adv, v_adv
real, dimension(:), allocatable, save :: theta_adv_a, theta_adv_b, qv_adv_a, qv_adv_b
real, dimension(:), allocatable, save :: u_adv_a, u_adv_b, v_adv_a, v_adv_b 
real, dimension(1) :: psurf_in
real, dimension(kl) :: tadv, qadv, um, vm, height_model
real, save :: tsurf_a, tsurf_b
character(len=*), intent(in) :: scm_mode, metforcing
logical, intent(in) :: fixtsurf

time_ktau = real(ktau)*dt

if ( scm_mode=="sublime" ) then
    
  if ( iarch_nudge==0 ) then
    iarch_nudge = 1 
    if ( allocated(ug) ) then
      deallocate( ug, vg, t_tend, q_tend )
      deallocate( uadv, vadv )
      deallocate( ug_file_b, vg_file_b, time_file )
      deallocate( ug_file_a, vg_file_a, ug_file, vg_file )
      deallocate( theta_adv, qv_adv, u_adv, v_adv )
      deallocate( theta_adv_a, theta_adv_b, qv_adv_a, qv_adv_b )
      deallocate( u_adv_a, u_adv_b, v_adv_a, v_adv_b )
      deallocate( height_file_a, height_file_b, height_file )
      call ccnf_close(ncid)
    end if
    time_a = -1.
    time_b = -1.

    allocate( ug(kl), vg(kl), t_tend(kl), q_tend(kl) )
    allocate( uadv(kl), vadv(kl) )
    
    call ccnf_open(metforcing,ncid,ncstatus)
    call ccnf_inq_dimlen(ncid,'force_layers',nlev)
    
    call ccnf_inq_dimlen(ncid,'time',ntimes)
    
    allocate( ug_file_b(nlev), vg_file_b(nlev), time_file(ntimes) )
    allocate( ug_file_a(nlev), vg_file_a(nlev), ug_file(nlev), vg_file(nlev) )
    allocate( theta_adv(nlev), qv_adv(nlev), u_adv(nlev), v_adv(nlev) )
    allocate( theta_adv_a(nlev), theta_adv_b(nlev), qv_adv_a(nlev), qv_adv_b(nlev) )
    allocate( u_adv_a(nlev), u_adv_b(nlev), v_adv_a(nlev), v_adv_b(nlev) )
    allocate( height_file_a(nlev), height_file_b(nlev), height_file(nlev) )
    
    do l = 1,ntimes
      time_file(l) = real(l-1)*3600.*(54./real(ntimes-1))
    end do

    spos(1:4) = (/ 1, 1, 1, iarch_nudge /)
    npos(1:4) = (/ 1, 1, nlev, 1 /)
    call ccnf_get_vara(ncid,'Forcing_height',spos,npos,height_file_b)
    
    ! geostropic winds
    call ccnf_get_vara(ncid,'Ug',spos,npos,ug_file_b)
    call ccnf_get_vara(ncid,'Vg',spos,npos,vg_file_b)

    ! theta advection
    call ccnf_get_vara(ncid,'Th_advection',spos,npos,theta_adv_b)
    
    ! water vapor mixing ratio advection
    call ccnf_get_vara(ncid,'Qv_advection',spos,npos,qv_adv_b)
    
    ! momentum advection
    call ccnf_get_vara(ncid,'U_advection',spos,npos,u_adv_b)
    call ccnf_get_vara(ncid,'V_advection',spos,npos,v_adv_b)
    
  end if
  
  if ( time_ktau>time_b ) then
    write(6,*) "Updating MET forcing"  
    iarch_nudge = iarch_nudge + 1
    time_a = time_file(iarch_nudge-1)
    time_b = time_file(iarch_nudge)
    
    height_file_a = height_file_b
    spos(1:4) = (/ 1, 1, 1, iarch_nudge /)
    npos(1:4) = (/ 1, 1, nlev, 1 /)
    call ccnf_get_vara(ncid,'Forcing_height',spos,npos,height_file_b)

    ug_file_a = ug_file_b
    call ccnf_get_vara(ncid,'Ug',spos,npos,ug_file_b)
    vg_file_a = vg_file_b
    call ccnf_get_vara(ncid,'Vg',spos,npos,vg_file_b)
    theta_adv_a = theta_adv_b
    call ccnf_get_vara(ncid,'Th_advection',spos,npos,theta_adv_b)
    qv_adv_a = qv_adv_b
    call ccnf_get_vara(ncid,'Qv_advection',spos,npos,qv_adv_b)    
    u_adv_a = u_adv_b    
    call ccnf_get_vara(ncid,'U_advection',spos,npos,u_adv_b)
    v_adv_a = v_adv_b
    call ccnf_get_vara(ncid,'V_advection',spos,npos,v_adv_b)
    
  end if

  height_model(1) = bet(1)*t(1,1)
  do k = 2,kl
    height_model(k) = height_model(k-1) + bet(k)*t(1,k) + betm(k)*t(1,k-1)
  end do
    
  x = (time_ktau - time_a)/(time_b-time_a)
  height_file(:) = (1.-x)*height_file_a(:) + x*height_file_b(:)
  ug_file(:) = (1.-x)*ug_file_a(:) + x*ug_file_b(:)
  vg_file(:) = (1.-x)*vg_file_a(:) + x*vg_file_b(:)
  theta_adv(:) = (1.-x)*theta_adv_a(:) + x*theta_adv_b(:)
  qv_adv(:) = (1.-x)*qv_adv_a(:) + x*qv_adv_b(:)
  u_adv(:) = (1.-x)*u_adv_a(:) + x*u_adv_b(:)
  v_adv(:) = (1.-x)*v_adv_a(:) + x*v_adv_b(:)
  
  call vinterp2m(height_file,height_model,ug_file,ug,nlev,kl)
  call vinterp2m(height_file,height_model,vg_file,vg,nlev,kl)
  call vinterp2m(height_file,height_model,u_adv,uadv,nlev,kl)
  call vinterp2m(height_file,height_model,v_adv,vadv,nlev,kl)
  call vinterp2m(height_file,height_model,theta_adv,tadv,nlev,kl)
  call vinterp2m(height_file,height_model,qv_adv,qadv,nlev,kl)
  do k = 1,kl
    tadv(k) = tadv(k)*(1.e5/(sig(k)*ps(1)))**(-rdry/cp) ! convert from potential temperature to temperature
  end do  
  
elseif ( scm_mode=="gabls4" ) then
    
  uadv(:) = 0.
  vadv(:) = 0.
    
  if ( iarch_nudge==0 ) then
    iarch_nudge = 1
    if ( allocated( t_force_a ) ) then
      deallocate( t_force_a, q_force_a, u_force_a, v_force_a )
      deallocate( t_force_b, q_force_b, u_force_b, v_force_b )
      deallocate( new_in, dat_in, sig_in )
      deallocate( ug, vg, t_tend, q_tend )
      call ccnf_close(ncid)  
    end if
    
    write(6,*) "Starting MET forcing"
    call ccnf_open(metforcing,ncid,ncstatus)
    call ccnf_inq_dimlen(ncid,'lev',nlev)
    spos(1:2) = (/ 1, iarch_nudge  /)
    npos(1:2) = (/ nlev, 1 /)  
    allocate( t_force_a(kl), q_force_a(kl), u_force_a(kl), v_force_a(kl) )
    allocate( t_force_b(kl), q_force_b(kl), u_force_b(kl), v_force_b(kl) )
    allocate( new_in(nlev), dat_in(nlev), sig_in(nlev) )
    allocate( ug(kl), vg(kl), t_tend(kl), q_tend(kl) )
    call ccnf_get_vara(ncid,'psurf',psurf_in(1:1))
    call ccnf_get_vara(ncid,'pf',spos(1:1),npos(1:1),dat_in)
    do k=1,nlev
      sig_in(nlev-k+1) = dat_in(k)/psurf_in(1)
    end do  
    call ccnf_get_vara(ncid,'time',iarch_nudge,time_b)  
    call ccnf_get_vara(ncid,'Tg',iarch_nudge,tsurf_b) 
    call ccnf_get_vara(ncid,'Ug',spos,npos,new_in)
    do k = 1,nlev
      dat_in(nlev-k+1) = new_in(k)
    end do
    call vinterpolate(dat_in,sig_in,u_force_b,sig,nlev,kl)
    call ccnf_get_vara(ncid,'Vg',spos,npos,new_in)
    do k = 1,nlev
      dat_in(nlev-k+1) = new_in(k)
    end do
    call vinterpolate(dat_in,sig_in,v_force_b,sig,nlev,kl)
    call ccnf_get_vara(ncid,'hadvT',spos,npos,new_in)
    do k = 1,nlev
      dat_in(nlev-k+1) = new_in(k)
    end do
    call vinterpolate(dat_in,sig_in,t_force_b,sig,nlev,kl)
    call ccnf_get_vara(ncid,'hadvQ',spos,npos,new_in)
    do k = 1,nlev
      dat_in(nlev-k+1) = new_in(k)
    end do
    call vinterpolate(dat_in,sig_in,q_force_b,sig,nlev,kl)
    ug(:) = u_force_b(:)
    vg(:) = v_force_b(:)
    t_tend(:) = t_force_b(:)
    q_tend(:) = q_force_b(:)
  end if

  if ( time_ktau>time_b ) then
    write(6,*) "Updating MET forcing"  
    iarch_nudge = iarch_nudge + 1
    spos(1:2) = (/ 1, iarch_nudge  /)
    npos(1:2) = (/ nlev, 1 /)  
    time_a = time_b
    tsurf_a = tsurf_b
    t_force_a(:) = t_force_b(:)
    q_force_a(:) = q_force_b(:)
    u_force_a(:) = u_force_b(:)
    v_force_a(:) = v_force_b(:)
    call ccnf_get_vara(ncid,'time',iarch_nudge,time_b)  
    call ccnf_get_vara(ncid,'Tg',iarch_nudge,tsurf_b) 
    call ccnf_get_vara(ncid,'Ug',spos(1:2),npos(1:2),new_in)
    do k = 1,nlev
      dat_in(nlev-k+1) = new_in(k)
    end do
    call vinterpolate(dat_in,sig_in,u_force_b,sig,nlev,kl)
    call ccnf_get_vara(ncid,'Vg',spos(1:2),npos(1:2),new_in)
    do k = 1,nlev
      dat_in(nlev-k+1) = new_in(k)
    end do
    call vinterpolate(dat_in,sig_in,v_force_b,sig,nlev,kl)
    call ccnf_get_vara(ncid,'hadvT',spos(1:2),npos(1:2),new_in)
    do k = 1,nlev
      dat_in(nlev-k+1) = new_in(k)
    end do
    call vinterpolate(dat_in,sig_in,t_force_b,sig,nlev,kl)
    call ccnf_get_vara(ncid,'hadvQ',spos(1:2),npos(1:2),new_in)
    do k = 1,nlev
      dat_in(nlev-k+1) = new_in(k)
    end do
    call vinterpolate(dat_in,sig_in,q_force_b,sig,nlev,kl)
  end if

  if ( ktau>0 ) then
    x = (time_ktau - time_a)/(time_b-time_a)
    ug(:) = (1.-x)*u_force_a(:) + x*u_force_b(:)
    vg(:) = (1.-x)*v_force_a(:) + x*v_force_b(:)
    tadv(:) = (1.-x)*t_force_a(:) + x*t_force_b(:)
    qadv(:) = (1.-x)*q_force_a(:) + x*q_force_b(:)

    if ( fixtsurf ) then
      tset(:) = (1.-x)*tsurf_a + x*tsurf_b
      tggsn(:,:) = tset(1)
      tss(:) = tset(1)
      tgg(:,:) = tset(1)
      if (nsib==6.or.nsib==7) then
        call cablesettemp(tset)
      end if
    end if
  end if

else
    
  write(6,*) "ERROR: Invalid option for scm_mode in nudging ",trim(scm_mode)
  stop
  
end if

if ( ktau>0 ) then
    
  um(:) = u(1,:)
  vm(:) = v(1,:)
  u(1,:) = (um(:)+dt*f(1)*(vm(:)-vg(:))+(dt*f(1))**2*ug(:))/(1.+(dt*f(1))**2)
  v(1,:) = (vm(:)-dt*f(1)*(um(:)-ug(:))+(dt*f(1))**2*vg(:))/(1.+(dt*f(1))**2)

  u(1,:) = u(1,:) + dt*uadv(:)
  v(1,:) = v(1,:) + dt*vadv(:)
  t(1,:) = t(1,:) + dt*tadv(:)
  qg(1,:) = qg(1,:) + dt*qadv(:)
  
  t_tend(:) = tadv(:)
  q_tend(:) = qadv(:)
    
end if


return
end subroutine nudgescm

subroutine replace_scm(lsmforcing)

use arrays_m                               ! Atmosphere dyamics prognostic arrays
use estab
use extraout_m                             ! Additional diagnostics
use infile                                 ! Input file routines
use morepbl_m                              ! Additional boundary layer diagnostics
use newmpar_m                              ! Grid parameters
use parm_m                                 ! Model configuration
use sigs_m                                 ! Atmosphere sigma levels
use soilsnow_m                             ! Soil, snow and surface data

implicit none

integer, save :: ncid
integer iarch, ncstatus
real t_lsm_a, t_lsm_b
real u_lsm_a, u_lsm_b
real v_lsm_a, v_lsm_b
real rh_lsm_a, rh_lsm_b
real sgdwn_lsm_a, sgdwn_lsm_b
real rgdwn_lsm_a, rgdwn_lsm_b
real pr_lsm_a, pr_lsm_b
real ps_lsm_a, ps_lsm_b
real, dimension(kl) :: qs, pf, ta
character(len=*), intent(in) :: lsmforcing

if ( dt/=1800 ) then
  write(6,*) "ERROR: Currently lsm_only mode requires dt=1800"
  stop
end if

if ( ktau==1 ) then
  call ccnf_open(lsmforcing,ncid,ncstatus)
end if    
 
! average variables over the integration period?
iarch = ktau
call ccnf_get_vara(ncid,'Temperature',iarch,t_lsm_a) 
call ccnf_get_vara(ncid,'Relative_humidity',iarch,rh_lsm_a)
call ccnf_get_vara(ncid,'Uwind',iarch,u_lsm_a)
call ccnf_get_vara(ncid,'Vwind',iarch,v_lsm_a)
call ccnf_get_vara(ncid,'SW_down',iarch,sgdwn_lsm_a)
call ccnf_get_vara(ncid,'LW_down',iarch,rgdwn_lsm_a)
call ccnf_get_vara(ncid,'Rain',iarch,pr_lsm_a)
call ccnf_get_vara(ncid,'PSFC',iarch,ps_lsm_a)

iarch = ktau + 1
call ccnf_get_vara(ncid,'Temperature',iarch,t_lsm_b) 
call ccnf_get_vara(ncid,'Relative_humidity',iarch,rh_lsm_b)
call ccnf_get_vara(ncid,'Uwind',iarch,u_lsm_b)
call ccnf_get_vara(ncid,'Vwind',iarch,v_lsm_b)
call ccnf_get_vara(ncid,'SW_down',iarch,sgdwn_lsm_b)
call ccnf_get_vara(ncid,'LW_down',iarch,rgdwn_lsm_b)
call ccnf_get_vara(ncid,'Rain',iarch,pr_lsm_b)
call ccnf_get_vara(ncid,'PSFC',iarch,ps_lsm_b)

t(1,1) = 0.5*(t_lsm_a+t_lsm_b)
u(1,1) = 0.5*(u_lsm_a+u_lsm_b)
v(1,1) = 0.5*(v_lsm_a+v_lsm_b)
sgsave(1) = 0.5*(sgdwn_lsm_a+sgdwn_lsm_b)*(1.-swrsave(1)*albvisnir(1,1)-swrsave(1)*albvisnir(1,2))
rgsave(1) = -0.5*(rgdwn_lsm_a+rgdwn_lsm_b)
condx(1) = 0.5*(pr_lsm_a+pr_lsm_b) 
ps(1) = 0.5*(ps_lsm_a+ps_lsm_b)
psl(1) = log(ps(1)/1.e5)
pf(1) = sig(1)*ps(1)
ta(1) = t(1,1)
qs(1) = qsat(pf(1),ta(1))
qg(1,1) = qs(1)*0.5*(rh_lsm_a+rh_lsm_b)/100.

return
end subroutine replace_scm
    
subroutine outputscm(scm_mode,timeoutput,profileoutput,lsmoutput)

use ateb
use arrays_m                               ! Atmosphere dynamics prognostic arrays
use const_phys                             ! Physical constants
use dates_m                                ! Date data
use estab
use extraout_m                             ! Additional diagnostics
use infile                                 ! Input file routines
use liqwpar_m                              ! Cloud water mixing ratios
use morepbl_m                              ! Additional boundary layer diagnostics
use newmpar_m                              ! Grid parameters
use nsibd_m                                ! Land-surface arrays
use parm_m                                 ! Model configuration
use pbl_m                                  ! Boundary layer arrays
use prec_m                                 ! Precipitation
use raddiag_m                              ! Radiation diagnostic
use sigs_m                                 ! Atmosphere sigma levels
use screen_m                               ! Screen level diagnostics
use soilsnow_m                             ! Soil, snow and surface data
use work2_m                                ! Diagnostic arrays

use scmarrays_m

implicit none

real, parameter :: nf90_fill_float = 9.9692099683868690e+36
real, dimension(1) :: aa      ! surface
real, dimension(1,kl) :: bb   ! full levels
real, dimension(1,kl+1) :: cc ! half levels
real, dimension(1,3) :: dd    ! snow
real, dimension(1,4) :: uu    ! urban
real, dimension(1,kl+1) :: wtflux
real, dimension(kl) :: zf, pf, rh
real, dimension(kl+1) :: zh
real, dimension(kl) :: qs, tmp
integer, save :: timencid, profilencid, lsmncid
integer, save :: iarch = 1
integer zfdim, zhdim, zsdim, tdim_time, tdim_prof, udim
integer, dimension(2) :: jdim
integer, dimension(2) :: spos, npos
integer icy, icm, icd, ich, icmi, ics
integer idnt
integer k
character(len=*), intent(in) :: timeoutput, profileoutput, lsmoutput, scm_mode
character(len=40) :: lname
character(len=33) :: grdtim
logical :: firstcall, lastcall

firstcall = ktau==0
lastcall = ktau==ntau

if ( scm_mode=="sublime" ) then

  if ( firstcall ) then
    write(6,*) "Creating time output file"  
    call ccnf_create(timeoutput,timencid)
    ! Turn off the data filling
    call ccnf_nofill(timencid)
    
    ! attributes
    call ccnf_put_attg(timencid,'model','CCAM+UCLEM')
    call ccnf_put_attg(timencid,'contact','Mathew Lipson <m.lipson@unsw.edu.au>')
    call ccnf_put_attg(timencid,'scmtype','Stretched grid climate model')
    call ccnf_put_attg(timencid,'ucmcomplexity','complex')
    call ccnf_put_attg(timencid,'ucmtiled','no')
    call ccnf_put_attg(timencid,'timestep',dt)
    
    if ( nvmix==7 ) then
      call ccnf_put_attg(timencid,'turbulencescheme','Jing')  
      call ccnf_put_attg(timencid,'eddydiff','?')
      call ccnf_put_attg(timencid,'massflux','?')
      call ccnf_put_attg(timencid,'lengthscale','?')
      call ccnf_put_attg(timencid,'kprofile','?')
    else if ( nvmix==6 ) then
      call ccnf_put_attg(timencid,'turbulencescheme','EDMF')
      call ccnf_put_attg(timencid,'eddydiff','TKE-EPS')
      call ccnf_put_attg(timencid,'massflux','updraft')
    else
      call ccnf_put_attg(timencid,'turbulencescheme','K-profile')  
      call ccnf_put_attg(timencid,'eddydiff','?')
      call ccnf_put_attg(timencid,'massflux','diagnosed')
      call ccnf_put_attg(timencid,'lengthscale','?')
      call ccnf_put_attg(timencid,'kprofile','Richardson')
    end if
    
    ! call ccnf_put_attg(timencid,'soillayers','6')
    ! call ccnf_put_attg(timencid,'snowlayers','1-3')
    ! call ccnf_put_attg(timencid,'surfaceprog','temp,moisture,ice,canopywater,snowtemp,snowmass,snowage')
    
    call ccnf_def_dim(timencid,'urban_layer',4,udim)
    call ccnf_def_dimu(timencid,'time',tdim_time)
    
    jdim(1) = udim
    call ccnf_def_var(timencid,'urban_layer','float',1,jdim(1:1),idnt)

    jdim(1) = tdim_time
    call ccnf_def_var(timencid,'time','float',1,jdim(1:1),idnt)
    call ccnf_put_att(timencid,idnt,'point_spacing','even')
    icy = kdate/10000
    icm = max(1, min(12, (kdate-icy*10000)/100))
    icd = max(1, min(31, (kdate-icy*10000-icm*100)))
    ich = ktime/100
    icmi = (ktime-ich*100)
    ics = 0
    write(grdtim,'("seconds since ",i4.4,"-",i2.2,"-",i2.2," ",2(i2.2,":"),i2.2)') icy,icm,icd,ich,icmi,ics
    call ccnf_put_att(timencid,idnt,'units',grdtim)    
    
    
    jdim(1) = tdim_time
    lname = 'longwave downward radiation at surface'
    call ccnf_def_var(timencid,'ldw','float',1,jdim(1:1),idnt)
    call ccnf_put_att(timencid,idnt,'long_name',lname)
    call ccnf_put_att(timencid,idnt,'units','W/m2')
    call ccnf_put_att(timencid,idnt,'missing_value',nf90_fill_float)
    lname = 'longwave upward radiation at surface'
    call ccnf_def_var(timencid,'lup','float',1,jdim(1:1),idnt)
    call ccnf_put_att(timencid,idnt,'long_name',lname)
    call ccnf_put_att(timencid,idnt,'units','W/m2')
    call ccnf_put_att(timencid,idnt,'missing_value',nf90_fill_float)
    lname = 'shortwave downward radiation at surface'
    call ccnf_def_var(timencid,'qdw','float',1,jdim(1:1),idnt)
    call ccnf_put_att(timencid,idnt,'long_name',lname)
    call ccnf_put_att(timencid,idnt,'units','W/m2')
    call ccnf_put_att(timencid,idnt,'missing_value',nf90_fill_float)
    lname = 'shorwave upward radiation at surface'
    call ccnf_def_var(timencid,'qup','float',1,jdim(1:1),idnt)
    call ccnf_put_att(timencid,idnt,'long_name',lname)
    call ccnf_put_att(timencid,idnt,'units','W/m2')
    call ccnf_put_att(timencid,idnt,'missing_value',nf90_fill_float)
    lname = 'sensible heat flux'
    call ccnf_def_var(timencid,'shf','float',1,jdim(1:1),idnt)
    call ccnf_put_att(timencid,idnt,'long_name',lname)
    call ccnf_put_att(timencid,idnt,'units','W/m2')
    call ccnf_put_att(timencid,idnt,'missing_value',nf90_fill_float)
    lname = 'latent heat flux (liq+sol)'
    call ccnf_def_var(timencid,'lhf','float',1,jdim(1:1),idnt)
    call ccnf_put_att(timencid,idnt,'long_name',lname)
    call ccnf_put_att(timencid,idnt,'units','W/m2')
    call ccnf_put_att(timencid,idnt,'missing_value',nf90_fill_float)
    lname = 'anthropogenic heat flux'
    call ccnf_def_var(timencid,'qf','float',1,jdim(1:1),idnt)
    call ccnf_put_att(timencid,idnt,'long_name',lname)
    call ccnf_put_att(timencid,idnt,'units','W/m2')
    call ccnf_put_att(timencid,idnt,'missing_value',nf90_fill_float)
    lname = 'urban storage heat flux'
    call ccnf_def_var(timencid,'g','float',1,jdim(1:1),idnt)
    call ccnf_put_att(timencid,idnt,'long_name',lname)
    call ccnf_put_att(timencid,idnt,'units','W/m2')
    call ccnf_put_att(timencid,idnt,'missing_value',nf90_fill_float)
    ! lname = 'Evaporation + sublimation flux'
    ! call ccnf_def_var(timencid,'evap','float',1,jdim(1:1),idnt)
    ! call ccnf_put_att(timencid,idnt,'long_name',lname)
    ! call ccnf_put_att(timencid,idnt,'units','mm/day')
    ! call ccnf_put_att(timencid,idnt,'missing_value',nf90_fill_float)
    lname = 'friction velocity'
    call ccnf_def_var(timencid,'ustar','float',1,jdim(1:1),idnt)
    call ccnf_put_att(timencid,idnt,'long_name',lname)
    call ccnf_put_att(timencid,idnt,'units','m/s')
    call ccnf_put_att(timencid,idnt,'missing_value',nf90_fill_float)
    ! lname = 'Precipitation (liq+sol) rate'
    ! call ccnf_def_var(timencid,'rain','float',1,jdim(1:1),idnt)
    ! call ccnf_put_att(timencid,idnt,'long_name',lname)
    ! call ccnf_put_att(timencid,idnt,'units','mm/day')
    ! call ccnf_put_att(timencid,idnt,'missing_value',nf90_fill_float)
    ! lname = 'Surface pressure'
    ! call ccnf_def_var(timencid,'psurf','float',1,jdim(1:1),idnt)
    ! call ccnf_put_att(timencid,idnt,'long_name',lname)
    ! call ccnf_put_att(timencid,idnt,'units','Pa')
    ! call ccnf_put_att(timencid,idnt,'missing_value',nf90_fill_float)
    lname = 'boundary layer height'
    call ccnf_def_var(timencid,'hpbl','float',1,jdim(1:1),idnt)
    call ccnf_put_att(timencid,idnt,'long_name',lname)
    call ccnf_put_att(timencid,idnt,'units','m')
    call ccnf_put_att(timencid,idnt,'missing_value',nf90_fill_float)
    lname = 'temperature skin layer'
    call ccnf_def_var(timencid,'tsk','float',1,jdim(1:1),idnt)
    call ccnf_put_att(timencid,idnt,'long_name',lname)
    call ccnf_put_att(timencid,idnt,'units','K')
    call ccnf_put_att(timencid,idnt,'missing_value',nf90_fill_float)
    ! lname = 'Radiative temperature if different from tsurf'
    ! call ccnf_def_var(timencid,'trad','float',1,jdim(1:1),idnt)
    ! call ccnf_put_att(timencid,idnt,'long_name',lname)
    ! call ccnf_put_att(timencid,idnt,'units','K')
    ! call ccnf_put_att(timencid,idnt,'missing_value',nf90_fill_float)
    lname = 'surface albedo'
    call ccnf_def_var(timencid,'alb','float',1,jdim(1:1),idnt)
    call ccnf_put_att(timencid,idnt,'long_name',lname)
    call ccnf_put_att(timencid,idnt,'units','0-1')
    call ccnf_put_att(timencid,idnt,'missing_value',nf90_fill_float)
    lname = 'roughness length momentum'
    call ccnf_def_var(timencid,'z0m','float',1,jdim(1:1),idnt)
    call ccnf_put_att(timencid,idnt,'long_name',lname)
    call ccnf_put_att(timencid,idnt,'units','m')
    call ccnf_put_att(timencid,idnt,'missing_value',nf90_fill_float)
    lname = 'roughness length heat'
    call ccnf_def_var(timencid,'z0h','float',1,jdim(1:1),idnt)
    call ccnf_put_att(timencid,idnt,'long_name',lname)
    call ccnf_put_att(timencid,idnt,'units','m')
    call ccnf_put_att(timencid,idnt,'missing_value',nf90_fill_float)
    lname = 'surface emissivity'
    call ccnf_def_var(timencid,'emis','float',1,jdim(1:1),idnt)
    call ccnf_put_att(timencid,idnt,'long_name',lname)
    call ccnf_put_att(timencid,idnt,'units','0-1')
    call ccnf_put_att(timencid,idnt,'missing_value',nf90_fill_float)
    lname = '2m temperature'
    call ccnf_def_var(timencid,'t2m','float',1,jdim(1:1),idnt)
    call ccnf_put_att(timencid,idnt,'long_name',lname)
    call ccnf_put_att(timencid,idnt,'units','K')
    call ccnf_put_att(timencid,idnt,'missing_value',nf90_fill_float)
    lname = '2m specific humidity'
    call ccnf_def_var(timencid,'q2m','float',1,jdim(1:1),idnt)
    call ccnf_put_att(timencid,idnt,'long_name',lname)
    call ccnf_put_att(timencid,idnt,'units','kg/kg')  
    call ccnf_put_att(timencid,idnt,'missing_value',nf90_fill_float)
    lname = '2m relative humidity'
    call ccnf_def_var(timencid,'rh2m','float',1,jdim(1:1),idnt)
    call ccnf_put_att(timencid,idnt,'long_name',lname)
    call ccnf_put_att(timencid,idnt,'units','0-100') 
    call ccnf_put_att(timencid,idnt,'missing_value',nf90_fill_float)
    lname = '10m u-component wind'
    call ccnf_def_var(timencid,'u10m','float',1,jdim(1:1),idnt)
    call ccnf_put_att(timencid,idnt,'long_name',lname)
    call ccnf_put_att(timencid,idnt,'units','m/s')
    call ccnf_put_att(timencid,idnt,'missing_value',nf90_fill_float)
    lname = '10m v-component wind'
    call ccnf_def_var(timencid,'v10m','float',1,jdim(1:1),idnt)
    call ccnf_put_att(timencid,idnt,'long_name',lname)
    call ccnf_put_att(timencid,idnt,'units','m/s')  
    call ccnf_put_att(timencid,idnt,'missing_value',nf90_fill_float)
    lname = '50m temperature'
    call ccnf_def_var(timencid,'t50m','float',1,jdim(1:1),idnt)
    call ccnf_put_att(timencid,idnt,'long_name',lname)
    call ccnf_put_att(timencid,idnt,'units','K')  
    call ccnf_put_att(timencid,idnt,'missing_value',nf90_fill_float)    
    lname = '50m specific humidity'
    call ccnf_def_var(timencid,'q50m','float',1,jdim(1:1),idnt)
    call ccnf_put_att(timencid,idnt,'long_name',lname)
    call ccnf_put_att(timencid,idnt,'units','kg/kg')  
    call ccnf_put_att(timencid,idnt,'missing_value',nf90_fill_float)
    lname = 'Relative humidity at 50 metre above the surface'
    call ccnf_def_var(timencid,'rh50m','float',1,jdim(1:1),idnt)
    call ccnf_put_att(timencid,idnt,'long_name',lname)
    call ccnf_put_att(timencid,idnt,'units','0-100')  
    call ccnf_put_att(timencid,idnt,'missing_value',nf90_fill_float)
    lname = '50m u-component wind'
    call ccnf_def_var(timencid,'u50m','float',1,jdim(1:1),idnt)
    call ccnf_put_att(timencid,idnt,'long_name',lname)
    call ccnf_put_att(timencid,idnt,'units','m/s')  
    call ccnf_put_att(timencid,idnt,'missing_value',nf90_fill_float)
    lname = '50m v-component wind'
    call ccnf_def_var(timencid,'v50m','float',1,jdim(1:1),idnt)
    call ccnf_put_att(timencid,idnt,'long_name',lname)
    call ccnf_put_att(timencid,idnt,'units','m/s')  
    call ccnf_put_att(timencid,idnt,'missing_value',nf90_fill_float)
    lname = 'cloud cover fraction'
    call ccnf_def_var(timencid,'cc','float',1,jdim(1:1),idnt)
    call ccnf_put_att(timencid,idnt,'long_name',lname)
    call ccnf_put_att(timencid,idnt,'units','0-1')  
    call ccnf_put_att(timencid,idnt,'missing_value',nf90_fill_float)

    lname = "sky view factor"
    call ccnf_def_var(timencid,'svf','float',1,jdim(1:1),idnt)
    call ccnf_put_att(timencid,idnt,'long_name',lname)
    call ccnf_put_att(timencid,idnt,'units','none')
    call ccnf_put_att(timencid,idnt,'missing_value',nf90_fill_float)
    lname = "aspect ratio"
    call ccnf_def_var(timencid,'asr','float',1,jdim(1:1),idnt)
    call ccnf_put_att(timencid,idnt,'long_name',lname)
    call ccnf_put_att(timencid,idnt,'units','none')
    call ccnf_put_att(timencid,idnt,'missing_value',nf90_fill_float)
    lname = "urban surface cover fraction"
    call ccnf_def_var(timencid,'urf','float',1,jdim(1:1),idnt)
    call ccnf_put_att(timencid,idnt,'long_name',lname)
    call ccnf_put_att(timencid,idnt,'units','none')
    call ccnf_put_att(timencid,idnt,'missing_value',nf90_fill_float)
    lname = "roof height"
    call ccnf_def_var(timencid,'rfh','float',1,jdim(1:1),idnt)
    call ccnf_put_att(timencid,idnt,'long_name',lname)
    call ccnf_put_att(timencid,idnt,'units','m')
    call ccnf_put_att(timencid,idnt,'missing_value',nf90_fill_float)
    lname = "standard deviation of roof heights"
    call ccnf_def_var(timencid,'sd_rfh','float',1,jdim(1:1),idnt)
    call ccnf_put_att(timencid,idnt,'long_name',lname)
    call ccnf_put_att(timencid,idnt,'units','m')
    call ccnf_put_att(timencid,idnt,'missing_value',nf90_fill_float)   
    lname = "emissivity wall"
    call ccnf_def_var(timencid,'emis_w','float',1,jdim(1:1),idnt)
    call ccnf_put_att(timencid,idnt,'long_name',lname)
    call ccnf_put_att(timencid,idnt,'units','none')
    call ccnf_put_att(timencid,idnt,'missing_value',nf90_fill_float)
    lname = "emissivity road"
    call ccnf_def_var(timencid,'emis_g','float',1,jdim(1:1),idnt)
    call ccnf_put_att(timencid,idnt,'long_name',lname)
    call ccnf_put_att(timencid,idnt,'units','none')
    call ccnf_put_att(timencid,idnt,'missing_value',nf90_fill_float)   
    lname = "emissivity roof"
    call ccnf_def_var(timencid,'emis_r','float',1,jdim(1:1),idnt)
    call ccnf_put_att(timencid,idnt,'long_name',lname)
    call ccnf_put_att(timencid,idnt,'units','none')
    call ccnf_put_att(timencid,idnt,'missing_value',nf90_fill_float)
    lname = "albedo wall"
    call ccnf_def_var(timencid,'alb_w','float',1,jdim(1:1),idnt)
    call ccnf_put_att(timencid,idnt,'long_name',lname)
    call ccnf_put_att(timencid,idnt,'units','none')
    call ccnf_put_att(timencid,idnt,'missing_value',nf90_fill_float)   
     lname = "albedo road"
    call ccnf_def_var(timencid,'alb_g','float',1,jdim(1:1),idnt)
    call ccnf_put_att(timencid,idnt,'long_name',lname)
    call ccnf_put_att(timencid,idnt,'units','none')
    call ccnf_put_att(timencid,idnt,'missing_value',nf90_fill_float)
    lname = "albedo roof"
    call ccnf_def_var(timencid,'alb_r','float',1,jdim(1:1),idnt)
    call ccnf_put_att(timencid,idnt,'long_name',lname)
    call ccnf_put_att(timencid,idnt,'units','none')
    call ccnf_put_att(timencid,idnt,'missing_value',nf90_fill_float)
    jdim(1) = udim
    jdim(2) = tdim_time
    lname = "thermal conductivity wall"
    call ccnf_def_var(timencid,'thc_w','float',2,jdim(1:2),idnt)
    call ccnf_put_att(timencid,idnt,'long_name',lname)
    call ccnf_put_att(timencid,idnt,'units','J/m/s/K')
    call ccnf_put_att(timencid,idnt,'missing_value',nf90_fill_float)
    lname = "thermal conductivity road"
    call ccnf_def_var(timencid,'thc_g','float',2,jdim(1:2),idnt)
    call ccnf_put_att(timencid,idnt,'long_name',lname)
    call ccnf_put_att(timencid,idnt,'units','J/m/s/K')
    call ccnf_put_att(timencid,idnt,'missing_value',nf90_fill_float)
    lname = "thermal conductivity roof"
    call ccnf_def_var(timencid,'thc_r','float',2,jdim(1:2),idnt)
    call ccnf_put_att(timencid,idnt,'long_name',lname)
    call ccnf_put_att(timencid,idnt,'units','J/m/s/K')
    call ccnf_put_att(timencid,idnt,'missing_value',nf90_fill_float)
    lname = "heat capacity wall"
    call ccnf_def_var(timencid,'hc_w','float',2,jdim(1:2),idnt)
    call ccnf_put_att(timencid,idnt,'long_name',lname)
    call ccnf_put_att(timencid,idnt,'units','J/m3/K')
    call ccnf_put_att(timencid,idnt,'missing_value',nf90_fill_float)
    lname = "heat capacity road"
    call ccnf_def_var(timencid,'hc_g','float',2,jdim(1:2),idnt)
    call ccnf_put_att(timencid,idnt,'long_name',lname)
    call ccnf_put_att(timencid,idnt,'units','J/m3/K')
    call ccnf_put_att(timencid,idnt,'missing_value',nf90_fill_float)
    lname = "heat capacity roof"
    call ccnf_def_var(timencid,'hc_r','float',2,jdim(1:2),idnt)
    call ccnf_put_att(timencid,idnt,'long_name',lname)
    call ccnf_put_att(timencid,idnt,'units','J/m3/K')
    call ccnf_put_att(timencid,idnt,'missing_value',nf90_fill_float)
    
    call ccnf_enddef(timencid)
    
    uu(1,:) = (/ 1., 2., 3., 4. /)
    spos(1) = 1
    npos(1) = 4
    call ccnf_put_vara(timencid,'urban_layer',spos(1),npos(1),uu(1,:))
    
    
    write(6,*) "Creating profile output file"  
    call ccnf_create(profileoutput,profilencid)
    ! Turn off the data filling
    call ccnf_nofill(profilencid)
    
    ! attributes
    call ccnf_put_attg(profilencid,'model','CCAM+UCLEM')
    call ccnf_put_attg(profilencid,'contact','Mathew Lipson <m.lipson@unsw.edu.au>')
    call ccnf_put_attg(profilencid,'scmtype','Stretched grid climate model')
    call ccnf_put_attg(profilencid,'ucmcomplexity','complex')
    call ccnf_put_attg(profilencid,'ucmtiled','no')
    call ccnf_put_attg(profilencid,'timestep',dt)
    
    if ( nvmix==7 ) then
      call ccnf_put_attg(profilencid,'turbulencescheme','Jing')  
      call ccnf_put_attg(profilencid,'eddydiff','?')
      call ccnf_put_attg(profilencid,'massflux','?')
      call ccnf_put_attg(profilencid,'lengthscale','?')
      call ccnf_put_attg(profilencid,'kprofile','?')
    else if ( nvmix==6 ) then
      call ccnf_put_attg(profilencid,'turbulencescheme','EDMF')
      call ccnf_put_attg(profilencid,'eddydiff','TKE-EPS')
      call ccnf_put_attg(profilencid,'massflux','updraft')
    else
      call ccnf_put_attg(profilencid,'turbulencescheme','K-profile')  
      call ccnf_put_attg(profilencid,'eddydiff','?')
      call ccnf_put_attg(profilencid,'massflux','diagnosed')
      call ccnf_put_attg(profilencid,'lengthscale','?')
      call ccnf_put_attg(profilencid,'kprofile','Richardson')
    end if
    
    ! call ccnf_put_attg(profilencid,'soillayers','6')
    ! call ccnf_put_attg(profilencid,'snowlayers','1-3')
    ! call ccnf_put_attg(profilencid,'surfaceprog','temp,moisture,ice,canopywater,snowtemp,snowmass,snowage')
    
    call ccnf_def_dim(profilencid,'levf',kl,zfdim)
    call ccnf_def_dim(profilencid,'levh',kl+1,zhdim)
    call ccnf_def_dim(profilencid,'levs',3,zsdim) ! snow levels
    call ccnf_def_dimu(profilencid,'time',tdim_prof)

    jdim(1) = zfdim
    call ccnf_def_var(profilencid,'levf','float',1,jdim(1:1),idnt)
    call ccnf_put_att(profilencid,idnt,'units','Pa')

    jdim(1) = zhdim
    call ccnf_def_var(profilencid,'levh','float',1,jdim(1:1),idnt)
    call ccnf_put_att(profilencid,idnt,'units','Pa')

    jdim(1) = zsdim
    call ccnf_def_var(profilencid,'levs','float',1,jdim(1:1),idnt)
    call ccnf_put_att(profilencid,idnt,'units','index')
    
    jdim(1) = tdim_prof
    call ccnf_def_var(profilencid,'time','float',1,jdim(1:1),idnt)
    call ccnf_put_att(profilencid,idnt,'point_spacing','even')
    icy = kdate/10000
    icm = max(1, min(12, (kdate-icy*10000)/100))
    icd = max(1, min(31, (kdate-icy*10000-icm*100)))
    ich = ktime/100
    icmi = (ktime-ich*100)
    ics = 0
    write(grdtim,'("seconds since ",i4.4,"-",i2.2,"-",i2.2," ",2(i2.2,":"),i2.2)') icy,icm,icd,ich,icmi,ics
    call ccnf_put_att(profilencid,idnt,'units',grdtim)    

    
    jdim(1) = zfdim
    jdim(2) = tdim_prof
    lname = 'height of full level'
    call ccnf_def_var(profilencid,'zf','float',2,jdim(1:2),idnt)
    call ccnf_put_att(profilencid,idnt,'long_name',lname)
    call ccnf_put_att(profilencid,idnt,'units','m')  
    call ccnf_put_att(profilencid,idnt,'missing_value',nf90_fill_float)
    lname = 'pressure at full level'
    call ccnf_def_var(profilencid,'pf','float',2,jdim(1:2),idnt)
    call ccnf_put_att(profilencid,idnt,'long_name',lname)
    call ccnf_put_att(profilencid,idnt,'units','Pa')  
    call ccnf_put_att(profilencid,idnt,'missing_value',nf90_fill_float)
    lname = 'temperature'
    call ccnf_def_var(profilencid,'t','float',2,jdim(1:2),idnt)
    call ccnf_put_att(profilencid,idnt,'long_name',lname)
    call ccnf_put_att(profilencid,idnt,'units','K')  
    call ccnf_put_att(profilencid,idnt,'missing_value',nf90_fill_float)
    ! lname = 'Temperature'
    ! call ccnf_def_var(profilencid,'temp','float',2,jdim(1:2),idnt)
    ! call ccnf_put_att(profilencid,idnt,'long_name',lname)
    ! call ccnf_put_att(profilencid,idnt,'units','K')  
    ! call ccnf_put_att(profilencid,idnt,'missing_value',nf90_fill_float)
    lname = 'potential temperature'
    call ccnf_def_var(profilencid,'th','float',2,jdim(1:2),idnt)
    call ccnf_put_att(profilencid,idnt,'long_name',lname)
    call ccnf_put_att(profilencid,idnt,'units','K')  
    call ccnf_put_att(profilencid,idnt,'missing_value',nf90_fill_float)
    lname = 'specific humidity'
    call ccnf_def_var(profilencid,'q','float',2,jdim(1:2),idnt)
    call ccnf_put_att(profilencid,idnt,'long_name',lname)
    call ccnf_put_att(profilencid,idnt,'units','kg/kg')  
    call ccnf_put_att(profilencid,idnt,'missing_value',nf90_fill_float)
    ! lname = 'Cloud water and ice'
    ! call ccnf_def_var(profilencid,'qc','float',2,jdim(1:2),idnt)
    ! call ccnf_put_att(profilencid,idnt,'long_name',lname)
    ! call ccnf_put_att(profilencid,idnt,'units','kg/kg')  
    ! call ccnf_put_att(profilencid,idnt,'missing_value',nf90_fill_float)
    lname = 'zonal component wind'
    call ccnf_def_var(profilencid,'u','float',2,jdim(1:2),idnt)
    call ccnf_put_att(profilencid,idnt,'long_name',lname)
    call ccnf_put_att(profilencid,idnt,'units','m/s')  
    call ccnf_put_att(profilencid,idnt,'missing_value',nf90_fill_float)
    lname = 'meridional component wind'
    call ccnf_def_var(profilencid,'v','float',2,jdim(1:2),idnt)
    call ccnf_put_att(profilencid,idnt,'long_name',lname)
    call ccnf_put_att(profilencid,idnt,'units','m/s')  
    call ccnf_put_att(profilencid,idnt,'missing_value',nf90_fill_float)

    jdim(1) = zfdim
    jdim(2) = tdim_prof
    lname = 'u-component geostrophic wind'
    call ccnf_def_var(profilencid,'ugeo','float',2,jdim(1:2),idnt)
    call ccnf_put_att(profilencid,idnt,'long_name',lname)
    call ccnf_put_att(profilencid,idnt,'units','m/s')  
    call ccnf_put_att(profilencid,idnt,'missing_value',nf90_fill_float)  
    lname = 'v-component geostrophic wind'
    call ccnf_def_var(profilencid,'vgeo','float',2,jdim(1:2),idnt)
    call ccnf_put_att(profilencid,idnt,'long_name',lname)
    call ccnf_put_att(profilencid,idnt,'units','m/s')  
    call ccnf_put_att(profilencid,idnt,'missing_value',nf90_fill_float)  
    lname = 'u-component momentum advection'
    call ccnf_def_var(profilencid,'dudt_ls','float',2,jdim(1:2),idnt)
    call ccnf_put_att(profilencid,idnt,'long_name',lname)
    call ccnf_put_att(profilencid,idnt,'units','m/s2')  
    call ccnf_put_att(profilencid,idnt,'missing_value',nf90_fill_float)  
    lname = 'v-component momentum advection'
    call ccnf_def_var(profilencid,'dvdt_ls','float',2,jdim(1:2),idnt)
    call ccnf_put_att(profilencid,idnt,'long_name',lname)
    call ccnf_put_att(profilencid,idnt,'units','m/s2')  
    call ccnf_put_att(profilencid,idnt,'missing_value',nf90_fill_float)  
    lname = 'temperature advection'
    call ccnf_def_var(profilencid,'dtdt_ls','float',2,jdim(1:2),idnt)
    call ccnf_put_att(profilencid,idnt,'long_name',lname)
    call ccnf_put_att(profilencid,idnt,'units','K/s')  
    call ccnf_put_att(profilencid,idnt,'missing_value',nf90_fill_float)  
    lname = 'moisture advection'
    call ccnf_def_var(profilencid,'dqdt_ls','float',2,jdim(1:2),idnt)
    call ccnf_put_att(profilencid,idnt,'long_name',lname)
    call ccnf_put_att(profilencid,idnt,'units','kg/kg/s')  
    call ccnf_put_att(profilencid,idnt,'missing_value',nf90_fill_float)  
    lname = 'vertical velocity'
    call ccnf_def_var(profilencid,'w','float',2,jdim(1:2),idnt)
    call ccnf_put_att(profilencid,idnt,'long_name',lname)
    call ccnf_put_att(profilencid,idnt,'units','m/s')  
    call ccnf_put_att(profilencid,idnt,'missing_value',nf90_fill_float)  
     
    jdim(1) = zhdim
    jdim(2) = tdim_prof
    lname = 'height of half level'
    call ccnf_def_var(profilencid,'zh','float',2,jdim(1:2),idnt)
    call ccnf_put_att(profilencid,idnt,'long_name',lname)
    call ccnf_put_att(profilencid,idnt,'units','m')  
    call ccnf_put_att(profilencid,idnt,'missing_value',nf90_fill_float)
    lname = 'pressure at half level'
    call ccnf_def_var(profilencid,'ph','float',2,jdim(1:2),idnt)
    call ccnf_put_att(profilencid,idnt,'long_name',lname)
    call ccnf_put_att(profilencid,idnt,'units','Pa')  
    call ccnf_put_att(profilencid,idnt,'missing_value',nf90_fill_float)
    lname = 'vertical temperature flux'
    call ccnf_def_var(profilencid,'wt','float',2,jdim(1:2),idnt)
    call ccnf_put_att(profilencid,idnt,'long_name',lname)
    call ccnf_put_att(profilencid,idnt,'units','K m/s')  
    call ccnf_put_att(profilencid,idnt,'missing_value',nf90_fill_float)
    lname = 'vertical moisture flux'
    call ccnf_def_var(profilencid,'wq','float',2,jdim(1:2),idnt)
    call ccnf_put_att(profilencid,idnt,'long_name',lname)
    call ccnf_put_att(profilencid,idnt,'units','kg/kg m/s')  
    call ccnf_put_att(profilencid,idnt,'missing_value',nf90_fill_float)
    lname = 'vertical flux u-component momentum'
    call ccnf_def_var(profilencid,'uw','float',2,jdim(1:2),idnt)
    call ccnf_put_att(profilencid,idnt,'long_name',lname)
    call ccnf_put_att(profilencid,idnt,'units','m2/s2')  
    call ccnf_put_att(profilencid,idnt,'missing_value',nf90_fill_float)
    lname = 'vertical flux v-component momentum'
    call ccnf_def_var(profilencid,'vw','float',2,jdim(1:2),idnt)
    call ccnf_put_att(profilencid,idnt,'long_name',lname)
    call ccnf_put_att(profilencid,idnt,'units','m/s')  
    call ccnf_put_att(profilencid,idnt,'missing_value',nf90_fill_float)
    ! lname = 'u-variance'
    ! call ccnf_def_var(profilencid,'uu','float',2,jdim(1:2),idnt)
    ! call ccnf_put_att(profilencid,idnt,'long_name',lname)
    ! call ccnf_put_att(profilencid,idnt,'units','m2/s2')  
    ! call ccnf_put_att(profilencid,idnt,'missing_value',nf90_fill_float)  
    ! lname = 'v-variance'
    ! call ccnf_def_var(profilencid,'vv','float',2,jdim(1:2),idnt)
    ! call ccnf_put_att(profilencid,idnt,'long_name',lname)
    ! call ccnf_put_att(profilencid,idnt,'units','m2/s2')  
    ! call ccnf_put_att(profilencid,idnt,'missing_value',nf90_fill_float)
    ! lname = 'w-variance'
    ! call ccnf_def_var(profilencid,'ww','float',2,jdim(1:2),idnt)
    ! call ccnf_put_att(profilencid,idnt,'long_name',lname)
    ! call ccnf_put_att(profilencid,idnt,'units','m2/s2')  
    ! call ccnf_put_att(profilencid,idnt,'missing_value',nf90_fill_float)  
    ! lname = 'Potential temperature variance'
    ! call ccnf_def_var(profilencid,'th^2','float',2,jdim(1:2),idnt)
    ! call ccnf_put_att(profilencid,idnt,'long_name',lname)
    ! call ccnf_put_att(profilencid,idnt,'units','K^2')  
    ! call ccnf_put_att(profilencid,idnt,'missing_value',nf90_fill_float)  
    lname = 'eddy diffusivity momentum'
    call ccnf_def_var(profilencid,'Km','float',2,jdim(1:2),idnt)
    call ccnf_put_att(profilencid,idnt,'long_name',lname)
    call ccnf_put_att(profilencid,idnt,'units','m2/s')  
    call ccnf_put_att(profilencid,idnt,'missing_value',nf90_fill_float)
    lname = 'eddy diffusivity heat'
    call ccnf_def_var(profilencid,'Kh','float',2,jdim(1:2),idnt)
    call ccnf_put_att(profilencid,idnt,'long_name',lname)
    call ccnf_put_att(profilencid,idnt,'units','m2/s')  
    call ccnf_put_att(profilencid,idnt,'missing_value',nf90_fill_float)
    ! lname = 'Massflux'
    ! call ccnf_def_var(profilencid,'mf','float',2,jdim(1:2),idnt)
    ! call ccnf_put_att(profilencid,idnt,'long_name',lname)
    ! call ccnf_put_att(profilencid,idnt,'units','kg/m2/s')  
    ! call ccnf_put_att(profilencid,idnt,'missing_value',nf90_fill_float)  
    jdim(1) = zfdim
    ! lname = 'Temperature tendency from radiation'
    ! call ccnf_def_var(profilencid,'dT_dt_rad','float',2,jdim(1:2),idnt)
    ! call ccnf_put_att(profilencid,idnt,'long_name',lname)
    ! call ccnf_put_att(profilencid,idnt,'units','K/s')  
    ! call ccnf_put_att(profilencid,idnt,'missing_value',nf90_fill_float)  
    ! lname = 'Temperature tendency from short wave radiation'
    ! call ccnf_def_var(profilencid,'dT_dt_swrad','float',2,jdim(1:2),idnt)
    ! call ccnf_put_att(profilencid,idnt,'long_name',lname)
    ! call ccnf_put_att(profilencid,idnt,'units','K/s')  
    ! call ccnf_put_att(profilencid,idnt,'missing_value',nf90_fill_float)  
    ! lname = 'Temperature tendency from long wave radiation'
    ! call ccnf_def_var(profilencid,'dT_dt_lwrad','float',2,jdim(1:2),idnt)
    ! call ccnf_put_att(profilencid,idnt,'long_name',lname)
    ! call ccnf_put_att(profilencid,idnt,'units','K/s')  
    ! call ccnf_put_att(profilencid,idnt,'missing_value',nf90_fill_float)  
    lname = 'turbulent kinetic energy'
    call ccnf_def_var(profilencid,'TKE','float',2,jdim(1:2),idnt)
    call ccnf_put_att(profilencid,idnt,'long_name',lname)
    call ccnf_put_att(profilencid,idnt,'units','m2/s2')  
    call ccnf_put_att(profilencid,idnt,'missing_value',nf90_fill_float)  
    lname = 'shear production'
    call ccnf_def_var(profilencid,'shear','float',2,jdim(1:2),idnt)
    call ccnf_put_att(profilencid,idnt,'long_name',lname)
    call ccnf_put_att(profilencid,idnt,'units','m2/s3')  
    call ccnf_put_att(profilencid,idnt,'missing_value',nf90_fill_float)  
    lname = 'buoyancy production'
    call ccnf_def_var(profilencid,'buoy','float',2,jdim(1:2),idnt)
    call ccnf_put_att(profilencid,idnt,'long_name',lname)
    call ccnf_put_att(profilencid,idnt,'units','m2/s3')  
    call ccnf_put_att(profilencid,idnt,'missing_value',nf90_fill_float)  
    lname = 'total transport'
    call ccnf_def_var(profilencid,'trans','float',2,jdim(1:2),idnt)
    call ccnf_put_att(profilencid,idnt,'long_name',lname)
    call ccnf_put_att(profilencid,idnt,'units','m2/s3')  
    call ccnf_put_att(profilencid,idnt,'missing_value',nf90_fill_float)  
    lname = 'dissipation'
    call ccnf_def_var(profilencid,'dissi','float',2,jdim(1:2),idnt)
    call ccnf_put_att(profilencid,idnt,'long_name',lname)
    call ccnf_put_att(profilencid,idnt,'units','m2/s3')  
    call ccnf_put_att(profilencid,idnt,'missing_value',nf90_fill_float)  
    jdim(1) = zhdim
    ! lname = 'Total turbulent energy'
    ! call ccnf_def_var(profilencid,'TTE','float',2,jdim(1:2),idnt)
    ! call ccnf_put_att(profilencid,idnt,'long_name',lname)
    ! call ccnf_put_att(profilencid,idnt,'units','m2/s2')  
    ! call ccnf_put_att(profilencid,idnt,'missing_value',nf90_fill_float)  

    jdim(1) = zsdim
    jdim(2) = tdim_prof
    ! lname = 'Snow temperature'
    ! call ccnf_def_var(profilencid,'tsn','float',2,jdim(1:2),idnt)
    ! call ccnf_put_att(profilencid,idnt,'long_name',lname)
    ! call ccnf_put_att(profilencid,idnt,'units','K')  
    ! call ccnf_put_att(profilencid,idnt,'missing_value',nf90_fill_float)  
      
    call ccnf_enddef(profilencid)
    
    bb(1,:) = sig(:)*ps(1)
    spos(1) = 1
    npos(1) = kl
    call ccnf_put_vara(profilencid,'levf',spos(1),npos(1),bb(1,:))

    cc(1,1:kl) = sigmh(1:kl)*ps(1)
    cc(1,kl+1) = 0.
    spos(1) = 1
    npos(1) = kl+1
    call ccnf_put_vara(profilencid,'levh',spos(1),npos(1),cc(1,:))
    
    dd(1,1) = sdepth(1,1)*0.5
    dd(1,2) = sdepth(1,1)+sdepth(1,2)*0.5
    dd(1,3) = sdepth(1,1)+sdepth(1,2)+sdepth(1,3)*0.5
    spos(1) = 1
    npos(1) = 3
    call ccnf_put_vara(profilencid,'levs',spos(1),npos(1),dd(1,:))
    
    iarch = 1  
  end if

  zf(1) = bet(1)*t(1,1)
  do k=2,kl
    zf(k) = zf(k-1) + bet(k)*t(1,k) + betm(k)*t(1,k-1)
  end do
  zf(:) = zf(:)/grav

  ! using t.  Maybe should use tv?
  zh(1) = 0.
  do k=2,kl
    zh(k) = zh(k-1) + t(1,k)*(-rdry/grav)*dsig(k)/sig(k)
  end do
  zh(kl+1) = 2.*zh(kl) - zh(kl-1)

  pf(1:kl) = sig(1:kl)*ps(1)

  qs(1:kl) = qsat(pf,t(1,:))
  rh(:) = qg(1,:)/qs(:)

  do k=1,kl
    wtflux(1,k) = wth_flux(1,k)*sigmh(k)**(rdry/cp)
  end do
  
  call ccnf_put_vara(timencid,'time',iarch,real(ktau)*dt)

  if ( ktau==0 ) then
    aa(:) = nf90_fill_float    
    call ccnf_put_vara(timencid,'ldw',iarch,aa(1))
    call ccnf_put_vara(timencid,'lup',iarch,aa(1))
    call ccnf_put_vara(timencid,'qdw',iarch,aa(1))
    call ccnf_put_vara(timencid,'qup',iarch,aa(1))
    call ccnf_put_vara(timencid,'shf',iarch,aa(1))
    call ccnf_put_vara(timencid,'lhf',iarch,aa(1))
    call ccnf_put_vara(timencid,'qf',iarch,aa(1))
    call ccnf_put_vara(timencid,'g',iarch,aa(1))
    ! call ccnf_put_vara(timencid,'evap',iarch,aa(1))
    call ccnf_put_vara(timencid,'ustar',iarch,aa(1))
    ! call ccnf_put_vara(timencid,'rain',iarch,aa(1))
    ! call ccnf_put_vara(timencid,'psurf',iarch,aa(1))
    call ccnf_put_vara(timencid,'hpbl',iarch,aa(1))
    call ccnf_put_vara(timencid,'tsk',iarch,aa(1))
    ! call ccnf_put_vara(timencid,'trad',iarch,aa(1))
    call ccnf_put_vara(timencid,'alb',iarch,aa(1))
    call ccnf_put_vara(timencid,'z0m',iarch,aa(1))
    call ccnf_put_vara(timencid,'z0h',iarch,aa(1))
    call ccnf_put_vara(timencid,'emis',iarch,aa(1))
    call ccnf_put_vara(timencid,'t2m',iarch,aa(1))
    call ccnf_put_vara(timencid,'q2m',iarch,aa(1))
    call ccnf_put_vara(timencid,'rh2m',iarch,aa(1))
    call ccnf_put_vara(timencid,'u10m',iarch,aa(1))
    call ccnf_put_vara(timencid,'v10m',iarch,aa(1))
    call ccnf_put_vara(timencid,'t50m',iarch,aa(1))
    call ccnf_put_vara(timencid,'q50m',iarch,aa(1))
    call ccnf_put_vara(timencid,'rh50m',iarch,aa(1))
    call ccnf_put_vara(timencid,'u50m',iarch,aa(1))
    call ccnf_put_vara(timencid,'v50m',iarch,aa(1))
    call ccnf_put_vara(timencid,'cc',iarch,aa(1))

    call ccnf_put_vara(timencid,'svf',iarch,aa(1))
    call ccnf_put_vara(timencid,'asr',iarch,aa(1))
    call ccnf_put_vara(timencid,'urf',iarch,aa(1))
    call ccnf_put_vara(timencid,'rfh',iarch,aa(1))
    call ccnf_put_vara(timencid,'sd_rfh',iarch,aa(1))
    call ccnf_put_vara(timencid,'emis_w',iarch,aa(1))
    call ccnf_put_vara(timencid,'emis_g',iarch,aa(1))
    call ccnf_put_vara(timencid,'emis_r',iarch,aa(1))
    call ccnf_put_vara(timencid,'alb_w',iarch,aa(1))
    call ccnf_put_vara(timencid,'alb_g',iarch,aa(1))
    call ccnf_put_vara(timencid,'alb_r',iarch,aa(1))
     
    spos(1) = 1
    spos(2) = iarch
    npos(1) = 4
    npos(2) = 1
    uu = nf90_fill_float
    call ccnf_put_vara(timencid,'thc_w',spos(1:2),npos(1:2),uu)
    call ccnf_put_vara(timencid,'thc_g',spos(1:2),npos(1:2),uu)
    call ccnf_put_vara(timencid,'thc_r',spos(1:2),npos(1:2),uu)
    call ccnf_put_vara(timencid,'hc_w',spos(1:2),npos(1:2),uu)
    call ccnf_put_vara(timencid,'hc_g',spos(1:2),npos(1:2),uu)
    call ccnf_put_vara(timencid,'hc_r',spos(1:2),npos(1:2),uu)
    
  else
      
    call ccnf_put_vara(timencid,'ldw',iarch,rgdn_ave(1))
    aa(:) = rgdn_ave(:) + rgn_ave(:)
    call ccnf_put_vara(timencid,'lup',iarch,aa(1))
    call ccnf_put_vara(timencid,'qdw',iarch,sgdn_ave(1))
    aa(:) = sgdn_ave(:) - sgn_ave(:)
    call ccnf_put_vara(timencid,'qup',iarch,aa(1))
    call ccnf_put_vara(timencid,'shf',iarch,fg(1))
    call ccnf_put_vara(timencid,'lhf',iarch,eg(1))
    call ccnf_put_vara(timencid,'qf',iarch,anthropogenic_flux(1))
    call ccnf_put_vara(timencid,'g',iarch,urban_storage_flux(1))
    !aa(:) = eg(:)*86400./hls ! this assumes sublimation
    ! call ccnf_put_vara(timencid,'evap',iarch,aa(1))
    call ccnf_put_vara(timencid,'ustar',iarch,ustar(1))
    ! call ccnf_put_vara(timencid,'rain',iarch,precip(1))
    ! call ccnf_put_vara(timencid,'psurf',iarch,ps(1))
    call ccnf_put_vara(timencid,'hpbl',iarch,pblh(1))
    call ccnf_put_vara(timencid,'tsk',iarch,tss(1))
    aa(:) = ((rgdn_ave(:)+rgn_ave(:))/(0.98*stefbo))**(0.25)
    ! call ccnf_put_vara(timencid,'trad',iarch,aa(1))
    aa(:) = swrsave*albvisnir(:,1) + (1.-swrsave)*albvisnir(:,2)
    call ccnf_put_vara(timencid,'alb',iarch,aa(1))
    call ccnf_put_vara(timencid,'z0m',iarch,zo(1))
    call ccnf_put_vara(timencid,'z0h',iarch,zoh(1))
    aa(1) = 1.*(1.-sigmu(1)) + urban_emiss(1)*sigmu(1)
    call ccnf_put_vara(timencid,'emis',iarch,aa(1))
    call ccnf_put_vara(timencid,'t2m',iarch,tscrn(1))
    aa(:) = qgscrn(:)/(1.+qgscrn(:))
    call ccnf_put_vara(timencid,'q2m',iarch,aa(1))
    call ccnf_put_vara(timencid,'rh2m',iarch,rhscrn(1))
    tmp(:)=u(1,:)
    call vout(tmp,aa(1),zf,10.,kl)
    call ccnf_put_vara(timencid,'u10m',iarch,aa(1))
    tmp(:)=v(1,:)
    call vout(tmp,aa(1),zf,10.,kl)
    call ccnf_put_vara(timencid,'v10m',iarch,aa(1))
    tmp=t(1,:)
    call vout(tmp,aa(1),zf,50.,kl)
    call ccnf_put_vara(timencid,'t50m',iarch,aa(1))
    tmp=qg(1,:)
    call vout(tmp,aa(1),zf,50.,kl)
    aa(1) = aa(1)/(1.+aa(1))
    call ccnf_put_vara(timencid,'q50m',iarch,aa(1))
    call vout(rh(:),aa(1),zf,50.,kl)
    call ccnf_put_vara(timencid,'rh50m',iarch,aa(1))
    tmp(:)=u(1,:)
    call vout(tmp,aa(1),zf,50.,kl)
    call ccnf_put_vara(timencid,'u50m',iarch,aa(1))
    tmp(:)=v(1,:)
    call vout(tmp,aa(1),zf,50.,kl)
    call ccnf_put_vara(timencid,'v50m',iarch,aa(1))
    aa(1) = cld_ave(1)
    call ccnf_put_vara(timencid,'cc',iarch,aa(1))

    call atebdeftype_export(aa(1:1),'hwratio',0)
    aa(1) = sqrt(aa(1)*aa(1)+1.) - aa(1) ! skyview factor for road
    call ccnf_put_vara(timencid,'svf',iarch,aa(1))
    call atebdeftype_export(aa(1:1),'hwratio',0)
    call ccnf_put_vara(timencid,'asr',iarch,aa(1))
    call atebdeftype_export(aa(1:1),'sigvegc',0)
    aa(1) = 1. - aa(1) ! urban cover after removing canyon vegetation
    call ccnf_put_vara(timencid,'urf',iarch,aa(1))
    call atebdeftype_export(aa(1:1),'bldheight',0)
    call ccnf_put_vara(timencid,'rfh',iarch,aa(1))
    aa(1) = 0.
    call ccnf_put_vara(timencid,'sd_rfh',iarch,aa(1))
    call atebdeftype_export(aa(1:1),'wallemiss',0)
    call ccnf_put_vara(timencid,'emis_w',iarch,aa(1))
    call atebdeftype_export(aa(1:1),'roademiss',0)
    call ccnf_put_vara(timencid,'emis_g',iarch,aa(1))
    call atebdeftype_export(aa(1:1),'roofemiss',0)
    call ccnf_put_vara(timencid,'emis_r',iarch,aa(1))
    call atebdeftype_export(aa(1:1),'wallalpha',0)
    call ccnf_put_vara(timencid,'alb_w',iarch,aa(1))
    call atebdeftype_export(aa(1:1),'roadalpha',0)
    call ccnf_put_vara(timencid,'alb_g',iarch,aa(1))
    call atebdeftype_export(aa(1:1),'roofalpha',0)
    call ccnf_put_vara(timencid,'alb_r',iarch,aa(1))
     
    spos(1) = 1
    spos(2) = iarch
    npos(1) = 4
    npos(2) = 1
    do k = 1,4
      write(lname,'("wallcond",(I1.1))') k
      call atebdeftype_export(uu(1:1,k),lname,0)  
    end do
    call ccnf_put_vara(timencid,'thc_w',spos(1:2),npos(1:2),uu)
    do k = 1,4
      write(lname,'("roadcond",(I1.1))') k
      call atebdeftype_export(uu(1:1,k),lname,0)  
    end do
    call ccnf_put_vara(timencid,'thc_g',spos(1:2),npos(1:2),uu)
    do k = 1,4
      write(lname,'("roofcond",(I1.1))') k
      call atebdeftype_export(uu(1:1,k),lname,0)  
    end do
    call ccnf_put_vara(timencid,'thc_r',spos(1:2),npos(1:2),uu)
    do k = 1,4
      write(lname,'("wallcp",(I1.1))') k
      call atebdeftype_export(uu(1:1,k),lname,0)  
    end do
    call ccnf_put_vara(timencid,'hc_w',spos(1:2),npos(1:2),uu)
    do k = 1,4
      write(lname,'("roadcp",(I1.1))') k
      call atebdeftype_export(uu(1:1,k),lname,0)  
    end do
    call ccnf_put_vara(timencid,'hc_g',spos(1:2),npos(1:2),uu)
    do k = 1,4
      write(lname,'("roofcp",(I1.1))') k
      call atebdeftype_export(uu(1:1,k),lname,0)  
    end do
    call ccnf_put_vara(timencid,'hc_r',spos(1:2),npos(1:2),uu)
  end if

  call ccnf_put_vara(profilencid,'time',iarch,real(ktau)*dt)

  spos(1:2) = (/ 1, iarch /)
  npos(1:2) = (/ kl, 1 /)
  bb(1,:) = zf(:)
  call ccnf_put_vara(profilencid,'zf',spos(1:2),npos(1:2),bb)
  bb(1,:) = pf(:)
  call ccnf_put_vara(profilencid,'pf',spos(1:2),npos(1:2),bb)
  call ccnf_put_vara(profilencid,'t',spos(1:2),npos(1:2),t)
  ! call ccnf_put_vara(profilencid,'temp',spos(1:2),npos(1:2),t)
  bb(1,:) = t(1,:)*(ps(1)*sig(:)/1.e5)**(-rdry/cp)
  call ccnf_put_vara(profilencid,'th',spos(1:2),npos(1:2),bb)
  bb(:,:) = qg(:,:)/(1.+qg(:,:))
  call ccnf_put_vara(profilencid,'q',spos(1:2),npos(1:2),bb)
  !bb(:,:) = qlg(:,:) + qfg(:,:)
  !bb(:,:) = bb(:,:)/(1.+bb(:,:))
  ! call ccnf_put_vara(profilencid,'qc',spos(1:2),npos(1:2),bb)
  call ccnf_put_vara(profilencid,'u',spos(1:2),npos(1:2),u)
  call ccnf_put_vara(profilencid,'v',spos(1:2),npos(1:2),v)

  spos(1) = 1
  spos(2) = iarch
  npos(1) = kl
  npos(2) = 1
  bb(1,:) = ug(:)
  call ccnf_put_vara(profilencid,'ugeo',spos(1:2),npos(1:2),bb)
  bb(1,:) = vg(:)
  call ccnf_put_vara(profilencid,'vgeo',spos(1:2),npos(1:2),bb)
  bb(1,:) = uadv(:)
  call ccnf_put_vara(profilencid,'dudt_ls',spos(1:2),npos(1:2),bb)
  bb(1,:) = vadv(:)
  call ccnf_put_vara(profilencid,'dvdt_ls',spos(1:2),npos(1:2),bb)
  bb(1,:) = t_tend(:)
  call ccnf_put_vara(profilencid,'dtdt_ls',spos(1:2),npos(1:2),bb)
  bb(1,:) = q_tend(:)
  call ccnf_put_vara(profilencid,'dqdt_ls',spos(1:2),npos(1:2),bb)
  bb(:,:) = nf90_fill_float
  call ccnf_put_vara(profilencid,'w',spos(1:2),npos(1:2),bb)

  spos(1) = 1
  spos(2) = iarch
  npos(1) = kl+1
  npos(2) = 1
  cc(1,:) = zh(:)
  call ccnf_put_vara(profilencid,'zh',spos(1:2),npos(1:2),cc)
  cc(1,1:kl) = sigmh(1:kl)*ps(1)
  cc(1,kl+1) = 0.
  call ccnf_put_vara(profilencid,'ph',spos(1:2),npos(1:2),cc)
  cc(1,1:kl) = wtflux(1,1:kl)
  cc(1,kl+1) = 0.
  call ccnf_put_vara(profilencid,'wt',spos(1:2),npos(1:2),cc)
  cc(1,1:kl) = wq_flux(1,1:kl)
  cc(1,kl+1) = 0.
  call ccnf_put_vara(profilencid,'wq',spos(1:2),npos(1:2),cc)
  cc(1,1:kl) = uw_flux(1,1:kl)
  cc(1,kl+1) = 0.
  call ccnf_put_vara(profilencid,'uw',spos(1:2),npos(1:2),cc)
  cc(1,1:kl) = vw_flux(1,1:kl)
  cc(1,kl+1) = 0.
  call ccnf_put_vara(profilencid,'vw',spos(1:2),npos(1:2),cc)
  ! cc(:,:) = nf90_fill_float
  ! call ccnf_put_vara(profilencid,'uu',spos(1:2),npos(1:2),cc)
  ! cc(:,:) = nf90_fill_float
  ! call ccnf_put_vara(profilencid,'vv',spos(1:2),npos(1:2),cc)
  ! cc(:,:) = nf90_fill_float
  ! call ccnf_put_vara(profilencid,'ww',spos(1:2),npos(1:2),cc)
  ! cc(:,:) = nf90_fill_float
  ! call ccnf_put_vara(profilencid,'th^2',spos(1:2),npos(1:2),cc)
  cc(1,1) = 0.
  cc(1,2:kl) = rkmsave(1,1:kl-1)
  cc(1,kl+1) = 0.
  call ccnf_put_vara(profilencid,'Km',spos(1:2),npos(1:2),cc)
  cc(1,1) = 0.
  cc(1,2:kl) = rkhsave(1,1:kl-1)
  cc(1,kl+1) = 0.
  call ccnf_put_vara(profilencid,'Kh',spos(1:2),npos(1:2),cc)
  if ( nvmix==6 ) then
    cc(1,1) = 0.
    cc(1,2:kl) = mfsave(1,1:kl-1)
    cc(1,kl+1) = 0.  
  else
    cc(1,:) = nf90_fill_float
  end if
  ! call ccnf_put_vara(profilencid,'mf',spos(1:2),npos(1:2),cc)
  npos(1) = kl
  ! bb(1,:) = sw_tend(1,1:kl) + lw_tend(1,1:kl)
  ! call ccnf_put_vara(profilencid,'dT_dt_rad',spos(1:2),npos(1:2),bb)
  ! bb(1,:) = sw_tend(1,1:kl)
  ! call ccnf_put_vara(profilencid,'dT_dt_swrad',spos(1:2),npos(1:2),bb)
  ! bb(1,:) = lw_tend(1,1:kl)
  ! call ccnf_put_vara(profilencid,'dT_dt_lwrad',spos(1:2),npos(1:2),bb)
  bb(1,:) = tkesave(1,1:kl)
  call ccnf_put_vara(profilencid,'TKE',spos(1:2),npos(1:2),bb)
  bb(1,:) = shearproduction(1,1:kl)
  call ccnf_put_vara(profilencid,'shear',spos(1:2),npos(1:2),bb)
  bb(1,:) = buoyproduction(1,1:kl)
  call ccnf_put_vara(profilencid,'buoy',spos(1:2),npos(1:2),bb)
  bb(1,:) = totaltransport(1,1:kl)
  call ccnf_put_vara(profilencid,'trans',spos(1:2),npos(1:2),bb)
  bb(1,:) = epssave(1,1:kl)
  call ccnf_put_vara(profilencid,'dissi',spos(1:2),npos(1:2),bb)
  npos(1) = kl+1
  ! cc(:,:) = nf90_fill_float
  ! call ccnf_put_vara(profilencid,'TTE',spos(1:2),npos(1:2),cc) 
  
  spos(1) = 1
  spos(2) = iarch
  npos(1) = 3
  npos(2) = 1
  ! call ccnf_put_vara(profilencid,'tsn',spos(1:2),npos(1:2),tggsn)

else if ( scm_mode=="gabls4" ) then

  ! if ( firstcall ) then
  !   write(6,*) "Creating time output file"  
  !   call ccnf_create(timeoutput,timencid)
  !   ! Turn off the data filling
  !   call ccnf_nofill(timencid)
    
  !   call ccnf_put_attg(timencid,'model','CCAM')
  !   call ccnf_put_attg(timencid,'contact','jing.duke@gmail.com')
  !   call ccnf_put_attg(timencid,'type','Stretched grid climate model')
  !   call ccnf_put_attg(timencid,'timestep',dt)
    
  !   if ( nvmix==7 ) then
  !     call ccnf_put_attg(timencid,'turbulencescheme','Jing')  
  !     call ccnf_put_attg(timencid,'eddydiff','?')
  !     call ccnf_put_attg(timencid,'massflux','?')
  !     call ccnf_put_attg(timencid,'lengthscale','?')
  !     call ccnf_put_attg(timencid,'kprofile','?')
  !   else if ( nvmix==6 ) then
  !     call ccnf_put_attg(timencid,'turbulencescheme','EDMF')
  !     call ccnf_put_attg(timencid,'eddydiff','TKE-EPS')
  !     call ccnf_put_attg(timencid,'massflux','updraft')
  !   else
  !     call ccnf_put_attg(timencid,'turbulencescheme','K-profile')  
  !     call ccnf_put_attg(timencid,'eddydiff','?')
  !     call ccnf_put_attg(timencid,'massflux','diagnosed')
  !     call ccnf_put_attg(timencid,'lengthscale','?')
  !     call ccnf_put_attg(timencid,'kprofile','Richardson')
  !   end if
    
  !   call ccnf_put_attg(timencid,'soillayers','6')
  !   call ccnf_put_attg(timencid,'snowlayers','1-3')
  !   call ccnf_put_attg(timencid,'surfaceprog','temp,moisture,ice,canopywater,snowtemp,snowmass,snowage')
    
  !   call ccnf_def_dimu(timencid,'time',tdim_time)

  !   jdim(1) = tdim_time
  !   call ccnf_def_var(timencid,'time','float',1,jdim(1:1),idnt)
  !   call ccnf_put_att(timencid,idnt,'point_spacing','even')
  !   icy = kdate/10000
  !   icm = max(1, min(12, (kdate-icy*10000)/100))
  !   icd = max(1, min(31, (kdate-icy*10000-icm*100)))
  !   ich = ktime/100
  !   icmi = (ktime-ich*100)
  !   ics = 0
  !   write(grdtim,'("seconds since ",i4.4,"-",i2.2,"-",i2.2," ",2(i2.2,":"),i2.2)') icy,icm,icd,ich,icmi,ics
  !   call ccnf_put_att(timencid,idnt,'units',grdtim)    
    
    
  !   jdim(1) = tdim_time
  !   lname = 'Long wave downward radiation at surface'
  !   call ccnf_def_var(timencid,'lwdw','float',1,jdim(1:1),idnt)
  !   call ccnf_put_att(timencid,idnt,'long_name',lname)
  !   call ccnf_put_att(timencid,idnt,'units','W/m2')
  !   call ccnf_put_att(timencid,idnt,'missing_value',nf90_fill_float)
  !   lname = 'Long wave upward radiation at surface'
  !   call ccnf_def_var(timencid,'lwup','float',1,jdim(1:1),idnt)
  !   call ccnf_put_att(timencid,idnt,'long_name',lname)
  !   call ccnf_put_att(timencid,idnt,'units','W/m2')
  !   call ccnf_put_att(timencid,idnt,'missing_value',nf90_fill_float)
  !   lname = 'Short wave downward radiation at surface'
  !   call ccnf_def_var(timencid,'swdw','float',1,jdim(1:1),idnt)
  !   call ccnf_put_att(timencid,idnt,'long_name',lname)
  !   call ccnf_put_att(timencid,idnt,'units','W/m2')
  !   call ccnf_put_att(timencid,idnt,'missing_value',nf90_fill_float)
  !   lname = 'Short wave upward radiation at surface'
  !   call ccnf_def_var(timencid,'swup','float',1,jdim(1:1),idnt)
  !   call ccnf_put_att(timencid,idnt,'long_name',lname)
  !   call ccnf_put_att(timencid,idnt,'units','W/m2')
  !   call ccnf_put_att(timencid,idnt,'missing_value',nf90_fill_float)
  !   lname = 'Sensible heat flux'
  !   call ccnf_def_var(timencid,'shf','float',1,jdim(1:1),idnt)
  !   call ccnf_put_att(timencid,idnt,'long_name',lname)
  !   call ccnf_put_att(timencid,idnt,'units','W/m2')
  !   call ccnf_put_att(timencid,idnt,'missing_value',nf90_fill_float)
  !   lname = 'Latent (Liq+sol) heat flux'
  !   call ccnf_def_var(timencid,'lhf','float',1,jdim(1:1),idnt)
  !   call ccnf_put_att(timencid,idnt,'long_name',lname)
  !   call ccnf_put_att(timencid,idnt,'units','W/m2')
  !   call ccnf_put_att(timencid,idnt,'missing_value',nf90_fill_float)
  !   lname = 'Evaporation + sublimation flux'
  !   call ccnf_def_var(timencid,'evap','float',1,jdim(1:1),idnt)
  !   call ccnf_put_att(timencid,idnt,'long_name',lname)
  !   call ccnf_put_att(timencid,idnt,'units','mm/day')
  !   call ccnf_put_att(timencid,idnt,'missing_value',nf90_fill_float)
  !   lname = 'Friction velocity'
  !   call ccnf_def_var(timencid,'ustar','float',1,jdim(1:1),idnt)
  !   call ccnf_put_att(timencid,idnt,'long_name',lname)
  !   call ccnf_put_att(timencid,idnt,'units','m/s')
  !   call ccnf_put_att(timencid,idnt,'missing_value',nf90_fill_float)
  !   lname = 'Precipitation (liq+sol) rate'
  !   call ccnf_def_var(timencid,'rain','float',1,jdim(1:1),idnt)
  !   call ccnf_put_att(timencid,idnt,'long_name',lname)
  !   call ccnf_put_att(timencid,idnt,'units','mm/day')
  !   call ccnf_put_att(timencid,idnt,'missing_value',nf90_fill_float)
  !   lname = 'Surface pressure'
  !   call ccnf_def_var(timencid,'psurf','float',1,jdim(1:1),idnt)
  !   call ccnf_put_att(timencid,idnt,'long_name',lname)
  !   call ccnf_put_att(timencid,idnt,'units','Pa')
  !   call ccnf_put_att(timencid,idnt,'missing_value',nf90_fill_float)
  !   lname = 'Boundary layer height'
  !   call ccnf_def_var(timencid,'hpbl','float',1,jdim(1:1),idnt)
  !   call ccnf_put_att(timencid,idnt,'long_name',lname)
  !   call ccnf_put_att(timencid,idnt,'units','m')
  !   call ccnf_put_att(timencid,idnt,'missing_value',nf90_fill_float)
  !   lname = 'Surface temperature'
  !   call ccnf_def_var(timencid,'tsurf','float',1,jdim(1:1),idnt)
  !   call ccnf_put_att(timencid,idnt,'long_name',lname)
  !   call ccnf_put_att(timencid,idnt,'units','K')
  !   call ccnf_put_att(timencid,idnt,'missing_value',nf90_fill_float)
  !   lname = 'Radiative temperature if different from tsurf'
  !   call ccnf_def_var(timencid,'trad','float',1,jdim(1:1),idnt)
  !   call ccnf_put_att(timencid,idnt,'long_name',lname)
  !   call ccnf_put_att(timencid,idnt,'units','K')
  !   call ccnf_put_att(timencid,idnt,'missing_value',nf90_fill_float)
  !   lname = 'Surface albedo'
  !   call ccnf_def_var(timencid,'alb','float',1,jdim(1:1),idnt)
  !   call ccnf_put_att(timencid,idnt,'long_name',lname)
  !   call ccnf_put_att(timencid,idnt,'units','0-1')
  !   call ccnf_put_att(timencid,idnt,'missing_value',nf90_fill_float)
  !   lname = 'Momentum roughness length'
  !   call ccnf_def_var(timencid,'z0m','float',1,jdim(1:1),idnt)
  !   call ccnf_put_att(timencid,idnt,'long_name',lname)
  !   call ccnf_put_att(timencid,idnt,'units','m')
  !   call ccnf_put_att(timencid,idnt,'missing_value',nf90_fill_float)
  !   lname = 'Heat roughness length'
  !   call ccnf_def_var(timencid,'z0h','float',1,jdim(1:1),idnt)
  !   call ccnf_put_att(timencid,idnt,'long_name',lname)
  !   call ccnf_put_att(timencid,idnt,'units','m')
  !   call ccnf_put_att(timencid,idnt,'missing_value',nf90_fill_float)
  !   lname = 'Surface emissivity'
  !   call ccnf_def_var(timencid,'emis','float',1,jdim(1:1),idnt)
  !   call ccnf_put_att(timencid,idnt,'long_name',lname)
  !   call ccnf_put_att(timencid,idnt,'units','0-1')
  !   call ccnf_put_att(timencid,idnt,'missing_value',nf90_fill_float)
  !   lname = '2m temperature'
  !   call ccnf_def_var(timencid,'t2m','float',1,jdim(1:1),idnt)
  !   call ccnf_put_att(timencid,idnt,'long_name',lname)
  !   call ccnf_put_att(timencid,idnt,'units','K')
  !   call ccnf_put_att(timencid,idnt,'missing_value',nf90_fill_float)
  !   lname = '2m specific humidity'
  !   call ccnf_def_var(timencid,'q2m','float',1,jdim(1:1),idnt)
  !   call ccnf_put_att(timencid,idnt,'long_name',lname)
  !   call ccnf_put_att(timencid,idnt,'units','kg/kg')  
  !   call ccnf_put_att(timencid,idnt,'missing_value',nf90_fill_float)
  !   lname = '2m relative humidity'
  !   call ccnf_def_var(timencid,'rh2m','float',1,jdim(1:1),idnt)
  !   call ccnf_put_att(timencid,idnt,'long_name',lname)
  !   call ccnf_put_att(timencid,idnt,'units','0-100') 
  !   call ccnf_put_att(timencid,idnt,'missing_value',nf90_fill_float)
  !   lname = '10m u-component wind'
  !   call ccnf_def_var(timencid,'u10m','float',1,jdim(1:1),idnt)
  !   call ccnf_put_att(timencid,idnt,'long_name',lname)
  !   call ccnf_put_att(timencid,idnt,'units','m/s')
  !   call ccnf_put_att(timencid,idnt,'missing_value',nf90_fill_float)
  !   lname = '10m v-component wind'
  !   call ccnf_def_var(timencid,'v10m','float',1,jdim(1:1),idnt)
  !   call ccnf_put_att(timencid,idnt,'long_name',lname)
  !   call ccnf_put_att(timencid,idnt,'units','m/s')  
  !   call ccnf_put_att(timencid,idnt,'missing_value',nf90_fill_float)
  !   lname = 'Temperature at 3.3 meter above the surface'
  !   call ccnf_def_var(timencid,'t3m','float',1,jdim(1:1),idnt)
  !   call ccnf_put_att(timencid,idnt,'long_name',lname)
  !   call ccnf_put_att(timencid,idnt,'units','K')  
  !   call ccnf_put_att(timencid,idnt,'missing_value',nf90_fill_float)    
  !   lname = 'Specific humdity at 3.3 meter above surface'
  !   call ccnf_def_var(timencid,'q3m','float',1,jdim(1:1),idnt)
  !   call ccnf_put_att(timencid,idnt,'long_name',lname)
  !   call ccnf_put_att(timencid,idnt,'units','kg/kg')  
  !   call ccnf_put_att(timencid,idnt,'missing_value',nf90_fill_float)
  !   lname = 'Relative humidity at 3.3 meter above the surface'
  !   call ccnf_def_var(timencid,'rh3m','float',1,jdim(1:1),idnt)
  !   call ccnf_put_att(timencid,idnt,'long_name',lname)
  !   call ccnf_put_att(timencid,idnt,'units','0-100')  
  !   call ccnf_put_att(timencid,idnt,'missing_value',nf90_fill_float)
  !   lname = 'u-component wind at 3.3 meter above surface'
  !   call ccnf_def_var(timencid,'u3m','float',1,jdim(1:1),idnt)
  !   call ccnf_put_att(timencid,idnt,'long_name',lname)
  !   call ccnf_put_att(timencid,idnt,'units','m/s')  
  !   call ccnf_put_att(timencid,idnt,'missing_value',nf90_fill_float)
  !   lname = 'v-component wind at 3.3 meter above surface'
  !   call ccnf_def_var(timencid,'v3m','float',1,jdim(1:1),idnt)
  !   call ccnf_put_att(timencid,idnt,'long_name',lname)
  !   call ccnf_put_att(timencid,idnt,'units','m/s')  
  !   call ccnf_put_att(timencid,idnt,'missing_value',nf90_fill_float)
  !   lname = 'Temperature at 8.8 meter above the surface'
  !   call ccnf_def_var(timencid,'t9m','float',1,jdim(1:1),idnt)
  !   call ccnf_put_att(timencid,idnt,'long_name',lname)
  !   call ccnf_put_att(timencid,idnt,'units','K')  
  !   call ccnf_put_att(timencid,idnt,'missing_value',nf90_fill_float)    
  !   lname = 'Specific humdity at 8.8 meter above surface'
  !   call ccnf_def_var(timencid,'q9m','float',1,jdim(1:1),idnt)
  !   call ccnf_put_att(timencid,idnt,'long_name',lname)
  !   call ccnf_put_att(timencid,idnt,'units','kg/kg')  
  !   call ccnf_put_att(timencid,idnt,'missing_value',nf90_fill_float)
  !   lname = 'Relative humidity at 8.8 meter above the surface'
  !   call ccnf_def_var(timencid,'rh9m','float',1,jdim(1:1),idnt)
  !   call ccnf_put_att(timencid,idnt,'long_name',lname)
  !   call ccnf_put_att(timencid,idnt,'units','0-100')  
  !   call ccnf_put_att(timencid,idnt,'missing_value',nf90_fill_float)
  !   lname = 'u-component wind at 8.8 meter above surface'
  !   call ccnf_def_var(timencid,'u9m','float',1,jdim(1:1),idnt)
  !   call ccnf_put_att(timencid,idnt,'long_name',lname)
  !   call ccnf_put_att(timencid,idnt,'units','m/s')  
  !   call ccnf_put_att(timencid,idnt,'missing_value',nf90_fill_float)
  !   lname = 'v-component wind at 8.8 meter above surface'
  !   call ccnf_def_var(timencid,'v9m','float',1,jdim(1:1),idnt)
  !   call ccnf_put_att(timencid,idnt,'long_name',lname)
  !   call ccnf_put_att(timencid,idnt,'units','m/s')  
  !   call ccnf_put_att(timencid,idnt,'missing_value',nf90_fill_float)
  !   lname = 'Temperature at 17.9 meter above the surface'
  !   call ccnf_def_var(timencid,'t18m','float',1,jdim(1:1),idnt)
  !   call ccnf_put_att(timencid,idnt,'long_name',lname)
  !   call ccnf_put_att(timencid,idnt,'units','K')  
  !   call ccnf_put_att(timencid,idnt,'missing_value',nf90_fill_float)    
  !   lname = 'Specific humdity at 17.9 meter above surface'
  !   call ccnf_def_var(timencid,'q18m','float',1,jdim(1:1),idnt)
  !   call ccnf_put_att(timencid,idnt,'long_name',lname)
  !   call ccnf_put_att(timencid,idnt,'units','kg/kg')  
  !   call ccnf_put_att(timencid,idnt,'missing_value',nf90_fill_float)
  !   lname = 'Relative humidity at 17.9 meter above the surface'
  !   call ccnf_def_var(timencid,'rh18m','float',1,jdim(1:1),idnt)
  !   call ccnf_put_att(timencid,idnt,'long_name',lname)
  !   call ccnf_put_att(timencid,idnt,'units','0-100')  
  !   call ccnf_put_att(timencid,idnt,'missing_value',nf90_fill_float)
  !   lname = 'u-component wind at 17.9 meter above surface'
  !   call ccnf_def_var(timencid,'u18m','float',1,jdim(1:1),idnt)
  !   call ccnf_put_att(timencid,idnt,'long_name',lname)
  !   call ccnf_put_att(timencid,idnt,'units','m/s')  
  !   call ccnf_put_att(timencid,idnt,'missing_value',nf90_fill_float)
  !   lname = 'v-component wind at 17.9 meter above surface'
  !   call ccnf_def_var(timencid,'v18m','float',1,jdim(1:1),idnt)
  !   call ccnf_put_att(timencid,idnt,'long_name',lname)
  !   call ccnf_put_att(timencid,idnt,'units','m/s')  
  !   call ccnf_put_att(timencid,idnt,'missing_value',nf90_fill_float)
  !   lname = 'Temperature at 25.3 meter above the surface'
  !   call ccnf_def_var(timencid,'t25m','float',1,jdim(1:1),idnt)
  !   call ccnf_put_att(timencid,idnt,'long_name',lname)
  !   call ccnf_put_att(timencid,idnt,'units','K')  
  !   call ccnf_put_att(timencid,idnt,'missing_value',nf90_fill_float)    
  !   lname = 'Specific humdity at 25.3 meter above surface'
  !   call ccnf_def_var(timencid,'q25m','float',1,jdim(1:1),idnt)
  !   call ccnf_put_att(timencid,idnt,'long_name',lname)
  !   call ccnf_put_att(timencid,idnt,'units','kg/kg')  
  !   call ccnf_put_att(timencid,idnt,'missing_value',nf90_fill_float)
  !   lname = 'Relative humidity at 25.3 meter above the surface'
  !   call ccnf_def_var(timencid,'rh25m','float',1,jdim(1:1),idnt)
  !   call ccnf_put_att(timencid,idnt,'long_name',lname)
  !   call ccnf_put_att(timencid,idnt,'units','0-100')  
  !   call ccnf_put_att(timencid,idnt,'missing_value',nf90_fill_float)
  !   lname = 'u-component wind at 25.3 meter above surface'
  !   call ccnf_def_var(timencid,'u25m','float',1,jdim(1:1),idnt)
  !   call ccnf_put_att(timencid,idnt,'long_name',lname)
  !   call ccnf_put_att(timencid,idnt,'units','m/s')  
  !   call ccnf_put_att(timencid,idnt,'missing_value',nf90_fill_float)
  !   lname = 'v-component wind at 25.3 meter above surface'
  !   call ccnf_def_var(timencid,'v25m','float',1,jdim(1:1),idnt)
  !   call ccnf_put_att(timencid,idnt,'long_name',lname)
  !   call ccnf_put_att(timencid,idnt,'units','m/s')  
  !   call ccnf_put_att(timencid,idnt,'missing_value',nf90_fill_float)
  !   lname = 'Temperature at 32.7 meter above the surface'
  !   call ccnf_def_var(timencid,'t33m','float',1,jdim(1:1),idnt)
  !   call ccnf_put_att(timencid,idnt,'long_name',lname)
  !   call ccnf_put_att(timencid,idnt,'units','K')  
  !   call ccnf_put_att(timencid,idnt,'missing_value',nf90_fill_float)    
  !   lname = 'Specific humdity at 32.7 meter above surface'
  !   call ccnf_def_var(timencid,'q33m','float',1,jdim(1:1),idnt)
  !   call ccnf_put_att(timencid,idnt,'long_name',lname)
  !   call ccnf_put_att(timencid,idnt,'units','kg/kg')  
  !   call ccnf_put_att(timencid,idnt,'missing_value',nf90_fill_float)
  !   lname = 'Relative humidity at 32.7 meter above the surface'
  !   call ccnf_def_var(timencid,'rh33m','float',1,jdim(1:1),idnt)
  !   call ccnf_put_att(timencid,idnt,'long_name',lname)
  !   call ccnf_put_att(timencid,idnt,'units','0-100')  
  !   call ccnf_put_att(timencid,idnt,'missing_value',nf90_fill_float)
  !   lname = 'u-component wind at 32.7 meter above surface'
  !   call ccnf_def_var(timencid,'u33m','float',1,jdim(1:1),idnt)
  !   call ccnf_put_att(timencid,idnt,'long_name',lname)
  !   call ccnf_put_att(timencid,idnt,'units','m/s')  
  !   call ccnf_put_att(timencid,idnt,'missing_value',nf90_fill_float)
  !   lname = 'v-component wind at 32.7 meter above surface'
  !   call ccnf_def_var(timencid,'v33m','float',1,jdim(1:1),idnt)
  !   call ccnf_put_att(timencid,idnt,'long_name',lname)
  !   call ccnf_put_att(timencid,idnt,'units','m/s')  
  !   call ccnf_put_att(timencid,idnt,'missing_value',nf90_fill_float)
  !   lname = 'Temperature at 41.9 meter above the surface'
  !   call ccnf_def_var(timencid,'t42m','float',1,jdim(1:1),idnt)
  !   call ccnf_put_att(timencid,idnt,'long_name',lname)
  !   call ccnf_put_att(timencid,idnt,'units','K')  
  !   call ccnf_put_att(timencid,idnt,'missing_value',nf90_fill_float)    
  !   lname = 'Specific humdity at 41.9 meter above surface'
  !   call ccnf_def_var(timencid,'q42m','float',1,jdim(1:1),idnt)
  !   call ccnf_put_att(timencid,idnt,'long_name',lname)
  !   call ccnf_put_att(timencid,idnt,'units','kg/kg')  
  !   call ccnf_put_att(timencid,idnt,'missing_value',nf90_fill_float)
  !   lname = 'Relative humidity at 41.9 meter above the surface'
  !   call ccnf_def_var(timencid,'rh42m','float',1,jdim(1:1),idnt)
  !   call ccnf_put_att(timencid,idnt,'long_name',lname)
  !   call ccnf_put_att(timencid,idnt,'units','0-100')  
  !   call ccnf_put_att(timencid,idnt,'missing_value',nf90_fill_float)
  !   lname = 'u-component wind at 41.9 meter above surface'
  !   call ccnf_def_var(timencid,'u42m','float',1,jdim(1:1),idnt)
  !   call ccnf_put_att(timencid,idnt,'long_name',lname)
  !   call ccnf_put_att(timencid,idnt,'units','m/s')  
  !   call ccnf_put_att(timencid,idnt,'missing_value',nf90_fill_float)
  !   lname = 'v-component wind at 41.9 meter above surface'
  !   call ccnf_def_var(timencid,'v42m','float',1,jdim(1:1),idnt)
  !   call ccnf_put_att(timencid,idnt,'long_name',lname)
  !   call ccnf_put_att(timencid,idnt,'units','m/s')  
  !   call ccnf_put_att(timencid,idnt,'missing_value',nf90_fill_float)
  !   lname = 'Total cloudcover fraction'
  !   call ccnf_def_var(timencid,'cc','float',1,jdim(1:1),idnt)
  !   call ccnf_put_att(timencid,idnt,'long_name',lname)
  !   call ccnf_put_att(timencid,idnt,'units','0-1')  
  !   call ccnf_put_att(timencid,idnt,'missing_value',nf90_fill_float)
  !   lname = 'Vertical flux u-component momentum at 3.3 meter above surface'
  !   call ccnf_def_var(timencid,'uw_3m','float',1,jdim(1:1),idnt)
  !   call ccnf_put_att(timencid,idnt,'long_name',lname)
  !   call ccnf_put_att(timencid,idnt,'units','m2/s2')  
  !   call ccnf_put_att(timencid,idnt,'missing_value',nf90_fill_float)
  !   lname = 'Vertical flux v-component momentum at 3.3 meter above surface'
  !   call ccnf_def_var(timencid,'vw_3m','float',1,jdim(1:1),idnt)
  !   call ccnf_put_att(timencid,idnt,'long_name',lname)
  !   call ccnf_put_att(timencid,idnt,'units','m2/s2')  
  !   call ccnf_put_att(timencid,idnt,'missing_value',nf90_fill_float)
  !   lname = 'Vertical temperature flux at 3.3 meter above surface'
  !   call ccnf_def_var(timencid,'wt_3m','float',1,jdim(1:1),idnt)
  !   call ccnf_put_att(timencid,idnt,'long_name',lname)
  !   call ccnf_put_att(timencid,idnt,'units','Km/s')  
  !   call ccnf_put_att(timencid,idnt,'missing_value',nf90_fill_float)
  !   lname = 'Turbulent kinetic energy at 3.3 meter above surface'
  !   call ccnf_def_var(timencid,'TKE_3m','float',1,jdim(1:1),idnt)
  !   call ccnf_put_att(timencid,idnt,'long_name',lname)
  !   call ccnf_put_att(timencid,idnt,'units','m^2/s^2')  
  !   call ccnf_put_att(timencid,idnt,'missing_value',nf90_fill_float)
  !   lname = 'Vertical flux u-component momentum at 7.03 meter above surface'
  !   call ccnf_def_var(timencid,'uw_7m','float',1,jdim(1:1),idnt)
  !   call ccnf_put_att(timencid,idnt,'long_name',lname)
  !   call ccnf_put_att(timencid,idnt,'units','m2/s2')  
  !   call ccnf_put_att(timencid,idnt,'missing_value',nf90_fill_float)
  !   lname = 'Vertical flux v-component momentum at 7.03 meter above surface'
  !   call ccnf_def_var(timencid,'vw_7m','float',1,jdim(1:1),idnt)
  !   call ccnf_put_att(timencid,idnt,'long_name',lname)
  !   call ccnf_put_att(timencid,idnt,'units','m2/s2')  
  !   call ccnf_put_att(timencid,idnt,'missing_value',nf90_fill_float)
  !   lname = 'Vertical temperature flux at 7.03 meter above surface'
  !   call ccnf_def_var(timencid,'wt_7m','float',1,jdim(1:1),idnt)
  !   call ccnf_put_att(timencid,idnt,'long_name',lname)
  !   call ccnf_put_att(timencid,idnt,'units','Km/s')  
  !   call ccnf_put_att(timencid,idnt,'missing_value',nf90_fill_float)
  !   lname = 'Turbulent kinetic energy at 7.03 meter above surface'
  !   call ccnf_def_var(timencid,'TKE_7m','float',1,jdim(1:1),idnt)
  !   call ccnf_put_att(timencid,idnt,'long_name',lname)
  !   call ccnf_put_att(timencid,idnt,'units','m^2/s^2')  
  !   call ccnf_put_att(timencid,idnt,'missing_value',nf90_fill_float)
  !   lname = 'Vertical flux u-component momentum at 15.43 meter above surface'
  !   call ccnf_def_var(timencid,'uw_15m','float',1,jdim(1:1),idnt)
  !   call ccnf_put_att(timencid,idnt,'long_name',lname)
  !   call ccnf_put_att(timencid,idnt,'units','m2/s2')  
  !   call ccnf_put_att(timencid,idnt,'missing_value',nf90_fill_float)
  !   lname = 'Vertical flux v-component momentum at 15.43 meter above surface'
  !   call ccnf_def_var(timencid,'vw_15m','float',1,jdim(1:1),idnt)
  !   call ccnf_put_att(timencid,idnt,'long_name',lname)
  !   call ccnf_put_att(timencid,idnt,'units','m2/s2')  
  !   call ccnf_put_att(timencid,idnt,'missing_value',nf90_fill_float)
  !   lname = 'Vertical temperature flux at 15.43 meter above surface'
  !   call ccnf_def_var(timencid,'wt_15m','float',1,jdim(1:1),idnt)
  !   call ccnf_put_att(timencid,idnt,'long_name',lname)
  !   call ccnf_put_att(timencid,idnt,'units','Km/s')  
  !   call ccnf_put_att(timencid,idnt,'missing_value',nf90_fill_float)
  !   lname = 'Turbulent kinetic energy at 15.43 meter above surface'
  !   call ccnf_def_var(timencid,'TKE_15m','float',1,jdim(1:1),idnt)
  !   call ccnf_put_att(timencid,idnt,'long_name',lname)
  !   call ccnf_put_att(timencid,idnt,'units','m^2/s^2')  
  !   call ccnf_put_att(timencid,idnt,'missing_value',nf90_fill_float)
  !   lname = 'Vertical flux u-component momentum at 22.79 meter above surface'
  !   call ccnf_def_var(timencid,'uw_23m','float',1,jdim(1:1),idnt)
  !   call ccnf_put_att(timencid,idnt,'long_name',lname)
  !   call ccnf_put_att(timencid,idnt,'units','m2/s2')  
  !   call ccnf_put_att(timencid,idnt,'missing_value',nf90_fill_float)
  !   lname = 'Vertical flux v-component momentum at 22.79 meter above surface'
  !   call ccnf_def_var(timencid,'vw_23m','float',1,jdim(1:1),idnt)
  !   call ccnf_put_att(timencid,idnt,'long_name',lname)
  !   call ccnf_put_att(timencid,idnt,'units','m2/s2')  
  !   call ccnf_put_att(timencid,idnt,'missing_value',nf90_fill_float)
  !   lname = 'Vertical temperature flux at 22.79 meter above surface'
  !   call ccnf_def_var(timencid,'wt_23m','float',1,jdim(1:1),idnt)
  !   call ccnf_put_att(timencid,idnt,'long_name',lname)
  !   call ccnf_put_att(timencid,idnt,'units','Km/s')  
  !   call ccnf_put_att(timencid,idnt,'missing_value',nf90_fill_float)
  !   lname = 'Turbulent kinetic energy at 22.79 meter above surface'
  !   call ccnf_def_var(timencid,'TKE_23m','float',1,jdim(1:1),idnt)
  !   call ccnf_put_att(timencid,idnt,'long_name',lname)
  !   call ccnf_put_att(timencid,idnt,'units','m^2/s^2')  
  !   call ccnf_put_att(timencid,idnt,'missing_value',nf90_fill_float)
  !   lname = 'Vertical flux u-component momentum at 30.15 meter above surface'
  !   call ccnf_def_var(timencid,'uw_30m','float',1,jdim(1:1),idnt)
  !   call ccnf_put_att(timencid,idnt,'long_name',lname)
  !   call ccnf_put_att(timencid,idnt,'units','m2/s2')  
  !   call ccnf_put_att(timencid,idnt,'missing_value',nf90_fill_float)
  !   lname = 'Vertical flux v-component momentum at 30.15 meter above surface'
  !   call ccnf_def_var(timencid,'vw_30m','float',1,jdim(1:1),idnt)
  !   call ccnf_put_att(timencid,idnt,'long_name',lname)
  !   call ccnf_put_att(timencid,idnt,'units','m2/s2')  
  !   call ccnf_put_att(timencid,idnt,'missing_value',nf90_fill_float)
  !   lname = 'Vertical temperature flux at 30.15 meter above surface'
  !   call ccnf_def_var(timencid,'wt_30m','float',1,jdim(1:1),idnt)
  !   call ccnf_put_att(timencid,idnt,'long_name',lname)
  !   call ccnf_put_att(timencid,idnt,'units','Km/s')  
  !   call ccnf_put_att(timencid,idnt,'missing_value',nf90_fill_float)
  !   lname = 'Turbulent kinetic energy at 30.15 meter above surface'
  !   call ccnf_def_var(timencid,'TKE_30m','float',1,jdim(1:1),idnt)
  !   call ccnf_put_att(timencid,idnt,'long_name',lname)
  !   call ccnf_put_att(timencid,idnt,'units','m^2/s^2')  
  !   call ccnf_put_att(timencid,idnt,'missing_value',nf90_fill_float)
  !   lname = 'Vertical flux u-component momentum at 37.51 meter above surface'
  !   call ccnf_def_var(timencid,'uw_38m','float',1,jdim(1:1),idnt)
  !   call ccnf_put_att(timencid,idnt,'long_name',lname)
  !   call ccnf_put_att(timencid,idnt,'units','m2/s2')  
  !   call ccnf_put_att(timencid,idnt,'missing_value',nf90_fill_float)
  !   lname = 'Vertical flux v-component momentum at 37.51 meter above surface'
  !   call ccnf_def_var(timencid,'vw_38m','float',1,jdim(1:1),idnt)
  !   call ccnf_put_att(timencid,idnt,'long_name',lname)
  !   call ccnf_put_att(timencid,idnt,'units','m2/s2')  
  !   call ccnf_put_att(timencid,idnt,'missing_value',nf90_fill_float)
  !   lname = 'Vertical temperature flux at 37.51 meter above surface'
  !   call ccnf_def_var(timencid,'wt_38m','float',1,jdim(1:1),idnt)
  !   call ccnf_put_att(timencid,idnt,'long_name',lname)
  !   call ccnf_put_att(timencid,idnt,'units','Km/s')  
  !   call ccnf_put_att(timencid,idnt,'missing_value',nf90_fill_float)
  !   lname = 'Turbulent kinetic energy at 37.51 meter above surface'
  !   call ccnf_def_var(timencid,'TKE_38m','float',1,jdim(1:1),idnt)
  !   call ccnf_put_att(timencid,idnt,'long_name',lname)
  !   call ccnf_put_att(timencid,idnt,'units','m^2/s^2')  
  !   call ccnf_put_att(timencid,idnt,'missing_value',nf90_fill_float)

  !   call ccnf_enddef(timencid)
    
  !   write(6,*) "Creating profile output file"  
  !   call ccnf_create(profileoutput,profilencid)
  !   ! Turn off the data filling
  !   call ccnf_nofill(profilencid)
    
  !   call ccnf_put_attg(profilencid,'model','CCAM')
  !   call ccnf_put_attg(profilencid,'contact','jing.duke@gmail.com')
  !   call ccnf_put_attg(profilencid,'type','Stretched grid climate model')
  !   call ccnf_put_attg(profilencid,'timestep',dt)
    
  !   if ( nvmix==7 ) then
  !     call ccnf_put_attg(profilencid,'turbulencescheme','Jing')  
  !     call ccnf_put_attg(profilencid,'eddydiff','?')
  !     call ccnf_put_attg(profilencid,'massflux','?')
  !     call ccnf_put_attg(profilencid,'lengthscale','?')
  !     call ccnf_put_attg(profilencid,'kprofile','?')
  !   else if ( nvmix==6 ) then
  !     call ccnf_put_attg(profilencid,'turbulencescheme','EDMF')
  !     call ccnf_put_attg(profilencid,'eddydiff','TKE-EPS')
  !     call ccnf_put_attg(profilencid,'massflux','updraft')
  !   else
  !     call ccnf_put_attg(profilencid,'turbulencescheme','K-profile')  
  !     call ccnf_put_attg(profilencid,'eddydiff','?')
  !     call ccnf_put_attg(profilencid,'massflux','diagnosed')
  !     call ccnf_put_attg(profilencid,'lengthscale','?')
  !     call ccnf_put_attg(profilencid,'kprofile','Richardson')
  !   end if
    
  !   call ccnf_put_attg(profilencid,'soillayers','6')
  !   call ccnf_put_attg(profilencid,'snowlayers','1-3')
  !   call ccnf_put_attg(profilencid,'surfaceprog','temp,moisture,ice,canopywater,snowtemp,snowmass,snowage')
    
  !   call ccnf_def_dim(profilencid,'levf',kl,zfdim)
  !   call ccnf_def_dim(profilencid,'levh',kl+1,zhdim)
  !   call ccnf_def_dim(profilencid,'levs',3,zsdim) ! snow levels
  !   call ccnf_def_dimu(profilencid,'time',tdim_prof)

  !   jdim(1) = zfdim
  !   call ccnf_def_var(profilencid,'levf','float',1,jdim(1:1),idnt)
  !   call ccnf_put_att(profilencid,idnt,'units','Pa')

  !   jdim(1) = zhdim
  !   call ccnf_def_var(profilencid,'levh','float',1,jdim(1:1),idnt)
  !   call ccnf_put_att(profilencid,idnt,'units','Pa')

  !   jdim(1) = zsdim
  !   call ccnf_def_var(profilencid,'levs','float',1,jdim(1:1),idnt)
  !   call ccnf_put_att(profilencid,idnt,'units','index')
    
  !   jdim(1) = tdim_prof
  !   call ccnf_def_var(profilencid,'time','float',1,jdim(1:1),idnt)
  !   call ccnf_put_att(profilencid,idnt,'point_spacing','even')
  !   icy = kdate/10000
  !   icm = max(1, min(12, (kdate-icy*10000)/100))
  !   icd = max(1, min(31, (kdate-icy*10000-icm*100)))
  !   ich = ktime/100
  !   icmi = (ktime-ich*100)
  !   ics = 0
  !   write(grdtim,'("seconds since ",i4.4,"-",i2.2,"-",i2.2," ",2(i2.2,":"),i2.2)') icy,icm,icd,ich,icmi,ics
  !   call ccnf_put_att(profilencid,idnt,'units',grdtim)    

    
  !   jdim(1) = zfdim
  !   jdim(2) = tdim_prof
  !   lname = 'Height of full level'
  !   call ccnf_def_var(profilencid,'zf','float',2,jdim(1:2),idnt)
  !   call ccnf_put_att(profilencid,idnt,'long_name',lname)
  !   call ccnf_put_att(profilencid,idnt,'units','m')  
  !   call ccnf_put_att(profilencid,idnt,'missing_value',nf90_fill_float)
  !   lname = 'Presure at full level'
  !   call ccnf_def_var(profilencid,'pf','float',2,jdim(1:2),idnt)
  !   call ccnf_put_att(profilencid,idnt,'long_name',lname)
  !   call ccnf_put_att(profilencid,idnt,'units','Pa')  
  !   call ccnf_put_att(profilencid,idnt,'missing_value',nf90_fill_float)
  !   lname = 'Temperature'
  !   call ccnf_def_var(profilencid,'t','float',2,jdim(1:2),idnt)
  !   call ccnf_put_att(profilencid,idnt,'long_name',lname)
  !   call ccnf_put_att(profilencid,idnt,'units','K')  
  !   call ccnf_put_att(profilencid,idnt,'missing_value',nf90_fill_float)
  !   lname = 'Temperature'
  !   call ccnf_def_var(profilencid,'temp','float',2,jdim(1:2),idnt)
  !   call ccnf_put_att(profilencid,idnt,'long_name',lname)
  !   call ccnf_put_att(profilencid,idnt,'units','K')  
  !   call ccnf_put_att(profilencid,idnt,'missing_value',nf90_fill_float)
  !   lname = 'Potential temperature'
  !   call ccnf_def_var(profilencid,'th','float',2,jdim(1:2),idnt)
  !   call ccnf_put_att(profilencid,idnt,'long_name',lname)
  !   call ccnf_put_att(profilencid,idnt,'units','K')  
  !   call ccnf_put_att(profilencid,idnt,'missing_value',nf90_fill_float)
  !   lname = 'Specific humidity'
  !   call ccnf_def_var(profilencid,'q','float',2,jdim(1:2),idnt)
  !   call ccnf_put_att(profilencid,idnt,'long_name',lname)
  !   call ccnf_put_att(profilencid,idnt,'units','kg/kg')  
  !   call ccnf_put_att(profilencid,idnt,'missing_value',nf90_fill_float)
  !   lname = 'Cloud water and ice'
  !   call ccnf_def_var(profilencid,'qc','float',2,jdim(1:2),idnt)
  !   call ccnf_put_att(profilencid,idnt,'long_name',lname)
  !   call ccnf_put_att(profilencid,idnt,'units','kg/kg')  
  !   call ccnf_put_att(profilencid,idnt,'missing_value',nf90_fill_float)
  !   lname = 'Zonal component wind'
  !   call ccnf_def_var(profilencid,'u','float',2,jdim(1:2),idnt)
  !   call ccnf_put_att(profilencid,idnt,'long_name',lname)
  !   call ccnf_put_att(profilencid,idnt,'units','m/s')  
  !   call ccnf_put_att(profilencid,idnt,'missing_value',nf90_fill_float)
  !   lname = 'Meridional component wind'
  !   call ccnf_def_var(profilencid,'v','float',2,jdim(1:2),idnt)
  !   call ccnf_put_att(profilencid,idnt,'long_name',lname)
  !   call ccnf_put_att(profilencid,idnt,'units','m/s')  
  !   call ccnf_put_att(profilencid,idnt,'missing_value',nf90_fill_float)

  !   jdim(1) = zfdim
  !   jdim(2) = tdim_prof
  !   lname = 'u-component geostrophic wind'
  !   call ccnf_def_var(profilencid,'ugeo','float',2,jdim(1:2),idnt)
  !   call ccnf_put_att(profilencid,idnt,'long_name',lname)
  !   call ccnf_put_att(profilencid,idnt,'units','m/s')  
  !   call ccnf_put_att(profilencid,idnt,'missing_value',nf90_fill_float)  
  !   lname = 'v-component geostrophic wind'
  !   call ccnf_def_var(profilencid,'vgeo','float',2,jdim(1:2),idnt)
  !   call ccnf_put_att(profilencid,idnt,'long_name',lname)
  !   call ccnf_put_att(profilencid,idnt,'units','m/s')  
  !   call ccnf_put_att(profilencid,idnt,'missing_value',nf90_fill_float)  
  !   lname = 'u-component momentum advection'
  !   call ccnf_def_var(profilencid,'dudt_ls','float',2,jdim(1:2),idnt)
  !   call ccnf_put_att(profilencid,idnt,'long_name',lname)
  !   call ccnf_put_att(profilencid,idnt,'units','m/s/s')  
  !   call ccnf_put_att(profilencid,idnt,'missing_value',nf90_fill_float)  
  !   lname = 'v-component momentum advection'
  !   call ccnf_def_var(profilencid,'dvdt_ls','float',2,jdim(1:2),idnt)
  !   call ccnf_put_att(profilencid,idnt,'long_name',lname)
  !   call ccnf_put_att(profilencid,idnt,'units','m/s/s')  
  !   call ccnf_put_att(profilencid,idnt,'missing_value',nf90_fill_float)  
  !   lname = 'Temperature advection'
  !   call ccnf_def_var(profilencid,'dtdt_ls','float',2,jdim(1:2),idnt)
  !   call ccnf_put_att(profilencid,idnt,'long_name',lname)
  !   call ccnf_put_att(profilencid,idnt,'units','K/s')  
  !   call ccnf_put_att(profilencid,idnt,'missing_value',nf90_fill_float)  
  !   lname = 'Moisture advection advection'
  !   call ccnf_def_var(profilencid,'dqdt_ls','float',2,jdim(1:2),idnt)
  !   call ccnf_put_att(profilencid,idnt,'long_name',lname)
  !   call ccnf_put_att(profilencid,idnt,'units','kg/kg/s')  
  !   call ccnf_put_att(profilencid,idnt,'missing_value',nf90_fill_float)  
  !   lname = 'vertical movement'
  !   call ccnf_def_var(profilencid,'w','float',2,jdim(1:2),idnt)
  !   call ccnf_put_att(profilencid,idnt,'long_name',lname)
  !   call ccnf_put_att(profilencid,idnt,'units','m/s')  
  !   call ccnf_put_att(profilencid,idnt,'missing_value',nf90_fill_float)  
     
  !   jdim(1) = zhdim
  !   jdim(2) = tdim_prof
  !   lname = 'Height of half level'
  !   call ccnf_def_var(profilencid,'zh','float',2,jdim(1:2),idnt)
  !   call ccnf_put_att(profilencid,idnt,'long_name',lname)
  !   call ccnf_put_att(profilencid,idnt,'units','m')  
  !   call ccnf_put_att(profilencid,idnt,'missing_value',nf90_fill_float)
  !   lname = 'Presure at half level'
  !   call ccnf_def_var(profilencid,'ph','float',2,jdim(1:2),idnt)
  !   call ccnf_put_att(profilencid,idnt,'long_name',lname)
  !   call ccnf_put_att(profilencid,idnt,'units','Pa')  
  !   call ccnf_put_att(profilencid,idnt,'missing_value',nf90_fill_float)
  !   lname = 'Vertical temperature flux'
  !   call ccnf_def_var(profilencid,'wt','float',2,jdim(1:2),idnt)
  !   call ccnf_put_att(profilencid,idnt,'long_name',lname)
  !   call ccnf_put_att(profilencid,idnt,'units','Km/s')  
  !   call ccnf_put_att(profilencid,idnt,'missing_value',nf90_fill_float)
  !   lname = 'Vertical moisture flux'
  !   call ccnf_def_var(profilencid,'wq','float',2,jdim(1:2),idnt)
  !   call ccnf_put_att(profilencid,idnt,'long_name',lname)
  !   call ccnf_put_att(profilencid,idnt,'units','kg/kg m/s')  
  !   call ccnf_put_att(profilencid,idnt,'missing_value',nf90_fill_float)
  !   lname = 'Vertical flux u-component momentum'
  !   call ccnf_def_var(profilencid,'uw','float',2,jdim(1:2),idnt)
  !   call ccnf_put_att(profilencid,idnt,'long_name',lname)
  !   call ccnf_put_att(profilencid,idnt,'units','m2/s2')  
  !   call ccnf_put_att(profilencid,idnt,'missing_value',nf90_fill_float)
  !   lname = 'Vertical flux v-component momentum'
  !   call ccnf_def_var(profilencid,'vw','float',2,jdim(1:2),idnt)
  !   call ccnf_put_att(profilencid,idnt,'long_name',lname)
  !   call ccnf_put_att(profilencid,idnt,'units','m/s')  
  !   call ccnf_put_att(profilencid,idnt,'missing_value',nf90_fill_float)
  !   lname = 'u-variance'
  !   call ccnf_def_var(profilencid,'uu','float',2,jdim(1:2),idnt)
  !   call ccnf_put_att(profilencid,idnt,'long_name',lname)
  !   call ccnf_put_att(profilencid,idnt,'units','m2/s2')  
  !   call ccnf_put_att(profilencid,idnt,'missing_value',nf90_fill_float)  
  !   lname = 'v-variance'
  !   call ccnf_def_var(profilencid,'vv','float',2,jdim(1:2),idnt)
  !   call ccnf_put_att(profilencid,idnt,'long_name',lname)
  !   call ccnf_put_att(profilencid,idnt,'units','m2/s2')  
  !   call ccnf_put_att(profilencid,idnt,'missing_value',nf90_fill_float)
  !   lname = 'w-variance'
  !   call ccnf_def_var(profilencid,'ww','float',2,jdim(1:2),idnt)
  !   call ccnf_put_att(profilencid,idnt,'long_name',lname)
  !   call ccnf_put_att(profilencid,idnt,'units','m2/s2')  
  !   call ccnf_put_att(profilencid,idnt,'missing_value',nf90_fill_float)  
  !   lname = 'Potential temperature variance'
  !   call ccnf_def_var(profilencid,'th^2','float',2,jdim(1:2),idnt)
  !   call ccnf_put_att(profilencid,idnt,'long_name',lname)
  !   call ccnf_put_att(profilencid,idnt,'units','K^2')  
  !   call ccnf_put_att(profilencid,idnt,'missing_value',nf90_fill_float)  
  !   lname = 'Eddy diffusivity momentum'
  !   call ccnf_def_var(profilencid,'Km','float',2,jdim(1:2),idnt)
  !   call ccnf_put_att(profilencid,idnt,'long_name',lname)
  !   call ccnf_put_att(profilencid,idnt,'units','m2/s')  
  !   call ccnf_put_att(profilencid,idnt,'missing_value',nf90_fill_float)
  !   lname = 'Eddy diffusivity heat'
  !   call ccnf_def_var(profilencid,'Kh','float',2,jdim(1:2),idnt)
  !   call ccnf_put_att(profilencid,idnt,'long_name',lname)
  !   call ccnf_put_att(profilencid,idnt,'units','m2/s')  
  !   call ccnf_put_att(profilencid,idnt,'missing_value',nf90_fill_float)
  !   lname = 'Massflux'
  !   call ccnf_def_var(profilencid,'mf','float',2,jdim(1:2),idnt)
  !   call ccnf_put_att(profilencid,idnt,'long_name',lname)
  !   call ccnf_put_att(profilencid,idnt,'units','kg/m2/s')  
  !   call ccnf_put_att(profilencid,idnt,'missing_value',nf90_fill_float)  
  !   jdim(1) = zfdim
  !   lname = 'Temperature tendency from radiation'
  !   call ccnf_def_var(profilencid,'dT_dt_rad','float',2,jdim(1:2),idnt)
  !   call ccnf_put_att(profilencid,idnt,'long_name',lname)
  !   call ccnf_put_att(profilencid,idnt,'units','K/s')  
  !   call ccnf_put_att(profilencid,idnt,'missing_value',nf90_fill_float)  
  !   lname = 'Temperature tendency from short wave radiation'
  !   call ccnf_def_var(profilencid,'dT_dt_swrad','float',2,jdim(1:2),idnt)
  !   call ccnf_put_att(profilencid,idnt,'long_name',lname)
  !   call ccnf_put_att(profilencid,idnt,'units','K/s')  
  !   call ccnf_put_att(profilencid,idnt,'missing_value',nf90_fill_float)  
  !   lname = 'Temperature tendency from long wave radiation'
  !   call ccnf_def_var(profilencid,'dT_dt_lwrad','float',2,jdim(1:2),idnt)
  !   call ccnf_put_att(profilencid,idnt,'long_name',lname)
  !   call ccnf_put_att(profilencid,idnt,'units','K/s')  
  !   call ccnf_put_att(profilencid,idnt,'missing_value',nf90_fill_float)  
  !   lname = 'Turbulent kinetic energy'
  !   call ccnf_def_var(profilencid,'TKE','float',2,jdim(1:2),idnt)
  !   call ccnf_put_att(profilencid,idnt,'long_name',lname)
  !   call ccnf_put_att(profilencid,idnt,'units','m^2/s^2')  
  !   call ccnf_put_att(profilencid,idnt,'missing_value',nf90_fill_float)  
  !   jdim(1) = zhdim
  !   lname = 'Total turbulent energy'
  !   call ccnf_def_var(profilencid,'TTE','float',2,jdim(1:2),idnt)
  !   call ccnf_put_att(profilencid,idnt,'long_name',lname)
  !   call ccnf_put_att(profilencid,idnt,'units','m^2/s^2')  
  !   call ccnf_put_att(profilencid,idnt,'missing_value',nf90_fill_float)  
  !   lname = 'Shear production'
  !   call ccnf_def_var(profilencid,'shear','float',2,jdim(1:2),idnt)
  !   call ccnf_put_att(profilencid,idnt,'long_name',lname)
  !   call ccnf_put_att(profilencid,idnt,'units','m2/s3')  
  !   call ccnf_put_att(profilencid,idnt,'missing_value',nf90_fill_float)  
  !   lname = 'Bouyancy production'
  !   call ccnf_def_var(profilencid,'buoy','float',2,jdim(1:2),idnt)
  !   call ccnf_put_att(profilencid,idnt,'long_name',lname)
  !   call ccnf_put_att(profilencid,idnt,'units','m2/s3')  
  !   call ccnf_put_att(profilencid,idnt,'missing_value',nf90_fill_float)  
  !   lname = 'Total transport'
  !   call ccnf_def_var(profilencid,'trans','float',2,jdim(1:2),idnt)
  !   call ccnf_put_att(profilencid,idnt,'long_name',lname)
  !   call ccnf_put_att(profilencid,idnt,'units','m2/s3')  
  !   call ccnf_put_att(profilencid,idnt,'missing_value',nf90_fill_float)  
  !   lname = 'Dissipation'
  !   call ccnf_def_var(profilencid,'dissi','float',2,jdim(1:2),idnt)
  !   call ccnf_put_att(profilencid,idnt,'long_name',lname)
  !   call ccnf_put_att(profilencid,idnt,'units','m2/s3')  
  !   call ccnf_put_att(profilencid,idnt,'missing_value',nf90_fill_float)  

  !   jdim(1) = zsdim
  !   jdim(2) = tdim_prof
  !   lname = 'Snow temperature'
  !   call ccnf_def_var(profilencid,'tsn','float',2,jdim(1:2),idnt)
  !   call ccnf_put_att(profilencid,idnt,'long_name',lname)
  !   call ccnf_put_att(profilencid,idnt,'units','K')  
  !   call ccnf_put_att(profilencid,idnt,'missing_value',nf90_fill_float)  
      
  !   call ccnf_enddef(profilencid)
    
  !   bb(1,:) = sig(:)*ps(1)
  !   spos(1) = 1
  !   npos(1) = kl
  !   call ccnf_put_vara(profilencid,'levf',spos(1),npos(1),bb(1,:))

  !   cc(1,1:kl) = sigmh(1:kl)*ps(1)
  !   cc(1,kl+1) = 0.
  !   spos(1) = 1
  !   npos(1) = kl+1
  !   call ccnf_put_vara(profilencid,'levh',spos(1),npos(1),cc(1,:))
    
  !   dd(1,1) = sdepth(1,1)*0.5
  !   dd(1,2) = sdepth(1,1)+sdepth(1,2)*0.5
  !   dd(1,3) = sdepth(1,1)+sdepth(1,2)+sdepth(1,3)*0.5
  !   spos(1) = 1
  !   npos(1) = 3
  !   call ccnf_put_vara(profilencid,'levs',spos(1),npos(1),dd(1,:))
    
  !   iarch = 1  
  ! end if

  ! zf(1) = bet(1)*t(1,1)
  ! do k=2,kl
  !   zf(k) = zf(k-1) + bet(k)*t(1,k) + betm(k)*t(1,k-1)
  ! end do
  ! zf(:) = zf(:)/grav

  ! ! using t.  Maybe should use tv?
  ! zh(1) = 0.
  ! do k=2,kl
  !   zh(k) = zh(k-1) + t(1,k)*(-rdry/grav)*dsig(k)/sig(k)
  ! end do
  ! zh(kl+1) = 2.*zh(kl) - zh(kl-1)

  ! pf(1:kl) = sig(1:kl)*ps(1)

  ! qs(1:kl) = qsat(pf,t(1,:))
  ! rh(:) = qg(1,:)/qs(:)

  ! do k=1,kl-1
  !   wtflux(1,k+1) = wth_flux(1,k)*sigmh(k)**(rdry/cp)
  ! end do
  ! wtflux(1,kl+1) = 0.

  ! call ccnf_put_vara(timencid,'time',iarch,real(ktau)*dt)

  ! if ( ktau==0 ) then
  !   aa(:) = nf90_fill_float    
  !   call ccnf_put_vara(timencid,'lwdw',iarch,aa(1))
  !   call ccnf_put_vara(timencid,'lwup',iarch,aa(1))
  !   call ccnf_put_vara(timencid,'swdw',iarch,aa(1))
  !   call ccnf_put_vara(timencid,'swup',iarch,aa(1))
  !   call ccnf_put_vara(timencid,'shf',iarch,aa(1))
  !   call ccnf_put_vara(timencid,'lhf',iarch,aa(1))
  !   call ccnf_put_vara(timencid,'evap',iarch,aa(1))
  !   call ccnf_put_vara(timencid,'ustar',iarch,aa(1))
  !   call ccnf_put_vara(timencid,'rain',iarch,aa(1))
  !   call ccnf_put_vara(timencid,'psurf',iarch,aa(1))
  !   call ccnf_put_vara(timencid,'hpbl',iarch,aa(1))
  !   call ccnf_put_vara(timencid,'tsurf',iarch,aa(1))
  !   call ccnf_put_vara(timencid,'trad',iarch,aa(1))
  !   call ccnf_put_vara(timencid,'alb',iarch,aa(1))
  !   call ccnf_put_vara(timencid,'z0m',iarch,aa(1))
  !   call ccnf_put_vara(timencid,'z0h',iarch,aa(1))
  !   call ccnf_put_vara(timencid,'emis',iarch,aa(1))
  !   call ccnf_put_vara(timencid,'t2m',iarch,aa(1))
  !   call ccnf_put_vara(timencid,'q2m',iarch,aa(1))
  !   call ccnf_put_vara(timencid,'rh2m',iarch,aa(1))
  !   call ccnf_put_vara(timencid,'u10m',iarch,aa(1))
  !   call ccnf_put_vara(timencid,'v10m',iarch,aa(1))
  !   call ccnf_put_vara(timencid,'t3m',iarch,aa(1))
  !   call ccnf_put_vara(timencid,'q3m',iarch,aa(1))
  !   call ccnf_put_vara(timencid,'rh3m',iarch,aa(1))
  !   call ccnf_put_vara(timencid,'u3m',iarch,aa(1))
  !   call ccnf_put_vara(timencid,'v3m',iarch,aa(1))
  !   call ccnf_put_vara(timencid,'t9m',iarch,aa(1))
  !   call ccnf_put_vara(timencid,'q9m',iarch,aa(1))
  !   call ccnf_put_vara(timencid,'rh9m',iarch,aa(1))
  !   call ccnf_put_vara(timencid,'u9m',iarch,aa(1))
  !   call ccnf_put_vara(timencid,'v9m',iarch,aa(1))
  !   call ccnf_put_vara(timencid,'t18m',iarch,aa(1))
  !   call ccnf_put_vara(timencid,'q18m',iarch,aa(1))
  !   call ccnf_put_vara(timencid,'rh18m',iarch,aa(1))
  !   call ccnf_put_vara(timencid,'u18m',iarch,aa(1))
  !   call ccnf_put_vara(timencid,'v18m',iarch,aa(1))
  !   call ccnf_put_vara(timencid,'t25m',iarch,aa(1))
  !   call ccnf_put_vara(timencid,'q25m',iarch,aa(1))
  !   call ccnf_put_vara(timencid,'rh25m',iarch,aa(1))
  !   call ccnf_put_vara(timencid,'u25m',iarch,aa(1))
  !   call ccnf_put_vara(timencid,'v25m',iarch,aa(1))
  !   call ccnf_put_vara(timencid,'t33m',iarch,aa(1))
  !   call ccnf_put_vara(timencid,'q33m',iarch,aa(1))
  !   call ccnf_put_vara(timencid,'rh33m',iarch,aa(1))
  !   call ccnf_put_vara(timencid,'u33m',iarch,aa(1))
  !   call ccnf_put_vara(timencid,'v33m',iarch,aa(1))
  !   call ccnf_put_vara(timencid,'t2m',iarch,aa(1))
  !   call ccnf_put_vara(timencid,'q42m',iarch,aa(1))
  !   call ccnf_put_vara(timencid,'rh42m',iarch,aa(1))
  !   call ccnf_put_vara(timencid,'u42m',iarch,aa(1))
  !   call ccnf_put_vara(timencid,'v42m',iarch,aa(1))
  !   call ccnf_put_vara(timencid,'cc',iarch,aa(1))
  !   call ccnf_put_vara(timencid,'uw_3m',iarch,aa(1))
  !   call ccnf_put_vara(timencid,'vw_3m',iarch,aa(1))
  !   call ccnf_put_vara(timencid,'wt_3m',iarch,aa(1))
  !   call ccnf_put_vara(timencid,'TKE_3m',iarch,aa(1))
  !   call ccnf_put_vara(timencid,'uw_7m',iarch,aa(1))
  !   call ccnf_put_vara(timencid,'vw_7m',iarch,aa(1))
  !   call ccnf_put_vara(timencid,'wt_7m',iarch,aa(1))
  !   call ccnf_put_vara(timencid,'TKE_7m',iarch,aa(1))
  !   call ccnf_put_vara(timencid,'uw_15m',iarch,aa(1))
  !   call ccnf_put_vara(timencid,'vw_15m',iarch,aa(1))
  !   call ccnf_put_vara(timencid,'wt_15m',iarch,aa(1))
  !   call ccnf_put_vara(timencid,'TKE_15m',iarch,aa(1))
  !   call ccnf_put_vara(timencid,'uw_23m',iarch,aa(1))
  !   call ccnf_put_vara(timencid,'vw_23m',iarch,aa(1))
  !   call ccnf_put_vara(timencid,'wt_23m',iarch,aa(1))
  !   call ccnf_put_vara(timencid,'TKE_23m',iarch,aa(1))
  !   call ccnf_put_vara(timencid,'uw_30m',iarch,aa(1))
  !   call ccnf_put_vara(timencid,'vw_30m',iarch,aa(1))
  !   call ccnf_put_vara(timencid,'wt_30m',iarch,aa(1))
  !   call ccnf_put_vara(timencid,'TKE_30m',iarch,aa(1))
  !   call ccnf_put_vara(timencid,'uw_38m',iarch,aa(1))
  !   call ccnf_put_vara(timencid,'vw_38m',iarch,aa(1))
  !   call ccnf_put_vara(timencid,'wt_38m',iarch,aa(1))
  !   call ccnf_put_vara(timencid,'TKE_38m',iarch,aa(1))
    
  ! else
      
  !   call ccnf_put_vara(timencid,'lwdw',iarch,rgdn_ave(1))
  !   aa(:) = rgdn_ave(:) + rgn_ave(:)
  !   call ccnf_put_vara(timencid,'lwup',iarch,aa(1))
  !   call ccnf_put_vara(timencid,'swdw',iarch,sgdn_ave(1))
  !   aa(:) = sgdn_ave(:) - sgn_ave(:)
  !   call ccnf_put_vara(timencid,'swup',iarch,aa(1))
  !   call ccnf_put_vara(timencid,'shf',iarch,fg(1))
  !   call ccnf_put_vara(timencid,'lhf',iarch,eg(1))
  !   aa(:) = eg(:)*86400./hls ! this assumes sublimation
  !   call ccnf_put_vara(timencid,'evap',iarch,aa(1))
  !   call ccnf_put_vara(timencid,'ustar',iarch,ustar(1))
  !   call ccnf_put_vara(timencid,'rain',iarch,precip(1))
  !   call ccnf_put_vara(timencid,'psurf',iarch,ps(1))
  !   call ccnf_put_vara(timencid,'hpbl',iarch,pblh(1))
  !   call ccnf_put_vara(timencid,'tsurf',iarch,tss(1))
  !   aa(:) = ((rgdn_ave(:)+rgn_ave(:))/(0.98*stefbo))**(0.25)
  !   call ccnf_put_vara(timencid,'trad',iarch,aa(1))
  !   aa(:) = swrsave*albvisnir(:,1) + (1.-swrsave)*albvisnir(:,2)
  !   call ccnf_put_vara(timencid,'alb',iarch,aa(1))
  !   call ccnf_put_vara(timencid,'z0m',iarch,zo(1))
  !   call ccnf_put_vara(timencid,'z0h',iarch,zoh(1))
  !   aa(:) = nf90_fill_float
  !   call ccnf_put_vara(timencid,'emis',iarch,aa(1))
  !   call ccnf_put_vara(timencid,'t2m',iarch,tscrn(1))
  !   aa(:) = qgscrn(:)/(1.+qgscrn(:))
  !   call ccnf_put_vara(timencid,'q2m',iarch,aa(1))
  !   call ccnf_put_vara(timencid,'rh2m',iarch,rhscrn(1))
  !   tmp(:)=u(1,:)
  !   call vout(tmp,aa(1),zf,10.,kl)
  !   call ccnf_put_vara(timencid,'u10m',iarch,aa(1))
  !   tmp(:)=v(1,:)
  !   call vout(tmp,aa(1),zf,10.,kl)
  !   call ccnf_put_vara(timencid,'v10m',iarch,aa(1))
  !   tmp=t(1,:)
  !   call vout(tmp,aa(1),zf,3.,kl)
  !   call ccnf_put_vara(timencid,'t3m',iarch,aa(1))
  !   tmp=qg(1,:)
  !   call vout(tmp,aa(1),zf,3.,kl)
  !   aa(1) = aa(1)/(1.+aa(1))
  !   call ccnf_put_vara(timencid,'q3m',iarch,aa(1))
  !   call vout(rh(:),aa(1),zf,3.,kl)
  !   call ccnf_put_vara(timencid,'rh3m',iarch,aa(1))
  !   tmp(:)=u(1,:)
  !   call vout(tmp,aa(1),zf,3.,kl)
  !   call ccnf_put_vara(timencid,'u3m',iarch,aa(1))
  !   tmp(:)=v(1,:)
  !   call vout(tmp,aa(1),zf,3.,kl)
  !   call ccnf_put_vara(timencid,'v3m',iarch,aa(1))
  !   tmp=t(1,:)
  !   call vout(tmp,aa(1),zf,9.,kl)
  !   call ccnf_put_vara(timencid,'t9m',iarch,aa(1))
  !   tmp=qg(1,:)
  !   call vout(tmp,aa(1),zf,9.,kl)
  !   aa(1) = aa(1)/(1.+aa(1))
  !   call ccnf_put_vara(timencid,'q9m',iarch,aa(1))
  !   call vout(rh(:),aa(1),zf,9.,kl)
  !   call ccnf_put_vara(timencid,'rh9m',iarch,aa(1))
  !   tmp(:)=u(1,:)
  !   call vout(tmp,aa(1),zf,9.,kl)
  !   call ccnf_put_vara(timencid,'u9m',iarch,aa(1))
  !   tmp(:)=v(1,:)
  !   call vout(tmp,aa(1),zf,9.,kl)
  !   call ccnf_put_vara(timencid,'v9m',iarch,aa(1))
  !   tmp=t(1,:)
  !   call vout(tmp,aa(1),zf,18.,kl)
  !   call ccnf_put_vara(timencid,'t18m',iarch,aa(1))
  !   tmp=qg(1,:)
  !   call vout(tmp,aa(1),zf,18.,kl)
  !   aa(1) = aa(1)/(1.+aa(1))
  !   call ccnf_put_vara(timencid,'q18m',iarch,aa(1))
  !   call vout(rh(:),aa(1),zf,18.,kl)
  !   call ccnf_put_vara(timencid,'rh18m',iarch,aa(1))
  !   tmp(:)=u(1,:)
  !   call vout(tmp,aa(1),zf,18.,kl)
  !   call ccnf_put_vara(timencid,'u18m',iarch,aa(1))
  !   tmp(:)=v(1,:)
  !   call vout(tmp,aa(1),zf,18.,kl)
  !   call ccnf_put_vara(timencid,'v18m',iarch,aa(1))
  !   tmp=t(1,:)
  !   call vout(tmp,aa(1),zf,25.,kl)
  !   call ccnf_put_vara(timencid,'t25m',iarch,aa(1))
  !   tmp=qg(1,:)
  !   call vout(tmp,aa(1),zf,25.,kl)
  !   aa(1) = aa(1)/(1.+aa(1))
  !   call ccnf_put_vara(timencid,'q25m',iarch,aa(1))
  !   call vout(rh(:),aa(1),zf,25.,kl)
  !   call ccnf_put_vara(timencid,'rh25m',iarch,aa(1))
  !   tmp(:)=u(1,:)
  !   call vout(tmp,aa(1),zf,25.,kl)
  !   call ccnf_put_vara(timencid,'u25m',iarch,aa(1))
  !   tmp(:)=v(1,:)
  !   call vout(tmp,aa(1),zf,25.,kl)
  !   call ccnf_put_vara(timencid,'v25m',iarch,aa(1))
  !   tmp=t(1,:)
  !   call vout(tmp,aa(1),zf,33.,kl)
  !   call ccnf_put_vara(timencid,'t33m',iarch,aa(1))
  !   tmp=qg(1,:)
  !   call vout(tmp,aa(1),zf,33.,kl)
  !   aa(1) = aa(1)/(1.+aa(1))
  !   call ccnf_put_vara(timencid,'q33m',iarch,aa(1))
  !   call vout(rh(:),aa(1),zf,33.,kl)
  !   call ccnf_put_vara(timencid,'rh33m',iarch,aa(1))
  !   tmp(:)=u(1,:)
  !   call vout(tmp,aa(1),zf,33.,kl)
  !   call ccnf_put_vara(timencid,'u33m',iarch,aa(1))
  !   tmp(:)=v(1,:)
  !   call vout(tmp,aa(1),zf,33.,kl)
  !   call ccnf_put_vara(timencid,'v33m',iarch,aa(1))
  !   tmp=t(1,:)
  !   call vout(tmp,aa(1),zf,42.,kl)
  !   call ccnf_put_vara(timencid,'t42m',iarch,aa(1))
  !   tmp=qg(1,:)
  !   call vout(tmp,aa(1),zf,42,kl)
  !   aa(1) = aa(1)/(1.+aa(1))
  !   call ccnf_put_vara(timencid,'q42m',iarch,aa(1))
  !   call vout(rh(:),aa(1),zf,42.,kl)
  !   call ccnf_put_vara(timencid,'rh42m',iarch,aa(1))
  !   tmp(:)=u(1,:)
  !   call vout(tmp,aa(1),zf,42.,kl)
  !   call ccnf_put_vara(timencid,'u42m',iarch,aa(1))
  !   tmp(:)=v(1,:)
  !   call vout(tmp,aa(1),zf,42.,kl)
  !   call ccnf_put_vara(timencid,'v42m',iarch,aa(1))
  !   aa(1) = cld_ave(1)
  !   call ccnf_put_vara(timencid,'cc',iarch,aa(1))
  !   tmp=uw_flux(1,:)
  !   call vout(tmp,aa(1),zh,3.,kl-1)
  !   call ccnf_put_vara(timencid,'uw_3m',iarch,aa(1))
  !   tmp=vw_flux(1,:)
  !   call vout(tmp,aa(1),zh,3.,kl-1)
  !   call ccnf_put_vara(timencid,'vw_3m',iarch,aa(1))
  !   tmp=wtflux(1,:)
  !   call vout(tmp,aa(1),zh,3.,kl-1)
  !   call ccnf_put_vara(timencid,'wt_3m',iarch,aa(1))
  !   tmp=tkesave(1,:)
  !   call vout(tmp,aa(1),zf,3.,kl)
  !   call ccnf_put_vara(timencid,'TKE_3m',iarch,aa(1))
  !   tmp=uw_flux(1,:)
  !   call vout(tmp,aa(1),zh,7.,kl-1)
  !   call ccnf_put_vara(timencid,'uw_7m',iarch,aa(1))
  !   tmp=vw_flux(1,:)
  !   call vout(tmp,aa(1),zh,7.,kl-1)
  !   call ccnf_put_vara(timencid,'vw_7m',iarch,aa(1))
  !   tmp=wtflux(1,:)
  !   call vout(tmp,aa(1),zh,7.,kl-1)
  !   call ccnf_put_vara(timencid,'wt_7m',iarch,aa(1))
  !   tmp=tkesave(1,:)
  !   call vout(tmp,aa(1),zf,7.,kl)
  !   call ccnf_put_vara(timencid,'TKE_7m',iarch,aa(1))
  !   tmp=uw_flux(1,:)
  !   call vout(tmp,aa(1),zh,15.,kl-1)
  !   call ccnf_put_vara(timencid,'uw_15m',iarch,aa(1))
  !   tmp=vw_flux(1,:)
  !   call vout(tmp,aa(1),zh,15.,kl-1)
  !   call ccnf_put_vara(timencid,'vw_15m',iarch,aa(1))
  !   tmp=wtflux(1,:)
  !   call vout(tmp,aa(1),zh,15.,kl-1)
  !   call ccnf_put_vara(timencid,'wt_15m',iarch,aa(1))
  !   tmp=tkesave(1,:)
  !   call vout(tmp,aa(1),zf,15.,kl)
  !   call ccnf_put_vara(timencid,'TKE_15m',iarch,aa(1))
  !   tmp=uw_flux(1,:)
  !   call vout(tmp,aa(1),zh,23.,kl-1)
  !   call ccnf_put_vara(timencid,'uw_23m',iarch,aa(1))
  !   tmp=vw_flux(1,:)
  !   call vout(tmp,aa(1),zh,23.,kl-1)
  !   call ccnf_put_vara(timencid,'vw_23m',iarch,aa(1))
  !   tmp=wtflux(1,:)
  !   call vout(tmp,aa(1),zh,23.,kl-1)
  !   call ccnf_put_vara(timencid,'wt_23m',iarch,aa(1))
  !   tmp=tkesave(1,:)
  !   call vout(tmp,aa(1),zf,23.,kl)
  !   call ccnf_put_vara(timencid,'TKE_23m',iarch,aa(1))
  !   tmp=uw_flux(1,:)
  !   call vout(tmp,aa(1),zh,30.,kl-1)
  !   call ccnf_put_vara(timencid,'uw_30m',iarch,aa(1))
  !   tmp=vw_flux(1,:)
  !   call vout(tmp,aa(1),zh,30.,kl-1)
  !   call ccnf_put_vara(timencid,'vw_30m',iarch,aa(1))
  !   tmp=wtflux(1,:)
  !   call vout(tmp,aa(1),zh,30.,kl-1)
  !   call ccnf_put_vara(timencid,'wt_30m',iarch,aa(1))
  !   tmp=tkesave(1,:)
  !   call vout(tmp,aa(1),zf,30.,kl)
  !   call ccnf_put_vara(timencid,'TKE_30m',iarch,aa(1))
  !   tmp=uw_flux(1,:)
  !   call vout(tmp,aa(1),zh,38.,kl-1)
  !   call ccnf_put_vara(timencid,'uw_38m',iarch,aa(1))
  !   tmp=vw_flux(1,:)
  !   call vout(tmp,aa(1),zh,38.,kl-1)
  !   call ccnf_put_vara(timencid,'vw_38m',iarch,aa(1))
  !   tmp=wtflux(1,:)
  !   call vout(tmp,aa(1),zh,38.,kl-1)
  !   call ccnf_put_vara(timencid,'wt_38m',iarch,aa(1))
  !   tmp=tkesave(1,:)
  !   call vout(tmp,aa(1),zf,38.,kl)
  !   call ccnf_put_vara(timencid,'TKE_38m',iarch,aa(1))
  ! end if

  ! call ccnf_put_vara(profilencid,'time',iarch,real(ktau)*dt)

  ! spos(1:2) = (/ 1, iarch /)
  ! npos(1:2) = (/ kl, 1 /)
  ! bb(1,:) = zf(:)
  ! call ccnf_put_vara(profilencid,'zf',spos(1:2),npos(1:2),bb)
  ! bb(1,:) = pf(:)
  ! call ccnf_put_vara(profilencid,'pf',spos(1:2),npos(1:2),bb)
  ! call ccnf_put_vara(profilencid,'t',spos(1:2),npos(1:2),t)
  ! call ccnf_put_vara(profilencid,'temp',spos(1:2),npos(1:2),t)
  ! bb(1,:) = t(1,:)*(ps(1)*sig(:)/1.e5)**(-rdry/cp)
  ! call ccnf_put_vara(profilencid,'th',spos(1:2),npos(1:2),bb)
  ! bb(:,:) = qg(:,:)/(1.+qg(:,:))
  ! call ccnf_put_vara(profilencid,'q',spos(1:2),npos(1:2),bb)
  ! bb(:,:) = qlg(:,:) + qfg(:,:)
  ! bb(:,:) = bb(:,:)/(1.+bb(:,:))
  ! call ccnf_put_vara(profilencid,'qc',spos(1:2),npos(1:2),bb)
  ! call ccnf_put_vara(profilencid,'u',spos(1:2),npos(1:2),u)
  ! call ccnf_put_vara(profilencid,'v',spos(1:2),npos(1:2),v)

  ! spos(1) = 1
  ! spos(2) = iarch
  ! npos(1) = kl
  ! npos(2) = 1
  ! bb(1,:) = ug(:)
  ! call ccnf_put_vara(profilencid,'ugeo',spos(1:2),npos(1:2),bb)
  ! bb(1,:) = vg(:)
  ! call ccnf_put_vara(profilencid,'vgeo',spos(1:2),npos(1:2),bb)
  ! bb(:,:) = nf90_fill_float
  ! call ccnf_put_vara(profilencid,'dudt_ls',spos(1:2),npos(1:2),bb)
  ! bb(:,:) = nf90_fill_float
  ! call ccnf_put_vara(profilencid,'dvdt_ls',spos(1:2),npos(1:2),bb)
  ! bb(1,:) = t_tend(:)
  ! call ccnf_put_vara(profilencid,'dtdt_ls',spos(1:2),npos(1:2),bb)
  ! bb(1,:) = q_tend(:)
  ! call ccnf_put_vara(profilencid,'dqdt_ls',spos(1:2),npos(1:2),bb)
  ! bb(:,:) = nf90_fill_float
  ! call ccnf_put_vara(profilencid,'w',spos(1:2),npos(1:2),bb)

  ! spos(1) = 1
  ! spos(2) = iarch
  ! npos(1) = kl+1
  ! npos(2) = 1
  ! cc(1,:) = zh(:)
  ! call ccnf_put_vara(profilencid,'zh',spos(1:2),npos(1:2),cc)
  ! cc(1,1:kl) = sigmh(1:kl)*ps(1)
  ! cc(1,kl+1) = 0.
  ! call ccnf_put_vara(profilencid,'ph',spos(1:2),npos(1:2),cc)
  ! cc(1,1:kl) = wtflux(1,1:kl)
  ! cc(1,kl+1) = 0.
  ! call ccnf_put_vara(profilencid,'wt',spos(1:2),npos(1:2),wtflux)
  ! cc(1,1:kl) = wq_flux(1,1:kl)
  ! cc(1,kl+1) = 0.
  ! call ccnf_put_vara(profilencid,'wq',spos(1:2),npos(1:2),wq_flux)
  ! cc(1,1:kl) = uw_flux(1,1:kl)
  ! cc(1,kl+1) = 0.
  ! call ccnf_put_vara(profilencid,'uw',spos(1:2),npos(1:2),uw_flux)
  ! cc(1,1:kl) = vw_flux(1,1:kl)
  ! cc(1,kl+1) = 0.
  ! call ccnf_put_vara(profilencid,'vw',spos(1:2),npos(1:2),vw_flux)
  ! cc(:,:) = nf90_fill_float
  ! call ccnf_put_vara(profilencid,'uu',spos(1:2),npos(1:2),cc)
  ! cc(:,:) = nf90_fill_float
  ! call ccnf_put_vara(profilencid,'vv',spos(1:2),npos(1:2),cc)
  ! cc(:,:) = nf90_fill_float
  ! call ccnf_put_vara(profilencid,'ww',spos(1:2),npos(1:2),cc)
  ! cc(:,:) = nf90_fill_float
  ! call ccnf_put_vara(profilencid,'th^2',spos(1:2),npos(1:2),cc)
  ! cc(1,1) = 0.
  ! cc(1,2:kl) = rkmsave(1,1:kl-1)
  ! cc(1,kl+1) = 0.
  ! call ccnf_put_vara(profilencid,'Km',spos(1:2),npos(1:2),cc)
  ! cc(1,1) = 0.
  ! cc(1,2:kl) = rkhsave(1,1:kl-1)
  ! cc(1,kl+1) = 0.
  ! call ccnf_put_vara(profilencid,'Kh',spos(1:2),npos(1:2),cc)
  ! if ( nvmix==6 ) then
  !   cc(1,1) = 0.
  !   cc(1,2:kl) = mfsave(1,1:kl-1)
  !   cc(1,kl+1) = 0.  
  ! else
  !   cc(1,:) = nf90_fill_float
  ! end if
  ! call ccnf_put_vara(profilencid,'mf',spos(1:2),npos(1:2),cc)
  ! npos(1) = kl
  ! bb(1,:) = sw_tend(1,1:kl) + lw_tend(1,1:kl)
  ! call ccnf_put_vara(profilencid,'dT_dt_rad',spos(1:2),npos(1:2),bb)
  ! bb(1,:) = sw_tend(1,1:kl)
  ! call ccnf_put_vara(profilencid,'dT_dt_swrad',spos(1:2),npos(1:2),bb)
  ! bb(1,:) = lw_tend(1,1:kl)
  ! call ccnf_put_vara(profilencid,'dT_dt_lwrad',spos(1:2),npos(1:2),bb)
  ! bb(1,:) = tkesave(1,1:kl)
  ! call ccnf_put_vara(profilencid,'TKE',spos(1:2),npos(1:2),bb)
  ! npos(1) = kl+1
  ! cc(:,:) = nf90_fill_float
  ! call ccnf_put_vara(profilencid,'TTE',spos(1:2),npos(1:2),cc)
  ! cc(:,:) = nf90_fill_float
  ! call ccnf_put_vara(profilencid,'shear',spos(1:2),npos(1:2),cc)
  ! cc(:,:) = nf90_fill_float
  ! call ccnf_put_vara(profilencid,'buoy',spos(1:2),npos(1:2),cc)
  ! cc(:,:) = nf90_fill_float
  ! call ccnf_put_vara(profilencid,'trans',spos(1:2),npos(1:2),cc)
  ! cc(:,:) = nf90_fill_float
  ! call ccnf_put_vara(profilencid,'dissi',spos(1:2),npos(1:2),cc)

  ! spos(1) = 1
  ! spos(2) = iarch
  ! npos(1) = 3
  ! npos(2) = 1
  ! call ccnf_put_vara(profilencid,'tsn',spos(1:2),npos(1:2),tggsn)

else
  write(6,*) "ERROR: Unknown option for scm_mode ",trim(scm_mode)
  stop
end if

if ( lastcall ) then
  write(6,*) "Closing output file"
  call ccnf_close(timencid)
  call ccnf_close(profilencid)
end if

iarch = iarch + 1

return
end subroutine outputscm

subroutine vout(dat_in,datout,zf_in,zfout,kl)

implicit none

integer kindex
integer, intent(in) :: kl
real x
real, dimension(kl) :: dat_in, zf_in
real, intent(in) :: zfout
real, intent(out) :: datout

do kindex=2,kl
  if ( zf_in(kindex)>=zfout ) then
    exit
  end if
end do

x = (zfout-zf_in(kindex-1))/(zf_in(kindex)-zf_in(kindex-1))
datout = (1.-x)*dat_in(kindex-1) + x*dat_in(kindex)

return
end subroutine vout
    
    
subroutine insoil
      
use newmpar_m                              ! Grid parameters
use soilv_m                                ! Soil parameters

implicit none
      
integer isoil, k

do isoil = 1,mxst
  cnsd(isoil)  = sand(isoil)*0.3+clay(isoil)*0.25+silt(isoil)*0.265
  hsbh(isoil)  = hyds(isoil)*abs(sucs(isoil))*bch(isoil) !difsat*etasat
  ibp2(isoil)  = nint(bch(isoil))+2
  i2bp3(isoil) = 2*nint(bch(isoil))+3
  write(6,"('isoil,ssat,sfc,swilt,hsbh ',i2,3f7.3,e11.4)") isoil,ssat(isoil),sfc(isoil),swilt(isoil),hsbh(isoil)
end do
cnsd(9) = 2.51

zshh(1)    = .5*zse(1)        ! not used (jlm)
zshh(ms+1) = .5*zse(ms)       ! not used (jlm)
ww(1) = 1.
do k = 2,ms
  zshh(k) = .5*(zse(k-1)+zse(k))  ! z(k)-z(k-1) (jlm)
  ww(k)   = zse(k)/(zse(k)+zse(k-1))
end do

return
end subroutine insoil    

subroutine vinterpolate(dat_in,sig_in,datout,sigout,nlev,kl)    
    
implicit none

integer, intent(in) :: nlev, kl
integer k, kin
real x
real, dimension(nlev), intent(in) :: dat_in, sig_in
real, dimension(kl), intent(in) :: sigout
real, dimension(kl), intent(out) :: datout

do k = 1,kl
  if ( sigout(k)>=sig_in(1) ) then
    ! extrapolate
    datout(k) = dat_in(1)
  elseif ( sigout(k)<=sig_in(nlev) ) then
    ! extrapolate
    datout(k) = dat_in(nlev)
  else
    ! interpolate
    do kin = 2,nlev
      if ( sigout(k)>sig_in(kin) ) exit
    end do      
    x = (sigout(k)-sig_in(kin-1))/(sig_in(kin)-sig_in(kin-1))
    datout(k) = (1.-x)*dat_in(kin-1) + x*dat_in(kin)
  end if
end do

return
end subroutine vinterpolate

subroutine vinterp2m(height_in,height_model,dat_in,datout,nlev,kl)

implicit none

integer, intent(in) :: nlev, kl
integer :: k, kin
real x
real, dimension(nlev), intent(in) :: height_in
real, dimension(nlev), intent(in) :: dat_in
real, dimension(kl), intent(in) :: height_model
real, dimension(kl), intent(out) :: datout

do k = 1,kl
  if ( height_model(k)<=height_in(1) ) then
    ! extrapolate
    datout(k) = dat_in(1)
  elseif ( height_model(k)>=height_in(nlev) ) then
    ! extrapolate
    datout(k) = dat_in(nlev)
  else
    ! interpolate
    do kin = 2,nlev
      if ( height_model(k)<height_in(kin) ) exit
    end do
    x = (height_model(k)-height_in(kin-1))/(height_in(kin)-height_in(kin-1))
    datout(k) = (1.-x)*dat_in(kin-1) + x*dat_in(kin)
  end if
end do

return
end subroutine vinterp2m
    
subroutine finishbanner

implicit none

! End banner
write(6,*) "=============================================================================="
write(6,*) "CCAM: Finished scm"
write(6,*) "=============================================================================="

return
end    

    
!--------------------------------------------------------------------
! Fix water vapour mixing ratio
subroutine fixqg(js,je)

use arrays_m                          ! Atmosphere dyamics prognostic arrays
use const_phys                        ! Physical constants
use liqwpar_m                         ! Cloud water mixing ratios
use newmpar_m                         ! Grid parameters
use parm_m                            ! Model configuration

implicit none

integer, intent(in) :: js, je
integer k
real, dimension(js:je) :: dumqtot, dumliq

if ( js<1 .or. je>ifull ) then
  write(6,*) "ERROR: Invalid index for fixqg"
  stop
end if

do k = 1,kl
  dumqtot(js:je) = qg(js:je,k) + qlg(js:je,k) + qfg(js:je,k) ! qtot
  dumqtot(js:je) = max( dumqtot(js:je), qgmin )
  dumliq(js:je)  = t(js:je,k) - hlcp*qlg(js:je,k) - hlscp*qfg(js:je,k)
  qfg(js:je,k)   = max( qfg(js:je,k), 0. ) 
  qlg(js:je,k)   = max( qlg(js:je,k), 0. )
  qrg(js:je,k)   = max( qrg(js:je,k), 0. )
  qsng(js:je,k)  = max( qsng(js:je,k), 0. )
  qgrg(js:je,k)  = max( qgrg(js:je,k), 0. )
  qg(js:je,k)    = dumqtot(js:je) - qlg(js:je,k) - qfg(js:je,k)
  qg(js:je,k)    = max( qg(js:je,k), 0. )
  t(js:je,k)     = dumliq(js:je) + hlcp*qlg(js:je,k) + hlscp*qfg(js:je,k)
end do

return
end subroutine fixqg   
    
!-------------------------------------------------------------------- 
! Check for NaN errors
subroutine nantest(message,js,je)

use aerosolldr, only : xtg,ssn,naero  ! LDR prognostic aerosols
use arrays_m                          ! Atmosphere dyamics prognostic arrays
use cc_mpi                            ! CC MPI routines
use cfrac_m                           ! Cloud fraction
use extraout_m                        ! Additional diagnostics
use liqwpar_m                         ! Cloud water mixing ratios
use morepbl_m                         ! Additional boundary layer diagnostics
use newmpar_m                         ! Grid parameters
use parm_m                            ! Model configuration
use pbl_m                             ! Boundary layer arrays
use work2_m                           ! Diagnostic arrays
use work3f_m                          ! Grid work arrays

implicit none

integer, intent(in) :: js, je
character(len=*), intent(in) :: message

if ( js<1 .or. je>ifull ) then
  write(6,*) "ERROR: Invalid index for nantest - ",trim(message)
  call ccmpi_abort(-1)
end if

if ( any(t(js:je,1:kl)/=t(js:je,1:kl)) ) then
  write(6,*) "ERROR: NaN detected in t on myid=",myid," at ",trim(message)
  call ccmpi_abort(-1)
end if

if ( any(t(js:je,1:kl)<75.) .or. any(t(js:je,1:kl)>425.) ) then
  write(6,*) "ERROR: Out-of-range detected in t on myid=",myid," at ",trim(message)
  write(6,*) "minval,maxval ",minval(t(js:je,1:kl)),maxval(t(js:je,1:kl))
  write(6,*) "minloc,maxloc ",minloc(t(js:je,1:kl)),maxloc(t(js:je,1:kl))
  call ccmpi_abort(-1)
end if

if ( any(u(js:je,1:kl)/=u(js:je,1:kl)) ) then
  write(6,*) "ERROR: NaN detected in u on myid=",myid," at ",trim(message)
  call ccmpi_abort(-1)
end if

if ( any(u(js:je,1:kl)<-400.) .or. any(u(js:je,1:kl)>400.) ) then
  write(6,*) "ERROR: Out-of-range detected in u on myid=",myid," at ",trim(message)
  write(6,*) "minval,maxval ",minval(u(js:je,1:kl)),maxval(u(js:je,1:kl))
  write(6,*) "minloc,maxloc ",minloc(u(js:je,1:kl)),maxloc(u(js:je,1:kl))
  call ccmpi_abort(-1) 
end if

if ( any(v(js:je,1:kl)/=v(js:je,1:kl)) ) then
  write(6,*) "ERROR: NaN detected in v on myid=",myid," at ",trim(message)
  call ccmpi_abort(-1)
end if

if ( any(v(js:je,1:kl)<-400.) .or. any(v(js:je,1:kl)>400.) ) then
  write(6,*) "ERROR: Out-of-range detected in v on myid=",myid," at ",trim(message)
  write(6,*) "minval,maxval ",minval(v(js:je,1:kl)),maxval(v(js:je,1:kl))
  write(6,*) "minloc,maxloc ",minloc(v(js:je,1:kl)),maxloc(v(js:je,1:kl))
  call ccmpi_abort(-1) 
end if

if ( any(qg(js:je,1:kl)/=qg(js:je,1:kl)) ) then
  write(6,*) "ERROR: NaN detected in qg on myid=",myid," at ",trim(message)
  call ccmpi_abort(-1)
end if

if ( any(qg(js:je,1:kl)<-1.e-8) .or. any(qg(js:je,1:kl)>7.e-2) ) then
  write(6,*) "ERROR: Out-of-range detected in qg on myid=",myid," at ",trim(message)
  write(6,*) "minval,maxval ",minval(qg(js:je,1:kl)),maxval(qg(js:je,1:kl))
  write(6,*) "minloc,maxloc ",minloc(qg(js:je,1:kl)),maxloc(qg(js:je,1:kl))
  call ccmpi_abort(-1) 
end if

if ( any(qlg(js:je,1:kl)/=qlg(js:je,1:kl)) ) then
  write(6,*) "ERROR: NaN detected in qlg on myid=",myid," at ",trim(message)
  call ccmpi_abort(-1)    
end if

if ( any(qlg(js:je,1:kl)<-1.e-8) .or. any(qlg(js:je,1:kl)>7.e-2) ) then
  write(6,*) "ERROR: Out-of-range detected in qlg on myid=",myid," at ",trim(message)
  write(6,*) "minval,maxval ",minval(qlg(js:je,1:kl)),maxval(qlg(js:je,1:kl))
  write(6,*) "minloc,maxloc ",minloc(qlg(js:je,1:kl)),maxloc(qlg(js:je,1:kl))
  call ccmpi_abort(-1) 
end if

if ( any(qfg(js:je,1:kl)/=qfg(js:je,1:kl)) ) then
  write(6,*) "ERROR: NaN detected in qfg on myid=",myid," at ",trim(message)
  call ccmpi_abort(-1)    
end if

if ( any(qfg(js:je,1:kl)<-1.e-8) .or. any(qfg(js:je,1:kl)>7.e-2) ) then
  write(6,*) "ERROR: Out-of-range detected in qfg on myid=",myid," at ",trim(message)
  write(6,*) "minval,maxval ",minval(qfg(js:je,1:kl)),maxval(qfg(js:je,1:kl))
  write(6,*) "minloc,maxloc ",minloc(qfg(js:je,1:kl)),maxloc(qfg(js:je,1:kl))
  call ccmpi_abort(-1) 
end if

if ( any(qrg(js:je,1:kl)/=qrg(js:je,1:kl)) ) then
  write(6,*) "ERROR: NaN detected in qrg on myid=",myid," at ",trim(message)
  call ccmpi_abort(-1)    
end if

if ( any(qrg(js:je,1:kl)<-1.e-8) .or. any(qrg(js:je,1:kl)>7.e-2) ) then
  write(6,*) "ERROR: Out-of-range detected in qrg on myid=",myid," at ",trim(message)
  write(6,*) "minval,maxval ",minval(qrg(js:je,1:kl)),maxval(qrg(js:je,1:kl))
  write(6,*) "minloc,maxloc ",minloc(qrg(js:je,1:kl)),maxloc(qrg(js:je,1:kl))
  call ccmpi_abort(-1) 
end if

if ( any(qsng(js:je,1:kl)/=qsng(js:je,1:kl)) ) then
  write(6,*) "ERROR: NaN detected in qsng on myid=",myid," at ",trim(message)
  call ccmpi_abort(-1)    
end if

if ( any(qsng(js:je,1:kl)<-1.e-8) .or. any(qsng(js:je,1:kl)>7.e-2) ) then
  write(6,*) "ERROR: Out-of-range detected in qsng on myid=",myid," at ",trim(message)
  write(6,*) "minval,maxval ",minval(qsng(js:je,1:kl)),maxval(qsng(js:je,1:kl))
  write(6,*) "minloc,maxloc ",minloc(qsng(js:je,1:kl)),maxloc(qsng(js:je,1:kl))
  call ccmpi_abort(-1) 
end if

if ( any(qgrg(js:je,1:kl)/=qgrg(js:je,1:kl)) ) then
  write(6,*) "ERROR: NaN detected in qgrg on myid=",myid," at ",trim(message)
  call ccmpi_abort(-1)    
end if

if ( any(qgrg(js:je,1:kl)<-1.e-8) .or. any(qgrg(js:je,1:kl)>7.e-2) ) then
  write(6,*) "ERROR: Out-of-range detected in qgrg on myid=",myid," at ",trim(message)
  write(6,*) "minval,maxval ",minval(qgrg(js:je,1:kl)),maxval(qgrg(js:je,1:kl))
  write(6,*) "minloc,maxloc ",minloc(qgrg(js:je,1:kl)),maxloc(qgrg(js:je,1:kl))
  call ccmpi_abort(-1) 
end if

if ( any(qlrad(js:je,1:kl)/=qlrad(js:je,1:kl)) ) then
  write(6,*) "ERROR: NaN detected in qlrad on myid=",myid," at ",trim(message)
  call ccmpi_abort(-1)    
end if

if ( any(qlrad(js:je,1:kl)<-1.e-8) .or. any(qlrad(js:je,1:kl)>7.e-2) ) then
  write(6,*) "ERROR: Out-of-range detected in qlrad on myid=",myid," at ",trim(message)
  write(6,*) "minval,maxval ",minval(qlrad(js:je,1:kl)),maxval(qlrad(js:je,1:kl))
  write(6,*) "minloc,maxloc ",minloc(qlrad(js:je,1:kl)),maxloc(qlrad(js:je,1:kl))
  call ccmpi_abort(-1) 
end if

if ( any(qfrad(js:je,1:kl)/=qfrad(js:je,1:kl)) ) then
  write(6,*) "ERROR: NaN detected in qfrad on myid=",myid," at ",trim(message)
  call ccmpi_abort(-1)    
end if

if ( any(qfrad(js:je,1:kl)<-1.e-8) .or. any(qfrad(js:je,1:kl)>7.e-2) ) then
  write(6,*) "ERROR: Out-of-range detected in qfrad on myid=",myid," at ",trim(message)
  write(6,*) "minval,maxval ",minval(qfrad(js:je,1:kl)),maxval(qfrad(js:je,1:kl))
  write(6,*) "minloc,maxloc ",minloc(qfrad(js:je,1:kl)),maxloc(qfrad(js:je,1:kl))
  call ccmpi_abort(-1) 
end if

if ( any(cfrac(js:je,1:kl)/=cfrac(js:je,1:kl)) ) then
  write(6,*) "ERROR: NaN detected in cfrac on myid=",myid," at ",trim(message)
  call ccmpi_abort(-1)    
end if

if ( any(cfrac(js:je,1:kl)<-1.e-8) .or. any(cfrac(js:je,1:kl)>1.) ) then
  write(6,*) "ERROR: Out-of-range detected in cfrac on myid=",myid," at ",trim(message)
  write(6,*) "minval,maxval ",minval(cfrac(js:je,1:kl)),maxval(cfrac(js:je,1:kl))
  write(6,*) "minloc,maxloc ",minloc(cfrac(js:je,1:kl)),maxloc(cfrac(js:je,1:kl))
  call ccmpi_abort(-1) 
end if

if ( any(rfrac(js:je,1:kl)/=rfrac(js:je,1:kl)) ) then
  write(6,*) "ERROR: NaN detected in rfrac on myid=",myid," at ",trim(message)
  call ccmpi_abort(-1)    
end if

if ( any(rfrac(js:je,1:kl)<-1.e-8) .or. any(rfrac(js:je,1:kl)>1.) ) then
  write(6,*) "ERROR: Out-of-range detected in rfrac on myid=",myid," at ",trim(message)
  write(6,*) "minval,maxval ",minval(rfrac(js:je,1:kl)),maxval(rfrac(js:je,1:kl))
  write(6,*) "minloc,maxloc ",minloc(rfrac(js:je,1:kl)),maxloc(rfrac(js:je,1:kl))
  call ccmpi_abort(-1) 
end if

if ( any(sfrac(js:je,1:kl)/=sfrac(js:je,1:kl)) ) then
  write(6,*) "ERROR: NaN detected in sfrac on myid=",myid," at ",trim(message)
  call ccmpi_abort(-1)    
end if

if ( any(sfrac(js:je,1:kl)<-1.e-8) .or. any(sfrac(js:je,1:kl)>1.) ) then
  write(6,*) "ERROR: Out-of-range detected in sfrac on myid=",myid," at ",trim(message)
  write(6,*) "minval,maxval ",minval(sfrac(js:je,1:kl)),maxval(sfrac(js:je,1:kl))
  write(6,*) "minloc,maxloc ",minloc(sfrac(js:je,1:kl)),maxloc(sfrac(js:je,1:kl))
  call ccmpi_abort(-1) 
end if

if ( any(gfrac(js:je,1:kl)/=gfrac(js:je,1:kl)) ) then
  write(6,*) "ERROR: NaN detected in gfrac on myid=",myid," at ",trim(message)
  call ccmpi_abort(-1)    
end if

if ( any(gfrac(js:je,1:kl)<-1.e-8) .or. any(gfrac(js:je,1:kl)>1.) ) then
  write(6,*) "ERROR: Out-of-range detected in gfrac on myid=",myid," at ",trim(message)
  write(6,*) "minval,maxval ",minval(gfrac(js:je,1:kl)),maxval(gfrac(js:je,1:kl))
  write(6,*) "minloc,maxloc ",minloc(gfrac(js:je,1:kl)),maxloc(gfrac(js:je,1:kl))
  call ccmpi_abort(-1) 
end if

if ( any(psl(js:je)/=psl(js:je)) ) then
  write(6,*) "ERROR: NaN detected in psl on myid=",myid," at ",trim(message)
  call ccmpi_abort(-1)
end if

if ( any(psl(js:je)<-1.4) .or. any(psl(js:je)>0.3) ) then
  write(6,*) "ERROR: Out-of-range detected in psl on myid=",myid," at ",trim(message)
  write(6,*) "minval,maxval ",minval(psl(js:je)),maxval(psl(js:je))
  write(6,*) "minloc,maxloc ",minloc(psl(js:je)),maxloc(psl(js:je))
  call ccmpi_abort(-1) 
end if

if ( any(ps(js:je)/=ps(js:je)) ) then
  write(6,*) "ERROR: NaN detected in ps on myid=",myid," at ",trim(message)
  call ccmpi_abort(-1)
end if

if ( any(tss(js:je)/=tss(js:je)) ) then
  write(6,*) "ERROR: NaN detected in tss on myid=",myid," at ",trim(message)
  call ccmpi_abort(-1)
end if

if ( any(tss(js:je)<75.) .or. any(tss(js:je)>425.) ) then
  write(6,*) "ERROR: Out-of-range detected in tss on myid=",myid," at ",trim(message)
  write(6,*) "minval,maxval ",minval(tss(js:je)),maxval(tss(js:je))
  write(6,*) "minloc,maxloc ",minloc(tss(js:je)),maxloc(tss(js:je))
  call ccmpi_abort(-1) 
end if

if ( abs(iaero)>=2 ) then
  if ( any(xtg(js:je,1:kl,1:naero)/=xtg(js:je,1:kl,1:naero)) ) then
    write(6,*) "ERROR: NaN detected in xtg on myid=",myid," at ",trim(message)
    call ccmpi_abort(-1)
  end if
  if ( any(xtg(js:je,1:kl,1:naero)<-1.e-8) .or. any(xtg(js:je,1:kl,1:naero)>6.5e-5) ) then
    write(6,*) "ERROR: Out-of-range detected in xtg on myid=",myid," at ",trim(message)
    write(6,*) "minval,maxval ",minval(xtg(js:je,1:kl,1:naero)),maxval(xtg(js:je,1:kl,1:naero))
    write(6,*) "minloc,maxloc ",minloc(xtg(js:je,1:kl,1:naero)),maxloc(xtg(js:je,1:kl,1:naero))
    call ccmpi_abort(-1) 
  end if  
  if ( any(ssn(js:je,1:kl,1:2)/=ssn(js:je,1:kl,1:2)) ) then
    write(6,*) "ERROR: NaN detected in ssn on myid=",myid," at ",trim(message)
    call ccmpi_abort(-1)
  end if
  if ( any(ssn(js:je,1:kl,1:2)<-1.e-8) .or. any(ssn(js:je,1:kl,1:2)>6.5e9) ) then
    write(6,*) "ERROR: Out-of-range detected in ssn on myid=",myid," at ",trim(message)
    write(6,*) "minval,maxval ",minval(ssn(js:je,1:kl,1:2)),maxval(ssn(js:je,1:kl,1:2))
    write(6,*) "minloc,maxloc ",minloc(ssn(js:je,1:kl,1:2)),maxloc(ssn(js:je,1:kl,1:2))
    call ccmpi_abort(-1) 
  end if    
end if

if ( any( fg(js:je)/=fg(js:je) ) ) then
  write(6,*) "ERROR: NaN detected in fg on myid=",myid," at ",trim(message)
  call ccmpi_abort(-1)
end if

if ( any( eg(js:je)/=eg(js:je) ) ) then
  write(6,*) "ERROR: NaN detected in eg on myid=",myid," at ",trim(message)
  call ccmpi_abort(-1)
end if

return
end subroutine nantest    
    
!subroutine gabls_flux(xtg,pps,pta,ppa,pu,pv,prhoa,zref,z0, &
!                      psfth,zustar,zvmod)
!
!implicit none
!
!real, parameter :: ZBM=4.8 ! constant used in the formulation
!real, parameter :: ZBH=7.8 ! constant used in the formulation
!real, parameter :: XKARMAN=0.4 ! Von-Karman constant
!real, parameter :: XP00=1.E5 ! reference pressure
!real, parameter :: XCPD=7.*XRD/2.
!real, parameter :: XRD=2./7.*XCPD
!real, parameter :: Z_Z0_O_Z0H = 10.0
!        
!real, dimension(ifull), intent(in) :: xtg, pps, pta, ppa, pu, pv, prhoa
!!        XTG= surface temperature
!!        PPS: surface pressure
!!        PTA: temperature of the 1st level
!!        PPA: pressure of the 1st level
!!        PU: zonal wind of the 1st level
!!        PV: meridional wind of the 1st level
!!        PRHOA= density of the 1st level
! 
!real, dimension(ifull), intent(in) :: zref, z0
!real, dimension(ifull), intent(out) :: psfth, zustar, zvmod
!
!real, dimension(ifull) :: zl
!real, dimension(ifull) :: psfu, psfv
!integer zi
!        
!!ZVMOD= WIND_THRESHOLD(SQRT(PU(:)**2. + PV(:)**2.),ZREF) !routine to avoid a null modulus of wind (cf below)
!ZVMOD= max(SQRT(PU(:)**2 + PV(:)**2),0.1*min(10.,Zref))        
!ZL(:) = -9999.
!
!DO ZI=1,50 ! iteration over 50 times
!        
!  ! computation of u*
!  ZUSTAR(:) = ZVMOD(:)/(LOG(ZREF(:)/Z0)/XKARMAN+ZBM/(ZL(:)*XKARMAN)*(ZREF(:)-Z0))
!        
!  ! computation of T*
!  ZTSTAR(:) = (PTA(:)*((XP00/PPA)**(XRD/XCPD))-XTG(:)*((XP00/PPS)**(XRD/XCPD))) &
!                /(LOG(ZREF(:)/Z0/Z_Z0_O_Z0H)/XKARMAN+ZBH/(ZL(:)*XKARMAN)*       &
!                (ZREF(:)-Z0/Z_Z0_O_Z0H))
!        
!  WHERE(ABS(ZTSTAR)<1.E-10)
!    ZTSTAR(:) = 1.E-10
!  END WHERE
!        
!  WHERE(ABS(ZUSTAR)<1.E-10)
!    ZUSTAR(:) = 1.E-10
!  END WHERE
!        
!  ZL(:)=(ZUSTAR(:)**2)/(XG/(PTA(:)*((XP00/PPA)**(XRD/XCPD)))*XKARMAN*ZTSTAR(:))
!        
!END DO
!        
!PSFU(:)  = -PU(:)/ZVMOD(:)*(ZUSTAR(:)**2)*PRHOA
!PSFV(:)  = -PV(:)/ZVMOD(:)*(ZUSTAR(:)**2)*PRHOA
!PSFTH(:) = -ZTSTAR(:)*ZUSTAR(:)*PRHOA*XCPD
!        
!!        Wind_threshold:
!!        = = = = = =  = =
!!        ;determine a threshold for the wind modulus
!!        _------------------------------------------
!!        Wind_threshold(Pwind,Zref) result(Pwind_new)
!!        Pwind_new=max(Pwind,0.1*min(10.,Zref))
! 
!return
!end subroutine gabls_flux
    
!--------------------------------------------------------------
! INTIAL PARAMETERS
blockdata main_blockdata

implicit none

include 'kuocom.h'           ! Convection parameters

! Vertical mixing options
data ncvmix/0/
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
data acon/.2/,bcon/.07/,rcm/.92e-5/
data rcrit_l/.75/,rcrit_s/.85/ 
! Cloud options
data ldr/1/,nclddia/1/,nstab_cld/0/,nrhcrit/10/,sigcll/.95/ 
data cldh_lnd/95./,cldm_lnd/85./,cldl_lnd/75./
data cldh_sea/95./,cldm_sea/90./,cldl_sea/80./
data ncloud/0/

end

