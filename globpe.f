      program globpe

!      PE model on conformal-cubic grid
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!      input files are :namelist (via file called "input")
!                       "nrun.dat"
!      data input and output file names are specified in namelist 'datafile'
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!     sign convention:
!                      u+ve eastwards  (on the panel)
!                      v+ve northwards (on the panel)

      use aerointerface                       ! Aerosol interface
      use arrays_m                            ! Atmosphere dyamics prognostic arrays
      use bigxy4_m                            ! Grid interpolation
      use carbpools_m, only : carbpools_init  ! Carbon pools
     &    ,fpn,frs,frp
      use cc_mpi                              ! CC MPI routines
      use cfrac_m                             ! Cloud fraction
      use dava_m                              ! Far-field nudging (weights)
      use davb_m                              ! Far-field nudging (host store)
      use diag_m                              ! Diagnostic routines
      use dpsdt_m                             ! Vertical velocity
      use epst_m                              ! Off-centre terms
      use extraout_m                          ! Additional diagnostics
      use gdrag_m, only : gdrag_init          ! Gravity wave drag
      use histave_m                           ! Time average arrays
      use indices_m                           ! Grid index arrays
      use infile                              ! Input file routines
      use kuocomb_m                           ! JLM convection
      use latlong_m                           ! Lat/lon coordinates
      use liqwpar_m                           ! Cloud water mixing ratios
      use map_m                               ! Grid map arrays
      use mlo, only : mlodiag,wlev,mxd        ! Ocean physics and prognostic arrays
     &   ,mindep,minwater
      use morepbl_m                           ! Additional boundary layer diagnostics
      use nharrs_m, only : nharrs_init        ! Non-hydrostatic atmosphere arrays
     &   ,lrestart
      use nlin_m                              ! Atmosphere non-linear dynamics
      use nsibd_m                             ! Land-surface arrays
      use parmhdff_m                          ! Horizontal diffusion parameters
      use pbl_m                               ! Boundary layer arrays
      use permsurf_m, only : permsurf_init    ! Fixed surface arrays
      use prec_m                              ! Precipitation
      use raddiag_m                           ! Radiation diagnostic
      use savuvt_m                            ! Saved dynamic arrays
      use savuv1_m                            ! Saved dynamic arrays
      use sbar_m                              ! Saved dynamic arrays
      use screen_m                            ! Screen level diagnostics
      use seaesfrad_m                         ! SEA-ESF radiation
      use sigs_m                              ! Atmosphere sigma levels
      use soil_m                              ! Soil and surface data
      use soilsnow_m                          ! Soil, snow and surface data
      use tbar2d_m, only : tbar2d_init        ! Atmosphere dynamics reference temperature
      use timeseries, only : write_ts         ! Tracer time series
      use tkeeps, only : tkeinit              ! TKE-EPS boundary layer
      use tracermodule, only :init_tracer     ! Tracer routines
     &   ,trfiles,tracer_mass,unit_trout
     &   ,interp_tracerflux,tracerlist
      use tracers_m                           ! Tracer data
      use unn_m                               ! Saved dynamic arrays
      use uvbar_m                             ! Saved dynamic arrays
      use vecs_m, only : vecs_init            ! Eigenvectors for atmosphere dynamics
      use vecsuv_m                            ! Map to cartesian coordinates
      use vegpar_m                            ! Vegetation arrays
      use vvel_m                              ! Additional vertical velocity
      use work2_m                             ! Diagnostic arrays
      use work3_m                             ! Mk3 land-surface diagnostic arrays
      use work3f_m                            ! Grid work arrays
      use work3sav_m                          ! Water and tracer saved arrays
      use workglob_m                          ! Additional grid interpolation
      use xarrs_m                             ! Saved dynamic arrays
      use xyzinfo_m                           ! Grid coordinate arrays

      implicit none

      include 'newmpar.h'                     ! Grid parameters
      include 'const_phys.h'                  ! Physical constants
      include 'darcdf.h'                      ! Netcdf data
      include 'dates.h'                       ! Date data
      include 'establ.h'                      ! Liquid saturation function
      include 'filnames.h'                    ! Filenames
      include 'kuocom.h'                      ! Convection parameters
      include 'mpif.h'                        ! MPI parameters
      include 'netcdf.inc'                    ! Netcdf parameters
      include 'parm.h'                        ! Model configuration
      include 'parmdyn.h'                     ! Dynamics parameters
      include 'parmgeom.h'                    ! Coordinate data
      include 'parmhor.h'                     ! Horizontal advection parameters
      include 'parmsurf.h'                    ! Surface parameters
      include 'parmvert.h'                    ! Vertical advection parameters
      include 'soilv.h'                       ! Soil parameters
      include 'stime.h'                       ! File date data
      include 'trcom2.h'                      ! Station data
      include 'version.h'                     ! Model version data
      
      integer leap
      common/leap_yr/leap                     ! Leap year (1 to allow leap years)
      integer nbarewet,nsigmf
      common/nsib/nbarewet,nsigmf             ! Land-surface options
      integer nnrad,idcld
      common/radnml/nnrad,idcld               ! Radiation options

      integer, dimension(8) :: tvals1, tvals2
      integer, dimension(8) :: nper3hr
      integer iaero, ier, igas, ilx, io_nest, iq, irest, isoil, itr1
      integer itr2, jalbfix, jlx, jtr1, jtr2, k,k2, kktau, mexrest
      integer mins_dt, mins_gmt, mspeca, mtimer_in, nalpha, newsnow
      integer ng, nlx, nmaxprsav, nmi, npa, npb, npc, n3hr, nsnowout
      integer nstagin, nstaguin, nwrite, nwtsav, mins_rad, mtimer_sav
      integer nn, i, j, mstn, nproc_in, ierr, nperhr, nscrn, nversion
      integer ierr2, kmax, isoth, nsig, lapsbot
      real, dimension(:,:), allocatable, save :: speed
      real, dimension(:), allocatable, save :: spare1,spare2
      real, dimension(:), allocatable, save :: spmean,div
      real, dimension(9) :: temparray, gtemparray
      real clhav, cllav, clmav, cltav, con, div_int, dsx, dtds, es
      real gke, hourst, hrs_dt, evapavge, precavge, preccavge, psavge
      real pslavge, pwater, rel_lat, rel_long, rlwup, spavge, pwatr
      real pwatr_l, qtot, aa, bb, cc, bb_2, cc_2, rat
      character(len=60) comm,comment
      character(len=47) header
      character(len=10) timeval
      character(len=8) rundate
      character(len=2) text
      logical odcalc

      namelist/defaults/nversion
      namelist/cardin/comment,dt,ntau,nwt,npa,npb,npc,nhorps,nperavg
     & ,ia,ib,ja,jb,itr1,jtr1,itr2,jtr2,id,jd,iaero,khdif,khor,nhorjlm
     & ,m,mex,mbd,nbd,ndi,ndi2,nem,nhor,nlv,nscrn
     & ,nmaxpr,nmi,nonl,nrot,nrad,ntaft,ntsea
     & ,ntsur,ntvdr,nvad,nvadh,nvmix,nxmap
     & ,restol,precon,kdate_s,ktime_s,leap,newtop,idcld,mup
     & ,lgwd,ngwd,rhsat
     & ,nextout,hdifmax,jalbfix
     & ,nalpha
     & ,nstag,nstagu,ntbar,nwrite
     & ,irest,nrun,nstn,rel_lat,rel_long,nrungcm,nsib
     & ,istn,jstn,iunp,nrotstn,slat,slon,zstn,name_stn
     & ,mexrest,mh_bs,ndept,nritch_t,nt_adv
     & ,mfix,mfix_qg,namip,amipo3,nh,npex,nhstest,nsemble
     & ,nspecial,panfg,panzo,nplens
     & ,newsnow,nsnowout,newrough,newsoilm,nglacier,newztsea
     & ,epsp,epsu,epsf
     & ,av_vmod,charnock,chn10,snmin,tss_sh,vmodmin,zobgin
     & ,rlong0,rlat0,schmidt
     & ,kbotdav,kbotu,nbox,nud_p,nud_q,nud_t,nud_uv,nud_hrs,nudu_hrs
     & ,nlocal,nvsplit,nbarewet,nsigmf,qgmin
     & ,io_clim ,io_in,io_nest,io_out,io_rest,io_spec,localhist   
     & ,m_fly,mstn,nqg,nurban,nmr,ktopdav,nud_sst,nud_sss
     & ,mfix_tr,mfix_ke,mfix_aero,kbotmlo,ktopmlo,mloalpha,nud_ouv
     & ,nud_sfh,bpyear,rescrn,helmmeth,nmlo,ol,mxd,mindep,minwater
      namelist/skyin/mins_rad,ndiur
      namelist/datafile/ifile,ofile,albfile,co2emfile,eigenv,
     &    hfile,icefile,mesonest,nmifile,o3file,radfile,restfile,
     &    rsmfile,scamfile,scrnfile,snowfile,so4tfile,soilfile,sstfile,
     &    surfile,tmaxfile,tminfile,topofile,trcfil,vegfile,zofile,
     &    smoistfile,soil2file,radonemfile,
     &    co2_00,radon_00,surf_00,co2_12,radon_12,surf_12,
     &    laifile,albnirfile,urbanfile,bathfile,vegprev,vegnext,
     &    cnsdir,salfile,oxidantfile
      namelist/kuonml/alflnd,alfsea
     &        ,cldh_lnd,cldm_lnd,cldl_lnd
     &        ,cldh_sea,cldm_sea,cldl_sea
     &        ,convfact,convtime,shaltime
     &        ,detrain,detrainx,dsig2,dsig4,entrain,fldown,iterconv
     &        ,ksc,kscmom,kscsea,ldr,mbase,mdelay,methdetr,methprec
     &        ,nbase,nclddia,ncvcloud
     &        ,ncvmix,nevapcc,nevapls,nkuo,nrhcrit,nstab_cld
     &        ,nuvconv,rhcv,rhmois,rhsat
     &        ,sigcb,sigcll,sig_ct,sigkscb,sigksct
     &        ,tied_con,tied_over,tied_rh,comm
     &        ,acon,bcon,rcm,rcrit_l,rcrit_s

      data nscrn/0/,nversion/0/,lapsbot/0/
      data npc/40/,nmi/0/,io_nest/1/,iaero/0/,newsnow/0/      
      data itr1/23/,jtr1/13/,itr2/25/,jtr2/11/
      data comment/' '/,comm/' '/,irest/1/,jalbfix/1/,nalpha/1/
      data mexrest/6/,mins_rad/60/,nwrite/0/,nsnowout/999999/

!#ifdef stacklimit
      ! For linux only
      call setstacklimit(-1)
!#endif

      !--------------------------------------------------------------
      ! INITALISE MPI ROUTINES
      call MPI_Init(ierr)
      call MPI_Comm_size(MPI_COMM_WORLD, nproc, ierr) ! Find number of processes
      call MPI_Comm_rank(MPI_COMM_WORLD, myid, ierr)  ! Find local processor id


      !--------------------------------------------------------------
      ! INITALISE LOGS
      call log_off()
      call log_setup()
#ifdef simple_timer
      call start_log(model_begin)
#endif


      !--------------------------------------------------------------
      ! READ NAMELISTS AND SET PARAMETER DEFAULTS
      ia=-1   ! diagnostic index
      ib=-1   ! diagnostic index
      ntbar=-1
      rel_lat=0.
      rel_long=0.
      ktau=0
      ol=20
      call initialparm

      ! All processors read the namelist, so no MPI comms are needed
      open(99,file="input",form="formatted",status="old")
      read (99, defaults)
      if (myid==0) then
        write(6,'(a10," running for nproc =",i7)')
     &                      version,nproc 
        write(6,*) 'Using defaults for nversion = ',nversion
      end if
      if(nversion.ne.0)call change_defaults(nversion)
      read (99, cardin)
      npc     =max(npc,1)
      nperday =nint(24.*3600./dt)
      nperhr  =nint(3600./dt)
      do n3hr=1,8
       nper3hr(n3hr)=nint(n3hr*3*3600/dt)
      enddo
      if (nwt==-99) nwt=nperday          ! set default nwt to 24 hours
      if (nperavg==-99) nperavg=nwt      ! set default nperavg to nwt
      if (nwrite==0) nwrite=nperday      ! only used for outfile IEEE
      if (nmlo.ne.0.and.abs(nmlo).le.9) then
        ol=max(ol,1)
      else
        ol=0
      end if
      wlev=ol
      mindep=max(0.,mindep)
      minwater=max(0.,minwater)
      read (99, skyin)
      kountr=nint(mins_rad*60./dt)  ! set default radiation to ~mins_rad m
      mins_rad=nint(kountr*dt/60.)  ! redefine to actual value
      read (99, datafile)
      read (99, kuonml)
      ngas=0
      read(99, trfiles, iostat=ierr)         ! try reading tracer namelist.  If no
      if (ierr.ne.0) rewind(99)              ! namelist is found, then disable
      if (tracerlist.ne.'') call init_tracer ! tracers and rewind namelist.

      !--------------------------------------------------------------
      ! READ TOPOGRAPHY FILE TO DEFINE CONFORMAL CUBIC GRID
      il_g=48
      rlong0=0.
      rlat0=0.
      schmidt=1.
      if (myid==0.and.io_in<=4) then
!       open new topo file and check its dimensions
!       here used to supply rlong0,rlat0,schmidt
!       Remanded of file is read in indata.f
        open(66,file=topofile,recl=2000,status='old')
        write(6,*) 'reading topofile header'
        read(66,*) ilx,jlx,rlong0,rlat0,schmidt,dsx,header
        write(6,*) 'ilx,jlx,rlong0,rlat0,schmidt,dsx ',
     &           ilx,jlx,rlong0,rlat0,schmidt,dsx,header
        il_g=ilx
      end if      ! (io_in<=4)
      call MPI_Bcast(il_g,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(rlong0,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(rlat0,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(schmidt,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)


      !--------------------------------------------------------------
      ! READ EIGENV FILE TO DEFINE VERTICAL LEVELS
      if (myid==0) then
!       Remanded of file is read in indata.f
        open(28,file=eigenv,status='old',form='formatted')
        read(28,*)kmax,lapsbot,isoth,nsig
        kl=kmax
        write(6,*)'kl,ol,lapsbot,isoth,nsig: ',
     &             kl,ol,lapsbot,isoth,nsig
      end if
      call MPI_Bcast(kl,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(lapsbot,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(isoth,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(lapsbot,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(nsig,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)


      !--------------------------------------------------------------
      ! DEFINE newmpar VARIABLES AND DEFAULTS
      jl_g=il_g+npanels*il_g      
#ifdef uniform_decomp
      if (myid==0) then
        write(6,*) "Using uniform grid decomposition"
      end if
      nxp=nint(sqrt(real(nproc)))
      nyp=nproc/nxp
      do while((mod(il_g,max(nxp,1)).ne.0.or.
     &         mod(nproc,max(nxp,1)).ne.0.or.
     &         mod(il_g,nyp).ne.0).and.nxp.gt.0)
        nxp=nxp-1
        nyp=nproc/max(nxp,1)
      end do
#else
      if (myid==0) then
        write(6,*) "Using face grid decomposition"
      end if
      if (mod(nproc,6).ne.0.and.mod(6,nproc).ne.0) then
        write(6,*) "ERROR: nproc must be a multiple of 6 or"
        write(6,*) "a factor of 6"
        stop
      end if
      nxp=max(1,nint(sqrt(real(nproc)/6.)))
      nyp=nproc/nxp
      do while((mod(il_g,max(nxp,1)).ne.0.or.
     &         mod(nproc/6,max(nxp,1)).ne.0.or.
     &         mod(jl_g,nyp).ne.0).and.nxp.gt.0)
        nxp=nxp-1
        nyp=nproc/max(nxp,1)
      end do
#endif
      if (nxp.eq.0) then
        write(6,*) "ERROR: Invalid number of processors for this grid"
        write(6,*) "Try increasing or decreasing nproc"
        stop
      end if
      ifull_g=il_g*jl_g
      ijk_g=il_g*jl_g*kl
      iquad=1+il_g*((8*npanels)/(npanels+4))
      il=il_g/nxp
      jl=jl_g/nyp
      ifull=il*jl
      ijk=il*jl*kl
!     The perimeter of the processor region has length 2*(il+jl).
!     The first row has 8 possible corner points per panel and the 
!     second has 16. In practice these are not all distinct so there could
!     be some optimisation.
#ifdef uniform_decomp
      npan=npanels+1
!     This should use jpan rather than jl. Will be far too big.
      iextra=(6*(il+jl)+24)*npan ! the *6 term is for extra vector row (e.g., inu,isu,iev,iwv)
#else      
      npan=max(1,(npanels+1)/nproc)
      iextra=6*(il+jl)+24*npan ! the *6 term is for extra vector row (e.g., inu,isu,iev,iwv)
#endif
      ! nrows_rad is a subgrid step for radiation routines
      nrows_rad=il_g/6
      do while(mod(jl,nrows_rad).ne.0)
        nrows_rad=nrows_rad-1
      end do
      if (myid.eq.0) then
        write(6,*) "il_g,jl_g,il,jl ",il_g,jl_g,il,jl
        write(6,*) "nxp,nyp         ",nxp,nyp
        write(6,*) "nrows_rad     = ",nrows_rad
      end if
      
      ! some default values for unspecified parameters
      if (ia.lt.0) ia=il/2
      if (ib.lt.0) ib=ia+3
      if (ntbar.lt.0) ntbar=kl/3
      if (ktopdav.lt.0) ktopdav=kl
      if (kbotmlo.lt.0) kbotmlo=ol
      if (ldr==0) mbase=0
      dsig4=max(dsig2+.01,dsig4)
      if(kbotu==0) kbotu=kbotdav
      if(mbd.ne.0.and.nbd.ne.0)then
        write(6,*) 'setting nbd=0 because mbd.ne.0'
        nbd=0
      endif
      nud_hrs=abs(nud_hrs)  ! just for people with old -ves in namelist
      if (nudu_hrs==0) nudu_hrs=nud_hrs
!     for 6-hourly output of sint_ave etc, want 6*60 = N*mins_rad      
      if ((nrad==4.or.nrad==5).and.mod(6*60,mins_rad).ne.0)
     &  stop 'prefer 6*60 = N*mins_rad '


      ! **** do namelist fixes above this ***
      

      !--------------------------------------------------------------
      ! DISPLAY NAMELIST
      if ( myid == 0 ) then   
        write(6,*)'Dynamics options A:'
        write(6,*)'   m    mex   mfix  mfix_qg   mup    nh    nonl',    
     &            '   npex  precon' 
        write(6,'(i5,9i7)')m,mex,mfix,mfix_qg,mup,nh,nonl,npex,precon
        write(6,*)'Dynamics options B:'
        write(6,*)'nritch_t nrot  ntbar  nxmap   epsp    epsu   epsf',
     &            '   restol'
        write (6,'(i5,3i7,1x,3f8.3,g9.2)')nritch_t,nrot,ntbar,nxmap,
     &                                    epsp,epsu,epsf,restol
        write(6,*)'Horizontal advection/interpolation options:'
        write(6,*)' ndept  nt_adv mh_bs  mhint '
        write (6,'(i5,11i7)') ndept,nt_adv,mh_bs,mhint
        write(6,*)'Horizontal wind staggering options:'
        write(6,*)'mstagpt nstag nstagu'
        write (6,'(i5,11i7)') mstagpt,nstag,nstagu
        write(6,*)'Vertical advection options:'
        write(6,*)'  nvad  nvadh  '
        write (6,'(i5,11i7)') nvad,nvadh
        if(nvad==4.or.nvad==-4)then
          write(6,*)'Vertical advection options for TVD:'
          write(6,*)' nimp   nthub  ntvd   ntvdr'
          write (6,'(i5,11i7)') nimp,nthub,ntvd,ntvdr
        endif
        write(6,*)'Horizontal mixing options:'
        write(6,*)' khdif  khor   nhor   nhorps nhorjlm'
        write (6,'(i5,11i7)') khdif,khor,nhor,nhorps,nhorjlm
        write(6,*)'Vertical mixing/physics options:'
        write(6,*)' nvmix nlocal nvsplit ncvmix lgwd   ngwd   '
        write (6,'(i5,6i7)') nvmix,nlocal,nvsplit,ncvmix,lgwd,ngwd
        write(6,*)'Cumulus convection options A:'
        write(6,*)' nkuo  sigcb sig_ct  rhcv  rhmois rhsat',
     &            ' convfact convtime shaltime'
        write (6,'(i5,6f7.2,9f8.2)')
     &    nkuo,sigcb,sig_ct,rhcv,rhmois,rhsat,
     &    convfact,convtime,shaltime
        write(6,*)'Cumulus convection options B:'
        write(6,*)' alflnd alfsea fldown iterconv',
     &            ' ncvcloud nevapcc nevapls nuvconv'
        write (6,'(3f7.2,i6,i10,4i8)') alflnd,alfsea,fldown,iterconv,
     &                              ncvcloud,nevapcc,nevapls,nuvconv
        write(6,*)'Cumulus convection options C:'
        write(6,*)' mbase mdelay methprec nbase detrain',
     &            ' entrain methdetr detrainx dsig2  dsig4'
        write (6,'(3i6,i9,f8.2,f9.2,i8,4f8.2)') mbase,mdelay,methprec,
     &              nbase,detrain,entrain,methdetr,detrainx,dsig2,dsig4
        write(6,*)'Shallow convection options:'
        write(6,*)'  ksc  kscsea kscmom sigkscb sigksct ',
     &            'tied_con tied_over tied_rh '
        write (6,'(i5,2i7,1x,3f8.3,2f10.3)') ksc,kscsea,kscmom,
     &              sigkscb,sigksct,tied_con,tied_over,tied_rh
        write(6,*)'Other moist physics options:'
        write(6,*)'  acon   bcon   qgmin      rcm    rcrit_l rcrit_s'
        write (6,'(2f7.2,2e10.2,2f7.2)') acon,bcon,qgmin,rcm,
     &                                   rcrit_l,rcrit_s
        write(6,*)'Radiation options:'
        write(6,*)' nrad  ndiur mins_rad kountr iaero  dt'
        write (6,'(i5,4i7,f10.2)') nrad,ndiur,mins_rad,kountr,iaero,dt
        write(6,*)'Diagnostic cloud options:'
        write(6,*)'  ldr nclddia nstab_cld nrhcrit sigcll '
        write (6,'(i5,i6,2i9,1x,f8.2)') ldr,nclddia,nstab_cld,nrhcrit,
     &                                  sigcll
        write(6,*)'Soil, canopy and PBL options A:'
        write(6,*)' jalbfix nalpha nbarewet newrough nglacier nrungcm',
     &            ' nsib  nsigmf'
        write (6,'(i5,9i8)')
     &          jalbfix,nalpha,nbarewet,newrough,nglacier,nrungcm,nsib,
     &          nsigmf
        write(6,*)'Soil, canopy and PBL options B:'
        write(6,*)' ntaft ntsea ntsur av_vmod tss_sh vmodmin  zobgin',
     &            ' charnock chn10'
        write (6,'(i5,2i6,4f8.2,f8.3,f9.5)') ntaft,ntsea,ntsur,
     &            av_vmod,tss_sh,vmodmin,zobgin,charnock,chn10    
        if(mbd.ne.0.or.nbd.ne.0)then
          write(6,*)'Nudging options:'
          write(6,*)' nbd    nud_p  nud_q  nud_t  nud_uv nud_hrs',
     &              ' nudu_hrs kbotdav  kbotu'
          write (6,'(i5,3i7,7i8)') 
     &      nbd,nud_p,nud_q,nud_t,nud_uv,nud_hrs,nudu_hrs,kbotdav,kbotu
        endif
        write(6,*)'Special and test options:'
        write(6,*)' mbd namip amipo3 newtop nhstest nplens nsemble',
     &            ' nspecial panfg panzo'
        write (6,'(2i5,L7,4i7,i8,f9.1,f8.4)') mbd,namip,amipo3,
     &          newtop,nhstest,nplens,nsemble,nspecial,panfg,panzo
        write(6,*)'I/O options:'
        write(6,*)' m_fly  io_in io_nest io_out io_rest  nwt',
     &            '  nperavg   nscrn'
        write (6,'(i5,4i7,3i8)') 
     &        m_fly,io_in,io_nest,io_out,io_rest,nwt,nperavg,nscrn
        if(ntrac.ne.0)then
          write(6,*)'Trace gas options:'
          write(6,*)' ngas   nllp   ntrac'
          write (6,'(i5,3i7)') ngas,nllp,ntrac
        endif

        write (6, cardin)
        if(nllp==0.and.nextout>=4)stop 'need nllp=3 for nextout>=4'
        write (6, skyin)
        write (6, datafile)
        write(6, kuonml)
      end if ! myid=0
      if (newtop>2) stop 'newtop>2 no longer allowed'
      if (mfix_qg>0.and.(nkuo==4.or.nvad==44))
     &        stop 'nkuo=4,nvad=44: mfix_qg>0 not allowed'
      if (mfix>3) stop 'mfix >3 not allowed now'


      !--------------------------------------------------------------
      ! INITIALISE ifull_g ALLOCATABLE ARRAYS
      call bigxy4_init(iquad)
      call xyzinfo_init(ifull_g,ifull,iextra)
      call indices_init(ifull_g,ifull,iextra,npanels,npan)
      call map_init(ifull_g,ifull,iextra)
      call latlong_init(ifull_g,ifull,iextra)      
      call vecsuv_init(ifull_g,ifull,iextra)


      !--------------------------------------------------------------
      ! SET UP CC GEOMETRY
      ! Only one processor calls setxyz to save memory with large grids
      if (myid==0) then
        call workglob_init(ifull_g)
        call setxyz(il_g,rlong0,rlat0,schmidt,x_g,y_g,z_g,wts_g,ax_g,
     &     ay_g,az_g,bx_g,by_g,bz_g,xx4,yy4,myid)
      end if
      ! Broadcast the following global arrays so that they can be
      ! decomposed into local arrays with ccmpi_setup
      call MPI_Bcast(ds,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(xx4,iquad*iquad,MPI_DOUBLE_PRECISION,0,
     &               MPI_COMM_WORLD,ierr)
      call MPI_Bcast(yy4,iquad*iquad,MPI_DOUBLE_PRECISION,0,
     &               MPI_COMM_WORLD,ierr)
      call MPI_Bcast(x_g,ifull_g,MPI_DOUBLE_PRECISION,0,
     &               MPI_COMM_WORLD,ierr)
      call MPI_Bcast(y_g,ifull_g,MPI_DOUBLE_PRECISION,0,
     &               MPI_COMM_WORLD,ierr)
      call MPI_Bcast(z_g,ifull_g,MPI_DOUBLE_PRECISION,0,
     &               MPI_COMM_WORLD,ierr)
      call MPI_Bcast(iw_g,ifull_g,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(is_g,ifull_g,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(ise_g,ifull_g,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(ie_g,ifull_g,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(ine_g,ifull_g,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(in_g,ifull_g,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(iwn_g,ifull_g,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(ien_g,ifull_g,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(inw_g,ifull_g,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(isw_g,ifull_g,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(ies_g,ifull_g,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(iws_g,ifull_g,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(inn_g,ifull_g,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(iss_g,ifull_g,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(iww_g,ifull_g,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(iee_g,ifull_g,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(iwu_g,ifull_g,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(isv_g,ifull_g,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(iwu2_g,ifull_g,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(isv2_g,ifull_g,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(ieu2_g,ifull_g,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(inv2_g,ifull_g,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(iev2_g,ifull_g,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(inu2_g,ifull_g,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(ieu_g,ifull_g,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(inv_g,ifull_g,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(iwwu2_g,ifull_g,MPI_INTEGER,0,MPI_COMM_WORLD,
     &               ierr)
      call MPI_Bcast(issv2_g,ifull_g,MPI_INTEGER,0,MPI_COMM_WORLD,
     &               ierr)
      call MPI_Bcast(ieeu2_g,ifull_g,MPI_INTEGER,0,MPI_COMM_WORLD,
     &               ierr)
      call MPI_Bcast(innv2_g,ifull_g,MPI_INTEGER,0,MPI_COMM_WORLD,
     &               ierr)
      call MPI_Bcast(lwws_g,npanels+1,MPI_INTEGER,0,MPI_COMM_WORLD,
     &               ierr)
      call MPI_Bcast(lwss_g,npanels+1,MPI_INTEGER,0,MPI_COMM_WORLD,
     &               ierr)
      call MPI_Bcast(lees_g,npanels+1,MPI_INTEGER,0,MPI_COMM_WORLD,
     &               ierr)
      call MPI_Bcast(less_g,npanels+1,MPI_INTEGER,0,MPI_COMM_WORLD,
     &               ierr)
      call MPI_Bcast(lwwn_g,npanels+1,MPI_INTEGER,0,MPI_COMM_WORLD,
     &               ierr)
      call MPI_Bcast(lwnn_g,npanels+1,MPI_INTEGER,0,MPI_COMM_WORLD,
     &               ierr)
      call MPI_Bcast(leen_g,npanels+1,MPI_INTEGER,0,MPI_COMM_WORLD,
     &               ierr)
      call MPI_Bcast(lenn_g,npanels+1,MPI_INTEGER,0,MPI_COMM_WORLD,
     &               ierr)
      call MPI_Bcast(lsww_g,npanels+1,MPI_INTEGER,0,MPI_COMM_WORLD,
     &               ierr)
      call MPI_Bcast(lssw_g,npanels+1,MPI_INTEGER,0,MPI_COMM_WORLD,
     &               ierr)
      call MPI_Bcast(lsee_g,npanels+1,MPI_INTEGER,0,MPI_COMM_WORLD,
     &               ierr)
      call MPI_Bcast(lsse_g,npanels+1,MPI_INTEGER,0,MPI_COMM_WORLD,
     &               ierr)
      call MPI_Bcast(lnww_g,npanels+1,MPI_INTEGER,0,MPI_COMM_WORLD,
     &               ierr)
      call MPI_Bcast(lnnw_g,npanels+1,MPI_INTEGER,0,MPI_COMM_WORLD,
     &               ierr)
      call MPI_Bcast(lnee_g,npanels+1,MPI_INTEGER,0,MPI_COMM_WORLD,
     &               ierr)
      call MPI_Bcast(lnne_g,npanels+1,MPI_INTEGER,0,MPI_COMM_WORLD,
     &               ierr)
      call MPI_Bcast(npann_g,14,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(npane_g,14,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(npanw_g,14,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(npans_g,14,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(em_g,ifull_g,MPI_REAL,0,MPI_COMM_WORLD,ierr)
      call ccmpi_setup()

      
      !--------------------------------------------------------------
      ! DEALLOCATE ifull_g ARRAYS WHERE POSSIBLE
      if (myid.ne.0) then
        deallocate(wts_g)
        deallocate(ax_g,bx_g,ay_g,by_g,az_g,bz_g)
        deallocate(emu_g,emv_g,f_g,fu_g,fv_g,dmdx_g,dmdy_g)
        deallocate(rlatt_g,rlongg_g)
        deallocate(iw_g,is_g,ise_g,ie_g,ine_g,in_g,iwn_g,ien_g)
        deallocate(inw_g,isw_g,ies_g,iws_g)
        deallocate(inn_g,iss_g,iww_g,iee_g,iwu_g,isv_g)
        deallocate(iwu2_g,isv2_g,ieu2_g,inv2_g,iev2_g,inu2_g,ieu_g)
        deallocate(inv_g,iwwu2_g,issv2_g,ieeu2_g,innv2_g)
        deallocate(lwws_g,lwss_g,lees_g,less_g,lwwn_g)
        deallocate(lwnn_g,leen_g,lenn_g,lsww_g)
        deallocate(lssw_g,lsee_g,lsse_g,lnww_g,lnnw_g)
        deallocate(lnee_g,lnne_g)
        deallocate(npann_g,npane_g,npanw_g,npans_g)
      end if

      
      !--------------------------------------------------------------
      ! INITIALISE LOCAL ARRAYS
      allocate(spare1(ifull),spare2(ifull))
      allocate(speed(ifull,kl))
      allocate(spmean(kl),div(kl))
      call arrays_init(ifull,iextra,kl)
      call carbpools_init(ifull,iextra,kl,nsib)
      call cfrac_init(ifull,iextra,kl)
      call dpsdt_init(ifull,iextra,kl)
      call epst_init(ifull,iextra,kl)
      call extraout_init(ifull,iextra,kl)
      call gdrag_init(ifull,iextra,kl)
      call histave_init(ifull,iextra,kl,ms)
      call kuocomb_init(ifull,iextra,kl)
      call liqwpar_init(ifull,iextra,kl)
      call nlin_init(ifull,iextra,kl)
      call morepbl_init(ifull,iextra,kl)
      call nharrs_init(ifull,iextra,kl)
      call nsibd_init(ifull,iextra,kl,nsib,nurban)
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
      call soil_init(ifull,iextra,kl,nrad,nsib)
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
      if (nvmix.eq.6) then
        call tkeinit(ifull,iextra,kl,0)
      end if
      if (ngas.gt.0) then
        call tracers_init(il,jl,kl,iextra)
      end if
      call work3sav_init(ifull,iextra,kl,ilt,jlt,klt,ngasmax) ! must occur after tracers_init
      if (nbd.ne.0.and.nud_hrs.ne.0) then
        call dava_init(ifull,iextra,kl)
        call davb_init(ifull,iextra,kl)
      end if
      ! Remaining arrays are allocated in indata.f, since their
      ! definition requires additional input data (e.g, land-surface)

      
      !--------------------------------------------------------------
      ! DISPLAY DIAGNOSTIC INDEX AND TIMER DATA
      con=180./pi
      if ( mydiag ) then
         write(6,*)'id,jd,rlongg,rlatt in degrees: ',
     &         id,jd,con*rlongg(idjd),con*rlatt(idjd)
      end if

      if (myid == 0) then
        call date_and_time(rundate)
        write(6,*)'RUNDATE IS ',rundate
        call date_and_time(time=timeval)
        write(6,*)'Starting time ',timeval
      end if


      !--------------------------------------------------------------
      ! READ INITIAL CONDITIONS
      if (myid==0) then
        ncid = ncopn ( ifile,0,ier ) ! 0 denotes read-only
        write(6,*)'ncid,ier,ifile ',ncid,ier,ifile
        write(6,*)'calling indata; will read from file ',ifile
      end if
      call indata(hourst,newsnow,jalbfix,iaero,lapsbot,isoth,nsig)
      ! do not close ncid as onthefly.f expects each file to have
      ! separate ncid numbers.

      ! max/min diagnostics      
      call maxmin(u,' u',ktau,1.,kl)
      call maxmin(v,' v',ktau,1.,kl)
      speed(:,:)=sqrt(u(1:ifull,:)**2+v(1:ifull,:)**2)  ! 3D 
      call maxmin(speed,'sp',ktau,1.,kl)
      call maxmin(t,' t',ktau,1.,kl)
      call maxmin(qg,'qg',ktau,1.e3,kl)
      call maxmin(qfg,'qf',ktau,1.e3,kl)
      call maxmin(qlg,'ql',ktau,1.e3,kl)
      call maxmin(wb,'wb',ktau,1.,ms)
      call maxmin(tggsn,'tggsn',ktau,1.,3)
      call maxmin(tgg,'tgg',ktau,1.,ms)
      pwatr_l=0.   ! in mm
      do k=1,kl
        do iq=1,ifull
         qtot=qg(iq,k)+qlg(iq,k)+qfg(iq,k)
         pwatr_l=pwatr_l-dsig(k)*wts(iq)*qtot*ps(iq)
        enddo
      enddo
      pwatr_l=pwatr_l/grav
      call MPI_Reduce ( pwatr_l, pwatr, 1, MPI_REAL, MPI_SUM, 0,
     &                  MPI_COMM_WORLD, ierr )
      if (myid==0) write (6,"('pwatr0 ',12f7.3)") pwatr
      if (nextout>=4) call setllp
      if (ntrac>0) then
        do ng=1,ntrac
          write (text,'("g",i1)')ng
          call maxmin(tr(:,:,ng),text,ktau,1.,kl)
        end do
      end if   ! (ntrac>0)


      !--------------------------------------------------------------
      ! OPEN MESONEST FILE
      if(mbd.ne.0.or.nbd.ne.0)then
         io_in=io_nest ! Needs to be seen by all processors
         if ( myid == 0 ) then
           ncid = ncopn(mesonest,0,ier )  ! 0 denotes read-only
           write(6,*)'ncid,ier,mesonest ',ncid,ier,mesonest
           call ncmsg("Reading infile",ier)
         endif ! myid == 0
      endif    ! (mbd.ne.0.or.nbd.ne.0)


      !--------------------------------------------------------------
      ! SETUP REMAINING PARAMETERS

      ! convection
!     sig(kuocb) occurs for level just BELOW sigcb
      kuocb=1
      do while(sig(kuocb+1)>=sigcb)
       kuocb=kuocb+1
      enddo
      if (myid==0)
     & write(6,*)'convective cumulus scheme: kuocb,sigcb = ',kuocb,sigcb

      ! horizontal diffusion 
      if(khdif==-99)then   ! set default khdif appropriate to resolution
        khdif=5
        if(myid==0)write(6,*)'Model has chosen khdif =',khdif
      endif
      do k=1,kl
       hdiff(k)=khdif*.1
      enddo
      if(khor>0)then
        do k=kl+1-khor,kl
         hdiff(k)=2.*hdiff(k-1)
        enddo
      elseif(khor<0)then
        do k=1,kl                    ! N.B. usually hdiff(k)=khdif*.1 
!!       increase hdiff between sigma=.15  and sigma=0., 0 to khor
         if(sig(k)<0.15)hdiff(k)=.1*max(1.,(1.-sig(k)/.15)*abs(khor))
        enddo
        if(myid==0)write(6,*)'khor,hdiff: ',khor,hdiff
      endif

!     *** be careful not to choose ia,ib,ja,jb to cover more than
!         one processor, or myid=0 may stop here!!!
      call printa('zs  ',zs,0,0,ia,ib,ja,jb,0.,.01)
      call printa('tss ',tss,0,0,ia,ib,ja,jb,200.,1.)
      if(mydiag)write(6,*)'wb(idjd) ',(wb(idjd,k),k=1,6)
      call printa('wb1   ',wb ,0,1,ia,ib,ja,jb,0.,100.)
      call printa('wb6  ',wb,0,ms,ia,ib,ja,jb,0.,100.)

      !--------------------------------------------------------------
      ! NRUN COUNTER
      if (myid==0) then
        open(11, file='nrun.dat',status='unknown')
        if(nrun==0)then
          read(11,*,iostat=ierr2) nrun
          nrun=nrun+1
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
      rndmax(:)=0.
      tmaxscr(:)=0.
      tminscr(:)=400.
      rhmaxscr(:)=0.
      rhminscr(:)=400.
      u10max(:)=0.
      v10max(:)=0.
      u1max(:)=0.
      v1max(:)=0.
      u2max(:)=0.
      v2max(:)=0.
      capemax(:)=0.
      u10mx(:)=0.
      tscr_ave(:)=0.
      qscrn_ave(:)=0.
      dew_ave(:)=0.
      epan_ave(:)=0.
      epot_ave(:)=0.
      eg_ave(:)=0.
      fg_ave(:)=0.
      rnet_ave(:)=0.
      sunhours(:)=0.
      ga_ave(:)=0.
      riwp_ave(:)=0.
      rlwp_ave(:)=0.
      evap(:)=0.
      precc(:)=0.
      precip(:)=0.
      rnd_3hr(:,8)=0. ! i.e. rnd24(:)=0.
      cbas_ave(:)=0.
      ctop_ave(:)=0.
      sno(:)=0.
      runoff(:)=0.
      wb_ave(:,:)=0.
      tsu_ave(:)=0.
      alb_ave(:)=0.
      fbeam_ave(:)=0.
      psl_ave(:)=0.
      mixdep_ave(:)=0.
      koundiag=0
      sint_ave(:) = 0.  ! solar_in_top
      sot_ave(:)  = 0.  ! solar_out_top
      soc_ave(:)  = 0.  ! solar_out_top (clear sky)
      sgdn_ave(:) = 0.  ! solar_ground (down-welling) +ve down
      sgn_ave(:)  = 0.  ! solar_ground (net) +ve down
      rtu_ave(:)  = 0.  ! LW_out_top 
      rtc_ave(:)  = 0.  ! LW_out_top (clear sky)
      rgdn_ave(:) = 0.  ! LW_ground (down-welling)  +ve down
      rgn_ave(:)  = 0.  ! LW_ground (net)  +ve up
      rgc_ave(:)  = 0.  ! LW_ground (clear sky)
      cld_ave(:)  = 0.
      cll_ave(:)  = 0.
      clm_ave(:)  = 0.
      clh_ave(:)  = 0.
      if (ngas>0) then
        traver=0.
      end if
      fpn_ave=0.
      frs_ave=0.
      frp_ave=0.


      !--------------------------------------------------------------
      ! OPEN OUTPUT FILES AND SAVE INITAL CONDITIONS
      if(nmi==0.and.nwt>0)then
!       write out the first ofile data set
        if (myid==0) write(6,*)'calling outfile'
        call outfile(20,rundate,nmi,nwrite,iaero)  ! which calls outcdf
        if(newtop<0) then
          if (myid==0) write(6,*) 'newtop<0 requires a stop here'
          stop  ! just for outcdf to plot zs  & write fort.22
        end if
      endif    ! (nmi==0.and.nwt.ne.0)


      !--------------------------------------------------------------
      ! INITIALISE DYNAMICS AND NUDGING
      dtin=dt
      n3hr=1   ! initial value at start of run
      if (myid==0) then
        write(6,*)'number of time steps per day = ',nperday
        write(6,*)'nper3hr,nper6hr .. ',nper3hr(:)
      end if
      mspeca=1
      if(mex.ne.1.and..not.lrestart)then
        mspeca=2
        dt=dtin*.5
      endif
      call gettin(0)             ! preserve initial mass & T fields; nmi too

      nmaxprsav=nmaxpr
      nwtsav=nwt
      hrs_dt = dtin/3600.      ! time step in hours
      mins_dt = nint(dtin/60.)  ! time step in minutes
      mtimer_in=mtimer
      nstagin=nstag    ! -ve nstagin gives swapping & its frequency
      nstaguin=nstagu  ! only the sign of nstaguin matters (chooses scheme)
      if(nstagin==5.or.nstagin<0)then
        nstag=4
        nstagu=4
        if(nstagin==5)then  ! for backward compatability
          nstagin=-1 
          nstaguin=5  
        endif
      endif
 
 
      !--------------------------------------------------------------
      ! BEGIN MAIN TIME LOOP
       if ( myid == 0 ) then
         call date_and_time(time=timeval,values=tvals1)
         write(6,*) "Start of loop time ", timeval
      end if
      call log_on()
#ifdef simple_timer
      call start_log(maincalc_begin)
#endif

      do 88 kktau=1,ntau   ! ****** start of main time loop
      ktau=kktau
      timer = timer + hrs_dt      ! timer now only used to give timeg
      timeg=mod(timer+hourst,24.)
      mtimer=mtimer_in+nint(ktau*dtin/60.)     ! 15/6/01 to allow dt < 1 minute
      mins_gmt=mod(mtimer+60*ktime/100,24*60)

      ! ***********************************************************************
      ! START DYNAMICS AND NUDGING
      ! ***********************************************************************

      ! NESTING ---------------------------------------------------------------
      call start_log(nestin_begin)
      if (nbd.ne.0) call nestin(iaero)
      call end_log(nestin_end)
      
      ! TRACERS ---------------------------------------------------------------
!     rml 17/02/06 interpolate tracer fluxes to current timestep
      if(ngas>0) call interp_tracerflux(kdate,hrs_dt)

      ! DYNAMICS --------------------------------------------------------------
      if(nstaguin>0.and.ktau>1)then   ! swapping here for nstaguin>0
        if(nstagin<0.and.mod(ktau,abs(nstagin))==0)then
          nstag=7-nstag  ! swap between 3 & 4
          nstagu=nstag
        endif
      endif

      do 79 mspec=mspeca,1,-1    ! start of introductory time loop
      dtds=dt/ds
      if(nvsplit<3.or.ktau==1)then
        un(1:ifull,:)=0. 
        vn(1:ifull,:)=0.
        tn(1:ifull,:)=0.
      elseif(nvsplit==3)then
        tn(1:ifull,:)=(t(1:ifull,:)-tx(1:ifull,:))/dt  ! tend. from phys. at end of previous step
        unn(1:ifull,:)=(u(1:ifull,:)-ux(1:ifull,:))/dt ! used in nonlin,upglobal
        vnn(1:ifull,:)=(v(1:ifull,:)-vx(1:ifull,:))/dt ! used in nonlin,upglobal
        t(1:ifull,:)=tx(1:ifull,:)   
        u(1:ifull,:)=ux(1:ifull,:)   
        v(1:ifull,:)=vx(1:ifull,:)   
      elseif(nvsplit==4)then
        unn(1:ifull,:)=(u(1:ifull,:)-ux(1:ifull,:))/dt ! used in nonlin,upglobal
        vnn(1:ifull,:)=(v(1:ifull,:)-vx(1:ifull,:))/dt ! used in nonlin,upglobal
        u(1:ifull,:)=ux(1:ifull,:)   
        v(1:ifull,:)=vx(1:ifull,:)   
        tn(1:ifull,:)=0.
      elseif(nvsplit==5)then
        tn(1:ifull,:)=(t(1:ifull,:)-tx(1:ifull,:))/dt ! t tend. from phys.
        t(1:ifull,:)=tx(1:ifull,:)   
        un(1:ifull,:)=0. 
        vn(1:ifull,:)=0.
      endif   ! (nvsplit<3.or.ktau==1) .. elseif ..

      call bounds(qg)
      call bounds(psl)
      if(mup.ne.1.or.(ktau==1.and.mspec==mspeca.and..not.lrestart))
     &    then
!       updps called first step or to permit clean restart option      
        call updps(0) 
      endif

!     set up tau +.5 velocities in ubar, vbar
      if(ktau<10.and.mydiag)then
        write(6,*) 'ktau,mex,mspec,mspeca:',ktau,mex,mspec,mspeca
      endif
      sbar(:,2:kl)=sdot(:,2:kl)
      if(ktau==1.and..not.lrestart)then
!       this sets (ubar,vbar) to ktau=1.5 values on 2nd time through
        ubar(:,:)=u(1:ifull,:)
        vbar(:,:)=v(1:ifull,:)
      elseif(mex==1)then
        ubar(:,:)=u(1:ifull,:)
        vbar(:,:)=v(1:ifull,:)
      elseif((ktau==2.and..not.lrestart).or.mex==2)then        
!       (tau+.5) from tau, tau-1
        ubar(:,:)=u(1:ifull,:)*1.5-savu(:,:)*.5
        vbar(:,:)=v(1:ifull,:)*1.5-savv(:,:)*.5
      elseif(mex==3)then
!       (tau+.5) from tau, tau-1, tau-2   ! ubar is savu1 here
        ubar(:,:)=u(1:ifull,:)+.5*(savu(:,:)-savu1(:,:))
        vbar(:,:)=v(1:ifull,:)+.5*(savv(:,:)-savv1(:,:))
       elseif(mex==30.and.(ktau>3.or.lrestart))then  ! using tau, tau-1, tau-2, tau-3
        do k=1,kl
         do iq=1,ifull
          bb=1.5*u(iq,k)-2.*savu(iq,k)+.5*savu1(iq,k)    ! simple b
          bb_2=(40.*u(iq,k)-35.*savu(iq,k)               ! cwqls b
     &          -16.*savu1(iq,k)+11.*savu2(iq,k))/34.
          cc=.5*u(iq,k)-savu(iq,k)+.5*savu1(iq,k)        ! simple c
          cc_2=(10.*u(iq,k)-13.*savu(iq,k)               ! cwqls c
     &          -4.*savu1(iq,k)+7.*savu2(iq,k))/34.
          aa=cc_2-cc
          rat=max(0.,min(1.,cc_2/(aa+sign(1.e-9,aa))))
          cc=rat*cc+(1.-rat)*cc_2 
          bb=rat*bb+(1.-rat)*bb_2 
          ubar(iq,k)=u(iq,k)+.5*bb+.25*cc
          bb=1.5*v(iq,k)-2.*savv(iq,k)+.5*savv1(iq,k)    ! simple b
          bb_2=(40.*v(iq,k)-35.*savv(iq,k)               ! cwqls b
     &          -16.*savv1(iq,k)+11.*savv2(iq,k))/34.
          cc=.5*v(iq,k)-savv(iq,k)+.5*savv1(iq,k)        ! simple c
          cc_2=(10.*v(iq,k)-13.*savv(iq,k)               ! cwqls c
     &          -4.*savv1(iq,k)+7.*savv2(iq,k))/34.
          aa=cc_2-cc
          rat=max(0.,min(1.,cc_2/(aa+sign(1.e-9,aa))))
          cc=rat*cc+(1.-rat)*cc_2 
          bb=rat*bb+(1.-rat)*bb_2 
          vbar(iq,k)=v(iq,k)+.5*bb+.25*cc
         enddo  ! iq loop
        enddo   ! k loop 
      else      ! i.e. mex >=4 and ktau>=3
!       (tau+.5) from tau, tau-1, tau-2   ! ubar is savu1 here
        ubar(:,:)=u(1:ifull,:)*15./8.-savu(:,:)*10./8.+savu1(:,:)*3./8.
        vbar(:,:)=v(1:ifull,:)*15./8.-savv(:,:)*10./8.+savv1(:,:)*3./8.
      endif    ! (ktau==1) .. else ..
      if (mod(ktau,nmaxpr)==0.and.mydiag)then
        nlx=max(2,nlv)  ! as savs not defined for k=1
        write (6,"(i4,' savu2,savu1,savu,u,ubar',5f8.2)") ktau,
     &      savu2(idjd,nlv),savu1(idjd,nlv),savu(idjd,nlv),
     &      u(idjd,nlv),ubar(idjd,nlv)
        write (6,"(i4,' savv2,savv1,savv,v,vbar',5f8.2)") ktau,
     &      savv2(idjd,nlv),savv1(idjd,nlv),savv(idjd,nlv),
     &      v(idjd,nlv),vbar(idjd,nlv)
      endif
      if(ktau>2.and.epsp>1..and.epsp<2.)then ! 
        if (mydiag.and.ktau==3)write(6,*)'using epsp= ',epsp
        do iq=1,ifull
         if((sign(1.,dpsdt(iq)).ne.sign(1.,dpsdtb(iq))).and.
     &      (sign(1.,dpsdtbb(iq)).ne.sign(1.,dpsdtb(iq))))then
           epst(iq)=epsp-1.
         else
           epst(iq)=0.
         endif
        enddo
      endif ! (ktau>2.and.epsp>1..and.epsp<2.)

      if (ktau<10.and.mydiag)then
        write(6,*)'savu,u,ubar ',ktau,savu(idjd,1),u(idjd,1),
     &                           ubar(idjd,1)
      endif
      if(ktau==1.and..not.lrestart.and.mspec==1.and.mex.ne.1)then
        u(1:ifull,:)=savu(1:ifull,:)  ! reset u,v to original values
        v(1:ifull,:)=savv(1:ifull,:)
      endif
      savu2(1:ifull,:)=savu1(1:ifull,:)  
      savv2(1:ifull,:)=savv1(1:ifull,:)
      savs1(1:ifull,:)=savs(1:ifull,:)  
      savu1(1:ifull,:)=savu(1:ifull,:)  
      savv1(1:ifull,:)=savv(1:ifull,:)
      savs(1:ifull,:) =sdot(1:ifull,2:kl)  
      savu(1:ifull,:) =u(1:ifull,:)  ! before any time-splitting occurs
      savv(1:ifull,:) =v(1:ifull,:)

      ! set diagnostic printout flag
      diag=(ktau>=abs(ndi).and.ktau<=ndi2)
      if(ndi<0)then
        if(ktau==ktau/ndi*ndi)diag=.true.
      endif

      ! save tracer arrays
      if(ngas>=1)then ! re-set trsav prior to vadv, hadv, hordif
!       N.B. nllp arrays are after gas arrays
        do igas=1,ngas
         do k=1,klt
          do iq=1,ilt*jlt 
           trsav(iq,k,igas)=tr(iq,k,igas) ! 4D for tr conservation in adjust5
          enddo
         enddo
        enddo
      endif      ! (ngas>=1)

      ! update non-linear dynamic terms
      call nonlin(iaero)
      if (diag)then
         if (mydiag) write(6,*) 'before hadv'
         call printa('tx  ',tx,ktau,nlv,ia,ib,ja,jb,0.,1.)
         if (mydiag) then
            nlx=min(nlv,kl-8)
            write(6,"('tx  ',9f8.2)") (tx(idjd,k),k=nlx,nlx+8)
            write(6,"('txe ',9f8.2)") (tx(ie(idjd),k),k=nlx,nlx+8)
            write(6,"('txw ',9f8.2)") (tx(iw(idjd),k),k=nlx,nlx+8)
            write(6,"('txn ',9f8.2)") (tx(in(idjd),k),k=nlx,nlx+8)
            write(6,"('txs ',9f8.2)") (tx(is(idjd),k),k=nlx,nlx+8)
            write(6,'(i2," qgv ",18f7.4)')ktau,(1000.*qg(idjd,k),k=1,kl)
         end if
         call printa('qgv ',qg,ktau,nlv,ia,ib,ja,jb,0.,1.e3)
      endif

!     evaluate horizontal advection for combined quantities
      call upglobal(iaero)
      if (diag)then
        if (mydiag) then
          write(6,*) 'after hadv'
          write (6,"('tx  ',9f8.2)") (tx(idjd,k),k=nlx,nlx+8)
        end if
        call printa('tx  ',tx,ktau,nlv,ia,ib,ja,jb,200.,1.)
        if (mydiag)write(6,'(i2," qgh ",18f7.4)')ktau,1000.*qg(idjd,:)
      end if

      if(nonl<0)then
        savt(1:ifull,:)=t(1:ifull,:)  ! can be used in nonlin during next step
      end if

      if(nstaguin<0.and.ktau>1)then  ! swapping here (lower down) for nstaguin<0
        if(nstagin<0.and.mod(ktau,abs(nstagin))==0)then
          nstag=7-nstag  ! swap between 3 & 4
          nstagu=nstag
        end if
      end if
      call adjust5(iaero)
      
      ! NESTING ---------------------------------------------------------------
      ! nesting now after mass fixers
      call start_log(nestin_begin)
      if (mspec==1) then
        if (mbd.ne.0) then
          call nestinb(iaero)
        else if (nbd.ne.0) then
          call davies
        end if
      end if
      call end_log(nestin_end)

      ! DYNAMICS --------------------------------------------------------------
      if(mspec==2)then     ! for very first step restore mass & T fields
        call gettin(1)
      endif    !  (mspec==2) 
      if(mfix_qg==0.or.mspec==2)then
        qfg(1:ifull,:)=max(qfg(1:ifull,:),0.) 
        qlg(1:ifull,:)=max(qlg(1:ifull,:),0.) 
        do k=1,kl
         do iq=1,ifull
          if(qfg(iq,k)+qlg(iq,k)<qgmin)qg(iq,k)=max(qg(iq,k),qgmin)
         enddo
        enddo
      endif  ! (mfix_qg==0.or.mspec==2)
79    dt=dtin                    ! ****** end of introductory time loop
      mspeca=1
      if(nvsplit>=3)then
        tx(1:ifull,:)=t(1:ifull,:)   ! saved for beginning of next step
        ux(1:ifull,:)=u(1:ifull,:)   
        vx(1:ifull,:)=v(1:ifull,:)   
      endif

      ! horizontal diffusion
      if (nhor<0) call hordifgt(iaero)  ! now not tendencies
      if (diag.and.mydiag) write(6,*) 'after hordifgt t ',t(idjd,:)

      ! ***********************************************************************
      ! START PHYSICS 
      ! ***********************************************************************
      call start_log(phys_begin)

      ! GWDRAG ----------------------------------------------------------------
      call start_log(gwdrag_begin)
      if (myid==0.and.nmaxpr==1) then
        write(6,*) "Before gwdrag"
      end if
      if (ngwd<0) call gwdrag  ! <0 for split - only one now allowed
      if (myid==0.and.nmaxpr==1) then
        write(6,*) "After gwdrag"
      end if
      call end_log(gwdrag_end)

      ! CONVECTION ------------------------------------------------------------
      call start_log(convection_begin)
      if (myid==0.and.nmaxpr==1) then
        write(6,*) "Before convection"
      end if
      select case(nkuo)
        case(5)
          call betts(t,qg,tn,land,ps) ! not called these days
        case(23,24)
          call convjlm(iaero)         ! split convjlm 
        case(46)
          call conjob(iaero)          ! split Arakawa-Gordon scheme
      end select
      if (nkuo/=0) then
        ! Not set in HS tests.
        cbas_ave(:)=cbas_ave(:)+condc(:)*(1.1-sig(kbsav(:))) ! diagnostic
        ctop_ave(:)=ctop_ave(:)+condc(:)*(1.1-sig(abs(ktsav(:)))) ! diagnostic
      end if
      if (myid==0.and.nmaxpr==1) then
        write(6,*) "After convection"
      end if
      call end_log(convection_end)

      ! CLOUD MICROPHYSICS ----------------------------------------------------
      call start_log(cloud_begin)
      if (myid==0.and.nmaxpr==1) then
        write(6,*) "Before cloud microphysics"
      end if
      select case(ldr)
        case(-2,-1,1,2)
          call leoncld(cfrac,cffall,iaero) ! LDR microphysics scheme
      end select
      do k=1,kl
       riwp_ave(:)=riwp_ave(:)-qfrad(:,k)*dsig(k)*ps(1:ifull)/grav ! ice water path
       rlwp_ave(:)=rlwp_ave(:)-qlrad(:,k)*dsig(k)*ps(1:ifull)/grav ! liq water path
      enddo
      if(nmaxpr==1.and.mydiag)then
        write (6,"('qfrad',3p9f8.3/5x,9f8.3)") qfrad(idjd,:)
        write (6,"('qlrad',3p9f8.3/5x,9f8.3)") qlrad(idjd,:)
        write (6,"('qf   ',3p9f8.3/5x,9f8.3)") qfg(idjd,:)
      endif
      rnd_3hr(:,8)=rnd_3hr(:,8)+condx(:)  ! i.e. rnd24(:)=rnd24(:)+condx(:)
      if (myid==0.and.nmaxpr==1) then
        write(6,*) "After cloud microphysics"
      end if
      call end_log(cloud_end)

      ! RADIATION -------------------------------------------------------------
      if (myid==0.and.nmaxpr==1) then
        write(6,*) "Before radiation"
      end if
      select case(nrad)
       case(4)
!       Fels-Schwarzkopf radiation
        odcalc=mod(ktau,kountr)==0 .or. ktau==1 ! ktau-1 better
        nnrad=kountr
        if (nhstest<0) then ! aquaplanet test -1 to -8  
         mtimer_sav=mtimer
         mtimer=mins_gmt     ! so radn scheme repeatedly works thru same day
        end if    ! (nhstest<0)
        call radrive (il*nrows_rad,odcalc,iaero)
        if (nhstest<0) then ! aquaplanet test -1 to -8  
          mtimer=mtimer_sav
        end if    ! (nhstest<0)
        if (nmaxpr==1) then
          ! Account for load bal explicitly rather than implicitly in
          ! the reduce in maxmin.
          call phys_loadbal
          call maxmin(slwa,'sl',ktau,.1,1)
        end if
       case(5)
        ! GFDL SEA-EFS radiation
        odcalc=mod(ktau,kountr)==0.or.ktau==1
        nnrad=kountr
        if (nhstest<0) then
         mtimer_sav=mtimer
         mtimer=mins_gmt
        end if
        call seaesfrad(il*nrows_rad,odcalc,iaero)
        if (nhstest<0) then
          mtimer=mtimer_sav
        end if
       case DEFAULT
!       use preset slwa array (use +ve nrad)
        slwa(:)=-10*nrad  
!       N.B. no rtt array for this nrad option
      end select
      if (myid==0.and.nmaxpr==1) then
        write(6,*) "After radiation"
      end if

      ! HELD & SUAREZ ---------------------------------------------------------
      if (ntsur<=1.or.nhstest==2) then ! Held & Suarez or no surf fluxes
       eg(:)=0.
       fg(:)=0.
       cdtq(:)=0.
       cduv(:)=0.
      end if     ! (ntsur<=1.or.nhstest==2) 
      if(nhstest==2) call hs_phys
      
      if (ntsur>1) then  ! should be better after convjlm
       if (diag) then
         call maxmin(u,'#u',ktau,1.,kl)
         call maxmin(v,'#v',ktau,1.,kl)
         call maxmin(t,'#t',ktau,1.,kl)
         call maxmin(qg,'qg',ktau,1.e3,kl)     
         call MPI_Barrier( MPI_COMM_WORLD, ierr ) ! stop others going past
       end if
         
       ! SURFACE FLUXES ---------------------------------------------
       ! (Includes ocean dynamics and mixing, as well as ice dynamics and thermodynamics)
       if (myid==0.and.nmaxpr==1) then
        write(6,*) "Before surface fluxes"
       end if
       call sflux(nalpha)
       epan_ave=epan_ave+epan  ! 2D 
       epot_ave=epot_ave+epot  ! 2D 
       ga_ave=ga_ave+ga        ! 2D 
       if (myid==0.and.nmaxpr==1) then
        write(6,*) "After surface fluxes"
       end if
       
       ! STATION OUTPUT ---------------------------------------------
       if (nstn>0.and.nrotstn(1)==0) call stationa ! write every time step
       
       ! DIAGNOSTICS ------------------------------------------------
       if (mod(ktau,nmaxpr)==0.and.mydiag) then
         write(6,*)
         write (6,
     .	  "('ktau =',i5,' gmt(h,m):',f6.2,i5,' runtime(h,m):',f7.2,i6)")
     .      ktau,timeg,mins_gmt,timer,mtimer
!        some surface (or point) diagnostics
         isoil = isoilm(idjd)
         write(6,*) 'land,isoil,ivegt,isflag ',
     &          land(idjd),isoil,ivegt(idjd),isflag(idjd)
         write (6,"('snage,snowd,alb   ',f8.4,2f8.2)")
     &      snage(idjd),snowd(idjd),albvisnir(idjd,1)
         write (6,"('sicedep,fracice,runoff ',3f8.2)")
     &            sicedep(idjd),fracice(idjd),runoff(idjd)
         write (6,"('tgg(1-6)   ',9f8.2)") (tgg(idjd,k),k=1,6)
         write (6,"('tggsn(1-3) ',9f8.2)") (tggsn(idjd,k),k=1,3)
         write (6,"('wb(1-6)    ',9f8.3)") (wb(idjd,k),k=1,6)
         write (6,"('wbice(1-6) ',9f8.3)") (wbice(idjd,k),k=1,6)
         write (6,"('smass(1-3) ',9f8.2)") (smass(idjd,k),k=1,3) ! as mm of water
         write (6,"('ssdn(1-3)  ',9f8.2)") (ssdn(idjd,k),k=1,3)
         pwater=0.   ! in mm
         iq=idjd
         div_int=0.
         do k=1,kl
          qtot=qg(iq,k)+qlg(iq,k)+qfg(iq,k)
           pwater=pwater-dsig(k)*qtot*ps(idjd)/grav
           div(k)=(u(ieu(iq),k)/emu(ieu(iq))
     &            -u(iwu(iq),k)/emu(iwu(iq))  
     &            +v(inv(iq),k)/emv(inv(iq))
     &            -v(isv(iq),k)/emv(isv(iq))) 
     &             *em(iq)**2/(2.*ds)  *1.e6
           div_int=div_int-div(k)*dsig(k)
         enddo
         write (6,"('pwater,condc,condx,rndmax,rmc',9f8.3)")
     &      pwater,condc(idjd),condx(idjd),rndmax(idjd),cansto(idjd)
         write (6,"('wetfac,sno,evap,precc,precip',
     &      6f8.2)") wetfac(idjd),sno(idjd),evap(idjd),precc(idjd),
     &      precip(idjd)
         write (6,"('tmin,tmax,tscr,tss,tpan',9f8.2)")
     &      tminscr(idjd),tmaxscr(idjd),tscrn(idjd),tss(idjd),
     &      tpan(idjd)
         write (6,"('u10,ustar,pblh',9f8.2)")
     &      u10(idjd),ustar(idjd),pblh(idjd)
         write (6,"('div_int,ps,qgscrn',5f8.2,f8.3)")
     &      div_int,.01*ps(idjd),1000.*qgscrn(idjd)
         write (6,"('dew_,eg_,epot,epan,eg,fg,ga',9f8.2)") 
     &      dew_ave(idjd),eg_ave(idjd),epot(idjd),epan(idjd),eg(idjd),
     &      fg(idjd),ga(idjd)
         write (6,"('zo,cduv', 
     &      2f8.5)") zo(idjd),cduv(idjd)/vmod(idjd)
         write (6,"('slwa,sint,sg,rt,rg    ',9f8.2)") 
     &      slwa(idjd),sintsave(idjd),sgsave(idjd),
     &      rtsave(idjd),rgsave(idjd)
         write (6,"('cll,clm,clh,clt ',9f8.2)") 
     &      cloudlo(idjd),cloudmi(idjd),cloudhi(idjd),cloudtot(idjd)
         write (6,"('u10max,v10max,rhmin,rhmax   ',9f8.2)")
     &               u10max(iq),v10max(iq),rhminscr(iq),rhmaxscr(iq)
         write (6,"('kbsav,ktsav,convpsav ',2i3,f8.4,9f8.2)")
     &               kbsav(idjd),ktsav(idjd),convpsav(idjd)
         write (6,"('t   ',9f8.3/4x,9f8.3)") t(idjd,:)
         write (6,"('u   ',9f8.3/4x,9f8.3)") u(idjd,:)
         write (6,"('v   ',9f8.3/4x,9f8.3)") v(idjd,:)
         write (6,"('qg  ',3p9f8.3/4x,9f8.3)") qg(idjd,:)
         write (6,"('qf  ',3p9f8.3/4x,9f8.3)") qfg(idjd,:)
         write (6,"('ql  ',3p9f8.3/4x,9f8.3)") qlg(idjd,:)
         write (6,"('cfrac',9f8.3/5x,9f8.3)") cfrac(idjd,:)
         do k=1,kl
          es=establ(t(idjd,k))
          spmean(k)=100.*qg(idjd,k)*
     &              max(ps(idjd)*sig(k)-es,1.)/(.622*es) ! max as for convjlm
         enddo
         write (6,"('rh  ',9f8.3/4x,9f8.3)") spmean(:)
         write (6,"('div ',9f8.3/4x,9f8.3)") div(:)
         write (6,"('omgf ',9f8.3/5x,9f8.3)")   ! in Pa/s
     &             ps(idjd)*dpsldt(idjd,:)
         write (6,"('sdot ',9f8.3/5x,9f8.3)") sdot(idjd,1:kl)
         if(nextout>=4)write (6,"('xlat,long,pres ',3f8.2)")
     &    tr(idjd,nlv,ngas+1),tr(idjd,nlv,ngas+2),tr(idjd,nlv,ngas+3)
       endif  ! (mod(ktau,nmaxpr)==0.and.mydiag)
      endif   ! (ntsur>1)

      ! VERTICAL MIXING ------------------------------------------------------
      if (myid==0.and.nmaxpr==1) then
       write(6,*) "Before PBL mixing"
      end if
      if (ntsur>=1) then ! calls vertmix but not sflux for ntsur=1
        call start_log(vertmix_begin)
        if(nmaxpr==1.and.mydiag)
     &    write (6,"('pre-vertmix t',9f8.3/13x,9f8.3)") t(idjd,:)
        call vertmix(iaero) 
        if(nmaxpr==1.and.mydiag)
     &    write (6,"('aft-vertmix t',9f8.3/13x,9f8.3)") t(idjd,:)
        call end_log(vertmix_end)
      endif  ! (ntsur>=1)
      if (myid==0.and.nmaxpr==1) then
       write(6,*) "After PBL mixing"
      end if

      ! AEROSOLS --------------------------------------------------------------
      if (abs(iaero).ge.2) then
        call start_log(aerosol_begin)
         if (myid==0.and.nmaxpr==1) then
         write(6,*) "Before aerosols"
        end if
        call aerocalc
        if (myid==0.and.nmaxpr==1) then
         write(6,*) "After aerosols"
        end if
        call end_log(aerosol_end)
      end if

      ! PHYSICS LOAD BALANCING ------------------------------------------------
!     This is the end of the physics. The next routine makes the load imbalance
!     overhead explicit rather than having it hidden in one of the diagnostic
!     calls.
      call phys_loadbal
      call end_log(phys_end)

      ! ***********************************************************************
      ! TRACERS
      ! ***********************************************************************
!     rml 16/02/06 call tracer_mass, write_ts
      if(ngas>0) then
        call tracer_mass(ktau,ntau) !also updates average tracer array
        call write_ts(ktau,ntau,dt)
!     rml 23/4/10 if final timestep write accumulated loss to tracer.stdout file
        if (ktau.eq.ntau) then
          if (myid == 0) then
            write(unit_trout,*) 'Methane loss accumulated over month'
            write(unit_trout,*) ktau,acloss_g(:)
          endif
        endif
      endif

      ! ***********************************************************************
      ! DIAGNOSTICS AND OUTPUT
      ! ***********************************************************************
      if(ndi==-ktau)then
        nmaxpr=1         ! diagnostic prints; reset 6 lines on
        if(ndi2==0)ndi2=ktau+40
      endif
      if(ktau==ndi2)then
         if(mydiag)write(6,*)'reset nmaxpr'
         nmaxpr=nmaxprsav
      endif
      if (mod(ktau,nmaxpr)==0.or.ktau==ntau)then
        call maxmin(u,' u',ktau,1.,kl)
        call maxmin(v,' v',ktau,1.,kl)
        speed(:,:)=u(1:ifull,:)**2+v(1:ifull,:)**2 ! 3D
        call average(speed,spmean,spavge)
        do k=1,kl
         spmean(k)=sqrt(spmean(k))
        enddo
        speed(:,:)=sqrt(speed(:,:)) ! 3D
        call maxmin(speed,'sp',ktau,1.,kl)
        call maxmin(t,' t',ktau,1.,kl)
        call maxmin(qg,'qg',ktau,1.e3,kl)
        call maxmin(qfg,'qf',ktau,1.e3,kl)
        call maxmin(qlg,'ql',ktau,1.e3,kl)
        call maxmin(sdot,'sd',ktau,1.,kl)  ! grid length units 
        if ( myid==0 ) then
           write(6,'("spmean ",9f8.3)') spmean
           write(6,'("spavge ",f8.3)') spavge
        end if
        call average(qg,spmean,spavge)
        if ( myid==0 ) then
           write(6,'("qgmean ",9f8.5)') spmean
           write(6,'("qgavge ",f8.5)') spavge
        end if
        if(ngas>0)then
          k2=min(2,klt)
          do ng=1,ngas
           write (text,'("g",i1)')ng
           call maxmin(tr(:,:,ng),text,ktau,1.,kl)
          enddo
          call average(tr(:,:,1),spmean,spavge)
          write(6,'("g1mean ",9f8.3)') spmean
          write(6,'("g1avge ",f8.3)') spavge
          call average(tr(:,:,k2),spmean,spavge)
          write(6,'("g2mean ",9f8.3)') spmean
          write(6,'("g2avge ",f8.3)') spavge
        endif   ! (ngas>0)
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
        psavge=0.
        pslavge=0.
        preccavge=0.
        precavge=0.
        do iq=1,ifull
         psavge=psavge+ps(iq)*wts(iq)
         pslavge=pslavge+psl(iq)*wts(iq)
         preccavge=preccavge+precc(iq)*wts(iq)
         precavge=precavge+precip(iq)*wts(iq)
        enddo
!       KE calculation, not taking into account pressure weighting
        gke = 0.
        do k=1,kl
          do iq=1,ifull
            gke = gke - 0.5 * wts(iq) * dsig(k) *
     &           ( u(iq,k)**2 + v(iq,k)**2 )
          end do
        end do
        cllav=0.
        clmav=0.
        clhav=0.
        cltav=0.
        do iq=1,ifull
          cllav=cllav+wts(iq)*cloudlo(iq)
          clmav=clmav+wts(iq)*cloudmi(iq)
          clhav=clhav+wts(iq)*cloudhi(iq)
          cltav=cltav+wts(iq)*cloudtot(iq)   
        end do

        ! All this combined into a single reduction
        temparray = (/ psavge, pslavge, preccavge, precavge, gke,
     &                 cllav, clmav,clhav, cltav /)
        call MPI_Reduce ( temparray, gtemparray, 9, MPI_REAL, MPI_SUM,0,
     &                  MPI_COMM_WORLD, ierr )
        if (myid==0) then
          write(6,97) gtemparray(1:5) ! psavge,pslavge,preccavge,precavge,gke
 97       format(' average ps, psl, precc, prec, gke: ',
     &           f10.2,f10.6,2f6.2,f7.2)
          write(6,971) gtemparray(6:9) ! cllav,clmav,clhav,cltav
 971      format(' global_average cll, clm, clh, clt: ',4f6.2)
        end if
        if (mydiag) then
          write(6,98) ktau,diagvals(ps)
 98       format(i7,' ps diag:',-2p9f7.1)
          if(t(idjd,kl)>258.)then
            write(6,*) 't(idjd,kl) > 258. for idjd = ',idjd
            write(6,91) ktau,(t(idjd,k),k=kl-8,kl)
 91         format(i7,'    t',9f7.2)
            write(6,92) ktau,(sdot(idjd,k),k=kl-8,kl)
 92         format(i7,' sdot',9f7.3)
          endif                ! (t(idjd,kl)>258.)
        end if                 ! myid==0
      endif      ! (mod(ktau,nmaxpr)==0)

!     update diag_averages and daily max and min screen temps 
!     N.B. runoff is accumulated in sflux
      tmaxscr  = max(tmaxscr,tscrn)
      tminscr  = min(tminscr,tscrn)
      rhmaxscr = max(rhmaxscr,rhscrn)
      rhminscr = min(rhminscr,rhscrn)
      rndmax   = max(rndmax,condx)
      capemax  = max(capemax,cape)
      u10mx    = max(u10mx,u10)  ! for hourly scrnfile
      dew_ave  = dew_ave-min(0.,eg)    
      eg_ave   = eg_ave+eg    
      fg_ave   = fg_ave+fg
      rnet_ave = rnet_ave+rnet
      tscr_ave = tscr_ave+tscrn 
      qscrn_ave = qscrn_ave+qgscrn 
      wb_ave   = wb_ave+wb
      tsu_ave  = tsu_ave+tss
      psl_ave  = psl_ave+psl
      call mlodiag(spare1,0)
      mixdep_ave=mixdep_ave+spare1
      spare1(:)=u(1:ifull,1)**2+v(1:ifull,1)**2
      spare2(:)=u(1:ifull,2)**2+v(1:ifull,2)**2
      do iq=1,ifull
       if (u10(iq)**2>u10max(iq)**2+v10max(iq)**2) then
         u10max(iq)=u10(iq)*u(iq,1)/max(.001,sqrt(spare1(iq)))
         v10max(iq)=u10(iq)*v(iq,1)/max(.001,sqrt(spare1(iq)))
       end if
       if (spare1(iq)>u1max(iq)**2+v1max(iq)**2) then
         u1max(iq)=u(iq,1)
         v1max(iq)=v(iq,1)
       end if
       if (spare2(iq)>u2max(iq)**2+v2max(iq)**2) then
         u2max(iq)=u(iq,2)
         v2max(iq)=v(iq,2)
       end if
      end do
      if (ngas>0) then
        do igas=1,ngas
          traver(:,:,igas)=traver(:,:,igas)+tr(1:ilt*jlt,:,igas)
        end do
      end if
      fpn_ave=fpn_ave+fpn
      frs_ave=frs_ave+frs
      frp_ave=frp_ave+frp

!     rnd03 to rnd21 are accumulated in mm     
      if (myid==0) then
        write(6,*) 'ktau,mod,nper3hr ',
     &            ktau,mod(ktau-1,nperday)+1,nper3hr(n3hr)
      end if
      if(mod(ktau-1,nperday)+1==nper3hr(n3hr))then
        rnd_3hr(:,n3hr)=rnd_3hr(:,8)
        if(nextout>=2)then
          spare1(:)=max(.001,sqrt(u(1:ifull,1)**2+v(1:ifull,1)**2))
          u10_3hr(:,n3hr)=u10(:)*u(1:ifull,1)/spare1(:)
          v10_3hr(:,n3hr)=u10(:)*v(1:ifull,1)/spare1(:)
          tscr_3hr(:,n3hr)=tscrn(:)
          do iq=1,ifull
           es=establ(t(iq,1))
           rh1_3hr(iq,n3hr)=100.*qg(iq,1)*
     &                        (ps(iq)*sig(1)-es)/(.622*es)
          enddo
        endif    ! (nextout==2)
        n3hr=n3hr+1
        if(n3hr>8)n3hr=1
      endif    ! (mod(ktau,nperday)==nper3hr(n3hr))

      if(ktau==ntau.or.mod(ktau,nperavg)==0)then
        dew_ave(:)   =   dew_ave(:)/min(ntau,nperavg)
        epan_ave(:)  =  epan_ave(:)/min(ntau,nperavg)
        epot_ave(:)  =  epot_ave(:)/min(ntau,nperavg)
        eg_ave(:)    =    eg_ave(:)/min(ntau,nperavg)
        fg_ave(:)    =    fg_ave(:)/min(ntau,nperavg)
        rnet_ave(:)  =  rnet_ave(:)/min(ntau,nperavg)
        sunhours(:)  =  sunhours(:)/min(ntau,nperavg)
        ga_ave(:)    =    ga_ave(:)/min(ntau,nperavg)
        riwp_ave(:)  =  riwp_ave(:)/min(ntau,nperavg)
        rlwp_ave(:)  =  rlwp_ave(:)/min(ntau,nperavg)
        tscr_ave(:)  =  tscr_ave(:)/min(ntau,nperavg)
        qscrn_ave(:) = qscrn_ave(:)/min(ntau,nperavg)
        do k=1,ms
          wb_ave(:,k)=wb_ave(:,k)/min(ntau,nperavg)
        end do
        tsu_ave(:)  = tsu_ave(:)/min(ntau,nperavg)
        psl_ave(:)  = psl_ave(:)/min(ntau,nperavg)
        mixdep_ave(:)=mixdep_ave(:)/min(ntau,nperavg)
        sgn_ave(:)  =  sgn_ave(:)/min(ntau,nperavg)  ! Dec07 because of solar fit
        if(myid==0)
     &    write(6,*) 'ktau,koundiag,nperavg =',ktau,koundiag,nperavg
        sint_ave(:) = sint_ave(:)/max(koundiag,1)
        sot_ave(:)  =  sot_ave(:)/max(koundiag,1)
        soc_ave(:)  =  soc_ave(:)/max(koundiag,1)
        sgdn_ave(:) =  sgdn_ave(:)/max(koundiag,1)
        rtu_ave(:)  =  rtu_ave(:)/max(koundiag,1)
        rtc_ave(:)  =  rtc_ave(:)/max(koundiag,1)
        rgdn_ave(:) =  rgdn_ave(:)/max(koundiag,1)
        rgn_ave(:)  =  rgn_ave(:)/max(koundiag,1)
        rgc_ave(:)  =  rgc_ave(:)/max(koundiag,1)
        cld_ave(:)  =  cld_ave(:)/max(koundiag,1)
        cll_ave(:)  =  cll_ave(:)/max(koundiag,1)
        clm_ave(:)  =  clm_ave(:)/max(koundiag,1)
        clh_ave(:)  =  clh_ave(:)/max(koundiag,1)
        alb_ave(:)  =  alb_ave(:)/max(koundiag,1)
        fbeam_ave(:)=  fbeam_ave(:)/max(koundiag,1)
        cbas_ave(:) = 1.1-cbas_ave(:)/max(1.e-4,precc(:))  ! 1.1 for no precc
        ctop_ave(:) = 1.1-ctop_ave(:)/max(1.e-4,precc(:))  ! 1.1 for no precc
        if (ngas>0) then
          traver=traver/min(ntau,nperavg)
        end if
        fpn_ave=fpn_ave/min(ntau,nperavg)
        frs_ave=frs_ave/min(ntau,nperavg)
        frp_ave=frp_ave/min(ntau,nperavg)
      end if    ! (ktau==ntau.or.mod(ktau,nperavg)==0)
      
      ! Update diagnostics for consistancy
      if (rescrn.gt.0) then
        call autoscrn
      end if
      
      if(ktau==ntau.or.mod(ktau,nwt)==0)then
        call log_off()
        call outfile(20,rundate,nmi,nwrite,iaero)  ! which calls outcdf
 
        if(ktau==ntau.and.irest==1) then
#ifdef simple_timer
          ! Don't include the time for writing the restart file
          call end_log(maincalc_end)
#endif
!         write restart file
          call outfile(19,rundate,nmi,nwrite,iaero)
          if(myid==0)
     &      write(6,*)'finished writing restart file in outfile'
#ifdef simple_timer
          call start_log(maincalc_begin)
#endif
        endif  ! (ktau==ntau.and.irest==1)
        call log_on()
      endif    ! (ktau==ntau.or.mod(ktau,nwt)==0)
      if(nstn>0.and.nrotstn(1)>0.and.mod(mtimer,60)==0)call stationb

      if(mod(ktau,nperavg)==0)then   
!       produce some diags & reset most averages once every nperavg      
        precavge=0.
        evapavge=0.
        do iq=1,ifull
         precavge=precavge+precip(iq)*wts(iq)
         evapavge=evapavge+evap(iq)*wts(iq)   ! in mm/day
        enddo
        pwatr=0.   ! in mm
        do k=1,kl
         do iq=1,ifull
          qtot=qg(iq,k)+qlg(iq,k)+qfg(iq,k)
          pwatr=pwatr-dsig(k)*wts(iq)*qtot*ps(iq)/grav
         enddo
        enddo
        temparray(1:3) = (/ precavge, evapavge, pwatr /)
        call MPI_Reduce ( temparray, gtemparray, 3, MPI_REAL, MPI_MAX,0,
     &                  MPI_COMM_WORLD, ierr )
        if ( myid == 0 ) then
           precavge = gtemparray(1)
           evapavge  = gtemparray(2)
           pwatr    = gtemparray(3)
           write(6,985) pwatr,precavge,evapavge ! MJT bug fix
        end if
985     format(' average pwatr,precc,prec,evap: ',4f7.3)
!       also zero most averaged fields every nperavg
        cbas_ave(:)=0.
        ctop_ave(:)=0.
        dew_ave(:)=0.
        epan_ave(:)=0.
        epot_ave(:)=0.
        eg_ave(:)=0.
        fg_ave(:)=0.
        rnet_ave(:)=0.
        sunhours(:)=0.
        riwp_ave(:)=0.
        rlwp_ave(:)=0.
        qscrn_ave(:) = 0.
        tscr_ave(:) = 0.
        wb_ave(:,:)=0.
        tsu_ave(:)=0.
        alb_ave(:)=0.
        fbeam_ave(:)=0.
        psl_ave(:)=0.
        mixdep_ave(:)=0.
        if (myid==0) write(6,*) 'resetting tscr_ave for ktau = ',ktau
        koundiag=0
        sint_ave(:) = 0.
        sot_ave(:)  = 0.
        soc_ave(:)  = 0.
        sgdn_ave(:) = 0.
        sgn_ave(:)  = 0.
        rtu_ave(:)  = 0.
        rtc_ave(:)  = 0.
        rgdn_ave(:) = 0.
        rgn_ave(:)  = 0.
        rgc_ave(:)  = 0.
        cld_ave(:)  = 0.
        cll_ave(:)  = 0.
        clm_ave(:)  = 0.
        clh_ave(:)  = 0.
!       zero evap, precip, precc, sno, runoff fields each nperavg (3/12/04) 
        evap(:)=0.  
        precip(:)=0.  ! converted to mm/day in outcdf
        precc(:)=0.   ! converted to mm/day in outcdf
        sno(:)=0.     ! converted to mm/day in outcdf
        runoff(:)=0.  ! converted to mm/day in outcdf
        if (ngas>0) traver=0.
        fpn_ave=0.
        frs_ave=0.
        frp_ave=0.
      endif  ! (mod(ktau,nperavg)==0)

      if(mod(ktau,nperday)==0)then   ! re-set at the end of each 24 hours
        if(ntau<10*nperday.and.nstn>0)then     ! print stn info
          do nn=1,nstn
           if ( .not. mystn(nn) ) cycle
           i=istn(nn)
           j=jstn(nn)
           iq=i+(j-1)*il
           write(6,956) ktau,iunp(nn),name_stn(nn),
     &      rnd_3hr(iq,4),rnd_3hr(iq,8),                  ! 12 hr & 24 hr
     &      tmaxscr(iq)-273.16+(zs(iq)/grav-zstn(nn))*stdlapse,
     &      tminscr(iq)-273.16+(zs(iq)/grav-zstn(nn))*stdlapse,
     &      tmaxscr(iq)-273.16,tminscr(iq)-273.16
956        format(i5,i3,a5,6f7.1)
          enddo
        endif  ! (ntau<10*nperday)
        rndmax (:) = 0.
        tmaxscr(:) = tscrn(:) 
        tminscr(:) = tscrn(:) 
        rhmaxscr(:) = rhscrn(:) 
        rhminscr(:) = rhscrn(:) 
        u10max(:)=0.
        v10max(:)=0.
        u1max(:)=0.
        v1max(:)=0.
        u2max(:)=0.
        v2max(:)=0.
        capemax(:)=0.
        rnd_3hr(:,8)=0.       ! i.e. rnd24(:)=0.
        if(nextout>=4)call setllp ! from Nov 11, reset once per day
        if(namip.ne.0)then ! not for last day, as day 1 of next month
          if (myid==0)
     &     write(6,*)
     &    'amipsst called at end of day for ktau,mtimer,namip ',
     &                                      ktau,mtimer,namip
          call amipsst
        endif ! (namip>0)
      elseif (namip.ne.0.and.nmlo.ne.0) then
        call amipsst
      endif   ! (mod(ktau,nperday)==0)

#ifdef vampir
      ! Flush vampir trace information to disk to save memory.
      call vtflush(ierr)
#endif
88    continue                   ! *** end of main time loop
#ifdef simple_timer
      call end_log(maincalc_end)
#endif
      call log_off()
      if (myid==0) then
         call date_and_time(time=timeval,values=tvals2)
         write(6,*) "End of time loop ", timeval
         write(6,*)'normal termination of run'
         call date_and_time(time=timeval)
         write(6,*) "End time ", timeval
         aa=3600*(tvals2(5)-tvals1(5)) + 
     &      60*(tvals2(6)-tvals1(6)) + (tvals2(7)-tvals1(7)) + 
     &      0.001 * (tvals2(8)-tvals1(8))
         if (aa.le.0.) aa=aa+86400.
         write(6,*) "Model time in main loop",aa
      end if
#ifdef simple_timer
      call end_log(model_end)
      call simple_timer_finalize
#endif

      call MPI_Finalize(ierr)

      end

      !--------------------------------------------------------------
      ! READ INTEGER TEXT FILES
      subroutine readint(filename,itss,ifullx)
      
      use cc_mpi            ! CC MPI routines
 
      implicit none
      
      include 'newmpar.h'   ! Grid parameters
      include 'parm.h'      ! Model configuration
      include 'parmgeom.h'  ! Coordinate data
            
      character *(*) filename
      character header*47
      integer ifullx,ilx,jlx,ierr
      integer itss(ifullx)
      integer glob2d(ifull_g)
      real rlong0x,rlat0x,schmidtx,dsx

      if( ifullx /= ifull ) then
        write(6,*) "Error, readint only works with ifull"
        write(6,*) "called with", ifullx
        stop
      end if
      if ( myid == 0 ) then
        write(6,*) 'reading data via readint from ',filename
        open(87,file=filename,status='old')
        read(87,*,iostat=ierr)
     &         ilx,jlx,rlong0x,rlat0x,schmidtx,dsx,header
        if ( ierr == 0 ) then
          write(6,*) ilx,jlx,rlong0x,rlat0x,schmidtx,dsx,header
          if(ilx.ne.il_g.or.jlx.ne.jl_g.or.rlong0x.ne.rlong0.
     &         or.rlat0x.ne.rlat0.or.schmidtx.ne.schmidt)
     &         stop 'wrong data file supplied'
          read(87,*) glob2d
          close(87)
        else if ( ierr < 0 ) then ! Error, so really unformatted file
          close(87)
          write(6,*) 'now doing unformatted read'
          open(87,file=filename,status='old',form='unformatted')
          read(87) glob2d
          close(87)
        else ! ierr > 0
          stop 'End of file occurred in readint'
        end if
        call ccmpi_distribute(itss, glob2d)
        write(6,*) trim(header), glob2d(id+(jd-1)*il_g)
      else
        call ccmpi_distribute(itss)
      end if
      end subroutine readint

      !--------------------------------------------------------------
      ! READ REAL TEXT FILES
      subroutine readreal(filename,tss,ifullx)
 
      use cc_mpi            ! CC MPI routines
 
      implicit none
      
      include 'newmpar.h'   ! Grid parameters
      include 'parm.h'      ! Model configuration
      include 'parmgeom.h'  ! Coordinate data

      character *(*) filename
      character header*47
      real tss(ifullx)
      real glob2d(ifull_g)
      real rlong0x,rlat0x,schmidtx,dsx
      integer ierr
      integer ilx,jlx,ifullx

      if( ifullx /= ifull ) then
        write(6,*) "Error, readreal only works with ifull"
        write(6,*) "called with", ifullx
        stop
      end if
      if ( myid == 0 ) then
        write(6,*) 'reading data via readreal from ',trim(filename)
        open(87,file=filename,status='old')
        read(87,*,iostat=ierr)
     &         ilx,jlx,rlong0x,rlat0x,schmidtx,dsx,header
        if ( ierr == 0 ) then
          write(6,*) ilx,jlx,rlong0x,rlat0x,schmidtx,dsx,header
          if(ilx.ne.il_g.or.jlx.ne.jl_g.or.rlong0x.ne.rlong0.
     &         or.rlat0x.ne.rlat0.or.schmidtx.ne.schmidt)
     &         stop 'wrong data file supplied'
          read(87,*) glob2d
          close(87)
        else if ( ierr < 0 ) then ! Error, so really unformatted file
          close(87)
          write(6,*) 'now doing unformatted read'
          open(87,file=filename,status='old',form='unformatted')
          read(87) glob2d
          close(87)
        else
          write(6,*) "error in readreal",trim(filename),ierr
          stop
        end if
        call ccmpi_distribute(tss, glob2d)
        write(6,*) trim(header), glob2d(id+(jd-1)*il_g)
      else
        call ccmpi_distribute(tss)
      end if
      end subroutine readreal


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
      
      integer iq,k
      
      do k=1,klt
       do iq=1,ilt*jlt        
        tr(iq,k,min(ntracmax,ngas+1))=rlatt(iq)*180./pi
        tr(iq,k,min(ntracmax,ngas+2))=rlongg(iq)*180./pi
        tr(iq,k,min(ntracmax,ngas+3))=.01*ps(iq)*sig(k)  ! in HPa
       enddo
      enddo
      if(nllp>=4)then   ! theta
        do k=1,klt
         do iq=1,ilt*jlt       
          tr(iq,k,min(ntracmax,ngas+4))=
     .	               t(iq,k)*(1.e-5*ps(iq)*sig(k))**(-rdry/cp)
         enddo
        enddo
      endif   ! (nllp>=4)
      if(nllp>=5)then   ! mixing_ratio (g/kg)
        do k=1,klt
         do iq=1,ilt*jlt       
          tr(iq,k,min(ntracmax,ngas+5))=1000.*qg(iq,k)
         enddo
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
      include 'parmvert.h'         ! Vertical advection parameters
      include 'soilv.h'            ! Soil parameters
      include 'stime.h'            ! File date data
      include 'trcom2.h'           ! Station data

      integer leap
      common/leap_yr/leap          ! Leap year (1 to allow leap years)
      integer nbarewet,nsigmf
      common/nsib/nbarewet,nsigmf  ! Land-surface options
      integer nnrad,idcld
      common/radnml/nnrad,idcld    ! Radiation options

!     for cardin
      data ia/1/,ib/3/,id/2/,ja/1/,jb/10/,jd/5/,nlv/1/,
     &     ndi/1/,ndi2/0/,nmaxpr/99/,          
     &     kdate_s/-1/,ktime_s/-1/,leap/0/,
     &     mbd/0/,nbd/0/,nbox/1/,kbotdav/4/,kbotu/0/,           
     &     nud_p/0/,nud_q/0/,nud_t/0/,nud_uv/1/,nud_hrs/24/,nudu_hrs/0/,
     &     ktopdav/-1/,nud_sst/0/,nud_sss/0/,kbotmlo/-1/,ktopmlo/1/,
     &     mloalpha/10/,nud_ouv/0/,nud_sfh/0/
      
!     Dynamics options A & B      
      data m/5/,mex/30/,mfix/3/,mfix_qg/1/,mup/1/,nh/0/,nonl/0/,npex/0/
      data nritch_t/300/,nrot/1/,nxmap/0/,
     &     epsp/-15./,epsu/0./,epsf/0./,precon/-2900/,restol/4.e-7/
      data helmmeth/0/
      data mfix_tr/0/,mfix_ke/0/,mfix_aero/0/
      data schmidt/1./,rlong0/0./,rlat0/90./,nrun/0/,nrunx/0/
!     Horiz advection options
      data ndept/1/,nt_adv/7/,mh_bs/4/
!     Horiz wind staggering options
c     data nstag/99/,nstagu/99/
      data nstag/-10/,nstagu/-1/
!     Vertical advection options
      data nvad/-4/,nvadh/2/,ntvdr/1/
!     Horizontal mixing options (now in initialparm)
      !data khdif/2/,khor/-8/,nhor/-157/,nhorps/-1/,nhorjlm/1/
!     Vertical mixing options
      data nvmix/3/,nlocal/6/,nvsplit/2/,ncvmix/0/,lgwd/0/,ngwd/-5/
!     Cumulus convection options
      data nkuo/23/,sigcb/1./,sig_ct/1./,rhcv/0./,rhmois/.1/,rhsat/1./,
     &     convfact/1.02/,convtime/.33/,shaltime/0./,      
     &     alflnd/1.1/,alfsea/1.1/,fldown/.6/,iterconv/3/,ncvcloud/0/,
     &     nevapcc/0/,nevapls/-4/,nuvconv/0/
     &     mbase/101/,mdelay/-1/,methprec/8/,nbase/-4/,detrain/.15/,
     &     entrain/.05/,methdetr/2/,detrainx/0./,dsig2/.15/,dsig4/.4/
!     Shallow convection options
      data ksc/-95/,kscsea/0/,kscmom/1/,sigkscb/.95/,sigksct/.8/,
     &     tied_con/2./,tied_over/0./,tied_rh/.75/
!     Other moist physics options
      data acon/.2/,bcon/.07/,qgmin/1.e-6/,rcm/.92e-5/,
     &     rcrit_l/.75/,rcrit_s/.85/ 
!     Radiation options
      data nrad/4/,ndiur/1/,idcld/1/
!     Diagnostic cloud options
      data ldr/1/,nclddia/1/,nstab_cld/0/,nrhcrit/10/,sigcll/.95/ 
      data cldh_lnd/95./,cldm_lnd/85./,cldl_lnd/75./
      data cldh_sea/95./,cldm_sea/90./,cldl_sea/80./
!     Soil, canopy, PBL options
      data nbarewet/0/,newrough/0/,nglacier/1/,
     &     nrungcm/-1/,nsib/3/,nsigmf/1/,
     &     ntaft/2/,ntsea/6/,ntsur/6/,av_vmod/.7/,tss_sh/1./,
     &     vmodmin/.2/,zobgin/.02/,charnock/.018/,chn10/.00125/
      data newsoilm/0/,newztsea/1/,newtop/1/,nem/2/                    
      data snmin/.11/  ! 1000. for 1-layer; ~.11 to turn on 3-layer snow
      data nurban/0/,nmr/0/,bpyear/0./,rescrn/0/
!     Special and test options
      data namip/0/,amipo3/.false./,nhstest/0/,nsemble/0/,nspecial/0/,
     &     panfg/4./,panzo/.001/,nplens/0/
!     I/O options
      data m_fly/4/,io_in/1/,io_out/1/,io_rest/1/
      data nperavg/-99/,nwt/-99/
      data io_clim/1/,io_spec/0/,nextout/3/,localhist/.false./
      data nstn/0/  
      data slat/nstnmax*-89./,slon/nstnmax*0./,iunp/nstnmax*6/,
     &     zstn/nstnmax*0./,name_stn/nstnmax*'   '/,nrotstn/nstnmax*0/  
!     Ocean options
      data nmlo/0/

c     initialize file names to something
      data albfile/' '/,icefile/' '/,maskfile/' '/
     &    ,snowfile/' '/,sstfile/' '/,topofile/' '/,zofile/' '/
     &    ,rsmfile/' '/,scamfile/' '/,soilfile/' '/,vegfile/' '/
     &    ,co2emfile/' '/,so4tfile/' '/
     &    ,smoistfile/' '/,soil2file/' '/,restfile/' '/
     &    ,radonemfile/' '/,surfile/' '/,surf_00/'s_00a '/
     &    ,surf_12/'s_12a '/,co2_00/' '/,co2_12/' '/,radon_00/' '/
     &    ,radon_12/' '/,ifile/' '/,ofile/' '/,nmifile/' '/
     &    ,eigenv/' '/,radfile/' '/,o3file/' '/,hfile/' '/,mesonest/' '/
     &          ,scrnfile/' '/,tmaxfile/' '/,tminfile/' '/,trcfil/' '/
     &    ,laifile/' '/,albnirfile/' '/,urbanfile/' '/,bathfile/' '/
     &    ,vegprev/' '/,vegnext/' '/,cnsdir/' '/,salfile/' '/
     &    ,oxidantfile/' '/
      data climcdf/'clim.cdf'/
      data monfil/'monthly.cdf'/,scrfcdf/'scrave.cdf'/
c     floating point:
      data timer/0./,mtimer/0/

c     stuff from insoil  for soilv.h
      data rlaim44/4.8, 6.3, 5., 3.75, 2.78, 2.5, 3.9, 2.77, 2.04, 2.6,  ! 1-10
     &          1.69, 1.9, 1.37, 1.5, 1.21, 1.58, 1.41, 2.3, 1.2, 1.71,  ! 11-20
     &          1.21, 2.3, 2.3, 1.2, 1.2, 1.87, 1., 3., .01, .01, 1.2,   ! 21-31
     &          6., 5.5, 5., 4.5, 5., 4., 3., 3.5, 1., 4., .5, 4., 0./   ! 32-44
      data rlais44/1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,              ! 1-10
     &           1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,                ! 11-20
     &           1., 1., 1., 1., .6, .6, .5, 1., 0., 0., 1.,            ! 21-31
     &           2., 2., 2., 2., 2., 1.5, 1.5, 1.5, 1., .5, .5, .5, 0./ ! 32-44
      data rsunc44/
     &      370., 330., 260., 200., 150., 130., 200., 150., 110., 160.,  ! 1-10
     &      100., 120.,  90.,  90.,  80.,  90.,  90., 150.,  80., 100.,  ! 11-20
     &       80.,  80.,  80.,  60.,  60., 120.,  80., 180., 2*995., 80., ! 21-31
     &      350., 4*300., 3*230., 150., 230., 995., 150., 9900./         ! 32-44
      data scveg44/0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,              ! 1-10
     &             0., 0., 0., 0., 0., .1, .1, .1, .1, .1,              ! 11-20
     &             .1, .2, .4, .2, .1, .1, .1, 0., 0., 0., 0.,          ! 21-31
     &         .05, 0., 0., 0., 0., .05, .05, .05, .1, 0., 0., .4, 0./  ! 32-44
      data slveg44/0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,           ! 1-10
     &             0., 0., 0., 0., 0., .1, .1, .1, .1, .1,           ! 11-20
     &             .1, .2, .4, .2, .1, .1, .1, 0., 0., 0., 0.,       ! 21-31
     &     1., 5.5, 3., 1., 3., 3., 3.5, 3., .5, 3.5, .1, 3.5, 0./   ! 32-44
      data froot/.05, .10, .35, .40, .10/       ! 10/02/99 veg. root distr.

      data silt/.08, .33, .17, .2, .06, .25, .15, .70, .33
     &          , .2, .33, .33, .17/                        ! with mxst=13
      data clay/.09, .3, .67, .2, .42, .48, .27, .17, .30
     &          , .2, .3, .3, .67/                          ! with mxst=13
      data sand/.83, .37, .17, .6, .52, .27, .58, .13, .37
     &          , .6, .37, .37, .17/                        ! with mxst=13
      data swilt/0., .072, .216, .286, .135, .219, .283, .175, .395  !eak
     &             , .216, .1142, .1547, .2864, .2498/
      data sfc/1.,  .143, .301, .367, .218, .31 , .37 , .255, .45 
     &            , .301, .22 , .25 , .367, .294/
      data ssat/2., .398, .479, .482, .443, .426, .482, .420, .451  !jlm
     &          , .479, .435, .451, .482, .476/
      data bch/4.2, 7.1, 11.4, 5.15, 10.4, 10.4, 7.12, 5.83, 7.1, 4.9
     &          , 5.39, 11.4, 8.52/      ! bch for gravity term
      data hyds/166.e-6, 4.e-6, 1.e-6, 21.e-6, 2.e-6, 1.e-6, 6.e-6,
     *             800.e-6, 1.e-6, 34.e-6, 7.e-6, 1.3e-6, 2.5e-6/  ! ksat
      data sucs/-.106, -.591, -.405, -.348, -.153, -.49, -.299,
     &          -.356, -.153, -.218, -.478, -.405, -.63/           ! phisat (m)
      data rhos/1600., 1595., 1381., 1373., 1476., 1521., 1373., 1537.,
     &          1455., 4*2600. /
      ! soil density changed to the line above by YP using the relationship
      ! rho = (1  - ssat) * 2650  ----- (3Nov2007)
      data  css/7* 850., 1920., 2100., 4*850./       ! heat capacity

      data zse/.022, .058, .154, .409, 1.085, 2.872/ ! layer thickness
!     so depths of centre of layers: .011, .051, .157, .4385, 1.1855, 3.164
!     with base at 4.6     

      end
      
      ! Initialise initial values that cannot be definined in main_blockdata
      subroutine initialparm
      
      use parmhdff_m ! Horizontal diffusion parameters
      
      implicit none
      
      ! Horizontal diffusion parameters
      nhor=-157
      nhorps=-1
      khor=-8
      khdif=2
      nhorjlm=1
      hdifmax=0.
      
      return
      end subroutine initialparm
      
      !--------------------------------------------------------------
      ! WRITE STATION DATA
      subroutine stationa

      use arrays_m           ! Atmosphere dyamics prognostic arrays
      use cc_mpi             ! CC MPI routines
      use diag_m             ! Diagnostic routines
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
      include 'establ.h'     ! Liquid saturation function
      include 'parm.h'       ! Model configuration
      include 'parmgeom.h'   ! Coordinate data
      include 'soilv.h'      ! Soil parameters
      include 'trcom2.h'     ! Station data

      integer leap
      common/leap_yr/leap    ! Leap year (1 to allow leap years)

      integer i,j,iq,iqt,isoil,k2,nn
      real coslong, sinlong, coslat, sinlat, polenx, poleny, polenz,
     &     zonx, zony, zonz, den, costh, sinth, uzon, vmer, rh1, rh2,
     &     es,wbav, rh_s

      coslong=cos(rlong0*pi/180.)   ! done here, where work2 has arrays
      sinlong=sin(rlong0*pi/180.)
      coslat=cos(rlat0*pi/180.)
      sinlat=sin(rlat0*pi/180.)
      polenx=-coslat
      poleny=0.
      polenz=sinlat
      do nn=1,nstn
!       Check if this station is in this processors region
        if ( .not. mystn(nn) ) cycle 
        if(ktau==1)write (iunp(nn),950) kdate,ktime,leap
950     format("#",i9,2i5)
        i=istn(nn)
        j=jstn(nn)
        iq=i+(j-1)*il
        zonx=            -polenz*y(iq)
        zony=polenz*x(iq)-polenx*z(iq)
        zonz=polenx*y(iq)
        den=sqrt( max(zonx**2+zony**2+zonz**2,1.e-7) ) 
        costh= (zonx*ax(iq)+zony*ay(iq)+zonz*az(iq))/den
        sinth=-(zonx*bx(iq)+zony*by(iq)+zonz*bz(iq))/den
        uzon= costh*u(iq,1)-sinth*v(iq,1)
        vmer= sinth*u(iq,1)+costh*v(iq,1)
        es=establ(t(iq,1))
        rh1=100.*qg(iq,1)*(ps(iq)*sig(1)-es)/(.622*es)
        es=establ(t(iq,2))
        rh2=100.*qg(iq,2)*(ps(iq)*sig(2)-es)/(.622*es)
        es=establ(tscrn(iq))
        rh_s=100.*qgscrn(iq)*(ps(iq)-es)/(.622*es)
        wbav=(zse(1)*wb(iq,1)+zse(2)*wb(iq,2)+zse(3)*wb(iq,3)
     &       +zse(4)*wb(iq,4))/(zse(1)+zse(2)+zse(3)+zse(4))
        iqt = min(iq,ilt*jlt) ! Avoid bounds problems if there are no tracers
        k2=min(2,klt)
        write (iunp(nn),951) ktau,tscrn(iq)-273.16,rnd_3hr(iq,8),
     &         tss(iq)-273.16,tgg(iq,1)-273.16,tgg(iq,2)-273.16,
     &         tgg(iq,3)-273.16,t(iq,1)-273.16,0.,
     &         wb(iq,1),wb(iq,2),
     &         cloudlo(iq),cloudmi(iq)+1.,cloudhi(iq)+2.,
     &         cloudtot(iq)+3.,
     &         fg(iq),eg(iq),0.,0.,rnet(iq),sgsave(iq),
     &         qg(iq,1)*1.e3,uzon,vmer,precc(iq),
     &         qg(iq,2)*1.e3,rh1,rh2,0.,0.,
     &         0.,0.,.01*ps(iq),wbav,
     &         epot(iq),qgscrn(iq)*1.e3,rh_s,u10(iq),uscrn(iq),
     &         condx(iq)              
!       N.B. qgscrn formula needs to be greatly improved
951     format(i4,6f7.2, 
     &         2f7.2, 2f6.3, 4f5.2,                ! t1 ... cld
     &         5f7.1,f6.1,f5.1,                    ! fg ... qg1
     &         2f6.1,f7.2, f5.1,2f6.1, 2(1x,f5.1), ! uu ... co2_2
     &         2(1x,f5.1) ,f7.1,f6.3,f7.1,5f6.1,   ! rad_1 ... rh_s
     &         f7.2)                               ! condx
        if(ktau==ntau)then
          write (iunp(nn),952)
952       format("#   2tscrn 3precip 4tss  5tgg1  6tgg2  7tgg3",
     &       "   8t1    9tgf  10wb1 11wb2 cldl cldm cldh  cld",
     &       "   16fg   17eg  18fgg  19egg  20rnet 21sg 22qg1",
     &       " 23uu   24vv 25precc qg2  rh1 28rh2 29co2_1 co2_2",
     &       " rad_1 rad_2  ps 34wbav 35epot qgscrn 37rh_s 38u10 uscrn",
     &       " 40condx")
          write (iunp(nn),953) land(iq),isoilm(iq),ivegt(iq),zo(iq),
     &                         zs(iq)/grav
953       format("# land,isoilm,ivegt,zo,zs/g: ",l2,2i3,2f9.3)
          isoil=max(1,isoilm(iq))
          write (iunp(nn),954) sigmf(iq),swilt(isoil),sfc(isoil),
     &                         ssat(isoil),0.5*sum(albvisnir(iq,:))
954       format("#sigmf,swilt,sfc,ssat,alb: ",5f7.3)
        endif
      enddo
      return	    
      end               

      !--------------------------------------------------------------
      ! WRITE STATION DATA
      subroutine stationb  ! primarily for ICTS

      use arrays_m            ! Atmosphere dyamics prognostic arrays
      use cc_mpi              ! CC MPI routines
      use cfrac_m             ! Cloud fraction
      use diag_m              ! Diagnostic routines
      use extraout_m          ! Additional diagnostics
      use histave_m           ! Time average arrays
      use infile              ! Input file routines
      use liqwpar_m           ! Cloud water mixing ratios
      use map_m               ! Grid map arrays
      use morepbl_m           ! Additional boundary layer diagnostics
      use nsibd_m             ! Land-surface arrays
      use pbl_m               ! Boundary layer arrays
      use prec_m              ! Precipitation
      use raddiag_m           ! Radiation diagnostic
      use screen_m            ! Screen level diagnostics
      use sigs_m              ! Atmosphere sigma levels
      use soil_m              ! Soil and surface data
      use soilsnow_m          ! Soil, snow and surface data
      use tracers_m           ! Tracer data
      use vecsuv_m            ! Map to cartesian coordinates
      use xyzinfo_m           ! Grid coordinate arrays

      implicit none

      include 'newmpar.h'     ! Grid parameters
      include 'const_phys.h'  ! Physical constants
      include 'dates.h'       ! Date data
      include 'establ.h'      ! Liquid saturation function
      include 'parm.h'        ! Model configuration
      include 'parmgeom.h'    ! Coordinate data
      include 'soilv.h'       ! Soil parameters
      include 'trcom2.h'      ! Station data

      real p(ifull,kl),tv(ifull,kl),energy(ifull,kl)
      real pmsl(ifull)
      integer i,j,iq,k,nn,kdateb,ktimeb,mtimerb,npres,iyr,imo,iday,itim
      integer niq,nr,iqq(9)
      real costh(9),sinth(9),energint(9),pwater(9),off,sca
      real uzon(9,kl),vmer(9,kl)
      real tt(9,4),uu(9,4),vv(9,4),pp(9,4),qgg(9,4)  ! in subr. stationb
      real press(kl),pres(4)
      data pres/850.,700.,500.,250./
      real coslong, sinlong, coslat, sinlat, polenx, poleny, polenz,
     &     zonx, zony, zonz, den, 
     &     qtot,div,div_inte,div_intq,fa,fb
      integer ix(15),jx(15)
      data ix/-1,0,1,1,1,0,-1,-1, -1,0,1,1,1,0,-1/  ! clockwise from NW
      data jx/1,1,1,0,-1,-1,-1,0, 1,1,1,0,-1,-1,-1/ ! clockwise from NW

      call mslp(pmsl,psl,zs,t(1:ifull,:))
      off=0.
      sca=1.
      tv(:,:)=t(1:ifull,:)+(.61*qg(1:ifull,:)-qfg(1:ifull,:)
     &                             -qlg(1:ifull,:))*t(1:ifull,:)  
      p(1:ifull,1)=zs(1:ifull)+bet(1)*tv(1:ifull,1)
      do k=2,kl
        p(1:ifull,k)=p(1:ifull,k-1)
     &            +bet(k)*tv(1:ifull,k)+betm(k)*tv(1:ifull,k-1)
      enddo    ! k  loop
      energy(:,:)=(cp-rdry)*t(1:ifull,:)+p(:,:)+
     &            .5*(u(1:ifull,:)**2+v(1:ifull,:)**2)
      coslong=cos(rlong0*pi/180.)   
      sinlong=sin(rlong0*pi/180.)
      coslat=cos(rlat0*pi/180.)
      sinlat=sin(rlat0*pi/180.)
      polenx=-coslat
      poleny=0.
      polenz=sinlat
      kdateb=kdate
      ktimeb=ktime
      mtimerb=mtimer
!     could alter the following to update kdateb more elegantly	    
      call datefix(kdateb,ktimeb,mtimerb)
      iyr=kdateb/10000
      imo=(kdateb-iyr*10000)/100
      iday=kdateb-iyr*10000-imo*100
      itim=ktimeb/100
      write(6,*) 'new kdateb,ktimeb,mtimerb ',kdateb,ktimeb,mtimerb
      write(6,*) 'new iyr,imo,iday,itim,off,sca ',iyr,imo,iday,itim,
     &                                            off,sca
      do nn=1,nstn
!      Check if this station is in this processors region
       if ( .not. mystn(nn) ) cycle 
         i=istn(nn)
         j=jstn(nn)
         nr=nrotstn(nn)
         iqq(1)=i+(j-1)*il
         iqq(2)=i+ix(nr)+(j+jx(nr)-1)*il
         iqq(3)=i+ix(nr+1)+(j+jx(nr+1)-1)*il
         iqq(4)=i+ix(nr+2)+(j+jx(nr+2)-1)*il
         iqq(5)=i+ix(nr+3)+(j+jx(nr+3)-1)*il
         iqq(6)=i+ix(nr+4)+(j+jx(nr+4)-1)*il
         iqq(7)=i+ix(nr+5)+(j+jx(nr+5)-1)*il
         iqq(8)=i+ix(nr+6)+(j+jx(nr+6)-1)*il
         iqq(9)=i+ix(nr+7)+(j+jx(nr+7)-1)*il
         do niq=1,9
          iq=iqq(niq)
          zonx=            -polenz*y(iq)
          zony=polenz*x(iq)-polenx*z(iq)
          zonz=polenx*y(iq)
          den=sqrt( max(zonx**2+zony**2+zonz**2,1.e-7) )  ! allow for poles
          costh(niq)= (zonx*ax(iq)+zony*ay(iq)+zonz*az(iq))/den
          sinth(niq)=-(zonx*bx(iq)+zony*by(iq)+zonz*bz(iq))/den
          do k=1,kl
           uzon(niq,k)=costh(niq)*u(iq,k)-sinth(niq)*v(iq,k)    
           vmer(niq,k)=sinth(niq)*u(iq,k)+costh(niq)*v(iq,k)    
          enddo
         enddo     ! niq=1,9
         do k=1,kl
          write (iunp(nn),"('T         ',i6,3i3,i4,f8.0,g8.1,9f8.3)") 
     &	         iyr,imo,iday,itim,k,off,sca,t(iqq(:),k)
          write (iunp(nn),"('U         ',i6,3i3,i4,f8.0,g8.1,9f8.3)") 
     &	         iyr,imo,iday,itim,k,off,sca,uzon(:,k)
          write (iunp(nn),"('V         ',i6,3i3,i4,f8.0,g8.1,9f8.3)") 
     &	         iyr,imo,iday,itim,k,off,sca,vmer(:,k)
          write (iunp(nn),"('GPH       ',i6,3i3,i4,f8.0,g8.1,9f8.0)") 
     &	         iyr,imo,iday,itim,k,off,sca,p(iqq(:),k)
          write (iunp(nn),"('QV        ',i6,3i3,i4,f8.0,g8.1,9f8.5)") 
     &	         iyr,imo,iday,itim,k,off,sca,qg(iqq(:),k)
          write (iunp(nn),"('QC        ',i6,3i3,i4,f8.0,g8.1,9f8.5)") 
     &	         iyr,imo,iday,itim,k,off,sca,qlg(iqq(:),k)
          write (iunp(nn),"('QI        ',i6,3i3,i4,f8.0,g8.1,9f8.5)") 
     &	         iyr,imo,iday,itim,k,off,sca,qfg(iqq(:),k)
          write (iunp(nn),"('CLC       ',i6,3i3,i4,f8.0,g8.1,9f8.3)") 
     &	         iyr,imo,iday,itim,k,off,sca,cfrac(iqq(:),k)
         enddo   ! k  loop
!        default values for 850, 700, 500, 250 (in case below ground)	  
         uu(:,:)=-99.
         vv(:,:)=-99.
         tt(:,:)=-99.
         pp(:,:)=-99.
         qgg(:,:)=-.99
         do niq=1,9
          iq=iqq(niq)
          do k=1,kl
           press(k)=.01*sig(k)*ps(iq)  ! in hPa
          enddo
          k=kl-1
          npres=4
22        k=k-1
          if(pres(npres)<press(k-1))then  ! e.g. 250>press(k)
!	    just use linear interpolation
            fa=(press(k-1)-pres(npres))/(press(k-1)-press(k))
            fb=1.-fa
            tt(niq,npres)=fa*t(iq,k)+fb*t(iq,k-1)
            uu(niq,npres)=fa*uzon(niq,k)+fb*uzon(niq,k-1)
            vv(niq,npres)=fa*vmer(niq,k)+fb*vmer(niq,k-1)
            pp(niq,npres)=fa*p(iq,k)+fb*p(iq,k-1)
            qgg(niq,npres)=fa*qg(iq,k)+fb*qg(iq,k-1)
            npres=npres-1
          endif
          if(k>2.and.npres>0)go to 22  
         enddo   ! niq loop
         write (iunp(nn),"('T         ',i6,3i3,' 850',f8.0,g8.1,9f8.3)") 
     &	  iyr,imo,iday,itim,off,sca,tt(:,1)
         write (iunp(nn),"('U         ',i6,3i3,' 850',f8.0,g8.1,9f8.3)") 
     &	  iyr,imo,iday,itim,off,sca,uu(:,1)
         write (iunp(nn),"('V         ',i6,3i3,' 850',f8.0,g8.1,9f8.3)") 
     &	  iyr,imo,iday,itim,off,sca,vv(:,1)
         write (iunp(nn),"('GPH       ',i6,3i3,' 850',f8.0,g8.1,9f8.0)")
     &	  iyr,imo,iday,itim,off,sca,pp(:,1)
         write (iunp(nn),"('QV        ',i6,3i3,' 850',f8.0,g8.1,9f8.5)")
     &	  iyr,imo,iday,itim,off,sca,qgg(:,1)
         write (iunp(nn),"('T         ',i6,3i3,' 700',f8.0,g8.1,9f8.3)") 
     &	  iyr,imo,iday,itim,off,sca,tt(:,2)
         write (iunp(nn),"('U         ',i6,3i3,' 700',f8.0,g8.1,9f8.3)") 
     &	  iyr,imo,iday,itim,off,sca,uu(:,2)
         write (iunp(nn),"('V         ',i6,3i3,' 700',f8.0,g8.1,9f8.3)") 
     &	  iyr,imo,iday,itim,off,sca,vv(:,2)
         write (iunp(nn),"('GPH       ',i6,3i3,' 700',f8.0,g8.1,9f8.0)") 
     &	  iyr,imo,iday,itim,off,sca,pp(:,2)
         write (iunp(nn),"('QV        ',i6,3i3,' 700',f8.0,g8.1,9f8.5)") 
     &	  iyr,imo,iday,itim,off,sca,qgg(:,2)
         write (iunp(nn),"('T         ',i6,3i3,' 500',f8.0,g8.1,9f8.3)") 
     &	  iyr,imo,iday,itim,off,sca,tt(:,3)
         write (iunp(nn),"('U         ',i6,3i3,' 500',f8.0,g8.1,9f8.3)") 
     &	  iyr,imo,iday,itim,off,sca,uu(:,3)
         write (iunp(nn),"('V         ',i6,3i3,' 500',f8.0,g8.1,9f8.3)") 
     &	  iyr,imo,iday,itim,off,sca,vv(:,3)
         write (iunp(nn),"('GPH       ',i6,3i3,' 500',f8.0,g8.1,9f8.0)")
     &	  iyr,imo,iday,itim,off,sca,pp(:,3)
         write (iunp(nn),"('QV        ',i6,3i3,' 500',f8.0,g8.1,9f8.5)")
     &	  iyr,imo,iday,itim,off,sca,qgg(:,3)
         write (iunp(nn),"('T         ',i6,3i3,' 250',f8.0,g8.1,9f8.3)") 
     &	  iyr,imo,iday,itim,off,sca,tt(:,4)
         write (iunp(nn),"('U         ',i6,3i3,' 250',f8.0,g8.1,9f8.3)") 
     &	  iyr,imo,iday,itim,off,sca,uu(:,4)
         write (iunp(nn),"('V         ',i6,3i3,' 250',f8.0,g8.1,9f8.3)") 
     &	  iyr,imo,iday,itim,off,sca,vv(:,4)
         write (iunp(nn),"('GPH       ',i6,3i3,' 250',f8.0,g8.1,9f8.0)")
     &	  iyr,imo,iday,itim,off,sca,pp(:,4)
         write (iunp(nn),"('QV        ',i6,3i3,' 250',f8.0,g8.1,9f8.5)") 
     &	  iyr,imo,iday,itim,off,sca,qgg(:,4)

         do niq=1,9
          iq=iqq(niq)
          pwater(niq)=0.   ! in mm
          energint(niq)=0.   ! in mm
          do k=1,kl
           qtot=qg(iq,k)+qlg(iq,k)+qfg(iq,k)
           pwater(niq)=pwater(niq)-dsig(k)*qtot*ps(iq)/grav
           energint(niq)=energint(niq)-dsig(k)*energy(iq,k)*ps(iq)/grav
          enddo  ! k   loop
         enddo   ! niq loop
!        integ div terms, just for niq=1	  
         iq=iqq(1)
         div_intq=0.
         div_inte=0.
         do k=1,kl
          div=(qg(iq+1,k)*u(iq+1,k)/emu(iq+1)
     &        -qg(iq-1,k)*u(iq-1,k)/emu(iq-1)  
     &        +qg(iq+il,k)*v(iq+il,k)/emv(iq+il)
     &        -qg(iq-il,k)*v(iq-il,k)/emv(iq-il)) 
     &              *em(iq)**2/(2.*ds)  *1.e6
          div_intq=div_intq-div*dsig(k)*ps(iq)/grav
          div=(energy(iq+1,k)*u(iq+1,k)/emu(iq+1)
     &        -energy(iq-1,k)*u(iq-1,k)/emu(iq-1)  
     &        +energy(iq+il,k)*v(iq+il,k)/emv(iq+il)
     &        -energy(iq-il,k)*v(iq-il,k)/emv(iq-il)) 
     &              *em(iq)**2/(2.*ds)  
          div_inte=div_inte-div*dsig(k)*ps(iq)/grav
        enddo  ! k   loop
         write (iunp(nn),"('CLCT      ',i6,3i3,'   1',f8.0,g8.1,9f8.3)") 
     &	         iyr,imo,iday,itim,off,sca,cloudtot(iqq(:)) 
         write (iunp(nn),"('DZ_PBL    ',i6,3i3,'   1',f8.0,g8.1,9f8.1)") 
     &	         iyr,imo,iday,itim,off,sca,pblh(iqq(:))
         write (iunp(nn),"('EVAP_S    ',i6,3i3,'   1',f8.0,g8.1,9f8.3)") 
     &	         iyr,imo,iday,itim,off,sca,evap(iqq(:))
         write(iunp(nn),"('IDIV_ENERG',i6,3i3,'   1',f8.0,2f8.0,8g8.1)") 
     &	         iyr,imo,iday,itim,off,sca,div_inte,(1.e20,k=2,9) ! niq=1 
         write(iunp(nn),"('IDIV_WATER',i6,3i3,'   1',f8.0,2f8.1,8g8.1)") 
     &	         iyr,imo,iday,itim,off,sca,div_intq,(1.e20,k=2,9) ! niq=1 
         sca=1.e9
         write(iunp(nn),"('IENERGY ',i8,3i3,'   1',f8.0,g8.1,-9p9f8.5)")
     &	         iyr,imo,iday,itim,off,sca,energint(:)
         sca=1.
         write (iunp(nn),"('ISOILW    ',i6,3i3,'   1',f8.0,g8.1,9f8.3)")
     &	   iyr,imo,iday,itim,off,sca,
     &                  wb(iqq(:),1)*zse(1)+wb(iqq(:),2)*zse(2)+
     &                  wb(iqq(:),3)*zse(3)+wb(iqq(:),4)*zse(4)+
     &                  wb(iqq(:),5)*zse(5)+wb(iqq(:),6)*zse(6)  !  check
         write (iunp(nn),"('IWATER    ',i6,3i3,'   1',f8.0,g8.1,9f8.3)")
     &	         iyr,imo,iday,itim,off,sca,pwater(:)
         write (iunp(nn),"('PMSL      ',i6,3i3,'   1',f8.0,g8.1,9f8.2)")
     &	         iyr,imo,iday,itim,off,sca,pmsl(iqq(:))*.01
         write (iunp(nn),"('PS        ',i6,3i3,'   1',f8.0,g8.1,9f8.2)") 
     &	         iyr,imo,iday,itim,off,sca,ps(iqq(:))*.01
         write (iunp(nn),"('QV_2M     ',i6,3i3,'   1',f8.0,g8.1,9f8.5)") 
     &	         iyr,imo,iday,itim,off,sca,qgscrn(iqq(:))
         write (iunp(nn),"('RELHUM_2M ',i6,3i3,'   1',f8.0,g8.1,9f8.2)")
     &	         iyr,imo,iday,itim,off,sca,rhscrn(iqq(:))
c     &	         rh_s(:)
         write (iunp(nn),"('RUNOFF    ',i6,3i3,'   1',f8.0,g8.1,9f8.3)")
     &	         iyr,imo,iday,itim,off,sca,runoff(iqq(:))
         write (iunp(nn),"('T_2M      ',i6,3i3,'   1',f8.0,g8.1,9f8.2)") 
     &	         iyr,imo,iday,itim,off,sca,tscrn(iqq(:))
         write (iunp(nn),"('T_SKIN    ',i6,3i3,'   1',f8.0,g8.1,9f8.2)") 
     &	         iyr,imo,iday,itim,off,sca,tss(iqq(:))
         write (iunp(nn),"('TMAX_2M   ',i6,3i3,'   1',f8.0,g8.1,9f8.2)") 
     &	         iyr,imo,iday,itim,off,sca,tmaxscr(iqq(:))
         write (iunp(nn),"('TMIN_2M   ',i6,3i3,'   1',f8.0,g8.1,9f8.2)") 
     &	         iyr,imo,iday,itim,off,sca,tminscr(iqq(:))
         write (iunp(nn),"('TOT_PREC  ',i6,3i3,'   1',f8.0,g8.1,9f8.2)")
     &	         iyr,imo,iday,itim,off,sca,precip(iqq(:))
         write (iunp(nn),"('VABS_10M  ',i6,3i3,'   1',f8.0,g8.1,9f8.3)")
     &	         iyr,imo,iday,itim,off,sca,u10(iqq(:))
         write (iunp(nn),"('W_SNOW    ',i6,3i3,'   1',f8.0,g8.1,9f8.3)") 
     &	         iyr,imo,iday,itim,off,sca,snowd(iqq(:))*.001  ! converts mm to m
         write (iunp(nn),"('SHFL_S    ',i6,3i3,'   1',f8.0,g8.1,9f8.2)")
     &	         iyr,imo,iday,itim,off,sca,fg(iqq(:))
         write (iunp(nn),"('LHFL_S    ',i6,3i3,'   1',f8.0,g8.1,9f8.2)")
     &	         iyr,imo,iday,itim,off,sca,eg(iqq(:))
         write (iunp(nn),"('SMHFL     ',i6,3i3,'   1',f8.0,g8.1,9f8.2)") 
     &	         iyr,imo,iday,itim,off,sca,snowflx(iqq(:))
      enddo   ! nn=1,nstn
      return    
      end
      
      !--------------------------------------------------------------
      ! PREVIOUS VERSION DEFAULT PARAMETERS
      subroutine change_defaults(nversion)

      use parmhdff_m              ! Horizontal diffusion parameters

      implicit none

      include 'newmpar.h'         ! Grid parameters
      include 'kuocom.h'          ! Convection parameters
      include 'parm.h'            ! Model configuration
      include 'parmdyn.h'         ! Dynamics parmaters
      include 'parmhor.h'         ! Horizontal advection parameters
      include 'parmsurf.h'        ! Surface parameters
      include 'parmvert.h'        ! Vertical advection parameters

      integer nbarewet,nsigmf
      common/nsib/nbarewet,nsigmf ! Land-surface options

      integer nversion

      if(nversion<907)then
        mfix=1         ! new is 3
        newrough=2     ! new is 0
        newtop=0       ! new is 1
        nvmix=5        ! new is 3
        ksc=0          ! new is -95
        sig_ct=.8      ! new is 1.
      endif
      if(nversion<904)then
        newtop=1       ! new is 0
        nvmix=3        ! new is 5
        ksc=-95        ! new is 0
        sig_ct=-.8     ! new is .8
      endif
      if(nversion<809)then
        nvmix=5        ! new is 3
        ksc=0          ! new is -95
        sig_ct=.8      ! new is -.8
      endif
      if(nversion<806)then
        nvmix=3        ! new is 5
        ksc=-95        ! new is 0
        nclddia=5      ! new is 1
      endif
      if(nversion==803)then
        restol=2.e-7   ! new is 4.e-7
      endif
      if(nversion<803)then
        restol=5.e-7   ! new is 2.e-7
        alflnd=1.15    ! new is 1.1
        alfsea=1.05    ! new is 1.1
        entrain=0.     ! new is .05
        ksc=0          ! new is -95
      endif
      if(nversion==709)ksc=99
      if(nversion<709)then
        precon=0       ! new is -2900
        restol=2.e-7   ! new is 5.e-7
        mbase=2000     ! new is 101
        mdelay=0       ! new is -1
        nbase=-2       ! new is -4
        sigkscb=-.2    ! new is .95
        sigksct=.75    ! new is .8
        tied_con=6.    ! new is 2.
        tied_over=2.   ! new is 0.
        tied_rh=.99    ! new is .75
      endif
      if(nversion<705)then
        nstag=5        ! new is -10
        nstagu=5       ! new is -1.
        detrain=.3     ! new is .15
      endif
      if(nversion<704)then
        m=6            ! new is 5
        mex=4          ! new is 30.
        npex=1         ! new is 0
        ntsur=2        ! new is 6
      endif
      if(nversion<703)then
        ntbar=4        ! new is 6
        ntsur=7        ! new is 2
        vmodmin=2.     ! new is .2
        nbase=1        ! new is -2
      endif
      if(nversion<702)then
        npex=3         ! new is 1
      endif
      if(nversion<701)then
        nbase=0        ! new is 1  new variable
      endif
      if(nversion<608)then
        npex=1         ! new is 3 
        epsp=-20.      ! new is -15.
      endif
      if(nversion<606)then
        epsp=1.1       ! new is -20.
        newrough=0     ! new is 2
        nstag=-10      ! new is 5
        nstagu=3       ! new is 5
        ntsur=6        ! new is 7
        nvad=4         ! new is -4
        mbase=10       ! new is 2000
      endif
      if(nversion<604)then
        mh_bs=3        ! new is 4
      endif
      if(nversion<602)then
        ntbar=9        ! new is 4
      endif
       if(nversion<601)then
        epsp=1.2       ! new is 1.1
        newrough=2     ! new is 0
        restol=1.e-6   ! new is 2.e-7
      endif
      if(nversion<511)then
        nstag=3        ! new is -10
        nvsplit=3      ! new is 2
!       mins_rad=120   ! new is 72   not in common block, so can't assign here
        detrain=.1     ! new is .3
        mbase=1        ! new is 10
        nuvconv=5      ! new is 0
        sigcb=.97      ! new is 1.
      endif
      if(nversion<510)then
        epsp=.1        ! new is 1.2
        epsu=.1        ! new is 0.
        khdif=5        ! new is 2
        khor=0         ! new is -8
        nbarewet=7     ! new is 0
        newrough=0     ! new is 2
        nhor=0         ! new is -157
        nhorps=1       ! new is -1
        nlocal=5       ! new is 6
        ntsur=7        ! new is 6
        nxmap=1        ! new is 0
!       jalbfix=0      ! new is 1  not in common block, so can't assign here
        tss_sh=0.      ! new is 1.
        zobgin=.05     ! new is .02
        detrain=.4     ! new is .1
        convtime=.3    ! new is .33
        iterconv=2     ! new is 3
        mbase=0        ! new is 1
        sigcb=1.       ! new is .97
        sigkscb=.98    ! new is -2.
        tied_rh=.75    ! new is .99
        ldr=2          ! new is 1
        rcm=1.e-5      ! new is .92e-5
      endif
      if(nversion<509)then
        ntsur=6        ! new is 7
      endif
      if(nversion<508)then
        mh_bs=1        ! new is 3
!       mhint=2        ! new is 0   can only be re-set in parameter statement
        nvmix=4        ! new is 3
        entrain=.3     ! new is 0.
      endif
       if(nversion<507)then
        nxmap=0        ! new is 1
      endif
       if(nversion<506)then
        mh_bs=4        ! new is 1
!       mhint=0        ! new is 2   can only be re-set in parameter statement
      endif
      if(nversion<503)then
        ntsur=5        ! new is 6
      endif
      if(nversion<411)then
        nstag=-3       ! new is 3
        nstagu=-3      ! new is 3
        nvadh=1        ! new is 2
        nhor=155       ! new is 0
        nlocal=1       ! new is 5
        nvsplit=2      ! new is 3
        ngwd=0         ! new is -5
        nevapls=5      ! new is -4
        nuvconv=0      ! new is 5
        detrain=.05    ! new is .4
        entrain=0.     ! new is .3
        detrainx=1.    ! new is 0.
        dsig2=.1       ! new is .15
        dsig4=.55      ! new is .4
        kscmom=0       ! new is 1
        ldr=1          ! new is 2
        nbarewet=2     ! new is 7
        av_vmod=1.     ! new is .7
        chn10=.00137   ! new is .00125
      endif
      return
      end
