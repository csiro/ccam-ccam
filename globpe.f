#ifdef scyld
      subroutine globpe 
#else
      Program globpe  
#endif
!      PE model on conformal-cubic grid
!      N.B. on a Cray, set ncray=1 in depts.f, latltoij
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!      input files are :namelist (via file called "input")
!                       "nrun.dat"
!      data input and output file names are specified in namelist 'datafile'
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!     sign convention:
!                      u+ve eastwards  (on the panel)
!                      v+ve northwards (on the panel)
      use cc_mpi
      use diag_m
!     rml 21/02/06 removed redundant tracer code (so2/o2 etc)
!     rml 16/02/06 use tracermodule, timeseries
      use tracermodule, only :init_tracer,trfiles,tracer_mass,unit_trout
     &                        ,interp_tracerflux
      use timeseries, only : write_ts
      implicit none
      include 'newmpar.h'
      include 'aalat.h'
      include 'arrays.h'   ! ts, t, u, v, psl, ps, zs
      include 'bigxy4.h' ! common/bigxy4/xx4(iquad,iquad),yy4(iquad,iquad)
      include 'const_phys.h'
      include 'darcdf.h'   ! idnc,ncid  - stuff for reading netcdf
      include 'dates.h'    ! dtin,mtimer
      include 'establ.h'
      include 'extraout.h'
      include 'filnames.h' ! list of files, read in once only
      include 'histave.h'
      include 'indices.h'
      include 'kuocom.h'
      include 'latlong.h'  ! rlatt,rlongg
      include 'liqwpar.h'  ! ifullw,qfg,qlg
      include 'map.h'      ! em, f, dpsldt, fu, fv, etc
      include 'morepbl.h'
      include 'mpif.h'
      include 'netcdf.inc' ! stuff for writing netcdf
      include 'nlin.h'
      include 'nsibd.h'    ! sib, tsigmf, isoilm
      include 'parm.h'
      include 'parmdyn.h'  
      include 'parmhor.h'  ! mhint, mh_bs, nt_adv, ndept
      include 'parmsurf.h' ! nplens
      include 'parmvert.h'
      include 'pbl.h'      ! cduv, cdtq, tss, qg
      include 'prec.h'
      include 'raddiag.h'
      include 'savuvt.h'
      include 'scamdim.h'  ! npmax
      include 'screen.h'   ! tscrn etc
      include 'sigs.h'
      include 'soil.h'     ! fracice
      include 'soilsnow.h'
      include 'soilv.h'
      include 'stime.h'    ! kdate_s,ktime_s  sought values for data read
      include 'tracers.h'  ! ngas, nllp, ntrac, tr
      include 'trcom2.h'   ! nstn,slat,slon,istn,jstn etc.
      include 'vecsuv.h'
      include 'vecsuva.h'
      include 'version.h'
      include 'vvel.h'     ! sdot
      include 'xarrs.h'
      include 'xyzinfo.h'  ! x,y,z,wts

      integer leap,nbarewet,nsigmf
      integer nnrad,idcld
      real alph_p,alph_pm,delneg,delpos,alph_q
      real dpsdt,dpsdtb,dpsdtbb
      common/dpsdt/dpsdt(ifull),dpsdtb(ifull),dpsdtbb(ifull) !globpe,adjust5,outcdf
      real epst
      common/epst/epst(ifull)
      common/leap_yr/leap  ! 1 to allow leap years
      common/mfixdiag/alph_p,alph_pm,delneg,delpos,alph_q
      common/nsib/nbarewet,nsigmf
      common/radnml/nnrad,idcld

      real savs1(ifull,2:kl),savu1(ifull,kl),savv1(ifull,kl)
      real savu2(ifull,kl),savv2(ifull,kl)
      real sbar
      common/sbar/sbar(ifull,2:kl)
      common/savuv1/savs1,savu1,savv1 
      real ubar(ifull,kl),vbar(ifull,kl)
      common/uvbar/ubar,vbar

      real taftfh(ifull),taftfhg(ifull)
      common/tafts/taftfh,taftfhg
      real unn,vnn
      common/unn/unn(ifull,kl),vnn(ifull,kl) ! nonlin,upglobal for nvsplit3,4
      real, dimension(ifull) :: dirad,dfgdt,degdt,wetfac,degdw,cie,
     &                          factch,qsttg,rho,zo,aft,fh,spare1,theta,
     &                          gamm,rg,vmod,dgdtg
      common/work2/dirad,dfgdt,degdt,wetfac,degdw,cie,factch,qsttg,rho,
     &     zo,aft,fh,spare1,theta,gamm,rg,vmod,dgdtg

      real, dimension(ifull) :: egg,evapxf,ewww,fgf,fgg,ggflux,rdg,rgg,
     &                          residf,ga,condxpr,fev,fes,fwtop,spare2
      integer, dimension(ifull) :: ism
      real, dimension(ifull,kl) :: dum3a,speed,spare
      real, dimension(2*ijk-16*ifull) :: dum3
      common/work3/egg,evapxf,ewww,fgf,fgg,ggflux,rdg,rgg,residf,ga,
     &     condxpr,fev,fes,ism,fwtop,spare2,dum3,dum3a,speed,spare

      real wblf,wbfice,sdepth,dum3b
      common/work3b/wblf(ifull,ms),wbfice(ifull,ms),sdepth(ifull,3),
     &              dum3b(ijk*2-2*ifull*ms-3*ifull)
      real rtt
      common/work3d/rtt(ifull,kl) ! just to pass between radriv90 & globpe
      real qccon, qlrad, qfrad
      common/work3f/qccon(ifull,kl),qlrad(ifull,kl),qfrad(ifull,kl) !leoncld etc
      real cfrac
      common/cfrac/cfrac(ifull,kl)     ! globpe,radriv90,vertmix,convjlm
      real qgsav, qfgsav, qlgsav, trsav
      common/work3sav/qgsav(ifull,kl),qfgsav(ifull,kl),qlgsav(ifull,kl)
     &             ,trsav(ilt*jlt,klt,ngasmax)  ! shared adjust5 & nonlin
      real spmean(kl),div(kl),omgf(ifull,kl),pmsl(ifull)
      equivalence (omgf,dpsldt),(pmsl,dum3a)
      integer, dimension(13), parameter :: mdays =
     &     (/31,28,31,30,31,30,31,31,30,31,30,31, 31/)
      logical odcalc
      character comm*60,comment*60,rundate*8,header*47,text*2
      character(len=10) :: timeval
      integer, dimension(8) :: tvals1, tvals2

!     Local variables
      integer iaero, ier, igas, ilx, io_nest, iq, irest, isoil, itr1,
     &     itr2, jalbfix, jlx, jtr1, jtr2, k,k2, kktau, 
     &     mexrest, mins_dt, mins_gmt, mspeca, mtimer_in,
     &     nalpha, newsnow, ng, nlx, nmaxprsav,
     &     nmi, npa, npb, npc, n3hr,      ! can remove npa,npb,npc
     &     nsnowout, nstagin, nstaguin, nwrite, nwtsav, mins_mbd,
     &     mins_rad, mtimer_sav, nsecs_diff, nn, i, j
     &     ,mstn  ! not used in 2006
      integer, dimension(8) :: nper3hr
      real clhav, cllav, clmav, cltav, con, 
     &     div_int, dsx, dtds, es, gke, hourst, hrs_dt,
     &     evapavge, precavge, preccavge, psavge,
     &     pslavge, pwater, rel_lat, rel_long, rlwup, spavge,
     &     pwatr, pwatr_l, qtot, aa,bb,cc,bb2,cc2,rat
      real, dimension(9) :: temparray, gtemparray ! For global sums
      integer :: nproc_in, ierr, nperhr, nscrn=0

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
     & ,rlong0,rlat0,schmidt   ! usually come in with topofile
     & ,kbotdav,kbotu,nbox,nud_p,nud_q,nud_t,nud_uv,nud_hrs,nudu_hrs
     & ,nlocal,nvsplit,nbarewet,nsigmf,qgmin
     & ,io_clim ,io_in,io_nest,io_out,io_rest,io_spec,nfly,localhist   
     & ,mstn,nqg   ! not used in 2006          
     & ,nurban ! MJT urban
      data npc/40/,nmi/0/,io_nest/1/,iaero/0/,newsnow/0/ 
      namelist/skyin/mins_rad,ndiur  ! kountr removed from here
      namelist/datafile/ifile,ofile,albfile,co2emfile,eigenv,
     &    hfile,icefile,mesonest,nmifile,o3file,radfile,restfile,
     &    rsmfile,scamfile,scrnfile,snowfile,so4tfile,soilfile,sstfile,
     &    surfile,tmaxfile,tminfile,topofile,trcfil,vegfile,zofile,
     &    smoistfile,soil2file,radonemfile,
     &    co2_00,radon_00,surf_00,co2_12,radon_12,surf_12
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
     &        ,acon,bcon,rcm,rcrit_l,rcrit_s   ! ldr stuff
      data itr1/23/,jtr1/13/,itr2/25/,jtr2/11/
      data comment/' '/,comm/' '/,irest/1/,jalbfix/1/,nalpha/1/
      data mexrest/6/,mins_rad/60/
      data nwrite/0/
      data nsnowout/999999/

      ! Check that declarations in include files match
      call check_dims()

#ifndef scyld
      call MPI_Init(ierr)       ! Start
#endif
      call MPI_Comm_size(MPI_COMM_WORLD, nproc_in, ierr) ! Find # of processes
      if ( nproc_in /= nproc ) then
         print*, "Error, model is compiled for ", nproc, " processors."
         print*, "Trying to run with", nproc_in
         call MPI_Abort(MPI_COMM_WORLD,ierr)
      end if
      call MPI_Comm_rank(MPI_COMM_WORLD, myid, ierr) ! Find my id

      call log_off()
      call log_setup()
#ifdef simple_timer
      call start_log(model_begin)
#endif

      ia=il/2
      ib=ia+3
!     ntbar=(kl+1)/4  ! just a default
      ntbar=kl/3      ! just a default 7/3/07

      ! Some extra default values
      rel_lat = 0.
      rel_long = 0.
      ! All processors read the namelist, standard input doesn't work properly
      open(99,file="input",form="formatted",status="old")
      if (myid==0)
     &  write(6,'(a10," compiled for il,jl,kl,nproc =",4i7)')
     &                       version,il_g,jl_g ,kl,nproc 
      read (99, cardin)
      nperday =nint(24.*3600./dt)
      nperhr  =nint(3600./dt)
      if(nmaxpr==99)nmaxpr=nperday/4
      do n3hr=1,8
       nper3hr(n3hr)=nint(n3hr*3*3600/dt)
      enddo
c     if(nstag==99)nstag=-nper3hr(2)   ! i.e. 6-hourly value
      if(nwt==-99)nwt=nperday          ! set default nwt to 24 hours
      if(nperavg==-99)nperavg=nwt      ! set default nperavg to nwt
      if(nwrite==0)nwrite=nperday  ! only used for outfile IEEE
      read (99, skyin)
      kountr=nint(mins_rad*60./dt)  ! set default radiation to ~mins_rad m
      mins_rad=nint(kountr*dt/60.)  ! redefine to actual value
      read (99, datafile)
      read (99, kuonml)
!     rml 16/02/06 read trfiles namelist and call init_tracer
      if(ngas>0) then
        read(99, trfiles)
        call init_tracer
      endif
     
      if(nudu_hrs==0)nudu_hrs=nud_hrs
      if(kbotu==0)kbotu=kbotdav
      dsig4=max(dsig2+.01,dsig4)
      if(mbd.ne.0.and.nbd.ne.0)then
        print *,'setting nbd=0 because mbd.ne.0'
        nbd=0
      endif
      if(nbd.ne.0)nud_hrs=abs(nud_hrs)  ! just for people with old -ves in namelist

      if ( myid == 0 ) then   ! **** do namelist fixes above this ***
      print *,'Dynamics options A:'
      print *,'   m    mex   mfix  mfix_qg   mup    nh    nonl',    
     &          '   npex  precon' 
      write(6,'(i5,9i7)')m,mex,mfix,mfix_qg,mup,nh,nonl,npex,precon
      print *,'Dynamics options B:'
      print *,'nritch_t nrot  ntbar  nxmap   epsp    epsu   epsf',
     &        '   restol'
      write (6,'(i5,3i7,1x,3f8.3,g9.2)')nritch_t,nrot,ntbar,nxmap,
     &          epsp,epsu,epsf,restol
      print *,'Horizontal advection/interpolation options:'
      print *,' ndept  nt_adv mh_bs  mhint '
      write (6,'(i5,11i7)') ndept,nt_adv,mh_bs,mhint
      print *,'Horizontal wind staggering options:'
      print *,'mstagpt nstag nstagu'
      write (6,'(i5,11i7)') mstagpt,nstag,nstagu
      if(nstag==0)stop 'need non-zero value for nstag'
      print *,'Vertical advection options:'
      print *,'  nvad  nvadh  '
      write (6,'(i5,11i7)') nvad,nvadh
      if(nvad==4.or.nvad==-4)then
        print *,'Vertical advection options for TVD:'
        print *,' nimp   nthub  ntvd   ntvdr'
        write (6,'(i5,11i7)') nimp,nthub,ntvd,ntvdr
      endif
      print *,'Horizontal mixing options:'
      print *,' khdif  khor   nhor   nhorps nhorjlm'
      write (6,'(i5,11i7)') khdif,khor,nhor,nhorps,nhorjlm
      print *,'Vertical mixing/physics options:'
      print *,' nvmix nlocal nvsplit ncvmix lgwd   ngwd   '
      write (6,'(i5,6i7)') nvmix,nlocal,nvsplit,ncvmix,lgwd,ngwd
      print *,'Cumulus convection options A:'
      print *,' nkuo  sigcb sig_ct  rhcv  rhmois rhsat',
     &        ' convfact convtime shaltime'
      write (6,'(i5,6f7.2,9f8.2)')
     &    nkuo,sigcb,sig_ct,rhcv,rhmois,rhsat,
     &    convfact,convtime,shaltime
      print *,'Cumulus convection options B:'
      print *,' alflnd alfsea fldown iterconv',
     &        ' ncvcloud nevapcc nevapls nuvconv'
      write (6,'(3f7.2,i6,i10,4i8)') alflnd,alfsea,fldown,iterconv,
     &    ncvcloud,nevapcc,nevapls,nuvconv
      print *,'Cumulus convection options C:'
      print *,' mbase mdelay methprec nbase detrain',
     &        ' entrain methdetr detrainx dsig2  dsig4'
      write (6,'(3i6,i9,f8.2,f9.2,i8,4f8.2)') mbase,mdelay,methprec,
     &              nbase,detrain,entrain,methdetr,detrainx,dsig2,dsig4
      print *,'Shallow convection options:'
      print *,'  ksc  kscsea kscmom sigkscb sigksct ',
     &        'tied_con tied_over tied_rh '
      write (6,'(i5,2i7,1x,3f8.3,2f10.3)') ksc,kscsea,kscmom,
     &    sigkscb,sigksct,tied_con,tied_over,tied_rh
      print *,'Other moist physics options:'
      print *,'  acon   bcon   qgmin      rcm    rcrit_l rcrit_s'
      write (6,'(2f7.2,2e10.2,2f7.2)') acon,bcon,qgmin,rcm,
     &                                 rcrit_l,rcrit_s
      print *,'Radiation options:'
      print *,' nrad  ndiur mins_rad kountr  dt'
      write (6,'(i5,3i7,f10.2)') nrad,ndiur,mins_rad,kountr,dt
!     for 6-hourly output of sint_ave etc, want 6*60 = N*mins_rad      
      if(nrad==4.and.mod(6*60,mins_rad).ne.0)
     &                              stop 'prefer 6*60 = N*mins_rad '
      print *,'Diagnostic cloud options:'
      print *,'  ldr nclddia nstab_cld nrhcrit sigcll '
      write (6,'(i5,i6,2i9,1x,f8.2)') ldr,nclddia,nstab_cld,nrhcrit,
     &                                sigcll
      print *,'Soil, canopy and PBL options A:'
      print *,' jalbfix nalpha nbarewet newrough nglacier nrungcm',
     &        ' nsib  nsigmf'
      write (6,'(i5,9i8)')
     &          jalbfix,nalpha,nbarewet,newrough,nglacier,nrungcm,nsib,
     &          nsigmf
      print *,'Soil, canopy and PBL options B:'
      print *,' ntaft ntsea ntsur av_vmod tss_sh vmodmin  zobgin',
     &        ' charnock chn10'
      write (6,'(i5,2i6,4f8.2,f8.3,f9.5)') ntaft,ntsea,ntsur,
     &          av_vmod,tss_sh,vmodmin,zobgin,charnock,chn10    
      if(mbd.ne.0.or.nbd.ne.0)then
        print *,'Nudging options:'
        print *,' nbd    nud_p  nud_q  nud_t  nud_uv nud_hrs nudu_hrs',
     &          ' kbotdav  kbotu'
        write (6,'(i5,3i7,7i8)') 
     &    nbd,nud_p,nud_q,nud_t,nud_uv,nud_hrs,nudu_hrs,kbotdav,kbotu
      endif
      print *,'Special and test options:'
      print *,' mbd namip amipo3 newtop nhstest nplens nsemble',
     &        ' nspecial panfg panzo'
      write (6,'(2i5,L7,4i7,i8,f9.1,f8.4)') mbd,namip,amipo3,
     &          newtop,nhstest,nplens,nsemble,nspecial,panfg,panzo
      print *,'I/O options:'
      print *,' nfly  io_in io_nest io_out io_rest  nwt  nperavg nscrn'
      write (6,'(i5,4i7,3i8)') 
     &          nfly,io_in,io_nest,io_out,io_rest,nwt,nperavg,nscrn
      if(ntrac.ne.0)then
        print *,'Trace gas options:'
        print *,' iradon  ico2  ngas   nllp   ntrac'
        write (6,'(i5,5i7)') iradon,ico2,ngas,nllp,ntrac
      endif

      write (6, cardin)
      if(nllp==0.and.nextout>=4)stop 'need nllp=3 for nextout>=4'
      write (6, skyin)
      write (6, datafile)
      write(6, kuonml)
      end if ! myid=0
      if(newtop>2)stop 'newtop>2 no longer allowed'
      if(mfix_qg>0.and.(nkuo==4.or.nvad==44))
     &        stop 'nkuo=4,nvad=44: mfix_qg>0 not allowed'
      if(mfix>3)stop 'mfix >3 not allowed now'

      ktau=0
      if(io_in<=4.and.myid==0)then
!       open new topo file and check its dimensions
!       here used to supply rlong0,rlat0,schmidt
        open(66,file=topofile,recl=2000,status='old')
        print *,'reading topofile header'
!       read(66,'(i3,i4,2f6.1,f5.2,f9.0,a47)')
c       read(66,'(i3,i4,2f6.1,f6.3,f8.0,a47)')
        read(66,*)
     &            ilx,jlx,rlong0,rlat0,schmidt,dsx,header
        print *,'ilx,jlx,rlong0,rlat0,schmidt,dsx ',
     &           ilx,jlx,rlong0,rlat0,schmidt,dsx,header
        if(ilx.ne.il_g.or.jlx.ne.jl_g)stop 'wrong topo file supplied'
      endif      ! (io_in<=4)
      call MPI_Bcast(rlong0,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(rlat0,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(schmidt,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
!     N.B. to display orog alone, run with io_in=4, nbd=0, nsib=0      

c     set up cc geometry
!     All processors call setxyz
      call setxyz(il_g,xx4,yy4,myid)
      call ccmpi_setup()
      call setaxu    ! globpea code
      call bounds(axu)
      call bounds(bxv)
      call bounds(ayu)
      call bounds(byv)
      call bounds(azu)
      call bounds(bzv)
      con=180./pi
      if ( mydiag ) then
         print *,'id,jd,rlongg,rlatt in degrees: ',
     &         id,jd,con*rlongg(idjd),con*rlatt(idjd)
      end if

      if ( myid == 0 ) then
         call date_and_time(rundate)
         print *,'RUNDATE IS ',rundate
         call date_and_time(time=timeval)
         print*, "Starting time ", timeval
      end if

!     open input files
!     Note that all processes read the eigenv file
      open(28,file=eigenv,status='old',form='formatted')
      if (myid==0) then
         if(io_in==2)
     &       open(10,file=ifile,form='formatted',status='old')
         if(abs(io_in)==3)
     &       open(10,file=ifile,form='unformatted',status='old')
         if(abs(io_in)==1)then ! netcdf input
            ncid = ncopn ( ifile,0,ier ) ! 0 denotes read-only
            print *,'ncid,ier,ifile ',ncid,ier,ifile
         endif
      end if

!     open input file for radiative data.
!     replaces block data bdata in new Fels-Schwarzkopf radiation code.
!     All processes read these
      if(nrad==4)then
        if (myid==0) print *,'Radiative data read from file ',radfile
        open(16,file=o3file,form='formatted',status='old')
        open(15,file=radfile,form='formatted',status='old')
      endif
      if(iaero.ne.0)then
         if (myid==0) print *,'so4total data read from file ',so4tfile
         call readreal(so4tfile,so4t,ifull)
      endif

!     set some 2D/3D default values
      wb(:,:)=.15
      snowd(:)=0.
      condx(:)=0.
      npc=max(npc,1)
c     if(ndi2>0)diag=.true.

!     indata for initial conditions
      if (myid==0) print *,'calling indata; will read from file ',ifile
!     N.B. first call to amipsst now done within indata
      call indata(hourst,newsnow,jalbfix)
      
      call maxmin(u,' u',ktau,1.,kl)
      call maxmin(v,' v',ktau,1.,kl)
      ! Note that u and v are extended.
      speed(:,:)=sqrt(u(1:ifull,:)**2+v(1:ifull,:)**2)  ! 3D 
      call maxmin(speed,'sp',ktau,1.,kl)
      call maxmin(t,' t',ktau,1.,kl)
      call maxmin(qg,'qg',ktau,1.e3,kl)
      call maxmin(qfg,'qf',ktau,1.e3,kl)
      call maxmin(qlg,'ql',ktau,1.e3,kl)
      call maxmin(wb,'wb',ktau,1.,ms)
      ! These were overlapped in a single call which doesn't work in the MPI
      ! version.
!      call maxmin(tggsn,'tgg',ktau,1.,ms+3)
      call maxmin(tggsn,'tggsn',ktau,1.,3)
      call maxmin(tgg,'tgg',ktau,1.,ms)
      pwatr_l=0.   ! in mm
      do k=1,kl
        do iq=1,ifull
         qtot=qg(iq,k)+qlg(iq,k)+qfg(iq,k)
         pwatr_l=pwatr_l-dsig(k)*wts(iq)*qtot*ps(iq)
        enddo
      enddo
      pwatr_l = pwatr_l/grav
      call MPI_Reduce ( pwatr_l, pwatr, 1, MPI_REAL, MPI_SUM, 0,
     &                  MPI_COMM_WORLD, ierr )
      if ( myid == 0 ) write (6,"('pwatr0 ',12f7.3)") pwatr
      if(nextout>=4)call setllp
      if(ntrac>0)then
        do ng=1,ntrac
         write (text,'("g",i1)')ng
         call maxmin(tr(:,:,ng),text,ktau,1.,kl)
        enddo
      endif   ! (ntrac>0)

      if (myid == 0 ) then
        close(10)
      end if
      if(mbd.ne.0.or.nbd.ne.0)then
         io_in=io_nest ! Needs to be seen by all processors
         if ( myid == 0 ) then
           if(abs(io_in)==1)then
             ncid = ncopn(mesonest,0,ier )  ! 0 denotes read-only
             print *,'ncid,ier,mesonest ',ncid,ier,mesonest
             if(ier.ne.0)then
               print *,'cannot open netcdf mesofile ',nf_strerror(ier)
               stop
             endif
           endif  ! (abs(io_in)==1)
           if(io_in==2)
     &       open(10,file=mesonest,form='FORMATTED',status='OLD')
           if(abs(io_in)==3)
     &       open(10,file=mesonest,form='UNFORMATTED',status='OLD')
         endif ! myid == 0
      endif    ! (mbd.ne.0.or.nbd.ne.0)
!     open output files; name is stored in namelist file
      if ( myid == 0 ) then
        if(io_out==2)
     &    open(20,file=ofile,form='formatted',status='unknown')
        if(io_out==3)
     &    open(20,file=ofile,form='unformatted',status='unknown')
c       if(ilt>1)open(37,file='tracers_latest',status='unknown')
      end if ! myid == 0

!     sig(kuocb) occurs for level just BELOW sigcb
      kuocb=1
      do while(sig(kuocb+1)>=sigcb)
       kuocb=kuocb+1
      enddo
      if (myid==0)
     &   print *,'convective cumulus scheme: kuocb,sigcb = ',kuocb,sigcb

      if(khdif==-99)then   ! set default khdif appropriate to resolution
        khdif=5
        if (myid==0) print *,'Model has chosen khdif =',khdif
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
        if (myid==0) print *,'khor,hdiff: ',khor,hdiff
      endif

!     *** be careful not to choose ia,ib,ja,jb to cover more than
!         one processor, or myid=0 may stop here!!!
      call printa('zs  ',zs,0,0,ia,ib,ja,jb,0.,.01)
      call printa('tss ',tss,0,0,ia,ib,ja,jb,200.,1.)
      if (mydiag) print *,'wb(idjd) ',(wb(idjd,k),k=1,6)
      call printa('wb1   ',wb ,0,1,ia,ib,ja,jb,0.,100.)
      call printa('wb6  ',wb,0,ms,ia,ib,ja,jb,0.,100.)

      if (myid==0) then
         open(11, file='nrun.dat',status='unknown')
         if(nrun==0)then
            read(11,*,end=227) nrun
 227        nrun=nrun+1
         endif                  ! nrun==0
         print *,'this is run ',nrun
         rewind 11
         write(11,*) nrun
         write(11,cardin)
         write(11,datafile)
         close(11)
      end if
      
!     Zero/set the diagnostic arrays
      rndmax(:)=0.
      tmaxscr(:)=0.
      tminscr(:)=400.
      rhmaxscr(:)=0.
      rhminscr(:)=400.
      u10max(:)=0.
      v10max(:)=0.
      u10mx(:)=0.
      tscr_ave(:)=0.
      qscrn_ave(:)=0.
      dew_ave(:)=0.
      epan_ave(:)=0.
      eg_ave(:)=0.
      fg_ave(:)=0.
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
      koundiag=0
      sint_ave(:) = 0.  ! solar_in_top
      sot_ave(:)  = 0.  ! solar_out_top
      soc_ave(:)  = 0.  ! solar_out_top (clear sky)
      sgdn_ave(:) = 0.  ! solar_ground (down-welling) +ve down
      sgn_ave(:)  = 0.  ! solar_ground (net) +ve down
      rtu_ave(:)  = 0.  ! LW_out_top 
      rtc_ave(:)  = 0.  ! LW_out_top (clear sky)
      rgdn_ave(:) = 0.  ! LW_ground (down-welling)  +ve up
      rgn_ave(:)  = 0.  ! LW_ground (net)  +ve up
      rgc_ave(:)  = 0.  ! LW_ground (clear sky)
      cld_ave(:)  = 0.
      cll_ave(:)  = 0.
      clm_ave(:)  = 0.
      clh_ave(:)  = 0.
 
      if(nmi==0.and.nwt>0)then
!       write out the first ofile data set 
        print *,'calling outfile myid= ',myid
        call outfile(20,rundate,nmi,nwrite)
        if(newtop<0)stop  ! just for outcdf to plot zs  & write fort.22
      endif    ! (nmi==0.and.nwt.ne.0)
      dtin=dt
      n3hr=1   ! initial value at start of run
      if (myid==0) then
         print *,'number of time steps per day = ',nperday
         print *,'nper3hr,nper6hr .. ',nper3hr(:)
      end if
      mspeca=1
      if(mex.ne.1)then
        mspeca=2
        dt=dtin*.5
      endif
      call gettin(0)             ! preserve initial mass & T fields; nmi too

      if(nbd.ne.0)call nestin
      if(mbd.ne.0)call nestinb(mins_mbd)  ! to get mins_mbd
      nmaxprsav=nmaxpr
      nwtsav=nwt
      hrs_dt = dtin/3600.      ! time step in hours
      mins_dt = nint(dtin/60.)  ! time step in minutes
      mtimer_in=mtimer
      nstagin=nstag
      nstaguin=nstagu
      if(nstagin==5.or.nstagin<0)then
        nstag=4
        nstagu=4
        if(nstagin==5)then  ! for backward compatability
          nstagin=-1
          nstaguin=5
        endif
      endif
 
      if ( myid == 0 ) then
         call date_and_time(time=timeval,values=tvals1)
         print*, "Start of loop time ", timeval
      end if
      call log_on()
#ifdef simple_timer
      call start_log(maincalc_begin)
#endif

      do 88 kktau=1,ntau   ! ****** start of main time loop
      ktau=kktau
      timer = timer + hrs_dt      ! timer now only used to give timeg
      timeg=mod(timer+hourst,24.)
!     mtimer=mtimer+mins_dt
      mtimer=mtimer_in+nint(ktau*dtin/60.)     ! 15/6/01 to allow dt < 1 minute
      mins_gmt=mod(mtimer+60*ktime/100,24*60)
      if(nbd.ne.0)call nestin
      if(nstaguin>0.and.ktau>1)then
        if(nstagin<0.and.mod(ktau,abs(nstagin))==0)then
          nstag=7-nstag  ! swap between 3 & 4
          nstagu=nstag
        endif
      endif

!     rml 17/02/06 interpolate tracer fluxes to current timestep
      if(ngas>0) call interp_tracerflux(kdate,hrs_dt)

      do 79 mspec=mspeca,1,-1    ! start of introductory time loop
      dtds=dt/ds
      if(nvsplit<3.or.ktau==1)then
        un(1:ifull,:)=0. 
        vn(1:ifull,:)=0.
        tn(1:ifull,:)=0.
      elseif(nvsplit==3)then
        tn(1:ifull,:)=(t(1:ifull,:)-tx(1:ifull,:))/dt ! tend. from phys. at end of previous step
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
      if(mup.ne.1.or.(ktau==1.and.mspec==mspeca))then
!       updps called first step or to permit clean restart option      
        call updps(0) 
      endif

!     set up tau +.5 velocities in ubar, vbar
      if(ktau<10.and.mydiag)then
       print*,'ktau,mex,mspec,mspeca:',ktau,mex,mspec,mspeca
      endif
      sbar(:,2:kl)=sdot(:,2:kl)
      if(ktau==1)then
!       this sets (ubar,vbar) to ktau=1.5 values on 2nd time through
        ubar(:,:)=u(1:ifull,:)
        vbar(:,:)=v(1:ifull,:)
      elseif(mex==1)then
        ubar(:,:)=u(1:ifull,:)
        vbar(:,:)=v(1:ifull,:)
      elseif(ktau==2.or.mex==2)then        
!       (tau+.5) from tau, tau-1
        ubar(:,:)=u(1:ifull,:)*1.5-savu(:,:)*.5
        vbar(:,:)=v(1:ifull,:)*1.5-savv(:,:)*.5
      elseif(mex==3)then
!       (tau+.5) from tau, tau-1, tau-2   ! ubar is savu1 here
c       sbar(:,:)=sdot(:,2:kl)+.5*(savs(:,:)-savs1(:,:))
        ubar(:,:)=u(1:ifull,:)+.5*(savu(:,:)-savu1(:,:))
        vbar(:,:)=v(1:ifull,:)+.5*(savv(:,:)-savv1(:,:))
       elseif(mex==30.and.ktau>3)then  ! using tau, tau-1, tau-2, tau-3
        do k=1,kl
         do iq=1,ifull
          bb=1.5*u(iq,k)-2.*savu(iq,k)+.5*savu1(iq,k)    ! simple b
          bb2=(40.*u(iq,k)-35.*savu(iq,k)                ! cwqls b
     &          -16.*savu1(iq,k)+11.*savu2(iq,k))/34.
          cc=.5*u(iq,k)-savu(iq,k)+.5*savu1(iq,k)        ! simple c
          cc2=(10.*u(iq,k)-13.*savu(iq,k)                ! cwqls c
     &          -4.*savu1(iq,k)+7.*savu2(iq,k))/34.
          aa=cc2-cc
          rat=max(0.,min(1.,cc2/(aa+sign(1.e-9,aa))))
          cc=rat*cc+(1.-rat)*cc2 
          bb=rat*bb+(1.-rat)*bb2 
          ubar(iq,k)=u(iq,k)+.5*bb+.25*cc
          bb=1.5*v(iq,k)-2.*savv(iq,k)+.5*savv1(iq,k)    ! simple b
          bb2=(40.*v(iq,k)-35.*savv(iq,k)                ! cwqls b
     &          -16.*savv1(iq,k)+11.*savv2(iq,k))/34.
          cc=.5*v(iq,k)-savv(iq,k)+.5*savv1(iq,k)        ! simple c
          cc2=(10.*v(iq,k)-13.*savv(iq,k)                ! cwqls c
     &          -4.*savv1(iq,k)+7.*savv2(iq,k))/34.
          aa=cc2-cc
          rat=max(0.,min(1.,cc2/(aa+sign(1.e-9,aa))))
          cc=rat*cc+(1.-rat)*cc2 
          bb=rat*bb+(1.-rat)*bb2 
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
c       write (6,"(i4,' savs2,savs1,savs,sdot,sbar',5f8.4)") ktau,
c    &      savs2(idjd,nlx),savs1(idjd,nlx),savs(idjd,nlx),
c    &      sdot(idjd,nlx),sbar(idjd,nlx)
        write (6,"(i4,' savu2,savu1,savu,u,ubar',5f8.2)") ktau,
     &      savu2(idjd,nlv),savu1(idjd,nlv),savu(idjd,nlv),
     &      u(idjd,nlv),ubar(idjd,nlv)
        write (6,"(i4,' savv2,savv1,savv,v,vbar',5f8.2)") ktau,
     &      savv2(idjd,nlv),savv1(idjd,nlv),savv(idjd,nlv),
     &      v(idjd,nlv),vbar(idjd,nlv)
      endif
      if(ktau>2.and.epsp>1..and.epsp<2.)then ! 
        if (mydiag.and.ktau==3)print *,'using epsp= ',epsp
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
       print *,'savu,u,ubar ',ktau,savu(idjd,1),u(idjd,1),ubar(idjd,1)
      endif
      if(ktau==1.and.mspec==1.and.mex.ne.1)then
        u(1:ifull,:)=savu(1:ifull,:)  ! reset u,v to original values
        v(1:ifull,:)=savv(1:ifull,:)
      endif
c     savs2(1:ifull,:)=savs1(1:ifull,:)  
      savu2(1:ifull,:)=savu1(1:ifull,:)  
      savv2(1:ifull,:)=savv1(1:ifull,:)
      savs1(1:ifull,:)=savs(1:ifull,:)  
      savu1(1:ifull,:)=savu(1:ifull,:)  
      savv1(1:ifull,:)=savv(1:ifull,:)
      savs(1:ifull,:) =sdot(1:ifull,2:kl)  
      savu(1:ifull,:) =u(1:ifull,:)  ! before any time-splitting occurs
      savv(1:ifull,:) =v(1:ifull,:)
c     if(mex.ne.4)sdot(:,2:kl)=sbar(:,:)   ! ready for vertical advection

      diag=.false.
      if(ktau>=abs(ndi).and.ktau<=ndi2)diag=.true.
      if(ndi<0)then
        if(ktau==ktau/ndi*ndi)diag=.true.
      endif
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

      call nonlin
      if (diag)then
         if (mydiag) print *,'before hadv'
         call printa('tx  ',tx,ktau,nlv,ia,ib,ja,jb,0.,1.)
         if (mydiag) then
            nlx=min(nlv,kl-8)
            write (6,"('tx  ',9f8.2)") (tx(idjd,k),k=nlx,nlx+8)
            write (6,"('txe ',9f8.2)") (tx(ie(idjd),k),k=nlx,nlx+8)
            write (6,"('txw ',9f8.2)") (tx(iw(idjd),k),k=nlx,nlx+8)
            write (6,"('txn ',9f8.2)") (tx(in(idjd),k),k=nlx,nlx+8)
            write (6,"('txs ',9f8.2)") (tx(is(idjd),k),k=nlx,nlx+8)
            write(6,'(i2," qgv ",18f7.4)')ktau,(1000.*qg(idjd,k),k=1,kl)
         end if
         call printa('qgv ',qg,ktau,nlv,ia,ib,ja,jb,0.,1.e3)
      endif

!     evaluate horizontal advection for combined quantities
      call upglobal
      if (diag)then
         if (mydiag) then
            print *,'after hadv'
            write (6,"('tx  ',9f8.2)") (tx(idjd,k),k=nlx,nlx+8)
         end if
         call printa('tx  ',tx,ktau,nlv,ia,ib,ja,jb,200.,1.)
         if (mydiag)write(6,'(i2," qgh ",18f7.4)')ktau,1000.*qg(idjd,:)
      endif

      if(nonl<0)then
        savt(1:ifull,:)=t(1:ifull,:)  ! can be used in nonlin during next step
      endif

!!      sbar(1:ifull,:) = savs1(1:ifull,:)  ! 3D really saving savs1 in sbar here 
!!      ubar(1:ifull,:) = savu1(1:ifull,:)  ! 3D really saving savu1 in ubar here 
!!      vbar(1:ifull,:) = savv1(1:ifull,:)  ! 3D really saving savv1 in vbar here 

      if(nstaguin<0.and.ktau>1)then
        if(nstagin<0.and.mod(ktau,abs(nstagin))==0)then
          nstag=7-nstag  ! swap between 3 & 4
          nstagu=nstag
        endif
      endif
      call adjust5
      if(mspec==1.and.nbd.ne.0)call davies  ! nesting now after mass fixers
      if(mbd.ne.0)then
       nsecs_diff=mod(nint(ktau*dt),60*mins_mbd)
       if(nsecs_diff<nint(dt)/2.or.  ! calc in secs, as dt may be < 1 m
     &    abs(nsecs_diff-60*mins_mbd)<=nint(dt)/2)then
          if (mydiag) print *,'bingo call nestinb: ktau '
     &                       ,ktau,mtimer,nsecs_diff 
          call nestinb(mins_mbd)
        endif
      endif

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

!      pwatr_l=0.   ! in mm
!      do k=1,kl
!        do iq=1,ifull
!	  qtot=qg(iq,k)+qlg(iq,k)+qfg(iq,k)
!          pwatr_l=pwatr_l-dsig(k)*wts(iq)*qtot*ps(iq)
!        enddo
!      enddo
!      pwatr_l = pwatr_l/grav
!      call MPI_Reduce ( pwatr_l, pwatr, 1, MPI_REAL, MPI_SUM, 0,
!     &                  MPI_COMM_WORLD, ierr )
!      if ( myid == 0 ) print*,'ktau,ave_pwatr ',ktau,pwatr

      if(nhor<0)call hordifgt  ! now not tendencies
      if (diag.and.mydiag)print *,'after hordifgt t ',t(idjd,:)
      call start_log(phys_begin)
      if(ngwd<0)call gwdrag  ! <0 for split - only one now allowed
      if(nkuo==23)call convjlm     ! split convjlm 
      if( nkuo /= 0 ) then
         ! Not set in HS tests.
         cbas_ave(:)=cbas_ave(:)+condc(:)*(1.1-sig(kbsav(:))) ! diagnostic
         ctop_ave(:)=ctop_ave(:)+condc(:)*(1.1-sig(abs(ktsav(:)))) ! diagnostic
      end if
      if(nkuo==46)call conjob    ! split Arakawa-Gordon scheme
      if(nkuo==5)call betts(t,qg,tn,land,ps) ! not called these days

      if(ldr.ne.0)then
c       print*,'Calling prognostic cloud scheme'
        call leoncld(cfrac)  !Output
        do k=1,kl
         riwp_ave(:)=riwp_ave(:)-qfrad(:,k)*dsig(k)*ps(1:ifull)/grav ! ice water path
         rlwp_ave(:)=rlwp_ave(:)-qlrad(:,k)*dsig(k)*ps(1:ifull)/grav ! liq water path
        enddo
        if(nmaxpr==1.and.mydiag)then
          write (6,"('qfrad',3p9f8.3/5x,9f8.3)") qfrad(idjd,:)
          write (6,"('qlrad',3p9f8.3/5x,9f8.3)") qlrad(idjd,:)
          write (6,"('qf   ',3p9f8.3/5x,9f8.3)") qfg(idjd,:)
        endif
      endif  ! (ldr.ne.0)
      rnd_3hr(:,8)=rnd_3hr(:,8)+condx(:)  ! i.e. rnd24(:)=rnd24(:)+condx(:)

!       put radiation here
        if(nrad==4) then
!         Fels-Schwarzkopf radiation
          odcalc=mod(ktau,kountr)==0 .or. ktau==1 ! ktau-1 better
          nnrad=kountr
          if(nhstest<0)then ! aquaplanet test -1 to -8  
           mtimer_sav=mtimer
           mtimer=mins_gmt     ! so radn scheme repeatedly works thru same day
          endif    ! (nhstest<0)
c         print *,'before radrive'
c         call maxmin(t,' t',ktau,1.,kl)
c         call maxmin(qg,'qg',ktau,1.e3,kl)
c         call maxmin(tgg,'tg',ktau,1.,ms)
          call radrive (odcalc,iaero)
          if(nhstest<0)then ! aquaplanet test -1 to -8  
            mtimer=mtimer_sav
          endif    ! (nhstest<0)
          t(1:ifull,:)=t(1:ifull,:)-dt*rtt(1:ifull,:) 
          if (nmaxpr==1) then
             ! Account for load bal explicitly rather than implicitly in
             ! the reduce in maxmin.
             call phys_loadbal
             call maxmin(rtt,'rt',ktau,1.e4,kl)
             call maxmin(slwa,'sl',ktau,.1,1)
          end if
        else
!         use preset slwa array (use +ve nrad)
          slwa(:)=-10*nrad  
!         N.B. no rtt array for this nrad option
        endif  !  (nrad==4)

        egg(:)=0.   ! reset for fort.60 files
        fgg(:)=0.   ! reset for fort.60 files
        if(ntsur<=1.or.nhstest==2)then ! Held & Suarez or no surf fluxes
         eg(:)=0.
         fg(:)=0.
         cdtq(:)=0.
         cduv(:)=0.
        endif     ! (ntsur<=1.or.nhstest==2) 
        if(nhstest==2)call hs_phys
        if(ntsur>1)then  ! should be better after convjlm
          if(diag)then
            call maxmin(u,'#u',ktau,1.,kl)
            call maxmin(v,'#v',ktau,1.,kl)
            call maxmin(t,'#t',ktau,1.,kl)
            call maxmin(qg,'qg',ktau,1.e3,kl)     
            call MPI_Barrier( MPI_COMM_WORLD, ierr ) ! stop others going past
          endif
          call sflux(nalpha)
          epan_ave = epan_ave+epan  ! 2D 
          ga_ave = ga_ave+ga        ! 2D   
          if(nstn>0.and.nrotstn(1)==0)call stationa ! write every time step
          if(mod(ktau,nmaxpr)==0.and.mydiag)then
          print *
          write (6,
     .	   "('ktau =',i5,' gmt(h,m):',f6.2,i5,' runtime(h,m):',f7.2,i6)")
     .	      ktau,timeg,mins_gmt,timer,mtimer
!         some surface (or point) diagnostics
          isoil = isoilm(idjd)
          print *,'land,isoil,ivegt,isflag ',
     &           land(idjd),isoil,ivegt(idjd),isflag(idjd)
          write (6,"('snage,snowd,osnowd,alb,tsigmf   ',f8.4,4f8.2)")
     &       snage(idjd),snowd(idjd),osnowd(idjd),alb(idjd),tsigmf(idjd)
          write (6,"('sicedep,fracice,runoff ',3f8.2)")
     &             sicedep(idjd),fracice(idjd),runoff(idjd)
          write (6,"('t1,otgsoil,theta,fev,fgf   ',9f8.2)") 
     &         t(idjd,1),otgsoil(idjd),theta(idjd),fev(idjd),fgf(idjd)
          write (6,"('tgg(1-6)   ',9f8.2)") (tgg(idjd,k),k=1,6)
          write (6,"('tggsn(1-3) ',9f8.2)") (tggsn(idjd,k),k=1,3)
          write (6,"('wb(1-6)    ',9f8.3)") (wb(idjd,k),k=1,6)
          write (6,"('wbice(1-6) ',9f8.3)") (wbice(idjd,k),k=1,6)
          write (6,"('wblf(1-6)  ',9f8.3)") (wblf(idjd,k),k=1,6)
          write (6,"('wbfice(1-6)',9f8.3)") (wbfice(idjd,k),k=1,6)
          write (6,"('smass(1-3) ',9f8.2)") (smass(idjd,k),k=1,3) ! as mm of water
          write (6,"('ssdn(1-3)  ',9f8.2)") (ssdn(idjd,k),k=1,3)
          write (6,"('sdepth(1-3)',9f8.2)") (sdepth(idjd,k),k=1,3) ! as m of snow
          pwater=0.   ! in mm
          iq=idjd
          div_int=0.
          do k=1,kl
           qtot=qg(iq,k)+qlg(iq,k)+qfg(iq,k)
            pwater=pwater-dsig(k)*qtot*ps(idjd)/grav
            div(k)=(u(ieu(iq),k)/emu(ieu(iq))
     &             -u(iwu(iq),k)/emu(iwu(iq))  
     &             +v(inv(iq),k)/emv(inv(iq))
     &             -v(isv(iq),k)/emv(isv(iq))) 
     &              *em(iq)**2/(2.*ds)  *1.e6
            div_int=div_int-div(k)*dsig(k)
          enddo
          write (6,"('pwater,condc,condx,rndmax,rmc',9f8.3)")
     &       pwater,condc(idjd),condx(idjd),rndmax(idjd),rmc(idjd)
          write (6,"('wetfac,sno,evap,precc,precip',
     &       6f8.2)") wetfac(idjd),sno(idjd),evap(idjd),precc(idjd),
     &       precip(idjd)
          write (6,"('tmin,tmax,tscr,tss,tgf,tpan',9f8.2)")
     &       tminscr(idjd),tmaxscr(idjd),tscrn(idjd),tss(idjd),
     &       tgf(idjd),tpan(idjd)
          write (6,"('u10,ustar,pblh',9f8.2)")
     &       u10(idjd),ustar(idjd),pblh(idjd)
          write (6,"('rgg,rdg,sgflux,div_int,ps,qgscrn',5f8.2,f8.3)")
     &       rgg(idjd),rdg(idjd),sgflux(idjd),
     &       div_int,.01*ps(idjd),1000.*qgscrn(idjd)
          write (6,"('dew_,eg_,epot,epan,eg,fg,ga',9f8.2)") 
     &       dew_ave(idjd),eg_ave(idjd),epot(idjd),epan(idjd),eg(idjd),
     &       fg(idjd),ga(idjd)
          write (6,"('taftfhg,degdt,gflux,dgdtg,zo,cduv', 
     &       f6.3,f7.2,2f8.2,2f8.5)") taftfhg(idjd),degdt(idjd),
     &       gflux(idjd),dgdtg(idjd),zo(idjd),cduv(idjd)/vmod(idjd)
          rlwup=(1.-tsigmf(idjd))*rgg(idjd)+tsigmf(idjd)*rdg(idjd)
          write (6,"('slwa,rlwup,sint,sg,rt,rg    ',9f8.2)") 
     &       slwa(idjd),rlwup,sintsave(idjd),sgsave(idjd),
     &       rtsave(idjd),rgsave(idjd)
          write (6,"('cll,clm,clh,clt ',9f8.2)") 
     &       cloudlo(idjd),cloudmi(idjd),cloudhi(idjd),cloudtot(idjd)
          write (6,"('u10max,v10max,rhmin,rhmax   ',9f8.2)")
     &                u10max(iq),v10max(iq),rhminscr(iq),rhmaxscr(iq)
          write (6,"('kbsav,ktsav,convpsav ',2i3,f8.4,9f8.2)")
     &                kbsav(idjd),ktsav(idjd),convpsav(idjd)
          write (6,"('t   ',9f8.2/4x,9f8.2)") t(idjd,:)
          write (6,"('u   ',9f8.2/4x,9f8.2)") u(idjd,:)
          write (6,"('v   ',9f8.2/4x,9f8.2)") v(idjd,:)
          write (6,"('qg  ',3p9f8.3/4x,9f8.3)") qg(idjd,:)
          write (6,"('qf  ',3p9f8.3/4x,9f8.3)") qfg(idjd,:)
          write (6,"('ql  ',3p9f8.3/4x,9f8.3)") qlg(idjd,:)
          write (6,"('cfrac',9f8.3/5x,9f8.3)") cfrac(idjd,:)
          do k=1,kl
           es=establ(t(idjd,k))
           spmean(k)=100.*qg(idjd,k)*(ps(idjd)*sig(k)-es)/(.622*es)
           spmean(k)=100.*qg(idjd,k)*(ps(idjd)*sig(k)-es)/(.622*es)
          enddo
          nlx=min(nlv,kl-8)
          write (6,"('rh(nlx+) ',9f8.2)") (spmean(k),k=nlx,nlx+8)
          write (6,"('div(nlx+)',9f8.2)") (div(k),k=nlx,nlx+8)
          write (6,"('omgf ',9f8.3/5x,9f8.3)")   ! in Pa/s
     &              ps(idjd)*omgf(idjd,:)
          write (6,"('sdot ',9f8.3/5x,9f8.3)") sdot(idjd,1:kl)
          if(nextout>=4)write (6,"('xlat,long,pres ',3f8.2)")
     &     tr(idjd,nlv,ngas+1),tr(idjd,nlv,ngas+2),tr(idjd,nlv,ngas+3)
          endif  ! (mod(ktau,nmaxpr)==0.and.mydiag)
        endif   ! (ntsur>1)
        if(ntsur>=1)then ! calls vertmix but not sflux for ntsur=1
          call vertmix 
        endif  ! (ntsur>=1)

!     This is the end of the physics. The next routine makes the load imbalance
!     overhead explicit rather than having it hidden in one of the diagnostic
!     calls.
      call phys_loadbal
      call end_log(phys_end)

!     rml 16/02/06 call tracer_mass, write_ts
      if(ngas>0) then
        call tracer_mass(ktau,ntau) !also updates average tracer array
        call write_ts(ktau,ntau,dt)
      endif


      if(ndi==-ktau)then
        nmaxpr=1         ! diagnostic prints; reset 6 lines on
        if(ndi2==0)ndi2=ktau+40
!       nwt=1
      endif
      if(ktau==ndi2)then
         if(mydiag)print *,'reset nmaxpr'
         nmaxpr=nmaxprsav
!        nwt=nwtsav
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
!       call maxmin(tggsn,'tgg',ktau,1.,ms+3)
        call maxmin(tggsn,'tggsn',ktau,1.,3)
        call maxmin(tgg,'tg',ktau,1.,ms)
        call maxmin(tgf,'tf',ktau,1.,ms)
        call maxmin(tss,'ts',ktau,1.,1)
        call maxmin(pblh,'pb',ktau,1.,1)
        call maxmin(precip,'pr',ktau,1.,1)
        call maxmin(precc,'pc',ktau,1.,1)
        call maxmin(convpsav,'co',ktau,1.,1)
        call maxmin(sno,'sn',ktau,1.,1)        ! as mm during timestep
        call maxmin(snowflx,'sm',ktau,1.,1)    ! snow melt flux
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
     &             ( u(iq,k)**2 + v(iq,k)**2 )
           enddo
        enddo
        cllav=0.
        clmav=0.
        clhav=0.
        cltav=0.
        do iq=1,ifull
         cllav=cllav+wts(iq)*cloudlo(iq)
         clmav=clmav+wts(iq)*cloudmi(iq)
         clhav=clhav+wts(iq)*cloudhi(iq)
         cltav=cltav+wts(iq)*cloudtot(iq)   
        enddo

        ! All this combined into a single reduction
        temparray = (/ psavge, pslavge, preccavge, precavge, gke,
     &                 cllav, clmav,clhav, cltav /)
        call MPI_Reduce ( temparray, gtemparray, 9, MPI_REAL, MPI_SUM,0,
     &                  MPI_COMM_WORLD, ierr )
        if ( myid == 0 ) then
           print 97, gtemparray(1:5) ! psavge,pslavge,preccavge,precavge,gke
 97        format(' average ps, psl, precc, prec, gke: ',
     &           f10.2,f10.6,2f6.2,f7.2)
           print 971, gtemparray(6:9) ! cllav,clmav,clhav,cltav
 971       format(' global_average cll, clm, clh, clt: ',4f6.2)
           print 972,alph_p,alph_pm,delneg,delpos,alph_q
 972       format(' alph_p,alph_pm,delneg,delpos,alph_q: ',5f8.4)
        end if
        if ( mydiag ) then
           print 98,ktau, diagvals(ps)
!     &          ((ps(ii+(jj-1)*il),ii=id-4,id+4,4),jj=jd-4,jd+4,4)
 98        format(i7,' ps diag:',-2p9f7.1)
           if(t(idjd,kl)>258.)then
              print *,'t(idjd,kl) > 258. for idjd = ',idjd
              print 91,ktau,(t(idjd,k),k=kl-8,kl)
 91           format(i7,'    t',9f7.2)
              print 92,ktau,(sdot(idjd,k),k=kl-8,kl)
 92           format(i7,' sdot',9f7.3)
           endif                ! (t(idjd,kl)>258.)
        end if                  ! myid==0
      endif      ! (mod(ktau,nmaxpr)==0)

!     update diag_averages and daily max and min screen temps 
!     N.B. runoff is accumulated in sflux
      tmaxscr = max(tmaxscr,tscrn)
      tminscr = min(tminscr,tscrn)
      rhmaxscr = max(rhmaxscr,rhscrn)
      rhminscr = min(rhminscr,rhscrn)
      rndmax   = max(rndmax,condx)
      u10mx    = max(u10mx,u10)  ! for hourly scrnfile
      dew_ave  = dew_ave-min(0.,eg)    
      eg_ave = eg_ave+eg    
      fg_ave = fg_ave+fg     
      tscr_ave = tscr_ave+tscrn 
      qscrn_ave = qscrn_ave+qgscrn 
      do iq=1,ifull
        if(u10(iq)**2>u10max(iq)**2+v10max(iq)**2)then
          spare1(iq)=max(.001,sqrt(u(iq,1)**2+v(iq,1)**2))  ! speed lev 1
          u10max(iq)=u10(iq)*u(iq,1)/spare1(iq)
          v10max(iq)=u10(iq)*v(iq,1)/spare1(iq)
       endif
      enddo

!     rnd03 to rnd21 are accumulated in mm     
      if (myid==0) then
         print *,'ktau,mod,nper3hr ',
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
        eg_ave(:)    =    eg_ave(:)/min(ntau,nperavg)
        fg_ave(:)    =    fg_ave(:)/min(ntau,nperavg)
        ga_ave(:)    =    ga_ave(:)/min(ntau,nperavg)
        riwp_ave(:)  =  riwp_ave(:)/min(ntau,nperavg)
        rlwp_ave(:)  =  rlwp_ave(:)/min(ntau,nperavg)
        tscr_ave(:)  =  tscr_ave(:)/min(ntau,nperavg)
        qscrn_ave(:) = qscrn_ave(:)/min(ntau,nperavg)
        if(myid==0)
     &       print *,'ktau,koundiag,nperavg =',ktau,koundiag,nperavg
        sint_ave(:) = sint_ave(:)/max(koundiag,1)
        sot_ave(:)  =  sot_ave(:)/max(koundiag,1)
        soc_ave(:)  =  soc_ave(:)/max(koundiag,1)
        sgdn_ave(:) =  sgdn_ave(:)/max(koundiag,1)
        sgn_ave(:)  =  sgn_ave(:)/min(ntau,nperday)  ! because of solar fit
        rtu_ave(:)  =  rtu_ave(:)/max(koundiag,1)
        rtc_ave(:)  =  rtc_ave(:)/max(koundiag,1)
        rgdn_ave(:) =  rgdn_ave(:)/max(koundiag,1)
        rgn_ave(:)  =  rgn_ave(:)/max(koundiag,1)
        rgc_ave(:)  =  rgc_ave(:)/max(koundiag,1)
        cld_ave(:)  =  cld_ave(:)/max(koundiag,1)
        cll_ave(:)  =  cll_ave(:)/max(koundiag,1)
        clm_ave(:)  =  clm_ave(:)/max(koundiag,1)
        clh_ave(:)  =  clh_ave(:)/max(koundiag,1)
        cbas_ave(:) = 1.1-cbas_ave(:)/max(1.e-4,precc(:))  ! 1.1 for no precc
        ctop_ave(:) = 1.1-ctop_ave(:)/max(1.e-4,precc(:))  ! 1.1 for no precc
      endif    ! (ktau==ntau.or.mod(ktau,nperavg)==0)
      
      if(nscrn==1.and.mod(ktau,nperhr)==0)then
        call outcdfs(rundate)  ! for scrnfile
        u10mx(:)=0.
      endif
      if(ktau==ntau.or.mod(ktau,nwt)==0)then
        call log_off()
        call outfile(20,rundate,nmi,nwrite)
 
        if(ktau==ntau.and.irest==1) then
#ifdef simple_timer
          ! Don't include the time for writing the restart file
          call end_log(maincalc_end)
#endif
!         write restart file
          if(io_rest==2)open
     &      (19,file=restfile,form='formatted',status='unknown')
          if(io_rest==3)open
     &      (19,file=restfile,form='unformatted',status='unknown')
          call outfile(19,rundate,nmi,nwrite)
          if(myid==0)print *,'finished writing restart file in outfile'
          close(19)
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
           print 985,pwatr,preccavge,precavge,evapavge
        end if
985     format(' average pwatr,precc,prec,evap: ',4f7.3)
!       also zero most averaged fields every nperavg
        cbas_ave(:)=0.
        ctop_ave(:)=0.
        dew_ave(:)=0.
        epan_ave(:)=0.
        eg_ave(:)=0.
        fg_ave(:)=0.
        riwp_ave(:)=0.
        rlwp_ave(:)=0.
        qscrn_ave(:) = 0.
        tscr_ave(:) = 0.
        if (myid==0) print *,'resetting tscr_ave for ktau = ',ktau
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
      endif  ! (mod(ktau,nperavg)==0)

      if(mod(ktau,nperday)==0)then   ! re-set at the end of each 24 hours
        if(ntau<10*nperday.and.nstn>0)then     ! print stn info
          do nn=1,nstn
           if ( .not. mystn(nn) ) cycle
           i=istn(nn)
           j=jstn(nn)
           iq=i+(j-1)*il
           print 956,ktau,iunp(nn),name_stn(nn),
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
        rnd_3hr(:,8)=0.       ! i.e. rnd24(:)=0.
        if(nextout>=4)call setllp ! from Nov 11, reset once per day
        if(namip>0)then
          if (myid==0)
     &    print *,'amipsst called at end of day for ktau,mtimer,namip ',
     &                                              ktau,mtimer,namip
          call amipsst
        endif ! (namip>0)
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
         print*, "End of time loop ", timeval
      end if

      if(myid==0) then
         print *,'normal termination of run'
         call date_and_time(time=timeval)
         print*, "End time ", timeval
         print*, "Model time in main loop", 3600*(tvals2(5)-tvals1(5)) + 
     &      60*(tvals2(6)-tvals1(6)) + (tvals2(7)-tvals1(7)) + 
     &      0.001 * (tvals2(8)-tvals1(8))
      end if
#ifdef simple_timer
      call end_log(model_end)
      call simple_timer_finalize
#endif

#ifndef scyld
      call MPI_Finalize(ierr)
#endif
 
      end

      subroutine readint(filename,itss,ifullx)
      use cc_mpi
      include 'newmpar.h'
      include 'parm.h'
      character *(*) filename
      character header*47
      integer itss(ifullx)
      integer glob2d(ifull_g)

      if( ifullx /= ifull ) then
         print*, "Error, readint only works with ifull"
         print*, "called with", ifullx
         stop
      end if
      if ( myid == 0 ) then
         print *,'reading data via readint from ',filename
         open(87,file=filename,status='old')
!        read(87,'(i3,i4,2f6.1,f6.3,f8.0,a47)',iostat=ierr)
         read(87,*,iostat=ierr)
     &          ilx,jlx,rlong0x,rlat0x,schmidtx,dsx,header
         if ( ierr == 0 ) then
            print *,ilx,jlx,rlong0x,rlat0x,schmidtx,dsx,header
            if(ilx.ne.il_g.or.jlx.ne.jl_g.or.rlong0x.ne.rlong0.
     &           or.rlat0x.ne.rlat0.or.schmidtx.ne.schmidt)
     &           stop 'wrong data file supplied'
            read(87,*) glob2d
            close(87)
         else if ( ierr < 0 ) then ! Error, so really unformatted file
            close(87)
            print *,'now doing unformatted read'
            open(87,file=filename,status='old',form='unformatted')
            read(87) glob2d
            close(87)
         else ! ierr > 0
            stop 'End of file occurred in readint'
         end if
         call ccmpi_distribute(itss, glob2d)
         print*, trim(header), glob2d(id+(jd-1)*il_g)
      else
         call ccmpi_distribute(itss)
      end if
      if (mydiag) print*, itss(idjd)
      end subroutine readint

      subroutine readreal(filename,tss,ifullx)
      use cc_mpi
      include 'newmpar.h'
      include 'parm.h'
      character *(*) filename
      character header*47
      real tss(ifullx)
      real glob2d(ifull_g)
      integer ierr

      if( ifullx /= ifull ) then
         print*, "Error, readreal only works with ifull"
         print*, "called with", ifullx
         stop
      end if
      if ( myid == 0 ) then
         print *,'reading data via readreal from ',trim(filename)
         open(87,file=filename,status='old')
         read(87,*,iostat=ierr)
     &          ilx,jlx,rlong0x,rlat0x,schmidtx,dsx,header
         if ( ierr == 0 ) then
            print *,ilx,jlx,rlong0x,rlat0x,schmidtx,dsx,header
            if(ilx.ne.il_g.or.jlx.ne.jl_g.or.rlong0x.ne.rlong0.
     &           or.rlat0x.ne.rlat0.or.schmidtx.ne.schmidt)
     &           stop 'wrong data file supplied'
            read(87,*) glob2d
            close(87)
         else if ( ierr < 0 ) then ! Error, so really unformatted file
            close(87)
            print *,'now doing unformatted read'
            open(87,file=filename,status='old',form='unformatted')
            read(87) glob2d
            close(87)
         else
           write(0,*) "error in readreal",trim(filename),ierr
           stop
         end if
         call ccmpi_distribute(tss, glob2d)
         print*, glob2d(id+(jd-1)*il_g)
      else
         call ccmpi_distribute(tss)
      end if
      if (mydiag) print*, tss(idjd)
      end subroutine readreal

      subroutine setllp
      use cc_mpi
!     sets tr arrays for lat, long, pressure if nextout>=4 &(nllp>=3)
      include 'newmpar.h'
      include 'aalat.h'
      include 'arrays.h'  ! ts, t, u, v, psl, ps, zs
      include 'const_phys.h'
      include 'sigs.h'
      include 'tracers.h'  ! ngas, nllp, ntrac
      do k=1,klt
       do iq=1,ilt*jlt        
        tr(iq,k,min(ntracmax,ngas+1))=alat(iq)
        tr(iq,k,min(ntracmax,ngas+2))=along(iq)
        tr(iq,k,min(ntracmax,ngas+3))=.01*ps(iq)*sig(k)  ! in HPa
       enddo
      enddo
c     if(nmaxpr==1.and.mydiag)then
c        print *,'in setllp; alat,along,tr1-3 ',alat(idjd),along(idjd),
c    &    	tr(idjd,k,ngas+1),tr(idjd,k,ngas+2),tr(idjd,k,ngas+3)
c     endif
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

      subroutine setaxu    ! globpea code
      use cc_mpi
      implicit none
      include 'newmpar.h'
      include 'indices.h'
      include 'vecsuv.h'   ! vecsuv info
      include 'vecsuva.h'  ! vecsuva info
      real wcua(ifull),wcva(ifull),
     &     wcud(ifull+iextra),wcvd(ifull+iextra)
      integer iq
      do iq=1,ifull   !  convert ax,bx to staggered positions
       wcua(iq)=.5*(ax(ieu(iq))+ax(iq))
       wcud(iq)=    ax(ieu(iq))-ax(iq)
       wcva(iq)=.5*(bx(inv(iq))+bx(iq))
       wcvd(iq)=    bx(inv(iq))-bx(iq)
      enddo   ! iq loop
      call boundsuv(wcud,wcvd)
      do iq=1,ifull   !  store ax,bx at staggered positions
       axu(iq)=wcua(iq)-(wcud(ieu(iq))-wcud(iwu(iq)))/16.
       bxv(iq)=wcva(iq)-(wcvd(inv(iq))-wcvd(isv(iq)))/16.
      enddo   ! iq loop
      do iq=1,ifull   !  convert ay,by to staggered positions
       wcua(iq)=.5*(ay(ieu(iq))+ay(iq))
       wcud(iq)=    ay(ieu(iq))-ay(iq)
       wcva(iq)=.5*(by(inv(iq))+by(iq))
       wcvd(iq)=    by(inv(iq))-by(iq)
      enddo   ! iq loop
      call boundsuv(wcud,wcvd)
      do iq=1,ifull   !  store ay,by at staggered positions
       ayu(iq)=wcua(iq)-(wcud(ieu(iq))-wcud(iwu(iq)))/16.
       byv(iq)=wcva(iq)-(wcvd(inv(iq))-wcvd(isv(iq)))/16.
      enddo   ! iq loop
      do iq=1,ifull   !  convert az,bz to staggered positions
       wcua(iq)=.5*(az(ieu(iq))+az(iq))
       wcud(iq)=    az(ieu(iq))-az(iq)
       wcva(iq)=.5*(bz(inv(iq))+bz(iq))
       wcvd(iq)=    bz(inv(iq))-bz(iq)
      enddo   ! iq loop
      call boundsuv(wcud,wcvd)
      do iq=1,ifull   !  store az,bz at staggered positions
       azu(iq)=wcua(iq)-(wcud(ieu(iq))-wcud(iwu(iq)))/16.
       bzv(iq)=wcva(iq)-(wcvd(inv(iq))-wcvd(isv(iq)))/16.
      enddo   ! iq loop
      return
      end subroutine setaxu

      blockdata main_blockdata
      include 'newmpar.h'
      include 'dates.h' ! ktime,kdate,timer,timeg,mtimer
      include 'filnames.h'  ! list of files, read in once only
      include 'kuocom.h'
      include 'indices.h'
      include 'mapproj.h'
      include 'parm.h'
      include 'parmdyn.h'  ! nstag,epsp,epsu
      include 'parmhor.h'  ! mhint, m_bs
      include 'parmsurf.h' ! nplens
      include 'parmvert.h'
      include 'scamdim.h'   ! passes npmax=ifull
      include 'soil.h'   
      include 'soilv.h'
      include 'stime.h'
      include 'trcom2.h'  ! nstn,slat,slon,istn,jstn etc.
      include 'vvel.h'    ! sdot
      real rlai,hc,disp,usuh,coexp,zruffs,trans,rt0us,rt1usa,rt1usb,
     &    xgsmax,xjmax0
      common /canopy/
     &    rlai,hc,disp,usuh,coexp,zruffs,trans,rt0us,rt1usa,rt1usb,
     &    xgsmax(0:44),xjmax0(0:44)
      integer, parameter :: ijkij=ijk+ifull
      data sdot/ijkij*0./   ! for vvel.h
      common/leap_yr/leap  ! 1 to allow leap years
      common/nsib/nbarewet,nsigmf
      common/radnml/nnrad,idcld   ! darlam, clddia

!     for cardin
      data ia/1/,ib/3/,id/2/,ja/1/,jb/10/,jd/5/,nlv/1/,
     &     ndi/1/,ndi2/0/,nmaxpr/99/,          
     &     kdate_s/-1/,ktime_s/-1/,leap/0/,
     &     mbd/0/,nbd/0/,nbox/1/,kbotdav/4/,kbotu/0/           
     &     nud_p/0/,nud_q/0/,nud_t/0/,nud_uv/1/,nud_hrs/24/,nudu_hrs/0/
      
!     Dynamics options A & B      
      data m/5/,mex/30/,mfix/1/,mfix_qg/1/,mup/1/,nh/0/,nonl/0/,npex/0/
      data nritch_t/300/,nrot/1/,nxmap/0/,
     &     epsp/-15./,epsu/0./,epsf/0./,precon/-2900/,restol/5.e-7/
      data schmidt/1./,rlong0/0./,rlat0/90./,nrun/0/,nrunx/0/
!     Horiz advection options
      data ndept/1/,nt_adv/7/,mh_bs/4/
!     Horiz wind staggering options
c     data nstag/99/,nstagu/99/
      data nstag/-10/,nstagu/-1/
!     Vertical advection options
      data nvad/-4/,nvadh/2/,ntvdr/1/
!     Horizontal mixing options
      data khdif/2/,khor/-8/,nhor/-157/,nhorps/-1/,nhorjlm/1/
!     Vertical mixing options
      data nvmix/3/,nlocal/6/,nvsplit/2/,ncvmix/0/,lgwd/0/,ngwd/-5/
!     Cumulus convection options
      data nkuo/23/,sigcb/1./,sig_ct/.8/,rhcv/0./,rhmois/.1/,rhsat/1./,
     &     convfact/1.02/,convtime/.33/,shaltime/0./,      
     &     alflnd/1.15/,alfsea/1.05/,fldown/.6/,iterconv/3/,ncvcloud/0/,
     &     nevapcc/0/,nevapls/-4/,nuvconv/0/
     &     mbase/101/,mdelay/-1/,methprec/8/,nbase/-4/,detrain/.15/,
     &     entrain/0./,methdetr/2/,detrainx/0./,dsig2/.15/,dsig4/.4/
!     Shallow convection options
      data ksc/99/,kscsea/0/,kscmom/1/,sigkscb/.95/,sigksct/.8/,
     &     tied_con/2./,tied_over/0./,tied_rh/.75/
!     Other moist physics options
      data acon/.2/,bcon/.07/,qgmin/1.e-6/,rcm/.92e-5/,
     &     rcrit_l/.75/,rcrit_s/.85/ 
!     Radiation options
      data nrad/4/,ndiur/1/,idcld/1/
!     Diagnostic cloud options
      data ldr/1/,nclddia/5/,nstab_cld/0/,nrhcrit/10/,sigcll/.95/ 
      data cldh_lnd/95./,cldm_lnd/85./,cldl_lnd/75./  ! not used for ldr
      data cldh_sea/95./,cldm_sea/90./,cldl_sea/80./  ! not used for ldr
!     Soil, canopy, PBL options
      data nbarewet/0/,newrough/2/,nglacier/1/,
     &     nrungcm/-1/,nsib/3/,nsigmf/1/,
     &     ntaft/2/,ntsea/6/,ntsur/6/,av_vmod/.7/,tss_sh/1./,
     &     vmodmin/.2/,zobgin/.02/,charnock/.018/,chn10/.00125/
      data nurban/0/ ! MJT urban
      data newsoilm/0/,newztsea/1/,newtop/1/,nem/2/                    
      data snmin/.11/  ! 1000. for 1-layer; ~.11 to turn on 3-layer snow
!     Special and test options
      data namip/0/,amipo3/.false./,nhstest/0/,nsemble/0/,nspecial/0/,
     &     panfg/4./,panzo/.001/,nplens/0/
!     I/O options
      data nfly/2/,io_in/1/,io_out/1/,io_rest/1/,nperavg/-99/,nwt/-99/
      data io_clim/1/,io_spec/0/,nextout/3/,localhist/.false./
      
      data nstn/0/  
      data slat/nstnmax*-89./,slon/nstnmax*0./,iunp/nstnmax*6/,
     &     zstn/nstnmax*0./,name_stn/nstnmax*'   '/,nrotstn/nstnmax*0/  

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
      data climcdf/'clim.cdf'/
      data monfil/'monthly.cdf'/,scrfcdf/'scrave.cdf'/
c     floating point:
      data timer/0./,mtimer/0/
      data du/19.76/,tanl/60.2/,rnml/130./,stl1/40./,stl2/10./

c     stuff from indata in soil.h
      data zoland/.16/

c     stuff from insoil  for soilv.h
      data rlaim44/4.8, 6.3, 5., 3.75, 2.78, 2.5, 3.9, 2.77, 2.04, 2.6,  ! 1-10
     &          1.69, 1.9, 1.37, 1.5, 1.21, 1.58, 1.41, 2.3, 1.2, 1.71,  ! 11-20
     &          1.21, 2.3, 2.3, 1.2, 1.2, 1.87, 1., 3., .01, .01, 1.2,   ! 21-31
     &          6., 5.5, 5., 4.5, 5., 4., 3., 3.5, 1., 4., .5, 4., 0./   ! 32-44
c    &          1.21, 2.3, 2.3, 1.2, 1.2, 1.87, 1., 3., .01, .01, .01,   ! 21-31
c    &          6., 6., 6., 6., 6., 4., 4., 4., 1., 4., .5, 4., 0./      ! 32-44
      data rlais44/1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,              ! 1-10
     &           1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,                ! 11-20
     &           1., 1., 1., 1., .6, .6, .5, 1., 0., 0., 1.,            ! 21-31
     &           2., 2., 2., 2., 2., 1.5, 1.5, 1.5, 1., .5, .5, .5, 0./ ! 32-44
c    &           1., 1., 1., 1., .6, .6, .5, 1., 0., 0., 0.,            ! 21-31
      data rsunc44/
     &      370., 330., 260., 200., 150., 130., 200., 150., 110., 160.,  ! 1-10
     &      100., 120.,  90.,  90.,  80.,  90.,  90., 150.,  80., 100.,  ! 11-20
     &       80.,  80.,  80.,  60.,  60., 120.,  80., 180., 2*995., 80., ! 21-31
     &      350., 4*300., 3*230., 150., 230., 995., 150., 9900./         ! 32-44
c    &       80.,  80.,  80.,  60.,  60., 120.,  80., 180., 3*995.,      ! 21-31
c    &      350., 4*300., 5*230., 995., 150., 9900./    ! eak: 8/12/98
      data scveg44/0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,              ! 1-10
     &             0., 0., 0., 0., 0., .1, .1, .1, .1, .1,              ! 11-20
     &             .1, .2, .4, .2, .1, .1, .1, 0., 0., 0., 0.,          ! 21-31
     &         .05, 0., 0., 0., 0., .05, .05, .05, .1, 0., 0., .4, 0./  ! 32-44
      data slveg44/0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,           ! 1-10
     &             0., 0., 0., 0., 0., .1, .1, .1, .1, .1,           ! 11-20
     &             .1, .2, .4, .2, .1, .1, .1, 0., 0., 0., 0.,       ! 21-31
     &     1., 5.5, 3., 1., 3., 3., 3.5, 3., .5, 3.5, .1, 3.5, 0./   ! 32-44
      data froot/.05, .10, .35, .40, .10/       ! 10/02/99 veg. root distr.
c     data froot/.20, .45, .20, .10, .05/

      data silt/.08, .33, .17, .2, .06, .25, .15, .70, .33
     &          , .2, .33, .33, .17/                        ! with mxst=13
      data clay/.09, .3, .67, .2, .42, .48, .27, .17, .30
     &          , .2, .3, .3, .67/                          ! with mxst=13
      data sand/.83, .37, .17, .6, .52, .27, .58, .13, .37
     &          , .6, .37, .37, .17/                        ! with mxst=13
      data swilt/0., .072, .216, .286, .135, .219, .283, .175, .395  !eak
!     data swilt/0., .010, .1  , .138, .135, .219, .283, .175, .395  !mm5
     &             , .216, .1142, .1547, .2864, .2498/
      data sfc/1.,  .143, .301, .367, .218, .31 , .37 , .255, .45 
     &            , .301, .22 , .25 , .367, .294/
      data ssat/2., .398, .479, .482, .443, .426, .482, .420, .451  !jlm
!     data ssat/2., .398, .479, .482, .443, .426, .482, .420, .450  !eak
!     data ssat/2., .339, .470, .468, .443, .426, .482, .420, .450  !mm5
     &          , .479, .435, .451, .482, .476/
      data bch/4.2, 7.1, 11.4, 5.15, 10.4, 10.4, 7.12, 5.83, 7.1, 4.9
     &          , 5.39, 11.4, 8.52/      ! bch for gravity term
      data hyds/166.e-6, 4.e-6, 1.e-6, 21.e-6, 2.e-6, 1.e-6, 6.e-6,
     *             800.e-6, 1.e-6, 34.e-6, 7.e-6, 1.3e-6, 2.5e-6/  ! ksat
      data sucs/-.106, -.591, -.405, -.348, -.153, -.49, -.299,
     &          -.356, -.153, -.218, -.478, -.405, -.63/           ! phisat (m)
      data rhos/7*2600., 1300.,  910., 4*2600./      ! soil density
      data  css/7* 850., 1920., 2100., 4*850./       ! heat capacity

      data zse/.022, .058, .154, .409, 1.085, 2.872/ ! layer thickness
!     data zse/.05, .15, .33, 1.05, 1.25, 1.35/      ! layer thickness <10/2/99
!     data zse/.05, .15, .30, 0.50, 1.0 , 1.5 /      ! was in indata

c!     following was set in setxyz -  now in blockdata at end of setxyz
c      data npann/  1, 2,107,  4,106,  6,  7,109,  9,112, 11, 12,102,101/
c      data npane/103, 3,  4,105,  5,110,108,  8, 10, 11,100,113, 13,  0/
c      data npanw/13,113,112,  1,  2,  4,104,102,  7,107,  8,  9,109, 12/
c      data npans/110, 0,  1,100,  3,103,  5,  6,106,  8,105, 10, 11,111/
c Leaf gsmax for forest (0.006), grass (0.008) and crop (0.012)
c littoral is regarded as forest, dense pasture between grass and crop
      data xgsmax  / 0.0,
     &     0.006,0.006,0.006,0.006,0.006,0.006,0.006,0.006,0.006,0.006,
     &     0.006,0.006,0.006,0.006,0.006,0.008,0.008,0.008,0.008,0.008,
     &     0.008,0.010,0.010,0.012,0.012,0.008,0.008,0.006,0.000,0.000,
     &     0.000,
     &  .006,.006,.006,.006,.006,.006,.008,.008,.006,.006,.0,0.010,0.0/
c littoral is regarded as forest, dense pasture between grass and crop
      data xjmax0 / 0.0,
     &     5e-5,5e-5,5e-5,5e-5,5e-5,5e-5,5e-5,5e-5,5e-5,5e-5,
     &     5e-5,5e-5,5e-5,5e-5,5e-5,10e-5,10e-5,10e-5,10e-5,10e-5,
     &     10e-5,15e-5,15e-5,15e-5,15e-5,10e-5,10e-5,5e-5,1e-5,1e-5,
     &     1e-5,
     &     5e-5,5e-5,5e-5,5e-5,5e-5,10e-5,10e-5,10e-5,5e-5,10e-5,
     &     1e-5,15e-5,1e-5/
c     &     25e-5,25e-5,20e-5,15e-5,25e-5,7e-5,7e-5,7e-5,15e-5,15e-5,
c     &     7e-5,25e-5,1e-5/ !Sellers 1996 J.Climate, I think they are too high

      end
      subroutine stationa
      use cc_mpi
      use diag_m
      implicit none
      include 'newmpar.h'
      include 'arrays.h'   ! ts, t, u, v, psl, ps, zs
      include 'const_phys.h'
      include 'dates.h'    ! dtin,mtimer
      include 'establ.h'
      include 'extraout.h'
      include 'map.h'      ! em, f, dpsldt, fu, fv, etc
      include 'nsibd.h'    ! sib, tsigmf, isoilm
      include 'morepbl.h'
      include 'parm.h'
      include 'pbl.h'      ! cduv, cdtq, tss, qg
      include 'prec.h'
      include 'screen.h'   ! tscrn etc
      include 'sigs.h'
      include 'soil.h'
      include 'soilsnow.h'
      include 'soilv.h'
      include 'tracers.h'  ! ngas, nllp, ntrac, tr
      include 'trcom2.h'   ! nstn,slat,slon,istn,jstn etc.
      include 'vecsuv.h'
      include 'xyzinfo.h'  ! x,y,z,wts
      common/leap_yr/leap  ! 1 to allow leap years
      real, dimension(ifull) :: dirad,dfgdt,degdt,wetfac,degdw,cie,
     &                          factch,qsttg,rho,zo,aft,fh,spare1,theta,
     &                          gamm,rg,vmod,dgdtg
      common/work2/dirad,dfgdt,degdt,wetfac,degdw,cie,factch,qsttg,rho,
     &     zo,aft,fh,spare1,theta,gamm,rg,vmod,dgdtg

      real, dimension(ifull) :: egg,evapxf,ewww,fgf,fgg,ggflux,rdg,rgg,
     &                          residf,ga,condxpr,fev,fes,fwtop,spare2
      integer, dimension(ifull) :: ism
      real, dimension(ifull,kl) :: dum3a,speed,spare
      real, dimension(2*ijk-16*ifull) :: dum3
      common/work3/egg,evapxf,ewww,fgf,fgg,ggflux,rdg,rgg,residf,ga,
     &     condxpr,fev,fes,ism,fwtop,spare2,dum3,dum3a,speed,spare

      real wblf,wbfice,sdepth,dum3b
      common/work3b/wblf(ifull,ms),wbfice(ifull,ms),sdepth(ifull,3),
     &              dum3b(ijk*2-2*ifull*ms-3*ifull)
      integer i,j,iq,ico2x,iqt,iradonx,isoil,k2,leap,nn
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
!            Check if this station is in this processors region
             if ( .not. mystn(nn) ) cycle 
             if(ktau==1)write (iunp(nn),950) kdate,ktime,leap
950          format("#",i9,2i5)
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
     &        +zse(4)*wb(iq,4))/(zse(1)+zse(2)+zse(3)+zse(4))
             ico2x=max(1,ico2)
             iradonx=max(1,iradon)
             iqt = min(iq,ilt*jlt) ! Avoid bounds problems if there are no tracers
             k2=min(2,klt)
             write (iunp(nn),951) ktau,tscrn(iq)-273.16,rnd_3hr(iq,8),
     &         tss(iq)-273.16,tgg(iq,1)-273.16,tgg(iq,2)-273.16,
     &         tgg(iq,3)-273.16,t(iq,1)-273.16,tgf(iq)-273.16,
     &         wb(iq,1),wb(iq,2),
     &         cloudlo(iq),cloudmi(iq)+1.,cloudhi(iq)+2.,
     &         cloudtot(iq)+3.,
     &         fg(iq),eg(iq),(1.-tsigmf(iq))*fgg(iq),
     &         (1.-tsigmf(iq))*egg(iq),rnet(iq),sgsave(iq),
     &         qg(iq,1)*1.e3,uzon,vmer,precc(iq),
     &         qg(iq,2)*1.e3,rh1,rh2,tr(iqt,1,ico2x),tr(iqt,k2,ico2x),
     &         tr(iqt,1,iradonx),tr(iqt,k2,iradonx) ,.01*ps(iq),wbav,
     &         epot(iq),qgscrn(iq)*1.e3,rh_s,u10(iq),uscrn(iq),
     &         condx(iq)
!              N.B. qgscrn formula needs to be greatly improved
951          format(i4,6f7.2, 
     &              2f7.2, 2f6.3, 4f5.2,                ! t1 ... cld
     &              5f7.1,f6.1,f5.1,                    ! fg ... qg1
     &              2f6.1,f7.2, f5.1,2f6.1, 2(1x,f5.1), ! uu ... co2_2
     &              2(1x,f5.1) ,f7.1,f6.3,f7.1,5f6.1,   ! rad_1 ... rh_s
     &              f7.2)                               ! condx
             if(ktau==ntau)then
               write (iunp(nn),952)
952            format("#   2tscrn 3precip 4tss  5tgg1  6tgg2  7tgg3",
     &       "   8t1    9tgf  10wb1 11wb2 cldl cldm cldh  cld",
     &       "   16fg   17eg  18fgg  19egg  20rnet 21sg 22qg1",
     &       " 23uu   24vv 25precc qg2  rh1 28rh2 29co2_1 co2_2",
     &       " rad_1 rad_2  ps 34wbav 35epot qgscrn 37rh_s 38u10 uscrn",
     &       " 40condx")
              write (iunp(nn),953) land(iq),isoilm(iq),ivegt(iq),zo(iq),
     &                              zs(iq)/grav
953           format("# land,isoilm,ivegt,zo,zs/g: ",l2,2i3,2f9.3)
              isoil=max(1,isoilm(iq))
               write (iunp(nn),954) sigmf(iq),swilt(isoil),sfc(isoil),
     &                              ssat(isoil),alb(iq)
954            format("#sigmf,swilt,sfc,ssat,alb: ",5f7.3)
!              rml 16/02/06 removed ico2em
               write (iunp(nn),955) i,j,radonem(iqt)
955            format("#i,j,radonem: ",2i4,f7.3)
             endif
           enddo
      return	    
      end               
      subroutine stationb  ! primarily for ICTS
      use cc_mpi
      use diag_m
      implicit none
      include 'newmpar.h'
      include 'arrays.h'   ! ts, t, u, v, psl, ps, zs
      include 'const_phys.h'
      include 'dates.h'    ! dtin,mtimer
      include 'establ.h'
      include 'extraout.h'
      include 'histave.h'
      include 'liqwpar.h'  ! ifullw,qfg,qlg
      include 'map.h'      ! em, f, dpsldt, fu, fv, etc
      include 'nsibd.h'    ! sib, tsigmf, isoilm
      include 'morepbl.h'
      include 'parm.h'
      include 'pbl.h'      ! cduv, cdtq, tss, qg
      include 'prec.h'
      include 'raddiag.h'
      include 'screen.h'   ! tscrn etc
      include 'sigs.h'
      include 'soil.h'
      include 'soilsnow.h'
      include 'soilv.h'
      include 'tracers.h'  ! ngas, nllp, ntrac, tr
      include 'trcom2.h'   ! nstn,slat,slon,istn,jstn etc.
      include 'vecsuv.h'
      include 'xyzinfo.h'  ! x,y,z,wts
      real cfrac
      common/cfrac/cfrac(ifull,kl)     ! globpe,radriv90,vertmix,convjlm
      real sgx,sgdnx,rgx,rgdnx,soutx,sintx,rtx
      common/radstuff/sgx(ifull),sgdnx(ifull),rgx(ifull),rgdnx(ifull),
     &             soutx(ifull),sintx(ifull),rtx(ifull)
      real p(ifull,kl),tv(ifull,kl),energy(ifull,kl)
      equivalence (tv,energy)
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
      print *,'new kdateb,ktimeb,mtimerb ',kdateb,ktimeb,mtimerb
      print *,'new iyr,imo,iday,itim,off,sca ',iyr,imo,iday,itim,off,sca
      do nn=1,nstn
c      print *,'nn,iunp,mystn ',nn,iunp(nn),mystn(nn)
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
c 	   print *,'nn,niq,costh,sinth ',nn,niq,costh(niq),sinth(niq)
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
c	    if(nn<3.and.niq==1)then
c             print *,'k,npres,pres,press ',k,npres,pres(npres),press(k)
c             print *,'fa,fb,pp,pk,pk-1 ',
c    &          fa,fb,pp(niq,npres),p(niq,k),p(niq,k-1)
c           endif
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
c        if(nn<6)print *,'energint ',energint
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
         write (iunp(nn),"('SWUP_S    ',i6,3i3,'   1',f8.0,g8.1,9f8.2)") 
     &	         iyr,imo,iday,itim,off,sca,sgdnx(iqq(:))-sgx(iqq(:)) 
         write (iunp(nn),"('SWDOWN_S  ',i6,3i3,'   1',f8.0,g8.1,9f8.2)") 
     &	         iyr,imo,iday,itim,off,sca,sgdnx(iqq(:))
         write (iunp(nn),"('LWUP_S    ',i6,3i3,'   1',f8.0,g8.1,9f8.2)") 
     &	         iyr,imo,iday,itim,off,sca,rgdnx(iqq(:))+rgx(iqq(:))
         write (iunp(nn),"('LWDOWN_S  ',i6,3i3,'   1',f8.0,g8.1,9f8.2)")
     &	         iyr,imo,iday,itim,off,sca,rgdnx(iqq(:))
         write (iunp(nn),"('SWUP_T    ',i6,3i3,'   1',f8.0,g8.1,9f8.2)") 
     &	         iyr,imo,iday,itim,off,sca,soutx(iqq(:))  
         write (iunp(nn),"('SWDOWN_T  ',i6,3i3,'   1',f8.0,g8.1,9f8.2)")
     &	         iyr,imo,iday,itim,off,sca,sintx(iqq(:))
         write (iunp(nn),"('LWUP_T    ',i6,3i3,'   1',f8.0,g8.1,9f8.2)")
     &	         iyr,imo,iday,itim,off,sca,rtx(iqq(:))  
         write (iunp(nn),"('SHFL_S    ',i6,3i3,'   1',f8.0,g8.1,9f8.2)")
     &	         iyr,imo,iday,itim,off,sca,fg(iqq(:))
         write (iunp(nn),"('LHFL_S    ',i6,3i3,'   1',f8.0,g8.1,9f8.2)")
     &	         iyr,imo,iday,itim,off,sca,eg(iqq(:))
         write (iunp(nn),"('SMHFL     ',i6,3i3,'   1',f8.0,g8.1,9f8.2)") 
     &	         iyr,imo,iday,itim,off,sca,snowflx(iqq(:))
      enddo   ! nn=1,nstn
      return    
      end               
