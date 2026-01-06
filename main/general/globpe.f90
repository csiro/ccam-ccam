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

!      PE model on conformal-cubic grid
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!      input files are :namelist (via file called "input")
!                       "nrun.dat"
!      data input and output file names are specified in namelist 'datafile'
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!     sign convention:
!                      u+ve eastwards  (on the panel)
!                      v+ve northwards (on the panel)

! Preprocessor directives:
!   CCAM         - support CCAM (required)
!   debug        - additional debugging checks, but runs slower
!   i8r8         - double precision mode
!   GPU          - target GPUs with OpenACC
!   GPUPHYSICS   - target GPUs for physical parameterisations with OpenACC.  Requires -DGPU.    
!   csircoupled  - CSIR coupled model
!   usempi3      - optimse communication with MPI shared memory (preferred)
!   share_ifullg - reduce shared memory with MPI, but requires usempi3 directive
!   vampir       - enable vampir profiling

program globpe

use aerointerface                          ! Aerosol interface
use aerosol_arrays                         ! Aerosol arrays
use amipsst_m                              ! AMIP SSTs
use arrays_m                               ! Atmosphere dyamics prognostic arrays
use bigxy4_m                               ! Grid interpolation
use cc_mpi                                 ! CC MPI routines
use cfrac_m                                ! Cloud fraction
use config_m                               ! CCAM parameter initialisation
use const_phys                             ! Physical constants
use dates_m                                ! Date data
use daviesnudge                            ! Far-field nudging
use diag_m                                 ! Diagnostic routines
use dpsdt_m                                ! Vertical velocity
use ensemble                               ! Ensemble
use epst_m                                 ! Off-centre terms
use estab                                  ! Liquid saturation function
use extraout_m                             ! Additional diagnostics
use filnames_m                             ! Filenames
use gdrag_m, only : gwdrag                 ! Gravity wave drag
use histave_m                              ! Time average arrays
use hordifg_m                              ! Horizontal diffusion
use hs_phys_m                              ! Held & Suarez
use indata                                 ! Data initialisation
use indices_m                              ! Grid index arrays
use infile                                 ! Input file routines
use kuocom_m                               ! JLM convection
use liqwpar_m                              ! Cloud water mixing ratios
use map_m                                  ! Grid map arrays
use mlodynamics                            ! Ocean dynamics
use module_ctrl_convection                 ! Interface for convection
use module_ctrl_microphysics               ! Interface for cloud microphysics
use module_ctrl_turbmix                    ! Boundary layer turbulent mixing
use morepbl_m                              ! Additional boundary layer diagnostics
use nesting                                ! Nesting and assimilation
use newmpar_m                              ! Grid parameters
use nharrs_m                               ! Non-hydrostatic atmosphere arrays
use nlin_m                                 ! Atmosphere non-linear dynamics
use outcdf                                 ! Output file routines
use parm_m                                 ! Model configuration
use parmdyn_m                              ! Dynamics parameters
use parmhdff_m                             ! Horizontal diffusion parameters
use pbl_m                                  ! Boundary layer arrays
use prec_m                                 ! Precipitation
use raddiag_m                              ! Radiation diagnostic
use river                                  ! River routing
use savuvt_m                               ! Saved dynamic arrays
use savuv1_m                               ! Saved dynamic arrays
use sbar_m                                 ! Saved dynamic arrays
use screen_m                               ! Screen level diagnostics
use scrnout_m                              ! Calculate diagnostics
use seaesfrad_m                            ! SEA-ESF radiation
use sflux_m                                ! Surface flux routines
use sigs_m                                 ! Atmosphere sigma levels
use soil_m                                 ! Soil and surface data
use soilsnow_m                             ! Soil, snow and surface data
use timeseries, only : write_ts            ! Tracer time series
use tracermodule, only : tracer_mass     & ! Tracer routines
   ,interp_tracerflux
use tracers_m                              ! Tracer data
use trvmix                                 ! Tracer mixing routines
use uvbar_m                                ! Saved dynamic arrays
use vvel_m                                 ! Additional vertical velocity
use work2_m                                ! Diagnostic arrays
use work3f_m                               ! Grid work arrays
use xarrs_m                                ! Saved dynamic arrays
use xyzinfo_m                              ! Grid coordinate arrays

#ifdef csircoupled
use vcom_ccam                              ! CSIR (SA) ocean model
#endif

implicit none
     
integer, dimension(8) :: tvals1, tvals2, nper3hr
integer, dimension(8) :: times_total_a, times_total_b
integer iq, k
integer mins_gmt, mspeca, mtimer_in
integer nlx, nmaxprsav, n3hr
integer nwtsav, mtimer_sav
integer jyear, jmonth, jday, jhour, jmin, mins
integer koundiag
real, dimension(:), allocatable, save :: spare1
real, dimension(3) :: temparray, gtemparray
real aa, bb, cc
real hourst, evapavge, precavge
real pwatr, bb_2, cc_2, rat
real(kind=8) :: tt_r8
logical oxidant_update, ltest
character(len=10) timeval


! Start model timer
call date_and_time(values=times_total_a)

! Compile options tests
#ifdef i8r8
if ( kind(iq)/=8 .or. kind(aa)/=8 ) then
  write(6,*) "ERROR: CCAM compiled for double precision, but single precision code was detected"
  stop
end if
#else
if ( kind(iq)/=4 .or. kind(aa)/=4 ) then
  write(6,*) "ERROR: CCAM compiled for single precision, but double precision code was detected"
  stop
end if
#endif

#ifdef share_ifullg
#ifndef usempi3
write(6,*) "ERROR: Compiling CCAM with -Dshare_ifullg requires -Dusempi3"
#endif
#endif


!--------------------------------------------------------------
! INITALISE MPI ROUTINES
call ccmpi_init

! Start banner
if ( myid==0 ) then
  write(6,*) "=============================================================================="
  write(6,*) "CCAM: Starting globpea"
  write(6,*) "=============================================================================="
end if


!----------------------------------------------------------------
! INITALISE TIMING LOGS, READ NAMELIST AND READ INITIAL CONDITIONS
call log_off
call log_setup
call globpe_init


!--------------------------------------------------------------
! OPEN OUTPUT FILES AND SAVE INITAL CONDITIONS
if ( nwt>0 ) then
  ! write out the first ofile data set
  if ( myid==0 ) write(6,*) "Calling outfile"
  call outfile(20,ofile,psl,u,v,t,qg)
end if    ! (nwt>0)
! CORDEX and 10min output do not currently write the zeroth time-step
!if ( surfile/=' ' ) then
!  call freqfile_cordex
!end if
!if ( freqfile/=' ' ) then
!  call freqfile_10
!end if
if ( newtop<0 ) then
  ! just for outcdf to plot zs  & write fort.22      
  if ( myid==0 ) write(6,*) "newtop<0 requires a stop here"
  call ccmpi_abort(-1)
end if


!-------------------------------------------------------------
! SETUP DIAGNOSTIC ARRAYS
allocate( spare1(ifull) )
do n3hr = 1,8
  nper3hr(n3hr) = nint(real(n3hr)*3.*3600./dt)
end do
n3hr = 1
nlx = 0                      ! diagnostic level
call zero_nperavg(koundiag)  ! reset average period diagnostics
call zero_nperhour           ! reset hourly period diagnostics
call zero_nperday            ! reset daily period diagnostics


!--------------------------------------------------------------
! INITIALISE DYNAMICS
if ( myid==0 ) then
  write(6,*) "number of time steps per day = ",nperday
end if
! use half time-step for initialisation
dtin = dt
mspeca = 1
if ( mex/=1 .and. ((.not.lrestart).or.always_mspeca) ) then
  mspeca = 2
  dt = 0.5*dtin
end if
call gettin(0) ! preserve initial mass & T fields


!--------------------------------------------------------------
! SET-UP TIMERS
mtimer_sav = 0      ! saved value for minute timer
nmaxprsav  = nmaxpr
nwtsav     = nwt
hourst     = real(nint(0.01*real(ktime))) + real(mod(ktime,100))/60. ! for tracers
mtimer_in  = mtimer
 

!--------------------------------------------------------------
! BEGIN MAIN TIME LOOP
if ( myid==0 ) then
  call date_and_time(time=timeval,values=tvals1)
  write(6,*) "Start of loop time ", timeval
end if
call log_on
call START_LOG(maincalc_begin)

do ktau = 1,ntau   ! ****** start of main time loop

  timer    = real(ktau)*dtin/3600.                     ! timer now only used to give timeg
  timeg    = mod( timer+hourst, 24. )                  ! UTC time for tracers
  mtimer   = mtimer_in + nint(real(ktau)*dtin/60.)     ! to allow dt < 1 minute
  mins_gmt = mod( mtimer+60*ktime/100, 1440 )          ! for radiation
  call getzinp(jyear,jmonth,jday,jhour,jmin,mins)      ! define mins as time since start of the year
  diag = ( ktau>=abs(ndi) .and. ktau<=ndi2 )           ! set diagnostic printout flag
  if ( ndi<0 ) then
    if ( mod(ktau,ndi)==0 ) then
      diag = .true.
    end if
  endif


  ! calculate MSE before dynamic (e.g. start of the advection)
  call calculate_dhdt_mse(1,ifull,mse_t1)

  
  ! ***********************************************************************
  ! ATMOSPHERE DYNAMICS
  ! ***********************************************************************

  call nantest("before atmosphere dynamics",1,ifull,"all")

  ! NESTING ---------------------------------------------------------------
  if ( nbd/=0 ) then
    ! Newtonian relaxiation
    call START_LOG(nestin_begin)
    call nestin
    call END_LOG(nestin_end)
    call nantest("after nesting",1,ifull,"nesting")
  end if


  ! DYNAMICS --------------------------------------------------------------
  if ( nstaguin>0 .and. ktau>=1 ) then   ! swapping here for nstaguin>0
    if ( nstagin<0 .and. mod(ktau-nstagoff,abs(nstagin))==0 ) then
      nstag  = 7 - nstag  ! swap between 3 & 4
      nstagu = nstag
    end if
  end if

  do mspec = mspeca,1,-1    ! start of introductory time loop

    un(1:ifull,:) = 0. 
    vn(1:ifull,:) = 0.
    tn(1:ifull,:) = 0.

    if ( mup/=1 .or. (ktau==1.and.mspec==mspeca.and.((.not.lrestart).or.always_mspeca)) ) then
      call bounds(psl)
      ! updps called first step or to permit clean restart option      
      call updps(0) 
    end if
    
    if ( ktau<10 .and. nmaxpr==1 ) then
      if ( myid==0 ) then
        write(6,*) 'ktau,mex,mspec,mspeca:',ktau,mex,mspec,mspeca
      end if
    end if
    
    ! set up tau +.5 velocities in ubar, vbar
    sbar(:,2:kl) = sdot(:,2:kl)
    if ( (ktau==1.and.((.not.lrestart).or.always_mspeca)) .or. mex==1 ) then
      ubar(1:ifull,1:kl) = u(1:ifull,1:kl)
      vbar(1:ifull,1:kl) = v(1:ifull,1:kl)
    else if ( (ktau==2.and.((.not.lrestart).or.always_mspeca)) .or. mex==2 ) then        
      ! (tau+.5) from tau, tau-1
      ubar(1:ifull,1:kl) = u(1:ifull,1:kl)*1.5 - savu(1:ifull,1:kl)*.5
      vbar(1:ifull,1:kl) = v(1:ifull,1:kl)*1.5 - savv(1:ifull,1:kl)*.5
    else if ( mex==3 )then
      ! (tau+.5) from tau, tau-1, tau-2   ! ubar is savu1 here
      ubar(1:ifull,1:kl) = u(1:ifull,1:kl)+.5*(savu(1:ifull,1:kl)-savu1(1:ifull,1:kl))
      vbar(1:ifull,1:kl) = v(1:ifull,1:kl)+.5*(savv(1:ifull,1:kl)-savv1(1:ifull,1:kl))
    else if ( mex==30 ) then  ! using tau, tau-1, tau-2, tau-3
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
          bb = 1.5*v(iq,k) - 2.*savv(iq,k) + .5*savv1(iq,k)                             ! simple b
          bb_2 = (40.*v(iq,k)-35.*savv(iq,k)-16.*savv1(iq,k)+11.*savv2(iq,k))/34.       ! cwqls b
          cc = .5*v(iq,k) - savv(iq,k) + .5*savv1(iq,k)                                 ! simple c
          cc_2 = (10.*v(iq,k)-13.*savv(iq,k)-4.*savv1(iq,k)+7.*savv2(iq,k))/34.         ! cwqls c
          aa = cc_2 - cc
          rat = max( 0., min( 1., cc_2/(aa+sign(1.e-9,aa)) ) )
          cc = rat*cc + (1.-rat)*cc_2 
          bb = rat*bb + (1.-rat)*bb_2 
          vbar(iq,k) = v(iq,k)+.5*bb+.25*cc
        end do ! iq loop
      end do   ! k loop 
    else       ! i.e. mex >=4 and ktau>=3
      ! (tau+.5) from tau, tau-1, tau-2   ! ubar is savu1 here
      ubar(1:ifull,1:kl) = (u(1:ifull,1:kl)*15.-savu(1:ifull,1:kl)*10.+savu1(1:ifull,1:kl)*3.)/8.
      vbar(1:ifull,1:kl) = (v(1:ifull,1:kl)*15.-savv(1:ifull,1:kl)*10.+savv1(1:ifull,1:kl)*3.)/8.
    end if     ! (ktau==1) .. else ..
      
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
          write(6,*) "using epsp= ",epsp
        end if
      end if
      where ( dpsdt(1:ifull)*dpsdtb(1:ifull)<0. .and. dpsdtbb(1:ifull)*dpsdtb(1:ifull)<0. )
        epst(1:ifull) = epsp - 1.
      elsewhere
        epst(1:ifull) = 0.
      end where
    end if ! (ktau>2.and.epsp>1..and.epsp<2.)

    if ( ktau<10 .and. mydiag ) then
      write(6,*)'savu,u,ubar ',ktau,savu(idjd,1),u(idjd,1),ubar(idjd,1)
    end if
    if ( ktau==1 .and. ((.not.lrestart).or.always_mspeca) .and. mspec==1 .and. mex/=1 ) then
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

    ! update non-linear dynamic terms
    call nonlin
      
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
    call upglobal
      
    if ( diag ) then
      if ( mydiag ) then
        write(6,*) 'after hadv'
        write (6,"('tx  ',9f8.2)") (tx(idjd,k),k=nlx,nlx+8)
      end if
      call printa('tx  ',tx,ktau,nlv,ia,ib,ja,jb,200.,1.)
      if ( mydiag ) then
        write(6,'(i2," qgh ",18f7.4)')ktau,1000.*qg(idjd,:)
      end if  
    end if

    if ( nstaguin<0 .and. ktau>=1 ) then  ! swapping here (lower down) for nstaguin<0
      if ( nstagin<0 .and. mod(ktau-nstagoff,abs(nstagin))==0 ) then
        nstag  = 7 - nstag  ! swap between 3 & 4
        nstagu = nstag
      end if
    end if
      
    ! Update the semi-implicit solution to the augumented geopotential
    call adjust5

    ! check for rounding errors and invalid data
    call fixqg(1,ifull)
    call nantest("after atmosphere dynamics",1,ifull,"dynamics")
      
      
    ! NESTING ---------------------------------------------------------------
    ! nesting now after mass fixers
    if ( mspec==1 ) then
      if ( mbd/=0 .or. (mbd_mlo/=0.and.namip==0) ) then
        ! scale-selective filter
        call START_LOG(nestin_begin)
        call nestinb
        call END_LOG(nestin_end)
        call nantest("after nesting",1,ifull,"nesting")        
      else if ( nbd/=0 ) then
        ! Newtonian relaxiation
        call START_LOG(nestin_begin)
        call davies
        call END_LOG(nestin_end)
        call nantest("after nesting",1,ifull,"nesting")        
      end if
    end if
    
      
    ! DYNAMICS --------------------------------------------------------------
    if ( mspec==2 ) then     ! for very first step restore mass & T fields
      call gettin(1)
    endif    !  (mspec==2) 
    dt = dtin
      
  end do ! ****** end of introductory time loop

  mspeca = 1
  
    
  ! HORIZONTAL DIFFUSION ----------------------------------------------------
  if ( nhor<0 ) then
    call START_LOG(hordifg_begin) 
    call hordifgt
    if ( diag .and. mydiag ) then
      write(6,*) 'after hordifgt t ',t(idjd,:)
    end if
    call END_LOG(hordifg_end)
    call nantest("after atm horizontal diffusion",1,ifull,"dynamics")        
  end if  


  ! calculate MSE after dynamic (advection)
  call calculate_dhdt_mse(1,ifull,mse_t2)
  dmsedt_adv = (mse_t2-mse_t1)/dt  
  

  ! ***********************************************************************
  ! RIVER ROUTING AND HYDROLOGY
  ! ***********************************************************************

  if ( nhstest>=0 ) then
    call START_LOG(river_begin)
    call rvrrouter
    call water_table_transport
    call END_LOG(river_end)
  end if  


    
  ! ***********************************************************************
  ! OCEAN DYNAMICS
  ! ***********************************************************************

  ! nmlo=0   Prescriped SSTs and sea-ice with JLM skin enhancement
  ! nmlo=1   1D mixed-layer-ocean model
  ! nmlo=2   nmlo=1 plus river-routing and horiontal diffusion
  ! nmlo=3   nmlo=2 plus 3D dynamics
  ! nmlo>9   Use external VCOM ocean model

  if ( mlomode(nmlo)=="ocn_dynamics" ) then
    ! DYNAMICS & DIFFUSION ------------------------------------------------
    call START_LOG(waterdynamics_begin)
    call mlohadv
    call END_LOG(waterdynamics_end)
  else if ( mlomode(nmlo)=="ocn_diffusion" ) then
    ! DIFFUSION -----------------------------------------------------------
    call START_LOG(waterdynamics_begin)
    call mlodiffusion
    call END_LOG(waterdynamics_end)
  end if
  if ( mlomode(nmlo)=="ocn_dynamics" .or. mlomode(nmlo)=="ocn_diffusion" ) then
    ! SEA ICE  ------------------------------------------------------------
    call START_LOG(waterdynamics_begin)
    call mlonewice
    call mloexpice("fracice",fracice,0)
    call mloexpice("thick",sicedep,0)
    call mloexpice("snowd",snowd,0)
    call mlosurf("sst",tss,0)
    call END_LOG(waterdynamics_end)
#ifdef csircoupled
  else
    ! ***********************************************************************
    ! VCOM ADVECTION
    ! ***********************************************************************
    call vcom_ccam_advect(fracice,sicedep,tss,tgg(:,1),tggsn(:,1))
    ! ***********************************************************************
    ! VCOM DIFFUSION
    ! ***********************************************************************
    call vcom_ccam_diffusion(fracice,sicedep,tss,tgg(:,1),tggsn(:,1))
#endif
  end if
    
      

  ! ***********************************************************************
  ! PHYSICS 
  ! ***********************************************************************
  call START_LOG(phys_begin)


  ! MISC ------------------------------------------------------------------

  ! radiation timer calculations
  if ( nrad==5 ) then
    if ( nhstest<0 ) then      ! aquaplanet test -1 to -8  
      mtimer_sav = mtimer
      mtimer     = mins_gmt    ! so radn scheme repeatedly works thru same day
      call seaesfrad_settime
      mtimer = mtimer_sav
    else
      call seaesfrad_settime 
    end if    ! (nhstest<0)      
  end if    
  ! aerosol timer (true indicates update oxidants, etc)
  oxidant_update = oxidant_timer<=mins-updateoxidant
  if ( oxidant_update ) then
    if ( nsib>0 ) then
      call START_LOG(aerosol_begin)
      if ( abs(iaero)>=2 ) then
        call aerocalc_init(mins)
      end if
      call END_LOG(aerosol_end)
    end if
  end if  
  ! interpolate tracer fluxes to current timestep
  if ( ngas>0 ) then
    call interp_tracerflux
  end if  

  ! initialse surface rainfall to zero (also initialised in convection)
  condc(1:ifull) = 0. ! default convective rainfall (assumed to be rain)
  condx(1:ifull) = 0. ! default total precip = rain + ice + snow + graupel (convection and large scale)
  conds(1:ifull) = 0. ! default total ice + snow (convection and large scale)
  condg(1:ifull) = 0. ! default total graupel (convection and large scale)
  ! Held & Suarez or no surf fluxes
  if ( ntsur<=1 .or. nhstest==2 ) then 
    eg(1:ifull)   = 0.
    fg(1:ifull)   = 0.
    cdtq(1:ifull) = 0.
    cduv(1:ifull) = 0.
  end if     ! (ntsur<=1.or.nhstest==2) 
  ! Save aerosol concentrations for outside convective fraction of grid box
  if ( abs(iaero)>=2 ) then
    xtosav(1:ifull,1:kl,1:naero) = xtg(1:ifull,1:kl,1:naero) ! Aerosol mixing ratio outside convective cloud
  end if
  call nantest("start of physics",1,ifull,"all")

  
  ! GWDRAG ----------------------------------------------------------------
  if ( nsib>0 ) then
    call START_LOG(gwdrag_begin)
    if ( ngwd<0 ) then
      call gwdrag  ! <0 for split - only one now allowed
    end if
    call END_LOG(gwdrag_end)
  end if
  call nantest("after gravity wave drag",1,ifull,"gwdrag")


  ! CONVECTION ------------------------------------------------------------
  do k = 1,kl  
    do iq = 1,ifull
      convh_ave(iq,k) = convh_ave(iq,k) - t(iq,k)*real(nperday)/real(nperavg)
    end do  
  end do
  if ( nsib>0 ) then
    call START_LOG(convection_begin)
    call ctrl_convection 
    call END_LOG(convection_end)
  end if
  call fixqg(1,ifull)
  call nantest("after convection",1,ifull,"conv")


  ! CLOUD MICROPHYSICS ----------------------------------------------------
  if ( nsib>0 ) then
    call START_LOG(cloud_begin)
    call ctrl_microphysics
    call END_LOG(cloud_end)
  end if
  do k = 1,kl  
    do iq = 1,ifull
      convh_ave(iq,k) = convh_ave(iq,k) + t(iq,k)*real(nperday)/real(nperavg)
    end do  
  end do
  call nantest("after cloud microphysics",1,ifull,"cloud") 
  
  
  ! RADIATION -------------------------------------------------------------
  ! calculate MSE before radiation
  call calculate_dhdt_mse(1,ifull,mse_t1)  
  if ( nsib>0 ) then
    call START_LOG(radnet_begin)
    rad_tend(1:ifull,1:kl) = rad_tend(1:ifull,1:kl) - t(1:ifull,1:kl)/dt
    select case ( nrad )
      case(4)
        ! Fels-Schwarzkopf radiation
        if ( nhstest<0 ) then    ! aquaplanet test -1 to -8  
          mtimer_sav = mtimer
          mtimer     = mins_gmt  ! so radn scheme repeatedly works thru same day
          call radrive(il*nrows_rad)
          mtimer = mtimer_sav
        else
          call radrive(il*nrows_rad)  
        end if    ! (nhstest<0)
      case(5)
        ! GFDL SEA-EFS radiation
        call seaesfrad(koundiag)
      case DEFAULT
        ! use preset slwa array (use +ve nrad)
        slwa(1:ifull) = -real(10*nrad)
    end select
    call END_LOG(radnet_end)
  end if ! ( nsib>0 )  
  do k = 1,kl
    t(1:ifull,k) = t(1:ifull,k) - dt*(sw_tend(1:ifull,k)+lw_tend(1:ifull,k))
    rad_tend(1:ifull,k) = rad_tend(1:ifull,k) + t(1:ifull,k)/dt
  end do
  ! calculate MSE after radiation
  call calculate_dhdt_mse(1,ifull,mse_t2)
  dmsedt_rad(:,:) = (mse_t2 (:,:)- mse_t1(:,:))/dt  
  call nantest("after radiation",1,ifull,"radiation")
    
  
  ! HELD & SUAREZ ---------------------------------------------------------
  if ( nhstest==2 ) then
    call hs_phys
  end if
  
  
  ! SURFACE FLUXES ---------------------------------------------
  ! Calculate sensible heat flux, latent heat flux, carbon cycle,
  ! sea-ice thermodynamocs and urban physics.  Ocean tubulent
  ! vertical mixing can be combined with vertical atmosphere
  ! turbulent mixing below.
  if ( nsib>0 ) then
    call START_LOG(sfluxnet_begin)
    if ( diag .and. ntiles==1 ) then
      call maxmin(u,'#u',ktau,1.,kl)
      call maxmin(v,'#v',ktau,1.,kl)
      call maxmin(t,'#t',ktau,1.,kl)
      call maxmin(qg,'qg',ktau,1.e3,kl)     
    end if
    if ( ntsur>1 ) then
      call sflux
    endif   ! (ntsur>1)    
    call END_LOG(sfluxnet_end)
  end if
  call nantest("after surface fluxes",1,ifull,"surface")

  
  ! AEROSOLS --------------------------------------------------------------
  ! Old time-split with aero_split=0.
  ! To be depreciated.
  if ( nsib>0 .and. aero_split==0 ) then
    call START_LOG(aerosol_begin)
    if ( abs(iaero)>=2 ) then
      call aerocalc(mins,0)
    end if
    call END_LOG(aerosol_end)
    call nantest("after aerosols",1,ifull,"aerosols (old)")
  end if  
  

  ! VERTICAL MIXING ------------------------------------------------------
  ! Turbulent mixing of atmosphere and optionally combined with ocean
  ! calculate MSE before vertical mixing
  call calculate_dhdt_mse(1,ifull,mse_t1)  
  if ( nsib>0 ) then
    call START_LOG(vertmix_begin)
    if ( nmaxpr==1 ) then
      if ( mydiag .and. ntiles==1 ) then
        write (6,"('pre-vertmix t',9f8.3/13x,9f8.3)") t(idjd,:)
      end if
    end if
    trb_tend(1:ifull,1:kl) = trb_tend(1:ifull,1:kl) - t(1:ifull,1:kl)/dt
    trb_qend(1:ifull,1:kl) = trb_qend(1:ifull,1:kl) - qg(1:ifull,1:kl)/dt - qlg(1:ifull,1:kl)/dt - qfg(1:ifull,1:kl)/dt
    if ( ntsur>=1 ) then
      call turbmix
    end if  ! (ntsur>=1)
    trb_tend(1:ifull,1:kl) = trb_tend(1:ifull,1:kl) + t(1:ifull,1:kl)/dt
    trb_qend(1:ifull,1:kl) = trb_qend(1:ifull,1:kl) + qg(1:ifull,1:kl)/dt + qlg(1:ifull,1:kl)/dt + qfg(1:ifull,1:kl)/dt
    if ( nmaxpr==1 ) then
      if ( mydiag .and. ntiles==1 ) then
        write (6,"('aft-vertmix t',9f8.3/13x,9f8.3)") t(idjd,:)
      end if
    end if
    call END_LOG(vertmix_end)
  end if
  call fixqg(1,ifull)
  ! calculate MSE after vertical mixing
  call calculate_dhdt_mse(1,ifull,mse_t2)
  dmsedt_pbl(:,:) = (mse_t2(:,:) - mse_t1(:,:))/dt  
  call nantest("after PBL mixing",1,ifull,"vmixing")

  
  ! AEROSOLS --------------------------------------------------------------
  ! New time-split with aero_split=1
  ! Emission driven prognostic aerosols.  Includes aerosol turbulent mixing.
  if ( nsib>0 ) then
    call START_LOG(aerosol_begin)
    if ( abs(iaero)>=2 ) then
      call aerocalc(mins,1)
    end if
    call END_LOG(aerosol_end)
    call nantest("after aerosols",1,ifull,"aerosols")
  end if

  
  ! TRACERS ---------------------------------------------------------------
  ! Turbulent mixing
  if ( nsib>0 ) then
    if ( ngas>0 ) then
      call tracervmix  
    end if
  end if

  
  ! MISC ------------------------------------------------------------------
  ! Update diagnostics for consistancy in history file
  call fixsat(1,ifull) ! if qg_fix>1, then removes supersaturated qg
  call nantest("after fixsat",1,ifull,"cloud")
  ! Convection diagnostic output
  cbas_ave(1:ifull) = cbas_ave(1:ifull) + condc(1:ifull)*(1.1-sig(kbsav(1:ifull)))      ! diagnostic
  ctop_ave(1:ifull) = ctop_ave(1:ifull) + condc(1:ifull)*(1.1-sig(abs(ktsav(1:ifull)))) ! diagnostic
  ! Microphysics diagnostic output
  rnd_3hr(1:ifull,8) = rnd_3hr(1:ifull,8) + real(condx(1:ifull),8)  ! i.e. rnd24(:)=rnd24(:)+condx(:)
  if ( rescrn>0 ) then
    call autoscrn(1,ifull)
  end if

  call END_LOG(phys_end)

  
  ! Diagnose CAPE, CIN and LI for cordex output  
  ! pcc2hist can calculate CAPE for standard output
  ltest = .false.
  if ( surfile/=' ' ) then
    ltest = ltest .or. mod(ktau,tbave)==0
  end if
  if ( freqfile/=' ' ) then
    ltest = ltest .or. mod(ktau,tbave10)==0  
  end if  
  if ( ltest ) then
    call START_LOG(cape_begin)  
    call capecalc
    call END_LOG(cape_end)
  end if  
  ! Calculate updraft helicity
  call uh_calc
  ! Update aerosol timer
  if ( oxidant_update ) then
    oxidant_timer = mins
  end if
  
  
  ! ***********************************************************************
  ! DIAGNOSTICS AND OUTPUT
  ! ***********************************************************************

  ! TIME AVERAGED OUTPUT ---------------------------------------
  ! update diag_averages and daily max and min screen temps 
  ! N.B. runoff is accumulated in sflux
  if ( ktau==1 ) then
    tmaxscr(:)  = tscrn(:) 
    tminscr(:)  = tscrn(:) 
    rhmaxscr(:) = rhscrn(:) 
    rhminscr(:) = rhscrn(:) 
  end if    
  call calculate_timeaverage(koundiag)


  ! TRACER OUTPUT ----------------------------------------------
  if ( ngas>0 ) then
    call tracer_mass !also updates average tracer array
    call write_ts(ktau,ntau,dt)
  endif

    
  ! STATION OUTPUT ---------------------------------------------
  if ( nstn>0 ) then
    call stationa ! write every time step
  end if
    
    
  ! DIAGNOSTICS ------------------------------------------------
  call write_diagnostics(mins_gmt,nmaxprsav)

    
  if ( myid==0 ) then
    write(6,*) 'ktau,mod,nper3hr ',ktau,mod(ktau-1,nperday)+1,nper3hr(n3hr)
  end if

  ! rnd03 to rnd21 are accumulated in mm     
  if ( mod(ktau-1,nperday)+1 == nper3hr(n3hr) ) then
    rnd_3hr(1:ifull,n3hr) = rnd_3hr(1:ifull,8)
    if ( nextout>=2 ) then
      spare1(:) = max( .001, sqrt(u(1:ifull,1)**2+v(1:ifull,1)**2) )
      u10_3hr(:,n3hr) = u10(:)*u(1:ifull,1)/spare1(:)
      v10_3hr(:,n3hr) = u10(:)*v(1:ifull,1)/spare1(:)
      tscr_3hr(:,n3hr) = tscrn(:)
      spare1(:) = establ(t(1:ifull,1)) ! spare1 = es
      rh1_3hr(1:ifull,n3hr) = 100.*qg(1:ifull,1)*(ps(1:ifull)*sig(1)-spare1(:))/(.622*spare1(:))
    end if    ! (nextout==2)
    n3hr = n3hr + 1
    if ( n3hr>8 ) n3hr = 1
  endif    ! (mod(ktau,nperday)==nper3hr(n3hr))
  
  ! Turn off log before writing output
  call log_off

  ! WRITE DATA TO HISTORY ---------------------------------
  if ( ktau==ntau .or. mod(ktau,nwt)==0 ) then
    call outfile(20,ofile,psl,u,v,t,qg)  ! which calls outcdf
  end if
  ! write high temporal frequency fields
  if ( surfile/=' ' ) then
    call freqfile_cordex
  end if
  if ( freqfile/=' ' ) then
    call freqfile_10
  end if

    
  ! ENSEMBLE --------------------------------------------------
  if ( ensemble_mode>0 ) then
    call START_LOG(ensemble_begin)
    call update_ensemble
    call END_LOG(ensemble_end)
    call fixqg(1,ifull)
  end if
  
  
  ! Turn on log after writing output
  call log_on
 
    
  ! TIME AVERAGED DIAGNOSTICS ---------------------------------
  if ( mod(ktau,nperavg)==0 ) then    
    ! produce some diags & reset most averages once every nperavg
    if ( nmaxpr==1 ) then
      precavge = sum(real(precip(1:ifull))*wts(1:ifull))
      evapavge = sum(real(evspsbl_ave(1:ifull))*wts(1:ifull))   ! in mm/day
      pwatr    = 0.   ! in mm
      do k = 1,kl
        pwatr = pwatr - sum(dsig(k)*wts(1:ifull)*(qg(1:ifull,k)+qlg(1:ifull,k)+qfg(1:ifull,k))*ps(1:ifull))/grav
      enddo
      temparray(1:3) = (/ precavge, evapavge, pwatr /)
      call ccmpi_reduce(temparray(1:3),gtemparray(1:3),"max",0,comm_world)
      if ( myid==0 ) then
        precavge = gtemparray(1)
        evapavge = gtemparray(2)
        pwatr    = gtemparray(3)
        write(6,985) pwatr,precavge,evapavge ! MJT bug fix
985     format(' average pwatr,precc,prec,evap: ',4f7.3)
      end if
    end if
    ! also zero most averaged fields every nperavg
    call zero_nperavg(koundiag)
  endif  ! (mod(ktau,nperavg)==0)
  
  ! HOURLY DIAGNOSTICS ---------------------------------------
  if ( mod(ktau,nperhr)==0 ) then
    call zero_nperhour
  end if

  ! DAILY DIAGNOSTICS ----------------------------------------
  if ( mod(ktau,nperday)==0 ) then   ! re-set at the end of each 24 hours
    call zero_nperday
  endif   ! (mod(ktau,nperday)==0)
  
    
  ! AMIP SSTs --------------------------------------------------
  if ( namip/=0 ) then
    call START_LOG(amipsst_begin)
    call amipsst
    call END_LOG(amipsst_end)
  end if

    
  ! Flush trace information to disk to save memory.
  call log_flush()
    
end do                  ! *** end of main time loop
  
call END_LOG(maincalc_end)

! WRITE DATA TO RESTART FILES ---------------------------
ktau = ntau
if ( irest==1 ) then
  ! write restart file
  call outfile(19,restfile,psl,u,v,t,qg)
  if ( myid==0 ) write(6,*) 'finished writing restart file in outfile'
endif

call log_off
  
  
!------------------------------------------------------------------
! SIMULATION COMPLETE
  
! Report timings of run
if ( myid==0 ) then
  call date_and_time(time=timeval,values=tvals2)
  write(6,*) "End of time loop ", timeval
  write(6,*) "Normal termination of run"
  write(6,*) "End time ", timeval
  call time_diff(tt_r8,tvals1,tvals2)
  aa = real( tt_r8 )
  write(6,*) "Model time in main loop",aa
end if
  
! close mesonest files
if ( mbd/=0 .or. nbd/=0 .or. (mbd_mlo/=0.and.namip==0) .or. ensemble_mode>0 ) then
  call histclose
end if

#ifdef csircoupled
! finalize VCOM
call vcom_finialize
#endif
  
call date_and_time(values=times_total_b)
call time_diff(total_time,times_total_a,times_total_b)
  
! report subroutine timings
call simple_timer_finalize

! Complete
if ( myid==0 ) then
  write(6,*) "------------------------------------------------------------------------------"
  write(6,*) "CCAM: globpea completed successfully"
  call finishbanner
end if

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
call nestin_exit
if ( mbd/=0 .and. nud_uv/=9 ) then
  call deallocateglobalpack
end if
nullify(xx4, yy4)
nullify(em_g)
nullify(x_g, y_g, z_g)

  
!****************************************************************

! finalize MPI comms
call ccmpi_finalize

end program
    
    
!--------------------------------------------------------------
! END OF CCAM LOG    
subroutine finishbanner

implicit none

! End banner
write(6,*) "=============================================================================="
write(6,*) "CCAM: Finished globpea"
write(6,*) "=============================================================================="

return
end subroutine finishbanner

!--------------------------------------------------------------
! SIMULATION TIME
! This version should be valid for at least 28 days
subroutine time_diff(t_diff,tcal1,tcal2)

use dates_m            ! Date data

implicit none

integer month, year, day, kdate_l
integer, dimension(12) :: mdays
integer, dimension(8), intent(in) :: tcal1, tcal2
real(kind=8), intent(out) :: t_diff
real(kind=8) day_diff

t_diff = sum( real(tcal2(5:8)-tcal1(5:8),8)*(/ 3600._8, 60._8, 1._8, 0.001_8 /) )

year = int( tcal1(1) )
month = int( tcal1(2) )
day = int( tcal1(3) )
kdate_l = year*10000 + month*100 + day
call calendar_function(mdays,kdate_l,cal_leap)

day_diff = real(tcal2(3) - tcal1(3),8)
if ( day_diff < 0._8 ) then
  day_diff = day_diff + real(mdays(month),8)
end if
t_diff = t_diff + day_diff*86400._8

return
end subroutine time_diff
    
!--------------------------------------------------------------
! PREPARE SPECIAL TRACER ARRAYS
! sets tr arrays for lat, long, pressure if nextout>=4 &(nllp>=3)
subroutine setllp
      
use arrays_m           ! Atmosphere dyamics prognostic arrays
use cc_mpi             ! CC MPI routines
use const_phys         ! Physical constants
use latlong_m          ! Lat/lon coordinates
use newmpar_m          ! Grid parameters
use sigs_m             ! Atmosphere sigma levels
use tracers_m          ! Tracer data
      
implicit none
      
integer k
      
if ( nllp<3 ) then
  write(6,*) "ERROR: Incorrect setting of nllp",nllp
  call ccmpi_abort(-1)
end if

if ( myid==0 ) then
  write(6,*) "Calling setllp with nllp = ",nllp
end if
      
do k = 1,kl
  tr(1:ifull,k,ngas+1) = rlatt(1:ifull)*180./pi
  tr(1:ifull,k,ngas+2) = rlongg(1:ifull)*180./pi
  tr(1:ifull,k,ngas+3) = .01*ps(1:ifull)*sig(k)  ! in HPa
enddo
if ( nllp >= 4 ) then   ! theta
  do k = 1,kl
    tr(1:ifull,k,ngas+4) = t(1:ifull,k)*(1.e-5*ps(1:ifull)*sig(k))**(-rdry/cp)
  enddo
endif   ! (nllp>=4)
if ( nllp >= 5 ) then   ! mixing_ratio (g/kg)
  do k = 1,kl
    tr(1:ifull,k,ngas+5) = 1000.*qg(1:ifull,k)
  enddo
endif   ! (nllp>=5)
      
return
end subroutine setllp

    
!--------------------------------------------------------------
! WRITE STATION DATA
subroutine stationa

use arrays_m           ! Atmosphere dyamics prognostic arrays
use cc_mpi             ! CC MPI routines
use const_phys         ! Physical constants
use dates_m            ! Date data
use diag_m             ! Diagnostic routines
use estab              ! Liquid saturation function
use extraout_m         ! Additional diagnostics
use indata             ! Data initialisation
use map_m              ! Grid map arrays
use morepbl_m          ! Additional boundary layer diagnostics
use newmpar_m          ! Grid parameters
use nsibd_m            ! Land-surface arrays
use parm_m             ! Model configuration
use parmgeom_m         ! Coordinate data
use pbl_m              ! Boundary layer arrays
use prec_m             ! Precipitation
use screen_m           ! Screen level diagnostics
use sigs_m             ! Atmosphere sigma levels
use soil_m             ! Soil and surface data
use soilsnow_m         ! Soil, snow and surface data
use soilv_m            ! Soil parameters
use tracers_m          ! Tracer data
use vecsuv_m           ! Map to cartesian coordinates
use vegpar_m           ! Vegetation arrays
use work2_m            ! Diagnostic arrays
use work3_m            ! Mk3 land-surface diagnostic arrays
use xyzinfo_m          ! Grid coordinate arrays

implicit none

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
  iqt = min( iq, il*jl ) ! Avoid bounds problems if there are no tracers
  k2  = min( 2, kl )
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
    
!--------------------------------------------------------------------
! Fix water vapour mixing ratio.  Needed for JLM convection.
subroutine fixqg(js,je)

use arrays_m                          ! Atmosphere dyamics prognostic arrays
use const_phys                        ! Physical constants
use cfrac_m                           ! Cloud fraction
use estab                             ! Liquid saturation function
use liqwpar_m                         ! Cloud water mixing ratios
use newmpar_m                         ! Grid parameters
use parm_m                            ! Model configuration
use sigs_m                            ! Atmosphere sigma levels

implicit none

integer, intent(in) :: js, je
integer k, iq
real qtot, tliq

! requires qg_fix>=1
if ( qg_fix>=1 ) then

  do k = 1,kl
    do iq = js,je  
      qtot = qg(iq,k) + qlg(iq,k) + qfg(iq,k)
      tliq = t(iq,k) - hlcp*qlg(iq,k) - hlscp*qfg(iq,k)
  
      qfg(iq,k)   = max( qfg(iq,k), 0. ) 
      qlg(iq,k)   = max( qlg(iq,k), 0. )
      qrg(iq,k)   = max( qrg(iq,k), 0. )
      qsng(iq,k)  = max( qsng(iq,k), 0. )
      qgrg(iq,k)  = max( qgrg(iq,k), 0. )
  
      qg(iq,k) = max( qtot - qlg(iq,k) - qfg(iq,k), 0. )
      t(iq,k)  = tliq + hlcp*qlg(iq,k) + hlscp*qfg(iq,k)
    end do  
  end do

end if ! qg_fix>=1  
  
return
end subroutine fixqg    
    
! To be depreciated - Limits saturation after turbulent mixing    
subroutine fixsat(js,je)

use arrays_m                          ! Atmosphere dyamics prognostic arrays
use const_phys                        ! Physical constants
use cfrac_m                           ! Cloud fraction
use estab                             ! Liquid saturation function
use liqwpar_m                         ! Cloud water mixing ratios
use morepbl_m                         ! Additional boundary layer diagnostics
use newmpar_m                         ! Grid parameters
use parm_m                            ! Model configuration
use sigs_m                            ! Atmosphere sigma levels

implicit none

integer, intent(in) :: js, je
integer k, iq
real, dimension(js:je,kl) :: zg
real, dimension(js:je) :: tliq, pk, qsi, deles, qsl
real qtot, qsw, fice
real dqsdt, hlrvap, al, qc
real, parameter :: tice = 233.16

! requires qg_fix>=2
if ( qg_fix<=1 ) return

! calculate approximate height above surface
zg(js:je,1) = bet(1)*t(js:je,1)/grav
do k = 2,kl
  zg(js:je,k) = zg(js:je,k-1) + (bet(k)*t(js:je,k) + betm(k)*t(js:je,k-1))/grav
end do

do k = 1,kl
    
  tliq(js:je) = t(js:je,k) - hlcp*qlg(js:je,k) - hlscp*qfg(js:je,k)
  pk(js:je) = ps(js:je)*sig(k)  
  qsi(js:je) = qsati(pk(js:je),tliq(js:je))  
  deles(js:je) = esdiffx(tliq(js:je))
  qsl(js:je) = qsi(js:je) + epsil*deles(js:je)/pk(js:je)
  
  do iq = js,je
    ! only apply below boundary layer height  
    if ( zg(iq,k) < pblh(iq) ) then
    
      qtot = max( qg(iq,k) + qlg(iq,k) + qfg(iq,k), 0. )
      qc   = max( qlg(iq,k) + qfg(iq,k), 0. )
      if ( qfg(iq,k)>1.e-8 ) then
        fice = min( qfg(iq,k)/(qfg(iq,k)+qlg(iq,k)), 1. )
      else
        fice = 0.
      end if
      qsw = fice*qsi(iq) + (1.-fice)*qsl(iq)
      hlrvap = (hl+fice*hlf)/rvap
      dqsdt = qsw*hlrvap/tliq(iq)**2
      al = 1./(1.+(hlcp+fice*hlfcp)*dqsdt)
      qc = max( al*(qtot - qsw), qc )
      if ( t(iq,k)>=tice ) then
        qfg(iq,k) = max( fice*qc, 0. )  
        qlg(iq,k) = max( qc-qfg(iq,k), 0. )
      end if

      qg(iq,k) = max( qtot - qlg(iq,k) - qfg(iq,k), 0. )
      t(iq,k)  = tliq(iq) + hlcp*qlg(iq,k) + hlscp*qfg(iq,k)
      
    end if ! zg<pblh  
  end do   ! iq 
  
end do     ! k

return
end subroutine fixsat    

!--------------------------------------------------------------
! Reset diagnostics for averaging period    
subroutine zero_nperavg(koundiag)

use aerosol_arrays, only :               & ! Aerosol arrays
     duste,dustwd,dustdd,dust_burden     &
    ,bce,bcwd,bcdd,bc_burden             &
    ,oce,ocwd,ocdd,oc_burden             &
    ,dmse,dms_burden                     &
    ,so2e,so2wd,so2dd,so2_burden         &
    ,so4e,so4wd,so4dd,so4_burden         &
    ,dmsso2o,so2so4o, salte,saltdd       &
    ,saltwd,salt_burden
use extraout_m                             ! Additional diagnostics
use histave_m                              ! Time average arrays
use liqwpar_m                              ! Cloud water mixing ratios
use morepbl_m                              ! Additional boundary layer diagnostics
use parm_m                                 ! Model configuration
use prec_m                                 ! Precipitation
use raddiag_m                              ! Radiation diagnostic
use sflux_m                                ! Surface flux routines
use soilsnow_m                             ! Soil, snow and surface data
use tracers_m                              ! Tracer data
use vvel_m                                 ! Additional vertical velocity

implicit none

integer, intent(inout) :: koundiag

convh_ave(:,:)       = 0.
cbas_ave(:)          = 0.
ctop_ave(:)          = 0.
dew_ave(:)           = 0._8
epan_ave(:)          = 0._8
epot_ave(:)          = 0._8
eg_ave(:)            = 0._8
fg_ave(:)            = 0._8
ga_ave(:)            = 0._8
anthropogenic_ave(:) = 0._8
anth_elecgas_ave(:)  = 0._8
anth_heating_ave(:)  = 0._8
anth_cooling_ave(:)  = 0._8
rnet_ave(:)          = 0._8
rhscr_ave(:)         = 0.
tscr_ave(:)          = 0.
wb_ave(:,:)          = 0.
wbice_ave(:,:)       = 0.
taux_ave(:)          = 0._8
tauy_ave(:)          = 0._8

! radiation
koundiag             = 0
sint_ave(:)          = 0._8
sot_ave(:)           = 0._8
soc_ave(:)           = 0._8
sgdn_ave(:)          = 0._8
sgdndir_ave(:)       = 0._8
sgn_ave(:)           = 0._8
rtu_ave(:)           = 0._8
rtc_ave(:)           = 0._8
rgdn_ave(:)          = 0._8
rgn_ave(:)           = 0._8
rgc_ave(:)           = 0._8
rgdc_ave(:)          = 0._8
sgc_ave(:)           = 0._8
sgdc_ave(:)          = 0._8
cld_ave(:)           = 0._8
cll_ave(:)           = 0._8
clm_ave(:)           = 0._8
clh_ave(:)           = 0._8
dni_ave(:)           = 0._8

! zero evap, precip, precc, sno, runoff fields each nperavg (3/12/04) 
evspsbl_ave(:)       = 0._8  ! converted to mm/day in outcdf
sbl_ave(:)           = 0._8  ! converted to mm/day in outcdf
precip(:)            = 0._8  ! converted to mm/day in outcdf
precc(:)             = 0._8  ! converted to mm/day in outcdf
sno(:)               = 0._8  ! converted to mm/day in outcdf
grpl(:)              = 0._8  ! converted to mm/day in outcdf
runoff_ave(:)        = 0._8  ! converted to mm/day in outcdf
runoff_surface_ave(:) = 0._8 ! converted to mm/day in outcdf
snowmelt_ave(:)      = 0._8  ! converted to mm/day in outcdf
cape_max(:)          = 0.
cape_ave(:)          = 0.
dhail1(:)            = 0.
dhail2(:)            = 0.
dhail3(:)            = 0.
dhail4(:)            = 0.
dhail5(:)            = 0.
hailrad_ave(:)       = 0.
hailrad_max(:)       = 0.
updraft_helicity_max(:) = -9.e9
updraft_helicity_min(:) = 9.e9

if ( ngas>0 ) then
  traver = 0.
end if

if ( ccycle/=0 ) then
  fnee_ave = 0.  
  fpn_ave  = 0.
  frd_ave  = 0.
  frp_ave  = 0.
  frpw_ave = 0.
  frpr_ave = 0.
  frs_ave  = 0.
  cnpp_ave = 0.
  cnbp_ave = 0.
  if ( diaglevel_carbon > 0 ) then
    fevc_ave = 0.
    plant_turnover_ave = 0.
    plant_turnover_wood_ave = 0.
  end if
end if

if ( abs(iaero)>=2 ) then
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
  salte         = 0.  ! Salt emissions
  saltdd        = 0.  ! Salt dry deposition
  saltwd        = 0.  ! Salt wet deposition
  salt_burden   = 0.  ! Salt burden
end if

if ( output_windmax/=0 ) then
  u_max = 0.
  v_max = 0.
  u10m_max = 0.
  v10m_max = 0.
end if

return
end subroutine zero_nperavg
    
!--------------------------------------------------------------
! Reset diagnostics for averaging period    
subroutine zero_nperhour

use histave_m                              ! Time average arrays

implicit none

prhour(:) = 0._8 ! for calculating prhmax

return
end subroutine zero_nperhour
    
    
!--------------------------------------------------------------
! Reset diagnostics for daily averages    
subroutine zero_nperday

use histave_m                              ! Time average arrays
use morepbl_m                              ! Additional boundary layer diagnostics
use parm_m                                 ! Model configuration
use prec_m                                 ! Precipitation
use raddiag_m                              ! Radiation diagnostic
use screen_m                               ! Screen level diagnostics

implicit none

rndmax(:)   = 0.
prhmax(:)   = 0.
tmaxscr(:)  = 100.
tminscr(:)  = 400.
rhmaxscr(:) = 0.
rhminscr(:) = 1000.
u10max(:)   = 0.
v10max(:)   = 0.
u1max(:)    = 0.
v1max(:)    = 0.
u2max(:)    = 0.
v2max(:)    = 0.
rnd_3hr(:,8)= 0._8       ! i.e. rnd24(:)=0.
tmaxurban(:)= 100.
tminurban(:)= 400.
wsgsmax(:)  = 0.
sunhours(:) = 0._8

if ( nextout>=4 ) then
  call setllp ! reset once per day
end if

return
end subroutine zero_nperday
    
!--------------------------------------------------------------
! Update diagnostics for averaging period    
subroutine calculate_timeaverage(koundiag)

use aerosol_arrays, only :               & ! Aerosol arrays
     duste,dustwd,dustdd,dust_burden     &
    ,bce,bcwd,bcdd,bc_burden             &
    ,oce,ocwd,ocdd,oc_burden             &
    ,dmse,dms_burden                     &
    ,so2e,so2wd,so2dd,so2_burden         &
    ,so4e,so4wd,so4dd,so4_burden         &
    ,dmsso2o,so2so4o,salte,saltdd        &
    ,saltwd,salt_burden
use arrays_m                               ! Atmosphere dyamics prognostic arrays
use const_phys                             ! Physical constants
use extraout_m                             ! Additional diagnostics
use histave_m                              ! Time average arrays
use liqwpar_m                              ! Cloud water mixing ratios
use mlo_ctrl                               ! Ocean physics control layer
use morepbl_m                              ! Additional boundary layer diagnostics
use newmpar_m                              ! Grid parameters
use nharrs_m                               ! Non-hydrostatic atmosphere arrays
use outcdf                                 ! Output file routines
use parm_m                                 ! Model configuration
use pbl_m                                  ! Boundary layer arrays
use prec_m                                 ! Precipitation
use raddiag_m                              ! Radiation diagnostic
use screen_m                               ! Screen level diagnostics
use sflux_m                                ! Surface flux routines
use soilsnow_m                             ! Soil, snow and surface data
use tracers_m                              ! Tracer data
use vvel_m                                 ! Additional vertical velocity
use work3_m                                ! Mk3 land-surface diagnostic arrays

implicit none

integer, intent(inout) :: koundiag
integer iq, k
real, dimension(ifull) :: spare1, spare2

precip(1:ifull)            = precip(1:ifull) + real(condx,8)
!precc(1:ifull)            = precc(1:ifull) + real(condc,8) ! currently inside convection parameterisations
sno(1:ifull)               = sno(1:ifull) + real(conds,8)
grpl(1:ifull)              = grpl(1:ifull) + real(condg,8)
rndmax(1:ifull)            = max( rndmax(1:ifull), condx )
prhour(1:ifull)            = prhour(1:ifull) + real(condx/3600.,8) ! condx/dt*dt/3600 to give kg/m2/hr
prhmax(1:ifull)            = max( prhmax(1:ifull), real(prhour) )

cape_max(1:ifull)          = max( cape_max(1:ifull), cape )
cape_ave(1:ifull)          = cape_ave(1:ifull) + cape

tmaxscr(1:ifull)           = max( tmaxscr(1:ifull), tscrn )
tminscr(1:ifull)           = min( tminscr(1:ifull), tscrn )
rhmaxscr(1:ifull)          = max( rhmaxscr(1:ifull), rhscrn )
rhminscr(1:ifull)          = min( rhminscr(1:ifull), rhscrn )
dew_ave(1:ifull)           = dew_ave(1:ifull) - real( min( 0., eg ), 8 )
epan_ave(1:ifull)          = epan_ave(1:ifull) + real(epan,8)
epot_ave(1:ifull)          = epot_ave(1:ifull) + real(epot,8) 
eg_ave(1:ifull)            = eg_ave(1:ifull) + real(eg,8)
fg_ave(1:ifull)            = fg_ave(1:ifull) + real(fg,8)
ga_ave(1:ifull)            = ga_ave(1:ifull) + real(ga,8)
anthropogenic_ave(1:ifull) = anthropogenic_ave(1:ifull) + real(anthropogenic_flux,8)
anth_elecgas_ave(1:ifull)  = anth_elecgas_ave(1:ifull) + real(urban_elecgas_flux,8)
anth_heating_ave(1:ifull)  = anth_heating_ave(1:ifull) + real(urban_heating_flux,8)
anth_cooling_ave(1:ifull)  = anth_cooling_ave(1:ifull) + real(urban_cooling_flux,8)
tmaxurban(1:ifull)         = max( tmaxurban(1:ifull), urban_tas )
tminurban(1:ifull)         = min( tminurban(1:ifull), urban_tas )
rnet_ave(1:ifull)          = rnet_ave(1:ifull) + real(rnet,8)
tscr_ave(1:ifull)          = tscr_ave(1:ifull) + tscrn 
rhscr_ave(1:ifull)         = rhscr_ave(1:ifull) + rhscrn 
wb_ave(1:ifull,1:ms)       = wb_ave(1:ifull,1:ms) + wb
wbice_ave(1:ifull,1:ms)    = wbice_ave(1:ifull,1:ms) + wbice
taux_ave(1:ifull)          = taux_ave(1:ifull) + real(taux,8)
tauy_ave(1:ifull)          = tauy_ave(1:ifull) + real(tauy,8)
wsgsmax(1:ifull)           = max( wsgsmax(1:ifull), wsgs )
runoff_ave(1:ifull)        = runoff_ave(1:ifull) + real(runoff,8)
runoff_surface_ave(1:ifull) = runoff_surface_ave(1:ifull) + real(runoff_surface,8)
snowmelt_ave(1:ifull)      = snowmelt_ave(1:ifull) + real(snowmelt,8)
evspsbl_ave(1:ifull)       = evspsbl_ave(1:ifull) + real(evspsbl,8)
sbl_ave(1:ifull)           = sbl_ave(1:ifull) + real(sbl,8)

spare1(:) = u(1:ifull,1)**2 + v(1:ifull,1)**2
spare2(:) = u(1:ifull,2)**2 + v(1:ifull,2)**2
do iq = 1,ifull
  if ( u10(iq)**2 > u10max(iq)**2 + v10max(iq)**2 ) then
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

if ( ngas>0 ) then
  traver(:,:,1:ngas) = traver(:,:,1:ngas) + tr(1:ifull,:,1:ngas)
end if

if ( ccycle/=0 ) then
  fnee_ave(1:ifull) = fnee_ave(1:ifull) + fnee  
  fpn_ave(1:ifull)  = fpn_ave(1:ifull) + fpn
  frd_ave(1:ifull)  = frd_ave(1:ifull) + frd
  frp_ave(1:ifull)  = frp_ave(1:ifull) + frp
  frpw_ave(1:ifull) = frpw_ave(1:ifull) + frpw
  frpr_ave(1:ifull) = frpr_ave(1:ifull) + frpr
  frs_ave(1:ifull)  = frs_ave(1:ifull) + frs
  cnpp_ave(1:ifull) = cnpp_ave(1:ifull) + cnpp
  cnbp_ave(1:ifull) = cnbp_ave(1:ifull) + cnbp
  if ( diaglevel_carbon > 0 ) then
    fevc_ave(1:ifull) = fevc_ave(1:ifull) + fevc
    plant_turnover_ave(1:ifull)      = plant_turnover_ave(1:ifull) + plant_turnover
    plant_turnover_wood_ave(1:ifull) = plant_turnover_wood_ave(1:ifull) + plant_turnover_wood
  end if
end if

sgn_ave(1:ifull)     = sgn_ave(1:ifull)  + real(sgsave(1:ifull),8)
sgdn_ave(1:ifull)    = sgdn_ave(1:ifull) + real(sgdn(1:ifull),8)
sgdndir_ave(1:ifull) = sgdndir_ave(1:ifull) + real(sgdndir(1:ifull),8)
if ( .not.always_mspeca ) then
  sint_ave(1:ifull)  = sint_ave(1:ifull) + real(sint(1:ifull),8)
  sot_ave(1:ifull)   = sot_ave(1:ifull)  + real(sout(1:ifull),8)
  soc_ave(1:ifull)   = soc_ave(1:ifull)  + real(soutclr(1:ifull),8)
  rtu_ave(1:ifull)   = rtu_ave(1:ifull)  + real(rt(1:ifull),8)
  rtc_ave(1:ifull)   = rtc_ave(1:ifull)  + real(rtclr(1:ifull),8)
  rgn_ave(1:ifull)   = rgn_ave(1:ifull)  + real(rgn(1:ifull),8)
  rgc_ave(1:ifull)   = rgc_ave(1:ifull)  + real(rgclr(1:ifull),8)
  rgdn_ave(1:ifull)  = rgdn_ave(1:ifull) + real(rgdn(1:ifull),8)
  rgdc_ave(1:ifull)  = rgdc_ave(1:ifull) + real(rgdclr(1:ifull),8)
  sgc_ave(1:ifull)   = sgc_ave(1:ifull)  + real(sgclr(1:ifull),8)
  sgdc_ave(1:ifull)  = sgdc_ave(1:ifull) + real(sgdclr(1:ifull),8)
  cld_ave(1:ifull)   = cld_ave(1:ifull)  + real(cloudtot(1:ifull),8)
  cll_ave(1:ifull)   = cll_ave(1:ifull)  + real(cloudlo(1:ifull),8)
  clm_ave(1:ifull)   = clm_ave(1:ifull)  + real(cloudmi(1:ifull),8)
  clh_ave(1:ifull)   = clh_ave(1:ifull)  + real(cloudhi(1:ifull),8)
end if
dni_ave(1:ifull)   = dni_ave(1:ifull)  + real(dni(1:ifull),8)
where ( sgdn(1:ifull)>120. )
  sunhours(1:ifull) = sunhours(1:ifull) + real(dt/3600.,8)
end where

if ( output_windmax/=0 ) then
  do k = 1,kl
    spare1 = u(1:ifull,k)**2 + v(1:ifull,k)**2
    spare2 = u_max(1:ifull,k)**2 + v_max(1:ifull,k)**2
    where ( spare2 > spare1 )
      u_max(1:ifull,k) = u(1:ifull,k)
      v_max(1:ifull,k) = v(1:ifull,k)
    end where
  end do
  spare1 = u(1:ifull,1)**2 + v(1:ifull,1)**2
  where ( u10(:)**2 > u10m_max(:)**2 + v10m_max(:)**2 )
    u10m_max(:) = u10(:)*u(1:ifull,1)/max(.001,sqrt(spare1(1:ifull)))
    v10m_max(:) = u10(:)*v(1:ifull,1)/max(.001,sqrt(spare1(1:ifull)))
  end where
end if

if ( ktau==ntau .or. mod(ktau,nperavg)==0 ) then
  cape_ave(1:ifull)          = cape_ave(1:ifull)/min(ntau,nperavg)
  dew_ave(1:ifull)           = dew_ave(1:ifull)/min(ntau,nperavg)
  epan_ave(1:ifull)          = epan_ave(1:ifull)/min(ntau,nperavg)
  epot_ave(1:ifull)          = epot_ave(1:ifull)/min(ntau,nperavg)
  eg_ave(1:ifull)            = eg_ave(1:ifull)/min(ntau,nperavg)
  fg_ave(1:ifull)            = fg_ave(1:ifull)/min(ntau,nperavg)
  ga_ave(1:ifull)            = ga_ave(1:ifull)/min(ntau,nperavg) 
  anthropogenic_ave(1:ifull) = anthropogenic_ave(1:ifull)/min(ntau,nperavg)
  anth_elecgas_ave(1:ifull)  = anth_elecgas_ave(1:ifull)/min(ntau,nperavg)
  anth_heating_ave(1:ifull)  = anth_heating_ave(1:ifull)/min(ntau,nperavg)
  anth_cooling_ave(1:ifull)  = anth_cooling_ave(1:ifull)/min(ntau,nperavg) 
  rnet_ave(1:ifull)          = rnet_ave(1:ifull)/min(ntau,nperavg)
  !riwp_ave(1:ifull)          = riwp_ave(1:ifull)/min(ntau,nperavg)
  !rlwp_ave(1:ifull)          = rlwp_ave(1:ifull)/min(ntau,nperavg)
  tscr_ave(1:ifull)          = tscr_ave(1:ifull)/min(ntau,nperavg)
  rhscr_ave(1:ifull)         = rhscr_ave(1:ifull)/min(ntau,nperavg)
  do k = 1,ms
    wb_ave(1:ifull,k)    = wb_ave(1:ifull,k)/min(ntau,nperavg)
    wbice_ave(1:ifull,k) = wbice_ave(1:ifull,k)/min(ntau,nperavg)
  end do
  taux_ave(1:ifull)    = taux_ave(1:ifull)/min(ntau,nperavg)
  tauy_ave(1:ifull)    = tauy_ave(1:ifull)/min(ntau,nperavg)
  sgn_ave(1:ifull)     = sgn_ave(1:ifull)/min(ntau,nperavg)  ! Dec07 because of solar fit
  sgdn_ave(1:ifull)    = sgdn_ave(1:ifull)/min(ntau,nperavg) ! because of solar fit
  sgdndir_ave(1:ifull) = sgdndir_ave(1:ifull)/min(ntau,nperavg) ! because of solar fit
  if ( always_mspeca ) then
    sint_ave(1:ifull)   = sint_ave(1:ifull)/max(koundiag,1)
    sot_ave(1:ifull)    = sot_ave(1:ifull)/max(koundiag,1)
    soc_ave(1:ifull)    = soc_ave(1:ifull)/max(koundiag,1)
    rtu_ave(1:ifull)    = rtu_ave(1:ifull)/max(koundiag,1)
    rtc_ave(1:ifull)    = rtc_ave(1:ifull)/max(koundiag,1)
    rgdn_ave(1:ifull)   = rgdn_ave(1:ifull)/max(koundiag,1)
    rgn_ave(1:ifull)    = rgn_ave(1:ifull)/max(koundiag,1)
    rgc_ave(1:ifull)    = rgc_ave(1:ifull)/max(koundiag,1)
    rgdc_ave(1:ifull)   = rgdc_ave(1:ifull)/max(koundiag,1)
    sgc_ave(1:ifull)    = sgc_ave(1:ifull)/max(koundiag,1)
    sgdc_ave(1:ifull)   = sgdc_ave(1:ifull)/max(koundiag,1)
    cld_ave(1:ifull)    = cld_ave(1:ifull)/max(koundiag,1)
    cll_ave(1:ifull)    = cll_ave(1:ifull)/max(koundiag,1)
    clm_ave(1:ifull)    = clm_ave(1:ifull)/max(koundiag,1)
    clh_ave(1:ifull)    = clh_ave(1:ifull)/max(koundiag,1)      
  else
    sint_ave(1:ifull)   = sint_ave(1:ifull)/min(ntau,nperavg)
    sot_ave(1:ifull)    = sot_ave(1:ifull)/min(ntau,nperavg)
    soc_ave(1:ifull)    = soc_ave(1:ifull)/min(ntau,nperavg)
    rtu_ave(1:ifull)    = rtu_ave(1:ifull)/min(ntau,nperavg)
    rtc_ave(1:ifull)    = rtc_ave(1:ifull)/min(ntau,nperavg)
    rgdn_ave(1:ifull)   = rgdn_ave(1:ifull)/min(ntau,nperavg)
    rgn_ave(1:ifull)    = rgn_ave(1:ifull)/min(ntau,nperavg)
    rgc_ave(1:ifull)    = rgc_ave(1:ifull)/min(ntau,nperavg)
    rgdc_ave(1:ifull)   = rgdc_ave(1:ifull)/min(ntau,nperavg)
    sgc_ave(1:ifull)    = sgc_ave(1:ifull)/min(ntau,nperavg)
    sgdc_ave(1:ifull)   = sgdc_ave(1:ifull)/min(ntau,nperavg)
    cld_ave(1:ifull)    = cld_ave(1:ifull)/min(ntau,nperavg)
    cll_ave(1:ifull)    = cll_ave(1:ifull)/min(ntau,nperavg)
    clm_ave(1:ifull)    = clm_ave(1:ifull)/min(ntau,nperavg)
    clh_ave(1:ifull)    = clh_ave(1:ifull)/min(ntau,nperavg)
  end if
  dni_ave(1:ifull)    = dni_ave(1:ifull)/min(ntau,nperavg) 
  cbas_ave(1:ifull)   = 1.1 - cbas_ave(1:ifull)/max(1.e-4,real(precc(:)))  ! 1.1 for no precc
  ctop_ave(1:ifull)   = 1.1 - ctop_ave(1:ifull)/max(1.e-4,real(precc(:)))  ! 1.1 for no precc
  hailrad_ave(1:ifull) = hailrad_ave(1:ifull)/min(ntau,nperavg)
  updraft_helicity_max(1:ifull) = max( updraft_helicity_max, updraft_helicity )
  updraft_helicity_min(1:ifull) = min( updraft_helicity_min, updraft_helicity )
 
  if ( ngas>0 ) then
    traver(1:ifull,1:kl,1:ngas) = traver(1:ifull,1:kl,1:ngas)/min(ntau,nperavg)
  end if

  if ( ccycle/=0 ) then
    fnee_ave(1:ifull)   = fnee_ave(1:ifull)/min(ntau,nperavg)  
    fpn_ave(1:ifull)    = fpn_ave(1:ifull)/min(ntau,nperavg)
    frd_ave(1:ifull)    = frd_ave(1:ifull)/min(ntau,nperavg)
    frp_ave(1:ifull)    = frp_ave(1:ifull)/min(ntau,nperavg)
    frpw_ave(1:ifull)   = frpw_ave(1:ifull)/min(ntau,nperavg)
    frpr_ave(1:ifull)   = frpr_ave(1:ifull)/min(ntau,nperavg)
    frs_ave(1:ifull)    = frs_ave(1:ifull)/min(ntau,nperavg)
    cnpp_ave(1:ifull)   = cnpp_ave(1:ifull)/min(ntau,nperavg)
    cnbp_ave(1:ifull)   = cnbp_ave(1:ifull)/min(ntau,nperavg)
    if ( diaglevel_carbon > 0 ) then
      fevc_ave(1:ifull)   = fevc_ave(1:ifull)/min(ntau,nperavg)
      plant_turnover_ave(1:ifull)      = plant_turnover_ave(1:ifull)/min(ntau,nperavg)
      plant_turnover_wood_ave(1:ifull) = plant_turnover_wood_ave(1:ifull)/min(ntau,nperavg)
    end if
  end if
   
  if ( abs(iaero)>=2 ) then
    duste         = duste/min(ntau,nperavg)        ! Dust emissions
    dustdd        = dustdd/min(ntau,nperavg)       ! Dust dry deposition
    dustwd        = dustwd/min(ntau,nperavg)       ! Dust wet deposition
    dust_burden   = dust_burden/min(ntau,nperavg)  ! Dust burden
    bce           = bce/min(ntau,nperavg)          ! Black carbon emissions
    bcdd          = bcdd/min(ntau,nperavg)         ! Black carbon dry deposition
    bcwd          = bcwd/min(ntau,nperavg)         ! Black carbon wet deposition
    bc_burden     = bc_burden/min(ntau,nperavg)    ! Black carbon burden
    oce           = oce/min(ntau,nperavg)          ! Organic carbon emissions
    ocdd          = ocdd/min(ntau,nperavg)         ! Organic carbon dry deposition
    ocwd          = ocwd/min(ntau,nperavg)         ! Organic carbon wet deposition
    oc_burden     = oc_burden/min(ntau,nperavg)    ! Organic carbon burden
    dmse          = dmse/min(ntau,nperavg)         ! DMS emissions
    dmsso2o       = dmsso2o/min(ntau,nperavg)      ! DMS -> SO2 oxidation
    so2e          = so2e/min(ntau,nperavg)         ! SO2 emissions
    so2so4o       = so2so4o/min(ntau,nperavg)      ! SO2 -> SO4 oxidation
    so2dd         = so2dd/min(ntau,nperavg)        ! SO2 dry deposition
    so2wd         = so2wd/min(ntau,nperavg)        ! SO2 wet deposiion
    so4e          = so4e/min(ntau,nperavg)         ! SO4 emissions
    so4dd         = so4dd/min(ntau,nperavg)        ! SO4 dry deposition
    so4wd         = so4wd/min(ntau,nperavg)        ! SO4 wet deposition
    dms_burden    = dms_burden/min(ntau,nperavg)   ! DMS burden
    so2_burden    = so2_burden/min(ntau,nperavg)   ! SO2 burden
    so4_burden    = so4_burden/min(ntau,nperavg)   ! SO4 burden
    salte         = salte/min(ntau,nperavg)        ! Salt emissions
    saltdd        = saltdd/min(ntau,nperavg)       ! Salt dry deposition
    saltwd        = saltwd/min(ntau,nperavg)       ! Salt wet deposition
    salt_burden   = salt_burden/min(ntau,nperavg)  ! Salt burden
  end if

end if    ! (ktau==ntau.or.mod(ktau,nperavg)==0)

return
end subroutine calculate_timeaverage

!--------------------------------------------------------------
! output diagnostics to log file    
subroutine write_diagnostics(mins_gmt,nmaxprsav)

use arrays_m                               ! Atmosphere dyamics prognostic arrays
use cc_mpi                                 ! CC MPI routines
use cfrac_m                                ! Cloud fraction
use const_phys                             ! Physical constants
use dates_m                                ! Date data
use diag_m                                 ! Diagnostic routines
use estab                                  ! Liquid saturation function
use extraout_m                             ! Additional diagnostics
use histave_m                              ! Time average arrays
use kuocom_m                               ! JLM convection
use liqwpar_m                              ! Cloud water mixing ratios
use morepbl_m                              ! Additional boundary layer diagnostics
use newmpar_m                              ! Grid parameters
use nsibd_m                                ! Land-surface arrays
use parm_m                                 ! Model configuration
use pbl_m                                  ! Boundary layer arrays
use prec_m                                 ! Precipitation
use raddiag_m                              ! Radiation diagnostic
use screen_m                               ! Screen level diagnostics
use sigs_m                                 ! Atmosphere sigma levels
use soil_m                                 ! Soil and surface data
use soilsnow_m                             ! Soil, snow and surface data
use sumdd_m                                ! High precision sum
use tracers_m                              ! Tracer data
use vegpar_m                               ! Vegetation arrays
use vvel_m                                 ! Additional vertical velocity
use work2_m                                ! Diagnostic arrays
use work3_m                                ! Mk3 land-surface diagnostic arrays
use xyzinfo_m                              ! Grid coordinate arrays

implicit none

integer, intent(in) :: mins_gmt, nmaxprsav
integer iq, k, isoil
real, dimension(ifull) :: duma
real, dimension(ifull,kl) :: dums
real, dimension(kl) :: spmean
real, dimension(ifull,9) :: tmparr
complex, dimension(9) :: local_sum, global_sum
real qtot, pwater, es, psavge, spavge, pslavge
real preccavge, precavge, gke, clhav, cllav
real clmav, cltav

if ( mod(ktau,nmaxpr)==0 .and. mydiag ) then
  write(6,*)
  write (6,"('ktau =',i5,' gmt(h,m):',f6.2,i5,' runtime(h,m):',f7.2,i6)") ktau,timeg,mins_gmt,timer,mtimer
  ! some surface (or point) diagnostics
  isoil = isoilm(idjd)
  write(6,*) 'land,isoil,ivegt,isflag ',land(idjd),isoil,ivegt(idjd),isflag(idjd)
  write (6,"('snage,snowd,alb   ',f8.4,2f8.2)") snage(idjd),snowd(idjd),albvisnir(idjd,1)
  write (6,"('sicedep,fracice,runoff ',3f8.2)") sicedep(idjd),fracice(idjd),runoff_ave(idjd)
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
  write (6,"('wetfac,sno,evap,precc,precip',6f8.2)") wetfac(idjd),sno(idjd),evspsbl_ave(idjd),precc(idjd),precip(idjd)
  write (6,"('tmin,tmax,tscr,tss,tpan',9f8.2)") tminscr(idjd),tmaxscr(idjd),tscrn(idjd),tss(idjd),tpan(idjd)
  write (6,"('u10,ustar,pblh',9f8.2)") u10(idjd),ustar(idjd),pblh(idjd)
  write (6,"('ps,qgscrn',5f8.2,f8.3)") .01*ps(idjd),1000.*qgscrn(idjd)
  write (6,"('dew_,eg_,epot,epan,eg,fg,ga',9f8.2)") dew_ave(idjd),eg_ave(idjd),epot(idjd),epan(idjd),eg(idjd),fg(idjd),ga(idjd)
  write (6,"('zo,cduv',2f8.5)") zo(idjd),cduv(idjd)/vmod(idjd)
  write (6,"('slwa,sint,sg,rt,rg    ',9f8.2)") slwa(idjd),sint(idjd),sgsave(idjd),rt(idjd),rgsave(idjd)
  write (6,"('cll,clm,clh,clt ',9f8.2)") cloudlo(idjd),cloudmi(idjd),cloudhi(idjd),cloudtot(idjd)
  write (6,"('u10max,v10max,rhmin,rhmax   ',9f8.2)") u10max(iq),v10max(iq),rhminscr(iq),rhmaxscr(iq)
  write (6,"('kbsav,ktsav,convpsav ',2i3,f8.4,9f8.2)") kbsav(idjd),ktsav(idjd),convpsav(idjd)
  spmean(:) = t(idjd,:)
  write (6,"('t   ',9f8.3/4x,9f8.3)") spmean(:)
  spmean(:) = u(idjd,:)
  write (6,"('u   ',9f8.3/4x,9f8.3)") spmean(:)
  spmean(:) = v(idjd,:)
  write (6,"('v   ',9f8.3/4x,9f8.3)") spmean(:)
  spmean(:) = qg(idjd,:)
  write (6,"('qg  ',9f8.3/4x,9f8.3)") spmean(:)
  spmean(:) = qfg(idjd,:)
  write (6,"('qf  ',9f8.3/4x,9f8.3)") spmean(:)
  spmean(:) = qlg(idjd,:)
  write (6,"('ql  ',9f8.3/4x,9f8.3)") spmean(:)
  spmean(:) = cfrac(idjd,:)
  write (6,"('cfrac',9f8.3/5x,9f8.3)") spmean(:)
  do k = 1,kl
    es        = establ(t(idjd,k))
    spmean(k) = 100.*qg(idjd,k)*max(ps(idjd)*sig(k)-es,1.)/(.622*es) ! max as for convjlm
  enddo
  write (6,"('rh  ',9f8.3/4x,9f8.3)") spmean(:)
  spmean(:) = ps(idjd)*dpsldt(idjd,:)
  write (6,"('omgf ',9f8.3/5x,9f8.3)") spmean(:) ! in Pa/s
  write (6,"('sdot ',9f8.3/5x,9f8.3)") sdot(idjd,1:kl)
  if ( nextout >= 4 ) then
    write (6,"('xlat,long,pres ',3f8.2)") tr(idjd,nlv,ngas+1),tr(idjd,nlv,ngas+2),tr(idjd,nlv,ngas+3)
  end if
endif  ! (mod(ktau,nmaxpr)==0.and.mydiag)
  
if ( ndi==-ktau ) then
  nmaxpr = 1         ! diagnostic prints; reset 6 lines on
  if ( ndi2==0 ) ndi2 = ktau + 40
endif
if ( ktau==ndi2 ) then
  if ( myid==0 ) write(6,*) 'reset nmaxpr'
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
  if ( myid==0 ) then
    write(6,'("spmean ",9f8.3)') spmean
    write(6,'("spavge ",f8.3)') spavge
  end if
  dums = qg(1:ifull,:)
  call average(dums,spmean,spavge)
  if ( myid==0 ) then
    write(6,'("qgmean ",9f8.5)') spmean
    write(6,'("qgavge ",f8.5)') spavge
  end if
  call maxmin(wb,'wb',ktau,1.,ms)
  call maxmin(tggsn,'tggsn',ktau,1.,3)
  call maxmin(tgg,'tg',ktau,1.,ms)
  call maxmin(tss,'ts',ktau,1.,1)
  call maxmin(pblh,'pb',ktau,1.,1)
  duma = real(precip)
  call maxmin(duma,'pr',ktau,1.,1)
  duma = real(precc)
  call maxmin(duma,'pc',ktau,1.,1)
  call maxmin(convpsav,'co',ktau,1.,1)
  duma = real(sno)
  call maxmin(duma,'sn',ktau,1.,1)        ! as mm during timestep
  call maxmin(rhscrn,'rh',ktau,1.,1)
  call maxmin(ps,'ps',ktau,.01,1)
  ! MJT notes, these lines need SUMDD
  local_sum(:) = cmplx(0.,0.)
  tmparr(1:ifull,1) = ps(1:ifull)*wts(1:ifull)
  tmparr(1:ifull,2) = psl(1:ifull)*wts(1:ifull)
  tmparr(1:ifull,3) = real(precc(1:ifull))*wts(1:ifull)
  tmparr(1:ifull,4) = real(precip(1:ifull))*wts(1:ifull)
  ! KE calculation, not taking into account pressure weighting  
  tmparr(1:ifull,5) = 0.
  do k = 1,kl
    tmparr(1:ifull,5) = tmparr(1:ifull,5) - 0.5 * wts(1:ifull) * dsig(k) * ( u(1:ifull,k)**2 + v(1:ifull,k)**2 )
  end do
  tmparr(1:ifull,6) = cloudlo(1:ifull)*wts(1:ifull)
  tmparr(1:ifull,7) = cloudmi(1:ifull)*wts(1:ifull)
  tmparr(1:ifull,8) = cloudhi(1:ifull)*wts(1:ifull)
  tmparr(1:ifull,9) = cloudtot(1:ifull)*wts(1:ifull)
  call drpdr_local_v(tmparr(:,1:9), local_sum(1:9))
  ! All this combined into a single reduction
  global_sum(:) = cmplx(0.,0.)
  call ccmpi_allreduce( local_sum, global_sum, "sumdr", comm_world )
  psavge    = real(global_sum(1))
  pslavge   = real(global_sum(2))
  preccavge = real(global_sum(3))
  precavge  = real(global_sum(4))
  gke       = real(global_sum(5))
  cllav     = real(global_sum(6))
  clmav     = real(global_sum(7))
  clhav     = real(global_sum(8))
  cltav     = real(global_sum(9))
  if ( myid==0 ) then
    write(6,97) psavge,pslavge,preccavge,precavge,gke
97  format(' average ps, psl, precc, prec, gke: ',f10.2,f10.6,2f6.2,f7.2)
    write(6,971) cllav,clmav,clhav,cltav
971 format(' global_average cll, clm, clh, clt: ',4f6.2)
  end if
  if ( mydiag ) then
    write(6,98) ktau,diagvals(ps)
98  format(i7,' ps diag:',9f9.1)
    if ( t(idjd,kl)>258. ) then
      write(6,*) 't(idjd,kl) > 258. for idjd = ',idjd
      write(6,91) ktau,(t(idjd,k),k=kl-8,kl)
91    format(i7,'    t',9f7.2)
      write(6,92) ktau,(sdot(idjd,k),k=kl-8,kl)
92    format(i7,' sdot',9f7.3)
    end if             ! (t(idjd,kl)>258.)
  end if               ! myid==0
endif                  ! (mod(ktau,nmaxpr)==0)

return
end subroutine write_diagnostics
    
!-------------------------------------------------------------------- 
! Check for NaN errors
subroutine nantest(message,js,je,mode)

use aerosol_arrays, only : xtg,naero  ! Aerosol arrays
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
integer, dimension(2) :: posmin, posmax
integer, dimension(3) :: posmin3, posmax3
character(len=*), intent(in) :: message, mode
logical t_check, uv_check, qv_check, qtot_check, qrad_check
logical ps_check, tss_check, aero_check, turb_check

if ( qg_fix<=-1 ) return

if ( js<1 .or. je>ifull ) then
  write(6,*) "ERROR: Invalid index for nantest - ",trim(message)
  call ccmpi_abort(-1)
end if

#ifdef debug
t_check = .true.
uv_check = .true.
qv_check = .true.
qtot_check = .true.
qrad_check = .true.
ps_check = .true.
tss_check = .true.
aero_check = .true.
turb_check = .true.
#else
select case(mode)
  case("gdrag")
    t_check = .false.
    uv_check = .true.
    qv_check = .false.
    qtot_check = .false.
    qrad_check = .false.
    ps_check = .false.
    tss_check = .false.
    aero_check = .false.
    turb_check = .false.
  case("conv")
    t_check = .true.
    uv_check = .true.
    qv_check = .true.
    qtot_check = .true.
    qrad_check = .false.
    ps_check = .false.
    tss_check = .false.
    aero_check = .false.
    turb_check = .false.   
  case("cloud")
    t_check = .true.
    uv_check = .false.
    qv_check = .true.
    qtot_check = .true.
    qrad_check = .true.
    ps_check = .false.
    tss_check = .false.
    aero_check = .false.
    turb_check = .false. 
  case("radiation")
    t_check = .true.
    uv_check = .false.
    qv_check = .false.
    qtot_check = .false.
    qrad_check = .false.
    ps_check = .false.
    tss_check = .false.
    aero_check = .false.
    turb_check = .false.     
  case("surface")
    t_check = .false.
    uv_check = .false.
    qv_check = .false.
    qtot_check = .false.
    qrad_check = .false.
    ps_check = .false.
    tss_check = .true.
    aero_check = .false.
    turb_check = .true.  
  case("vmixing")
    t_check = .true.
    uv_check = .true.
    qv_check = .true.
    qtot_check = .true.
    qrad_check = .false.
    ps_check = .false.
    tss_check = .false.
    aero_check = .false.
    turb_check = .false. 
  case("aerosol")
    t_check = .false.
    uv_check = .false.
    qv_check = .false.
    qtot_check = .false.
    qrad_check = .false.
    ps_check = .false.
    tss_check = .false.
    aero_check = .true.
    turb_check = .false.
  case("nesting")
    t_check = .true.
    uv_check = .true.
    qv_check = .true.
    qtot_check = .false.
    qrad_check = .false.
    ps_check = .true.
    tss_check = .false.
    aero_check = .true.
    turb_check = .false.
  case("dynamics")
    t_check = .true.
    uv_check = .true.
    qv_check = .true.
    qtot_check = .true.
    qrad_check = .false.
    ps_check = .true.
    tss_check = .true.
    aero_check = .true.
    turb_check = .false.
  case default ! all
    t_check = .true.
    uv_check = .true.
    qv_check = .true.
    qtot_check = .true.
    qrad_check = .true.
    ps_check = .true.
    tss_check = .true.
    aero_check = .true.
    turb_check = .true.
end select
#endif  

if ( t_check ) then  
  if ( any(t(js:je,1:kl)/=t(js:je,1:kl)) ) then
    write(6,*) "ERROR: NaN detected in t on myid=",myid," at ",trim(message)
    call ccmpi_abort(-1)
  end if
  if ( any(t(js:je,1:kl)<75.) .or. any(t(js:je,1:kl)>450.) ) then
    write(6,*) "ERROR: Out-of-range detected in t on myid=",myid," at ",trim(message)
    write(6,*) "minval,maxval ",minval(t(js:je,1:kl)),maxval(t(js:je,1:kl))
    posmin = minloc(t(js:je,1:kl))
    posmax = maxloc(t(js:je,1:kl))
    posmin(1) = posmin(1) + js - 1
    posmax(1) = posmax(1) + js - 1
    posmin(1) = iq2iqg(posmin(1))
    posmax(1) = iq2iqg(posmax(1))
    write(6,*) "minloc,maxloc ",posmin,posmax
    call ccmpi_abort(-1)
  end if
end if  

if ( uv_check ) then
  if ( any(u(js:je,1:kl)/=u(js:je,1:kl)) ) then
    write(6,*) "ERROR: NaN detected in u on myid=",myid," at ",trim(message)
    call ccmpi_abort(-1)
  end if
  if ( any(u(js:je,1:kl)<-400.) .or. any(u(js:je,1:kl)>400.) ) then
    write(6,*) "ERROR: Out-of-range detected in u on myid=",myid," at ",trim(message)
    write(6,*) "minval,maxval ",minval(u(js:je,1:kl)),maxval(u(js:je,1:kl))
    posmin = minloc(u(js:je,1:kl))
    posmax = maxloc(u(js:je,1:kl))
    posmin(1) = posmin(1) + js - 1
    posmax(1) = posmax(1) + js - 1
    posmin(1) = iq2iqg(posmin(1))
    posmax(1) = iq2iqg(posmax(1))
    write(6,*) "minloc,maxloc ",posmin,posmax
    call ccmpi_abort(-1) 
  end if
  if ( any(v(js:je,1:kl)/=v(js:je,1:kl)) ) then
    write(6,*) "ERROR: NaN detected in v on myid=",myid," at ",trim(message)
    call ccmpi_abort(-1)
  end if
  if ( any(v(js:je,1:kl)<-400.) .or. any(v(js:je,1:kl)>400.) ) then
    write(6,*) "ERROR: Out-of-range detected in v on myid=",myid," at ",trim(message)
    write(6,*) "minval,maxval ",minval(v(js:je,1:kl)),maxval(v(js:je,1:kl))
    posmin = minloc(v(js:je,1:kl))
    posmax = maxloc(v(js:je,1:kl))
    posmin(1) = posmin(1) + js - 1
    posmax(1) = posmax(1) + js - 1
    posmin(1) = iq2iqg(posmin(1))
    posmax(1) = iq2iqg(posmax(1))
    write(6,*) "minloc,maxloc ",posmin,posmax
    call ccmpi_abort(-1) 
  end if
end if  

if ( qv_check ) then
  if ( any(qg(js:je,1:kl)/=qg(js:je,1:kl)) ) then
    write(6,*) "ERROR: NaN detected in qg on myid=",myid," at ",trim(message)
    call ccmpi_abort(-1)
  end if
  if ( any(qg(js:je,1:kl)<-1.e-8) .or. any(qg(js:je,1:kl)>9.e-2) ) then
    write(6,*) "ERROR: Out-of-range detected in qg on myid=",myid," at ",trim(message)
    write(6,*) "minval,maxval ",minval(qg(js:je,1:kl)),maxval(qg(js:je,1:kl))
    posmin = minloc(qg(js:je,1:kl))
    posmax = maxloc(qg(js:je,1:kl))
    posmin(1) = posmin(1) + js - 1
    posmax(1) = posmax(1) + js - 1
    posmin(1) = iq2iqg(posmin(1))
    posmax(1) = iq2iqg(posmax(1))
    write(6,*) "minloc,maxloc ",posmin,posmax
    call ccmpi_abort(-1) 
  end if
end if  

if ( qtot_check ) then
  if ( any(qlg(js:je,1:kl)/=qlg(js:je,1:kl)) ) then
    write(6,*) "ERROR: NaN detected in qlg on myid=",myid," at ",trim(message)
    call ccmpi_abort(-1)    
  end if
  if ( any(qlg(js:je,1:kl)<-1.e-8) .or. any(qlg(js:je,1:kl)>2.e-1) ) then
    write(6,*) "ERROR: Out-of-range detected in qlg on myid=",myid," at ",trim(message)
    write(6,*) "minval,maxval ",minval(qlg(js:je,1:kl)),maxval(qlg(js:je,1:kl))
    posmin = minloc(qlg(js:je,1:kl))
    posmax = maxloc(qlg(js:je,1:kl))
    posmin(1) = posmin(1) + js - 1
    posmax(1) = posmax(1) + js - 1
    posmin(1) = iq2iqg(posmin(1))
    posmax(1) = iq2iqg(posmax(1))
    write(6,*) "minloc,maxloc ",posmin,posmax
    call ccmpi_abort(-1) 
  end if
  if ( any(qfg(js:je,1:kl)/=qfg(js:je,1:kl)) ) then
    write(6,*) "ERROR: NaN detected in qfg on myid=",myid," at ",trim(message)
    call ccmpi_abort(-1)    
  end if
  if ( any(qfg(js:je,1:kl)<-1.e-8) .or. any(qfg(js:je,1:kl)>2.e-1) ) then
    write(6,*) "ERROR: Out-of-range detected in qfg on myid=",myid," at ",trim(message)
    write(6,*) "minval,maxval ",minval(qfg(js:je,1:kl)),maxval(qfg(js:je,1:kl))
    posmin = minloc(qfg(js:je,1:kl))
    posmax = maxloc(qfg(js:je,1:kl))
    posmin(1) = posmin(1) + js - 1
    posmax(1) = posmax(1) + js - 1
    posmin(1) = iq2iqg(posmin(1))
    posmax(1) = iq2iqg(posmax(1))
    write(6,*) "minloc,maxloc ",posmin,posmax
    call ccmpi_abort(-1) 
  end if
  if ( any(qrg(js:je,1:kl)/=qrg(js:je,1:kl)) ) then
    write(6,*) "ERROR: NaN detected in qrg on myid=",myid," at ",trim(message)
    call ccmpi_abort(-1)    
  end if
  if ( any(qrg(js:je,1:kl)<-1.e-8) ) then
    write(6,*) "ERROR: Out-of-range detected in qrg on myid=",myid," at ",trim(message)
    write(6,*) "minval ",minval(qrg(js:je,1:kl))
    posmin = minloc(qrg(js:je,1:kl))
    posmax = maxloc(qrg(js:je,1:kl))
    posmin(1) = posmin(1) + js - 1
    posmax(1) = posmax(1) + js - 1
    posmin(1) = iq2iqg(posmin(1))
    posmax(1) = iq2iqg(posmax(1))
    write(6,*) "minloc,maxloc ",posmin,posmax
    call ccmpi_abort(-1) 
  end if
  if ( any(qsng(js:je,1:kl)/=qsng(js:je,1:kl)) ) then
    write(6,*) "ERROR: NaN detected in qsng on myid=",myid," at ",trim(message)
    call ccmpi_abort(-1)    
  end if
  if ( any(qsng(js:je,1:kl)<-1.e-8) ) then
    write(6,*) "ERROR: Out-of-range detected in qsng on myid=",myid," at ",trim(message)
    write(6,*) "minval ",minval(qsng(js:je,1:kl))
    posmin = minloc(qsng(js:je,1:kl))
    posmax = maxloc(qsng(js:je,1:kl))
    posmin(1) = posmin(1) + js - 1
    posmax(1) = posmax(1) + js - 1
    posmin(1) = iq2iqg(posmin(1))
    posmax(1) = iq2iqg(posmax(1))
    write(6,*) "minloc,maxloc ",posmin,posmax
    call ccmpi_abort(-1) 
  end if
  if ( any(qgrg(js:je,1:kl)/=qgrg(js:je,1:kl)) ) then
    write(6,*) "ERROR: NaN detected in qgrg on myid=",myid," at ",trim(message)
    call ccmpi_abort(-1)    
  end if
  if ( any(qgrg(js:je,1:kl)<-1.e-8) ) then
    write(6,*) "ERROR: Out-of-range detected in qgrg on myid=",myid," at ",trim(message)
    write(6,*) "minval ",minval(qgrg(js:je,1:kl))
    posmin = minloc(qgrg(js:je,1:kl))
    posmax = maxloc(qgrg(js:je,1:kl))
    posmin(1) = posmin(1) + js - 1
    posmax(1) = posmax(1) + js - 1
    posmin(1) = iq2iqg(posmin(1))
    posmax(1) = iq2iqg(posmax(1))
    write(6,*) "minloc,maxloc ",posmin,posmax
    call ccmpi_abort(-1) 
  end if
  if ( any(rfrac(js:je,1:kl)/=rfrac(js:je,1:kl)) ) then
    write(6,*) "ERROR: NaN detected in rfrac on myid=",myid," at ",trim(message)
    call ccmpi_abort(-1)    
  end if
  if ( any(rfrac(js:je,1:kl)<-1.e-8) .or. any(rfrac(js:je,1:kl)>1.001) ) then
    write(6,*) "ERROR: Out-of-range detected in rfrac on myid=",myid," at ",trim(message)
    write(6,*) "minval,maxval ",minval(rfrac(js:je,1:kl)),maxval(rfrac(js:je,1:kl))
    posmin = minloc(rfrac(js:je,1:kl))
    posmax = maxloc(rfrac(js:je,1:kl))
    posmin(1) = posmin(1) + js - 1
    posmax(1) = posmax(1) + js - 1
    posmin(1) = iq2iqg(posmin(1))
    posmax(1) = iq2iqg(posmax(1))
    write(6,*) "minloc,maxloc ",posmin,posmax
    call ccmpi_abort(-1) 
  end if
  if ( any(sfrac(js:je,1:kl)/=sfrac(js:je,1:kl)) ) then
    write(6,*) "ERROR: NaN detected in sfrac on myid=",myid," at ",trim(message)
    call ccmpi_abort(-1)    
  end if
  if ( any(sfrac(js:je,1:kl)<-1.e-8) .or. any(sfrac(js:je,1:kl)>1.001) ) then
    write(6,*) "ERROR: Out-of-range detected in sfrac on myid=",myid," at ",trim(message)
    write(6,*) "minval,maxval ",minval(sfrac(js:je,1:kl)),maxval(sfrac(js:je,1:kl))
    posmin = minloc(sfrac(js:je,1:kl))
    posmax = maxloc(sfrac(js:je,1:kl))
    posmin(1) = posmin(1) + js - 1
    posmax(1) = posmax(1) + js - 1
    posmin(1) = iq2iqg(posmin(1))
    posmax(1) = iq2iqg(posmax(1))
    write(6,*) "minloc,maxloc ",posmin,posmax
    call ccmpi_abort(-1) 
  end if
  if ( any(gfrac(js:je,1:kl)/=gfrac(js:je,1:kl)) ) then
    write(6,*) "ERROR: NaN detected in gfrac on myid=",myid," at ",trim(message)
    call ccmpi_abort(-1)    
  end if
  if ( any(gfrac(js:je,1:kl)<-1.e-8) .or. any(gfrac(js:je,1:kl)>1.001) ) then
    write(6,*) "ERROR: Out-of-range detected in gfrac on myid=",myid," at ",trim(message)
    write(6,*) "minval,maxval ",minval(gfrac(js:je,1:kl)),maxval(gfrac(js:je,1:kl))
    posmin = minloc(gfrac(js:je,1:kl))
    posmax = maxloc(gfrac(js:je,1:kl))
    posmin(1) = posmin(1) + js - 1
    posmax(1) = posmax(1) + js - 1
    posmin(1) = iq2iqg(posmin(1))
    posmax(1) = iq2iqg(posmax(1))
    write(6,*) "minloc,maxloc ",posmin,posmax
    call ccmpi_abort(-1) 
  end if
end if
  
if ( qrad_check ) then
  if ( any(qlrad(js:je,1:kl)/=qlrad(js:je,1:kl)) ) then
    write(6,*) "ERROR: NaN detected in qlrad on myid=",myid," at ",trim(message)
    call ccmpi_abort(-1)    
  end if
  if ( any(qlrad(js:je,1:kl)<-1.e-8) .or. any(qlrad(js:je,1:kl)>8.e-2) ) then
    write(6,*) "ERROR: Out-of-range detected in qlrad on myid=",myid," at ",trim(message)
    write(6,*) "minval,maxval ",minval(qlrad(js:je,1:kl)),maxval(qlrad(js:je,1:kl))
    posmin = minloc(qlrad(js:je,1:kl))
    posmax = maxloc(qlrad(js:je,1:kl))
    posmin(1) = posmin(1) + js - 1
    posmax(1) = posmax(1) + js - 1
    posmin(1) = iq2iqg(posmin(1))
    posmax(1) = iq2iqg(posmax(1))
    write(6,*) "minloc,maxloc ",posmin,posmax
    call ccmpi_abort(-1) 
  end if
  if ( any(qfrad(js:je,1:kl)/=qfrad(js:je,1:kl)) ) then
    write(6,*) "ERROR: NaN detected in qfrad on myid=",myid," at ",trim(message)
    call ccmpi_abort(-1)    
  end if
  if ( any(qfrad(js:je,1:kl)<-1.e-8) .or. any(qfrad(js:je,1:kl)>8.e-2) ) then
    write(6,*) "ERROR: Out-of-range detected in qfrad on myid=",myid," at ",trim(message)
    write(6,*) "minval,maxval ",minval(qfrad(js:je,1:kl)),maxval(qfrad(js:je,1:kl))
    posmin = minloc(qfrad(js:je,1:kl))
    posmax = maxloc(qfrad(js:je,1:kl))
    posmin(1) = posmin(1) + js - 1
    posmax(1) = posmax(1) + js - 1
    posmin(1) = iq2iqg(posmin(1))
    posmax(1) = iq2iqg(posmax(1))
    write(6,*) "minloc,maxloc ",posmin,posmax
    call ccmpi_abort(-1) 
  end if
  if ( any(cfrac(js:je,1:kl)/=cfrac(js:je,1:kl)) ) then
    write(6,*) "ERROR: NaN detected in cfrac on myid=",myid," at ",trim(message)
    call ccmpi_abort(-1)    
  end if
  if ( any(cfrac(js:je,1:kl)<-1.e-8) .or. any(cfrac(js:je,1:kl)>1.001) ) then
    write(6,*) "ERROR: Out-of-range detected in cfrac on myid=",myid," at ",trim(message)
    write(6,*) "minval,maxval ",minval(cfrac(js:je,1:kl)),maxval(cfrac(js:je,1:kl))
    posmin = minloc(cfrac(js:je,1:kl))
    posmax = maxloc(cfrac(js:je,1:kl))
    posmin(1) = posmin(1) + js - 1
    posmax(1) = posmax(1) + js - 1
    posmin(1) = iq2iqg(posmin(1))
    posmax(1) = iq2iqg(posmax(1))
    write(6,*) "minloc,maxloc ",posmin,posmax
    call ccmpi_abort(-1) 
  end if
end if  

if ( ps_check ) then
  if ( any(psl(js:je)/=psl(js:je)) ) then
    write(6,*) "ERROR: NaN detected in psl on myid=",myid," at ",trim(message)
    call ccmpi_abort(-1)
  end if
  if ( any(psl(js:je)<-1.6) .or. any(psl(js:je)>0.6) ) then
    write(6,*) "ERROR: Out-of-range detected in psl on myid=",myid," at ",trim(message)
    write(6,*) "minval,maxval ",minval(psl(js:je)),maxval(psl(js:je))
    posmin(1:1) = minloc(psl(js:je))
    posmax(1:1) = maxloc(psl(js:je))
    posmin(1) = posmin(1) + js - 1
    posmax(1) = posmax(1) + js - 1
    posmin(1) = iq2iqg(posmin(1))
    posmax(1) = iq2iqg(posmax(1))
    write(6,*) "minloc,maxloc ",posmin(1:1),posmax(1:1)
    call ccmpi_abort(-1) 
  end if
  if ( any(ps(js:je)/=ps(js:je)) ) then
    write(6,*) "ERROR: NaN detected in ps on myid=",myid," at ",trim(message)
    call ccmpi_abort(-1)
  end if
end if

if ( tss_check ) then
  if ( any(tss(js:je)/=tss(js:je)) ) then
    write(6,*) "ERROR: NaN detected in tss on myid=",myid," at ",trim(message)
    call ccmpi_abort(-1)
  end if
  if ( any(tss(js:je)<75.) .or. any(tss(js:je)>425.) ) then
    write(6,*) "ERROR: Out-of-range detected in tss on myid=",myid," at ",trim(message)
    write(6,*) "minval,maxval ",minval(tss(js:je)),maxval(tss(js:je))
    posmin(1:1) = minloc(tss(js:je))
    posmax(1:1) = maxloc(tss(js:je))
    posmin(1) = posmin(1) + js - 1
    posmax(1) = posmax(1) + js - 1
    posmin(1) = iq2iqg(posmin(1))
    posmax(1) = iq2iqg(posmax(1))
    write(6,*) "minloc,maxloc ",posmin(1:1),posmax(1:1)
    call ccmpi_abort(-1) 
  end if
end if

if ( aero_check .and. abs(iaero)>=2 ) then
  if ( any(xtg(js:je,1:kl,1:naero)/=xtg(js:je,1:kl,1:naero)) ) then
    write(6,*) "ERROR: NaN detected in xtg on myid=",myid," at ",trim(message)
    call ccmpi_abort(-1)
  end if
  if ( any(xtg(js:je,1:kl,1:naero)<-1.e-8) .or. any(xtg(js:je,1:kl,1:naero)>2.e-3) ) then
    write(6,*) "ERROR: Out-of-range detected in xtg on myid=",myid," at ",trim(message)
    write(6,*) "minval,maxval ",minval(xtg(js:je,1:kl,1:naero)),maxval(xtg(js:je,1:kl,1:naero))
    posmin3 = minloc(xtg(js:je,1:kl,1:naero))
    posmax3 = maxloc(xtg(js:je,1:kl,1:naero))
    posmin3(1) = posmin3(1) + js - 1
    posmax3(1) = posmax3(1) + js - 1
    posmin3(1) = iq2iqg(posmin3(1))
    posmax3(1) = iq2iqg(posmax3(1))
    write(6,*) "minloc,maxloc ",posmin3,posmax3
    call ccmpi_abort(-1) 
  end if  
end if

if ( turb_check ) then
  if ( any( fg(js:je)/=fg(js:je) ) ) then
    write(6,*) "ERROR: NaN detected in fg on myid=",myid," at ",trim(message)
    call ccmpi_abort(-1)
  end if
  if ( any( eg(js:je)/=eg(js:je) ) ) then
    write(6,*) "ERROR: NaN detected in eg on myid=",myid," at ",trim(message)
    call ccmpi_abort(-1)
  end if
end if

return
end subroutine nantest

!-------------------------------------------------------------------- 
! Calculate change in moist static energy
subroutine calculate_dhdt_mse(js,je,mse_t1)

use arrays_m                               ! Atmosphere dyamics prognostic arrays
use const_phys                             ! Physical constants
use newmpar_m                              ! Grid parameters
use sigs_m                                 ! Atmosphere sigma levels
  
implicit none 
  
integer, intent(in) :: js, je
integer             :: iq, k
real, dimension(js:je, kl), intent(out) :: mse_t1
real, dimension(js:je, kl) :: zo

zo(js:je,1) = bet(1)*t(js:je,1)/grav ! heights above surface
do k = 2,kl
  zo(js:je,k) = zo(js:je,k-1) + (bet(k)*t(js:je,k)+betm(k)*t(js:je,k-1))/grav ! heights above surface
end do

do k = 1,kl
  do iq = js,je
    mse_t1(iq,k)=grav*zo(iq,k)+cp*t(iq,k)+hl*qg(iq,k)
  end do
end do

return
end subroutine calculate_dhdt_mse
