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
    
! This subroutine initialises CCAM prognostic variables and
! surface forcings.  Currently, we assume that disk io is
! slower than MPI, so surface forcings are loaded on a
! single processor and then distributed using MPI_distrubute.
! Note that nested files are split over processors (see
! onthefly.f).

module indata

private
public indataf

interface datacheck
  module procedure rdatacheck, idatacheck
end interface datacheck

contains
    
subroutine indataf(hourst,jalbfix,lapsbot,isoth,nsig,io_nest)
     
use aerointerface                                ! Aerosol interface
use aerosolldr, only : xtg,naero                 ! LDR prognostic aerosols
use arrays_m                                     ! Atmosphere dyamics prognostic arrays
use ateb                                         ! Urban
use bigxy4_m                                     ! Grid interpolation
use cable_ccam, only : loadcbmparm,loadtile      ! CABLE interface
use cc_mpi                                       ! CC MPI routines
use daviesnudge                                  ! Far-field nudging
use diag_m                                       ! Diagnostic routines
use epst_m                                       ! Off-centre terms
use extraout_m                                   ! Additional diagnostics
use gdrag_m                                      ! Gravity wave drag
use indices_m                                    ! Grid index arrays
use infile                                       ! Input file routines
use latlong_m                                    ! Lat/lon coordinates
use liqwpar_m                                    ! Cloud water mixing ratios
use map_m                                        ! Grid map arrays
use mlo                                          ! Ocean physics and prognostic arrays
use mlodynamics                                  ! Ocean dynamics
use morepbl_m                                    ! Additional boundary layer diagnostics
use nharrs_m, only : lrestart                    ! Non-hydrostatic atmosphere arrays
use nsibd_m                                      ! Land-surface arrays
use onthefly_m                                   ! Input interpolation routines
use pbl_m                                        ! Boundary layer arrays
use permsurf_m                                   ! Fixed surface arrays
use river                                        ! River routing
use sigs_m                                       ! Atmosphere sigma levels
use soil_m                                       ! Soil and surface data
use soilsnow_m                                   ! Soil, snow and surface data
use timeseries, only : init_ts                   ! Tracer time series
use tracermodule, only : readtracerflux          ! Tracer routines
use tracers_m                                    ! Tracer data
use vecs_m                                       ! Eigenvectors for atmosphere dynamics
use vecsuv_m                                     ! Map to cartesian coordinates
use vegpar_m                                     ! Vegetation arrays
use xyzinfo_m                                    ! Grid coordinate arrays
      
implicit none
      
include 'newmpar.h'                              ! Grid parameters
include 'const_phys.h'                           ! Physical constants
include 'darcdf.h'                               ! Netcdf data
include 'dates.h'                                ! Date data
include 'filnames.h'                             ! Filenames
include 'parm.h'                                 ! Model configuration
include 'parmdyn.h'                              ! Dynamics parmaters
include 'parmgeom.h'                             ! Coordinate data
include 'soilv.h'                                ! Soil parameters
include 'stime.h'                                ! File date data
include 'trcom2.h'                               ! Station data

integer, parameter :: jlmsigmf=1  ! 1 for jlm fixes to dean's data
integer, parameter :: nfixwb=2    ! 0, 1 or 2; wb fixes with nrungcm=1
integer, parameter :: ntest=0
integer, parameter :: klmax=100   ! Maximum vertical levels

!     for the held-suarez test
real, parameter :: delty = 60.    ! pole to equator variation in equal temperature
real, parameter :: deltheta = 10. ! vertical variation
real, parameter :: rkappa = 2./7.

integer, intent(in) :: jalbfix
integer, intent(inout) :: io_nest
integer ii, imo, indexi, indexl, indexs, ip, iq, isoil, isoth
integer iveg, iyr, jj, k, kdate_sav, ktime_sav, l
integer nface, nn, nsig, i, j, n
integer ierr, ic, jc, iqg, ig, jg
integer isav, jsav, ier, lapsbot

character(len=160) :: co2in,radonin,surfin
character(len=80) :: header

real, intent(out) :: hourst
real, dimension(ifull) :: zss, aa, zsmask
real, dimension(ifull) :: dep, depth, rlai
real, dimension(ifull,5) :: duma
real, dimension(ifull,2) :: ocndwn
real, dimension(ifull,wlev,4) :: mlodwn
real, dimension(ifull,kl,naero) :: xtgdwn
real, dimension(ifull,kl,9) :: dumb
real, dimension(:,:), allocatable :: glob2d
real, dimension(:), allocatable :: davt_g
real, dimension(3*kl+1) :: dumc
real, dimension(1:9) :: swilt_diag, sfc_diag
real, dimension(1:ms) :: wb_tmpry
real rlonx,rlatx,alf
real c, cent
real coslat, coslong, costh, den, diffb, diffg, dist
real epsmax, fracs, fracwet, ftsoil, gwdfac, hefact
real polenx, poleny, polenz, pslavge
real rad, radu, radv, ri, rj, rlat_d, rlon_d
real rlatd, rlongd
real sinlat, sinlong, sinth, snalb,sumdsig, thet, tsoil
real uzon, vmer, wet3, zonx, zony, zonz, zsdiff, tstom
real xbub, ybub, xc, yc, zc, xt, yt, zt, tbubb, emcent
real deli, delj, centi, distnew, distx, rhs, ril2
real newzo,visalb,niralb

! The following look-up tables are for the Mk3 land-surface scheme
real, dimension(44), parameter :: vegpmin = (/                           &
                    .98,.85,.85,.5,.2,.1 ,.85,.5,.2,.5,                  & ! 1-10
                    .2,.1 ,.5,.2,.1 ,.1,.1 ,.85,.5,.2,                   & ! 11-20
                    .1 ,.85,.60,.50,.5 ,.2,.1 ,.5, .0, .0, .4,           & ! 21-31
                    .98,.75,.75,.75,.5,.86,.65,.79,.3, .42,.02,.54,0./)    ! 32-44
real, dimension(44), parameter :: vegpmax = (/                           &
                    .98,.85,.85,.5,.7,.60,.85,.5,.5,.5,                  & ! 1-10
                    .5,.50,.5,.6,.60,.4,.40,.85,.5,.8,                   & ! 11-20
                    .20,.85,.85,.50,.80,.7,.40,.5, .0, .0, .6,           & ! 21-31
                    .98,.75,.75,.75,.5,.86,.65,.79,.3, .42,.02,.54,0./)    ! 32-44
real, dimension(12), parameter :: fracsum =                              &
                     (/-.5,-.5,-.3,-.1,.1, .3, .5, .5, .3, .1,-.1,-.3/)
real, dimension(44), parameter :: fracwets = (/                          &
                    .5, .5, .5, .5, .5, .5, .5, .5, .5, .5,              & !  1-10 summer
                    .5, .5, .5, .5, .5, .5, .5, .5, .5, .5,              & ! 11-20 summer
                    .5, .5, .5, .5, .5, .5, .5, .5, .5, .5, .5,          & ! 21-31 summer
       .5,.5, .3, .3, .3, .15, .15, .15, .1, .15, .02, .35, .5/)           ! 32-44 summer
real, dimension(44), parameter :: fracwetw = (/                          &
                    .5, .5, .5, .5, .5, .5, .5, .5, .5, .5,              & !  1-10 winter
                    .5, .5, .5, .5, .5, .5, .5, .5, .5, .5,              & ! 11-20 winter
                    .5, .5, .5, .5, .5, .5, .5, .5, .5, .5, .5,          & ! 21-31 winter
       .5,.5, .6, .6, .6, .25, .3 , .25, .2, .25, .05, .6, .5 /)           ! 32-44 winter
real vegpsig(44)
data vegpsig/ .98,.85,.85,.5,.2,.05,.85,.5,.2,.5,                        & ! 1-10
              .2,.05,.5,.2,.05,.2,.05,.85,.5,.2,                         & ! 11-20
              .05,.85,.85,.55,.65,.2,.05,.5, .0, .0, .5,                 & ! 21-31
              .98,.75,.75,.75,.5,.86,.65,.79,.3, .42,.02,.54,0./           ! 32-44

real, dimension(klmax) :: tbarr,qgin,zbub
real :: pmsl=1.010e5, thlapse=3.e-3, tsea=290., gauss=2.
real :: heightin=2000., hfact=0.1, uin=0., vin=0.

namelist/tin/gauss,heightin,hfact,pmsl,qgin,tbarr,tsea,uin,vin,thlapse

call START_LOG(indata_begin)

!--------------------------------------------------------------
! SET DEFAULT VALUES
tgg(:,:)=280.
tggsn(:,:)=280.
wb(:,:)=.15
snowd(:)=0. 
condx(:)=0.
zolnd(:)=zobgin
eg(:)=0.
fg(:)=0.
cduv(:)=0.
cdtq(:)=0.
swrsave(:)=0.5
hourst=0.
albvissav(:)=-1.
albvisnir(:,:)=0.3
vlai(:)=0.
ivegt(:)=1
isoilm(:)=1
zs(:)=0.
zsmask(:)=0.
he(:)=0.         
land(:)=.false.
kdate=kdate_s
ktime=ktime_s

!--------------------------------------------------------------
! READ AND PROCESS ATMOSPHERE SIGMA LEVELS
if (myid==0) then
  read(28,*)(dumc(k),k=1,kl),(dumc(2*kl+k),k=1,kl),    &
        (bam(k),k=1,kl),((emat(k,l),k=1,kl),l=1,kl),   &
        ((einv(k,l),k=1,kl),l=1,kl),(qvec(k),k=1,kl),  &
        ((tmat(k,l),k=1,kl),l=1,kl)
  write(6,*) 'kl,lapsbot,sig from eigenv file: ',kl,lapsbot,dumc(1:kl)
  ! File has an sigmh(kl+1) which is not required. Causes bounds violation
  ! to read this.
  read(28,*)(dumc(kl+k),k=1,kl)
  close(28) 
  write(6,*) 'tbar: ',dumc(2*kl+1:3*kl)
  write(6,*) 'bam:  ',bam
       
  ! test netcdf for CABLE input
  if (nsib>=6) then
    call ccnf_open(vegfile,ncidveg,ierr)
    if (ierr==0) then
      dumc(3*kl+1)=1.
    else
      dumc(3*kl+1)=0.  
    end if
  else if (nsib==5) then
    call ccnf_open(vegfile,ncidveg,ierr)
    if (ierr==0) then
      dumc(3*kl+1)=1.
    else
      dumc(3*kl+1)=0.  
    end if
  else
    dumc(3*kl+1)=0.
  end if
endif ! (myid==0)

! distribute vertical and vegfile data to all processors
! dumc(1:kl)   = sig,   dumc(kl+1:2*kl) = sigmh, dumc(2*kl+1:3*kl) = tbar
! dumc(3*kl+1) = lncveg
call ccmpi_bcast(dumc(1:3*kl+1),0,comm_world)
sig   =dumc(1:kl)
sigmh =dumc(kl+1:2*kl)
tbar  =dumc(2*kl+1:3*kl)
lncveg=nint(dumc(3*kl+1))

dsig(1:kl-1)=sigmh(2:kl)-sigmh(1:kl-1)
dsig(kl)=-sigmh(kl)
sumdsig=0.
do k=1,kl
  sumdsig=sumdsig-dsig(k)
  tbardsig(k)=0.
enddo
if (myid==0) write(6,*)'dsig,sumdsig ',dsig,sumdsig
if (isoth>=0) then
  dtmax=1./(sig(1)*log(sig(1)/sig(2)))
  tbardsig(1)=dtmax*(tbar(1)-tbar(2))
  do k=2,kl-1
    tbardsig(k)=(tbar(k+1)-tbar(k-1))/(2.*dsig(k))
  enddo
endif
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
if (myid==0) then
  write(6,*)'rata ',rata
  write(6,*)'ratb ',ratb
  write(6,*)'ratha ',ratha
  write(6,*)'rathb ',rathb
end if
   
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

if ( myid==0 ) then
  write(6,*) 'bet  ',bet
  write(6,*) 'betm ',betm
end if
!if (nh/=0) then
! Non-hydrostatic case
if ( nh==2 .and. lapsbot/=3 ) stop 'nh=2 needs lapsbot=3'
if ( abs(epsp)<=1. ) then
  ! exact treatment when epsp is constant
  call eig(sig,sigmh,tbar,lapsbot,isoth,dt,epsp,epsh,nsig,bet,betm,nh)
else
  call eig(sig,sigmh,tbar,lapsbot,isoth,dt,0.,0.,nsig,bet,betm,nh)
end if
!else
!  ! MJT notes - The hydrostatic case could have called
!  ! eig and avoided a ccmpi_bcast.  However, since eig
!  ! does not always exactly reproduce the input file 
!  ! emat, einv and bam, then we keep the bcast for 
!  ! backwards compatibility
!  call ccmpi_bcast(bam,0,comm_world)
!  call ccmpi_bcast(emat,0,comm_world)
!  call ccmpi_bcast(einv,0,comm_world)
!endif  ! (nh/=0)


! zmin here is approx height of the lowest level in the model
zmin = -rdry*280.*log(sig(1))/grav
if ( myid==0 ) write(6,*) 'zmin = ',zmin


!--------------------------------------------------------------
! READ OROGRAPHY (io_in and nhstest)
!     read in fresh zs, land-sea mask (land where +ve), variances
if ( io_in<=4 .and. nhstest>=0 ) then
  if ( myid==0 ) then
    allocate( glob2d(ifull_g,3) )
    if ( lnctopo==1 ) then
      write(6,*) 'read zs from topofile'
      call surfread(glob2d(:,1),'zs',netcdfid=ncidtopo)
      glob2d(:,1) = grav*glob2d(:,1)
      write(6,*) 'read land-sea fraction'
      call surfread(glob2d(:,2),'lsm',netcdfid=ncidtopo)
      write(6,*) 'read he'
      call surfread(glob2d(:,3),'tsd',netcdfid=ncidtopo)
      call ccnf_close(ncidtopo)
    else
      write(6,*) 'read zs from topofile'
      read(66,*,iostat=ierr) glob2d(:,1)
      if ( ierr/=0 ) stop 'end-of-file reached on topofile'
      write(6,*) 'read land-sea fraction'
      read(66,*,iostat=ierr) glob2d(:,2)
      if ( ierr/=0 ) stop 'end-of-file reached on topofile'
      write(6,*) 'read he'
      read(66,*,iostat=ierr) glob2d(:,3)
      if ( ierr/=0 ) stop 'end-of-file reached on topofile'
      close(66)
    end if
    call ccmpi_distribute(duma(:,1:3),glob2d(:,1:3))
    deallocate(glob2d)
  else
    call ccmpi_distribute(duma(:,1:3))
  end if
  zs(1:ifull)     = duma(:,1)
  zsmask(1:ifull) = duma(:,2)
  he(1:ifull)     = duma(:,3)
  if ( mydiag ) write(6,*) 'he read in from topofile',he(idjd)

  ! special options for orography         
  if ( nspecial==2 ) then  ! to flood Madagascar, or similar 
    do iq=1,ifull
      if(rlatt(iq)*180./pi>-28.and.rlatt(iq)*180./pi<-11.and. &
         rlongg(iq)*180./pi>42.and.rlongg(iq)*180./pi<51.)then
        if(zs(iq)>0.)then   
          write(6,*) 'zeroing iq,zs ',iq,zs(iq)
          zs(iq)=-.6   ! usual sea points come in as -1.  
          zsmask(iq)=0.
        endif  ! (zs(iq)>0.)
      endif    ! (rlatt(iq)*180./pi  ....)
    enddo
  endif       ! (nspecial==2)
  if(nspecial==31)then
    do iq=1,ifull
      if(rlongg(iq)*180./pi>60.and.rlongg(iq)*180./pi<240..and.  &
         rlatt(iq)*180./pi>20.and.rlatt(iq)*180./pi<60.)         &
         zs(iq)=.8*zs(iq)
    enddo
  endif
  if(nspecial==32)then
    do iq=1,ifull
      if(rlongg(iq)*180./pi>60.and.rlongg(iq)*180./pi<240..and.  &
         rlatt(iq)*180./pi>20.and.rlatt(iq)*180./pi<60.)         &
         zs(iq)=.5*zs(iq)
    enddo
  endif
  if(nspecial==33)then
    do iq=1,ifull
      if(rlongg(iq)*180./pi>60.and.rlongg(iq)*180./pi<240..and.  &
         rlatt(iq)*180./pi>20.and.rlatt(iq)*180./pi<60.)         &
         zs(iq)=.2*zs(iq)
    enddo
  endif

  land(1:ifull)=zsmask(1:ifull)>=0.5
 
else                   ! aquaplanet test -1 to -8 or -22
  zs(:)=0.             ! or pgb from June 2003
  zsmask(:)=0.
  he(:)=0.         
  land(:)=.false.
endif  ! (io_in<=4.and.nhstest>=0)  ..else..

if ( mydiag ) then
  write(6,"('zs#_topof ',9f8.1)") diagvals(zs)
  write(6,"('he#_topof ',9f8.1)") diagvals(he)
  write(6,"('zs#_mask ',9f8.2)") diagvals(zsmask)
end if


!-----------------------------------------------------------------
! The following parameterisations require special input data.
! Therefore, these routines are initialised in indata.f instead of
! initialised in globpe.f


!--------------------------------------------------------------
! READ SURFACE DATA (nsib and nspecial)
! nsib=3 (original land surface scheme with original 1deg+Dean's datasets)
! nsib=5 (original land surface scheme with MODIS datasets)
! nsib=6 (CABLE land surface scheme with internal screen diagnostics)
! nsib=7 (CABLE land surface scheme with CCAM screen diagnostics)
if (nsib>=1) then
  call insoil
  call rdnsib
  if (nsib==6.or.nsib==7) then
    ! albvisnir at this point holds soil albedo for cable initialisation
    call loadcbmparm(vegfile,vegprev,vegnext,phenfile,casafile)
    ! albvisnir at this point holds net albedo
  elseif (nsib==3) then
    ! special options for standard land surface scheme
    if(nspecial==35)then      ! test for Andy Cottrill
      do iq=1,ifull
        rlongd=rlongg(iq)*180./pi
        rlatd=rlatt(iq)*180./pi
        if(rlatd>-32..and.rlatd<-23.5)then
          if(rlongd>145..and.rlongd<=150.)ivegt(iq)=4
          if(rlongd>150..and.rlongd<154.)ivegt(iq)=2
        endif
      enddo
    endif  ! (nspecial==35)
    ! zap vegetation over SEQ for Andy
    if(nspecial==41)then
      do iq=1,ifull
        rlongd=rlongg(iq)*180./pi
        rlatd=rlatt(iq)*180./pi
        if(rlatd>-32. .and. rlatd<-23.5)then
          if(rlongd>145. .and. rlongd<=152.)ivegt(iq)=4 
          if(rlongd>152. .and. rlongd< 154.)ivegt(iq)=2 
        endif
      enddo
    endif  ! (nspecial==41)
    do iq=1,ifull
      ! check for littoral veg over Oz      
      rlongd=rlongg(iq)*180./pi
      rlatd=rlatt(iq)*180./pi
      if(rlongd>110.and.rlongd<155.and.rlatd>-45.and.rlatd<-10)then
        if(ivegt(iq)==28)then
          write(6,*)'littoral vegt ',iq,rlongd,rlatd
          if(rlongd>150.and.rlongd<152.and.rlatd>-28.and.rlatd<-26)then
            ivegt(iq)=24   ! fix-up of Graetz data for Andy from July '07
          endif
        endif
      endif
    enddo
    do iq=1,ifull
      if(land(iq))then  
        ! following line from 16/3/06 avoids sand on Himalayas        
        if(zs(iq)>2000.*grav.and.isoilm(iq)<3)isoilm(iq)=3
      endif
    enddo
    ! put in Antarctica ice-shelf fixes 5/3/07
    do iq=1,ifull
      if(zs(iq)<=0.)then
        rlongd=rlongg(iq)*180./pi
        rlatd=rlatt(iq)*180./pi
        if((rlongd>165..and.rlongd<195..and.rlatd<-77.2-(rlongd-165.)/30.).or.      & ! Ross shelf
           (rlongd>300..and.rlongd<330..and.rlatd<-75.2-(rlongd-300.)*2.8/30.).or.  & ! Ronne shelf
           (rlongd>68..and.rlongd<75..and.rlatd<-64.-(rlongd-60.)*5.8/15.))then       ! Amery shelf
             zs(iq)=1.
             land(iq)=.true.
             isoilm(iq)=9
             ivegt(iq)=42
             if(mydiag)write(6,*)'setting sea to ice sheet for iq = ',iq
        endif
      endif  ! (zs(iq)<=0.)
    enddo
  end if ! (nsib==6.or.nsib==7) ..else..
  if (nsib/=6.and.nsib/=7) then
    ! JJK special option for adjusting surface albedo and roughness
    if(nspecial<-10)then
      do iq=1,ifull
        rlongd=rlongg(iq)*180./pi
        rlatd=rlatt(iq)*180./pi
        if((rlatdn<rlatd .and. rlatd<rlatdx).and.(rlongdn<rlongd .and. rlongd<rlongdx))then
          ! assume nspecial = -vvir where vv = % vis alb and ir = % nir alb
          newzo=real(abs(nspecial)/10000)
          visalb=real((abs(nspecial)-nint(newzo)*10000)/100)
          niralb=real(abs(nspecial)-nint(newzo)*10000-nint(visalb)*100)
          if ( newzo  > 1. ) zolnd(iq)=newzo/1000. ! input mm, output m
          if ( visalb > 1. ) albvisnir(iq,1)=visalb/100. ! (to make 0-1)
          if ( niralb > 1. ) albvisnir(iq,2)=niralb/100. ! (to make 0-1)
          if ( .not. land(iq) ) then
            land(iq)=.true.
            ivegt(iq)=42
            isoilm(iq)=7
          end if
        endif
      enddo
    endif  ! (nspecial<-10)	
  end if ! (nsib/=6.and.nsib/=7)
end if   ! nsib>=1


!-----------------------------------------------------------------
! INITIALISE URBAN SCHEME (nurban)
! nurban=0 no urban
! nurban=1 urban (save in restart file)
! nurban=-1 urban (save in history and restart files)
if (nurban/=0) then
  if (myid==0) write(6,*) 'Initialising ateb urban scheme'
  if (lncveg==1) then
    call surfread(sigmu,'urban',netcdfid=ncidveg)
  else
    call surfread(sigmu,'urban',filename=urbanfile)
    sigmu=sigmu*0.01
  end if
  where (.not.land(1:ifull).or.sigmu<0.01)
    sigmu(:)=0.
  end where
  call atebinit(ifull,sigmu(:),0)
else
  sigmu(:)=0.
  call atebdisable(0) ! disable urban
end if

if (myid==0.and.lncveg==1) then
  call ccnf_close(ncidveg)
end if


! fixes for Dean's land-use for CTM
if ( nsib==5 ) then
  where ( sigmu(:)>0.5.and.rlatt(:)*180./pi>-45.and.rlatt(:)*180./pi<-10..and. &
          rlongg(:)*180./pi>112..and.rlongg(:)*180./pi<154.4 )
    ivegt(:) = 31 ! urban
  end where
  where ( isoilm_in(:)==-1.and.rlatt(:)*180./pi>-45.and.rlatt(:)*180./pi<-10..and. &
          rlongg(:)*180./pi>112..and.rlongg(:)*180./pi<154.4 )
    ivegt(:) = 29 ! lake
  end where
end if


!**************************************************************
!**************************************************************
! No changes to land, isoilm or ivegt arrays after this point
!**************************************************************
!**************************************************************


!--------------------------------------------------------------
! LAND SURFACE ERROR CHECKING
if(nsib>=1)then   !  check here for soil & veg mismatches
  if (mydiag) write(6,*)'idjd,land,isoil,ivegt ',idjd,land(idjd),isoilm(idjd),ivegt(idjd)
  do iq=1,ifull
    if(land(iq))then
      if(ivegt(iq)==0)then
        write(6,*)'stopping because nsib>=1 and veg type not defined for iq = ',iq
        write(6,*)'lat,long ',rlatt(iq)*180./pi,rlongg(iq)*180./pi
        call ccmpi_abort(-1)
      endif  ! (ivegt(iq)==0)
      if(isoilm(iq)==0)then
        write(6,*)'stopping because nsib>=1 and soil type not defined for iq = ',iq
        call ccmpi_abort(-1)
      endif  ! (isoilm(iq)==0)
    endif    ! (land(iq))
  enddo    !  iq loop
endif      ! (nsib>=1)


!-----------------------------------------------------------------
! INTIALISE MIXED LAYER OCEAN (nmlo)
! nmlo<0 same as below, but save all data in history file
! nmlo=0 no mixed layer ocean
! nmlo=1 mixed layer ocean (KPP)
! nmlo=2 same as 1, but with Smag horz diffusion and river routing
! nmlo=3 same as 2, but with horizontal and vertical advection
if (nmlo/=0.and.abs(nmlo)<=9) then
  if (myid==0) write(6,*) 'Initialising MLO'
  call surfread(dep,'depth',filename=bathfile)
  where (land)
    dep=0.
  elsewhere
    dep=max(dep,2.*minwater)
  end where
  call mloinit(ifull,dep,0)
  call mlodyninit
end if
if ( abs(nmlo)>=2 .or. nriver==1 ) then
  if (myid==0) write(6,*) 'Initialising river routing'
  call rvrinit
end if


!-----------------------------------------------------------------
! INITIALISE PROGNOSTIC AEROSOLS (iaero)
! iaero=1 prescribed aerosol direct effect
! iaero=2 LDR prognostic aerosol (direct only with nrad=4, direct+indirect with nrad=5)
select case(abs(iaero))
  case(1)
    if(myid==0)write(6,*)'so4total data read from file ',so4tfile
    call readreal(so4tfile,so4t,ifull)
  case(2)
    call load_aerosolldr(so4tfile,oxidantfile,kdate)
end select

  
!--------------------------------------------------------------
! INITIALISE TRACERS (ngas)
if ( myid==0 ) then
  write(6,*)'nllp,ngas,ntrac,ilt,jlt,klt ',nllp,ngas,ntrac,ilt,jlt,klt
end if
if ( ngas>0 ) then
  if ( myid==0 ) write(6,*)'Initialising tracers'
  ! tracer initialisation (if start of run)
  call init_ts(ngas,dt)
  call readtracerflux(kdate)
endif   
  

!--------------------------------------------------------------
! DEFINE FIXED SURFACE ARRAYS
! Note that now only land and water points are allowed, as
! sea-ice can change during the run
indexl=0
do iq=1,ifull
  if(land(iq))then  ! land
    indexl=indexl+1
    iperm(indexl)=iq
  endif ! (land(iq))
enddo   ! iq loop
ipland=indexl
indexi=ipland
ipsea=ifull
indexs=ipsea+1
do iq=1,ifull
  if(.not.land(iq))then
    indexs=indexs-1     ! sea point
    iperm(indexs)=iq    ! sea point
  endif  ! (sicedep(iq)>0.)
enddo   ! iq loop
ipsice=indexs-1
if (mydiag) write(6,*)'ipland,ipsea: ',ipland,ipsea

!-----------------------------------------------------------------
! READ INITIAL CONDITIONS FROM IFILE (io_in)
ncid=-1  ! initialise nc handle with no files open
if ( io_in<4 ) then
  if ( myid==0 ) then
    write(6,*) 'Read initial conditions from ifile'
  end if
  call histopen(ncid,ifile,ier) ! open parallel initial condition files (onthefly will close ncid)
  call ncmsg("ifile",ier)       ! report error messages
  if (myid==0) then
    write(6,*) 'ncid,ifile ',ncid,trim(ifile)
  end if
  kdate_sav=kdate_s
  ktime_sav=ktime_s
  zss=zs(1:ifull)
  if (abs(io_in)==1) then
    call onthefly(0,kdate,ktime,psl,zss,tss,sicedep,fracice,t,u,v, &
                  qg,tgg,wb,wbice,snowd,qfg,qlg,qrg,qsng,qgrg,     &
                  tggsn,smass,ssdn,ssdnn,snage,isflag,mlodwn,      &
                  ocndwn,xtgdwn)
  endif   ! (abs(io_in)==1)
  call histclose
  if(mydiag)then
    write(6,*)'ds,zss',ds,zss(idjd)
    write(6,*)'kdate_sav,ktime_sav ',kdate_sav,ktime_sav
    write(6,*)'kdate_s,ktime_s >= ',kdate_s,ktime_s
    write(6,*)'kdate,ktime ',kdate,ktime
    write(6,"(' wbice(1-ms)',9f7.3)")(wbice(idjd,k),k=1,ms)
  endif
  if(kdate/=kdate_sav.or.ktime/=ktime_sav)then
    write(6,*) 'stopping in indata, not finding correct kdate/ktime'
    write(6,*) "kdate,    ktime     ",kdate,ktime
    write(6,*) "kdate_sav,ktime_sav ",kdate_sav,ktime_sav
    call ccmpi_abort(-1)
  endif

  ! adjust input for differences in orography
  if(newtop==2)then
    ! reduce sea tss to mslp      e.g. for qcca in ncep gcm
    do iq=1,ifull
      if(tss(iq)<0.)then
        if (abs(zss(iq))>1000.) then
          write(6,*)'zss,tss_sea in, out',iq,zss(iq),tss(iq),tss(iq)-zss(iq)*stdlapse/grav
        end if
        tss(iq)=tss(iq)-zss(iq)*stdlapse/grav
      endif
    enddo
  endif                  ! (newtop==2)
  tss(1:ifull)=abs(tss(1:ifull)) ! not done in infile because -ve needed for onthefly
  hourst=.01*ktime
  if ( myid==0 ) then
    write(6,*)'rlongg(1),rlongg(ifull) ',rlongg(1),rlongg(ifull)
    write(6,*)'using em: ',(em(ii),ii=1,10)
    write(6,*)'using  f: ',(f(ii),ii=1,10)
    write(6,*)'in indata hourst = ',hourst
    write(6,*)'sigmas: ',sig
    write(6,*)'sigmh: ',sigmh
  end if
  if ( mydiag ) then
    write(6,*)'newtop, zsold, zs,tss_in,land ',newtop,zss(idjd),zs(idjd),tss(idjd),land(idjd)
  end if
  if (newtop>=1.and..not.lrestart) then    
    if (nproc==1) then
      pslavge=sum(psl(1:ifull)*wts(1:ifull))
      write (6,"(' initial pslavge ',f10.6)") pslavge
    endif 
    do iq=1,ifull
      if (land(iq)) then
        tss(iq)=tss(iq)+(zss(iq)-zs(iq))*stdlapse/grav
        do k=1,ms
          tgg(iq,k)=tgg(iq,k)+(zss(iq)-zs(iq))*stdlapse/grav
        enddo
      endif     ! (land(iq))
    enddo        ! iq loop
    if ( mydiag ) then
      write(6,*)'newtop>=1 new_land_tss,zsold,zs: ',tss(idjd),zss(idjd),zs(idjd)
      ! compensate psl, t(,,1), qg as read in from infile
      write(6,"(' zs#  in     ',9f8.1)") diagvals(zs)
      write(6,"(' zss# in     ',9f8.1)") diagvals(zss)
      write(6,"(' 100*psl#  in',9f8.2)") 100.*diagvals(psl)
      write(6,*) 'now call retopo from indata'
    end if ! ( mydiag )
    call retopo(psl,zss,zs,t,qg)
    if(nmaxpr==1.and.mydiag)then
      write(6,"(' 100*psl# out',9f8.2)") 100.*diagvals(psl)
    endif
    if (nproc==1) then
      pslavge=sum(psl(1:ifull)*wts(1:ifull))
      write (6,"(' after retopo pslavge ',f10.6)") pslavge
    endif 
  endif   ! (newtop>=1.and..not.lrestart)

  qg(1:ifull,:)=max(qg(1:ifull,:),0.)
  ps(1:ifull)=1.e5*exp(psl(1:ifull))
        
else

  ! read in namelist for uin,vin,tbarr etc. for special runs
  if (myid==0) write(6,*)'Read IC from namelist tinit'
  read (99, tin)
  if (myid==0) write(6, tin)

endif   ! (io_in<4)


!-----------------------------------------------------------------
! SPECIAL OPTIONS FOR INITIAL CONDITIONS (nspecial, io_in and nhstest)
if(nspecial>100)then ! increase ps globally by nspecial Pa
  ps(1:ifull)=ps(1:ifull)+nspecial
  psl(:)=log(1.e-5*ps(:))
endif  ! (nspecial>100)       

if (nsib==3) then
  ! put in Antarctica ice-shelf fixes 5/3/07
  do iq=1,ifull
    if(zs(iq)<=0.)then
      rlongd=rlongg(iq)*180./pi
      rlatd=rlatt(iq)*180./pi
      if((rlongd>165..and.rlongd<195..and.rlatd<-77.2-(rlongd-165.)/30.).or.     & ! Ross shelf
         (rlongd>300..and.rlongd<330..and.rlatd<-75.2-(rlongd-300.)*2.8/30.).or. & ! Ronne shelf
         (rlongd>68..and.rlongd<75..and.rlatd<-64.-(rlongd-60.)*5.8/15.))then      ! Amery shelf
           sicedep(iq)=0.
           snowd(iq)=max(snowd(iq),100.)  ! max from Dec 07
           if(mydiag)write(6,*)'setting sea to ice sheet for iq = ',iq
      endif
    endif  ! (zs(iq)<=0.)
  enddo
endif

if(io_in>=5)then
  if (nsib/=0) then
    stop 'ERROR: io_in>=5 requiers nsib=0'
  end if
  ! for rotated coordinate version, see jmcg's notes
  coslong=cos(rlong0*pi/180.)
  sinlong=sin(rlong0*pi/180.)
  coslat=cos(rlat0*pi/180.)
  sinlat=sin(rlat0*pi/180.)
  polenx=-coslat
  poleny=0.
  polenz=sinlat
  write(6,*) 'polenx,poleny,polenz ',polenx,poleny,polenz
  cent=.5*(il_g+1)  ! True center of face
  do k=1,kl
    do iq=1,ifull
      t(iq,k)=tbarr(k)
      qg(iq,k)=qgin(k)
      psl(iq)=.01
    enddo               ! iq loop
    do j=1,jpan
      do i=1,ipan
        ! Need to add offsets to get proper face indices
        do n=1,npan
          rad=sqrt((i+ioff-cent)**2+(j+joff-cent)**2)
          radu=sqrt((i+ioff+.5-cent)**2+(j+joff-cent)**2)
          radv=sqrt((i+ioff-cent)**2+(j+joff+.5-cent)**2)
          iq=indp(i,j,n)
          u(iq,k)=uin*max(1.-radu/(.5*il_g),0.)
          v(iq,k)=vin*max(1.-radv/(.5*il_g),0.)
          if(io_in>=7.and.k==kl)then
            ps(iq)=1.e5*(1.-log(1. + thlapse*zs(iq)/(grav*tsea))  *grav/(cp*thlapse)) **(cp/rdry)
            psl(iq)= log(1.e-5*ps(iq))
          endif
        enddo         ! n loop
      enddo           ! i loop
    enddo             ! j loop
  enddo               ! k loop
endif                 ! io_in>=5

if(io_in==8)then
  ! assign u and v from zonal and meridional uin and vin (no schmidt here)
  ! with zero at poles
  do iq=1,ifull
    psl(iq)=.01
    uzon=uin * abs(cos(rlatt(iq)))
    vmer=vin * abs(cos(rlatt(iq)))
    ! den=sqrt( max(x(iq)**2 + y(iq)**2,1.e-7) )  ! allow for poles
    ! costh=(-y(iq)*ax(iq) + x(iq)*ay(iq))/den
    ! sinth=az(iq)/den
    ! set up unit zonal vector components
    zonx=real(            -polenz*y(iq))
    zony=real(polenz*x(iq)-polenx*z(iq))
    zonz=real(polenx*y(iq)             )
    den=sqrt( max(zonx**2 + zony**2 + zonz**2,1.e-7) ) ! allow for poles
    costh= (zonx*ax(iq)+zony*ay(iq)+zonz*az(iq))/den
    sinth=-(zonx*bx(iq)+zony*by(iq)+zonz*bz(iq))/den
    do k=1,kl
      ! calculate u and v relative to the cc grid,
      u(iq,k)= costh*uzon+sinth*vmer
      v(iq,k)=-sinth*uzon+costh*vmer
    enddo  ! k loop
  enddo      ! iq loop
  ! special option to display panels
  if(uin<0.)then
    do n=1,npan
      do j=1,jpan
        do i=1,ipan
          iq=indp(i,j,n)
          u(iq,:) = n - noff
          t(iq,:) = 0.0001 + n - noff
        enddo
      enddo
    enddo
  endif
endif

! for the held-suarez hs test (also aquaplanet initial)
if (io_in==10) then
  vin=0.
  do k=1,kl
    do iq=1,ifull
      ! set up unit zonal vector components
      zonx=real(            -polenz*y(iq))
      zony=real(polenz*x(iq)-polenx*z(iq))
      zonz=real(polenx*y(iq)             )
      den=sqrt( max(zonx**2 + zony**2 + zonz**2,1.e-7) ) ! allow for poles
      costh= (zonx*ax(iq)+zony*ay(iq)+zonz*az(iq))/den
      sinth=-(zonx*bx(iq)+zony*by(iq)+zonz*bz(iq))/den
      ! set the temperature to the equilibrium zonal mean
      t(iq,k) = max ( 200., (315. - delty*sin(rlatt(iq))**2 - deltheta*log(sig(k))*cos(rlatt(iq))**2)*sig(k)**rkappa )
      ! set zonal wind to an approximate equilibrium
      uin = 125. * sin(2.*rlatt(iq))**2 * sig(k)*(1.-sig(k))/(1. + 10.*(sig(k)-0.25)**2 )
      u(iq,k)= costh*uin+sinth*vin
      v(iq,k)=-sinth*uin+costh*vin
      if(iq==idjd.and.k==nlv.and.mydiag)then
        write(6,*)'indata setting u,v for h-s'
        write(6,*)'iq,k,ax,ay,az',iq,k,ax(iq),ay(iq),az(iq)
        write(6,*)'costh,sinth,x,y,z',costh,sinth,x(iq),y(iq),z(iq)
        write(6,*)'uin,vin,u,v',uin,vin,u(iq,k),v(iq,k)
      endif
      qg(iq,k) = 0.
      ps(iq) = 1.e5
      psl(iq) = .01
      zs(iq) = 0.
    enddo               ! iq loop
  enddo ! k loop
endif   ! (io_in==10) held-suarez test (also aquaplanet initial)

if(io_in==11)then
  ! advection test, once around globe per 10 days
  ! only non-rotated set up so far
  vmer=0.
  ! assign u and v from zonal and meridional winds
  do iq=1,ifull
    den=real(sqrt( max(x(iq)**2 + y(iq)**2,1.e-7_8) )) ! allow for poles
    costh=real(-y(iq)*ax(iq) + x(iq)*ay(iq))/den
    sinth=az(iq)/den
    uzon=2.*pi*rearth/(10.*86400) * abs(cos(rlatt(iq)))
    psl(iq)=.01
    ps(iq)=1.e5*exp(psl(iq))
    f(iq)=0.
    fu(iq)=0.
    fv(iq)=0.
    do k=1,kl
      ! calculate u and v relative to the cc grid,
      ! using components of gaussian grid u and v with theta
      u(iq,k)= costh*uzon+sinth*vmer
      v(iq,k)=-sinth*uzon+costh*vmer
      t(iq,k)=tbarr(k)
      qg(iq,k)=1.e-6
      if(rlongg(iq)>0.and.rlongg(iq)<10.*pi/180.) qg(iq,k)=10.e-3
    enddo  ! k loop
  enddo      ! iq loop
endif  ! io_in==11

if(io_in==27)then
  ! cold bubble test 
  ! constants: 
  Tbubb=300. ! potential temperature
  tss(:)=Tbubb
  do k=1,kl
    ! Height of sigma-levels at t=0
    zbub(k)=(cp*Tbubb/grav)*(1.-sig(k)**(rdry/cp)) !Hydrostatic lapse rate
    ! phi(k)=g*zbub(k)  !geopotential
  enddo
  ! u and v on the cc grid,
  u(:,:)=0.
  v(:,:)=0.
  ps(:) = 1.e5        !initial surface pressure
  psl(:) = .01        !initial lnps
  zs(:)=0.            !zero topography 
  qg(:,:)=1.e-6       !Moisture? Dry atmosphere for bubble  
  f(:)=0.
  fu(:)=0.
  fv(:)=0.
  do k=1,kl
    ! Height of sigma-levels at t=0
    t(:,k)=Tbubb-zbub(k)*grav/cp  !environmental temperature
    if(sig(k)<.4)t(:,k)=t(:,k-1)
  enddo
  ! Inserting the bubble
  ! radius of bubble in the horizontal:
  xt=4000.
  yt=4000.
  zt=2000.
  ic=il_g/2    ! Indices on the global grid
  jc=3*il_g/2
  emcent=em_g(3*il_g*il_g/2)
  if(myid==0) write(6,*) 'emcent, ds/emcent ',emcent,ds/emcent
  xc=ic*ds/emcent
  yc=jc*ds/emcent
  zc=3000. !center in the vertical
  do n=1,npan
    do j=1,jpan
      do i=1,ipan
        iq=indp(i,j,n)          ! Index on this processor
        iqg=indg(i,j,n)         ! Global index
        jg = 1 + (iqg-1)/il_g   ! Global indices
        ig = iqg - (jg-1)*il_g
        xbub=ig*ds/emcent  !gridpoint horizontal distance
        ybub=jg*ds/emcent  !gridpoint horizontal distance
        do k=1,kl
          dist=sqrt(((xbub-xc)/xt)**2+((ybub-yc)/yt)**2+((zbub(k)-zc)/zt)**2)
          if(dist<=1)then
            t(iq,k)=t(iq,k)-15.*((cos(dist*pi/2.))**2)
          endif
        enddo
      enddo
    enddo
  enddo
  call printa('t   ',t,0,nlv,ia,ib,ja,jb,200.,1.)
endif  ! io_in==27, cold bubble test

if(io_in>=4)then   ! i.e. for special test runs without infile
  ! set default tgg etc to level 2 temperatures, if not read in above
  do iq=1,ifull
    tgg(iq,ms) =t(iq,2)   ! just for io_in>=4
    tgg(iq,2) =t(iq,2)    ! just for io_in>=4
    tss(iq)=t(iq,2)
    if(.not.land(iq))then
      tss(iq)=tsea
    endif
  enddo   ! iq loop
endif

! aquaplanet (APE) test
if(nhstest<0)then  
  fracice(:)=0.
  sicedep(:)=0.
  snowd(:)=0.
  tss(:)=273.16
  do iq=1,ifull
    if((nhstest==-1.or.nhstest<=-6).and.abs(rlatt(iq))<pi/3.) tss(iq)=273.16 +27.*(1.-sin(1.5*rlatt(iq))**2)   ! Expt 1 control
    if(nhstest==-2.and.abs(rlatt(iq))<pi/3.) tss(iq)=273.16 +27.*(1.-3.*abs(rlatt(iq))/pi)                     ! Expt 2 peaked
    if(nhstest==-3.and.abs(rlatt(iq))<pi/3.) tss(iq)=273.16 +27.*(1.-sin(1.5*rlatt(iq))**4)                    ! Expt 3 flat
    if(nhstest==-4.and.abs(rlatt(iq))<pi/3.) tss(iq)=273.16 +13.5*(1.-sin(1.5*rlatt(iq))**2)                 & ! Expt 4 qobs
                                                            +13.5*(1.-sin(1.5*rlatt(iq))**4)                   ! Expt 4 control5n
    if(nhstest==-5.and.rlatt(iq)>pi/36..and.rlatt(iq)<pi/3.)                                                 &
                                                   tss(iq)=273.16 +27.*(1.-sin((rlatt(iq)-pi/36.)*90./55.)**2) ! Expt 5
    if(nhstest==-5.and.rlatt(iq)>-pi/3..and.rlatt(iq)<=pi/36.)                                               &
                                                   tss(iq)=273.16 +27.*(1.-sin((rlatt(iq)-pi/36.)*90./65.)**2) ! Expt 5
    if(nhstest==-6.and.abs(rlatt(iq))<pi/12.)then                                                              ! Expt6 1keq
      if(rlongg(iq)<pi/6..or.rlongg(iq)>11.*pi/6.) tss(iq)=tss(iq)+1.*(cos(3.*rlongg(iq))*cos(6.*rlatt(iq)))**2 
    endif  ! (nhstest==-6....)                      
    if(nhstest==-7.and.abs(rlatt(iq))<pi/12.)then                                                              ! Expt7 3keq
      if(rlongg(iq)<pi/6..or.rlongg(iq)>11.*pi/6.) tss(iq)=tss(iq)+3.*(cos(3.*rlongg(iq))*cos(6.*rlatt(iq)))**2  
    endif  ! (nhstest==-7....)
    if(nhstest==-8.and.abs(rlatt(iq))<pi/6.) tss(iq)=tss(iq)+3.*cos(rlongg(iq))*cos(3.*rlatt(iq))**2           ! Expt 8  3kw1
  enddo   ! iq loop
  do k=1,ms
    tgg(:,k)=tss(:)
    wb(:,k)=0.
  enddo
  if(io_in>4)then    ! not reading initial input file
    do k=1,kl
      do iq=1,ifull
        qg(iq,k)=.01*(1.-abs(rlatt(iq))*2./pi)**3*sig(k)**3  ! Nov 2004
      enddo
    enddo
  endif  ! (io_in>4)
  if(nhstest==-13)then ! a jlm test
    f(:)=0.
    u(:,:)=0.
    v(:,:)=0.
    t(:,:)=300.
    do iq=1,ifull
      if(rlatt(iq)*180./pi>10.and.rlatt(iq)*180./pi<13.)t(iq,1)=305.
    end do
  endif
  if(io_in==5)then ! for chess-board colouring of grid, via tss
                   ! but get some blurring from cc2hist interps
    qg(:,:)=0.
    ps(:)=1.e5
    psl(:)=.01
    zs(:)=0.
    uin=300.
    if(nhstest==-2)uin=300.5
    vin=1.
    vin=3.
    if(nhstest==-3)vin=0. ! just for 6 colours of panels
    if(myid==0)then
      allocate(davt_g(ifull_g))
      do n=0,5
        do j=1,il_g
          do i=1,il_g
            davt_g(indglobal(i,j,n))=uin+mod(n,3)+vin*mod(mod(i,2)+mod(j,2),2)
          enddo ! i loop
        enddo  ! j loop
      enddo   ! n loop
      call ccmpi_distribute(tss,davt_g)
      deallocate(davt_g)
    else
      call ccmpi_distribute(tss)
    endif ! myid==0
  endif   ! io_in==5
endif    ! (nhstest<0)


!--------------------------------------------------------------
! DEFINE SURFACE DATA PRESETS (nrungcm)

!     nrungcm<0 controls presets for snowd, wb, tgg and other soil variables
!     they can be: preset/read_in_from_previous_run
!                  written_out/not_written_out    after 24 h as follows:
!          nrungcm = -1  preset           | not written to separate file
!                    -2  preset           |     written  
!                    -3  read_in          |     written  (usual for NWP)
!                    -4  read_in          | not written  (usual for netCDF input)
!                    -5  read_in (not wb) |     written  (should be good)
!                    -6 same as -1 bit tapered wb over dry interio of Aust
!                    >5 like -1 but sets most wb percentages

! preset soil data
if ( .not.lrestart ) then
  if((nrungcm.le.-1.and.nrungcm.ge.-2).or.nrungcm==-6.or.nrungcm>5)then
    ! presetting wb when no soil moisture available initially
    iyr=kdate/10000
    imo=(kdate-10000*iyr)/100
    do iq=1,ifull
      if(land(iq))then
        iveg=ivegt(iq)
        if (nsib==6.or.nsib==7) iveg=1
        isoil=isoilm(iq)
        rlonx=rlongg(iq)*180./pi
        rlatx=rlatt(iq)*180./pi
        ! fracsum(imo) is .5 for nh summer value, -.5 for nh winter value
        fracs=sign(1.,rlatt(iq))*fracsum(imo)  ! +ve for local summer
        if(nrungcm>5)then
          fracwet=.01*nrungcm   ! e.g. 50 gives .5
        else
          fracwet=(.5+fracs)*fracwets(iveg)+(.5-fracs)*fracwetw(iveg)
          ! N.B. for all Dean's points, fracwet=fracwets=fracwetw=.5           
        endif
        wb(iq,ms)= (1.-fracwet)*swilt(isoilm(iq))+fracwet*sfc(isoilm(iq)) 
        if(abs(rlatx)<18.)wb(iq,ms)=sfc(isoilm(iq)) ! tropics
        ! following jlm subtropics from Aug 2003 (.1/.9), (.6, .4)
        if(rlatx<20..and.rlatx>8.) then
          wb(iq,ms)=(.35-.5*fracsum(imo))*swilt(isoilm(iq))+(.65+.5*fracsum(imo))*sfc(isoilm(iq)) ! NH
        end if
        if(rlatx>-16..and.rlatx<-8.) then
          wb(iq,ms)=(.35+.5*fracsum(imo))*swilt(isoilm(iq))+(.65-.5*fracsum(imo))*sfc(isoilm(iq)) ! SH
        end if
        if(rlatx>-32..and.rlatx<-22..and.rlonx>117..and.rlonx<146.) then
          if(nrungcm==-6)then
            ! following tapers it over 2 degrees lat/long
            alf=.5*min(abs(rlonx-117.),abs(rlonx-146.),abs(rlatx+22.),abs(rlatx+32.),2.)
            wb(iq,ms)=alf*swilt(isoilm(iq))+(1.-alf)*wb(iq,ms)
          else
            wb(iq,ms)=swilt(isoilm(iq)) ! dry interior of Australia
          endif
        endif
      endif    ! (land(iq))
    enddo     ! iq loop
    do k=1,ms-1
      wb(:,k)=wb(:,ms)
    enddo    !  k loop

    do iq=1,ifull
      ! fix for antarctic snow
      if(land(iq).and.rlatt(iq)*180./pi<-60.)snowd(iq)=max(snowd(iq),400.)
    enddo   ! iq loop

    if ( mydiag ) then
      iveg=ivegt(idjd)
      if (nsib==6.or.nsib==7) iveg=1
      isoil=isoilm(idjd)
      if (isoil>0) then
        write(6,*)'isoil,iveg,month,fracsum,rlatt: ',isoil,iveg,imo,fracsum(imo),rlatt(idjd)
        fracs=sign(1.,rlatt(idjd))*fracsum(imo) ! +ve for local summer
        fracwet=(.5+fracs)*fracwets(iveg)+(.5-fracs)*fracwetw(iveg)
        write(6,*)'fracs,fracwet,initial_wb: ',fracs,fracwet,wb(idjd,ms)
      end if
    end if
  endif       !  ((nrungcm==-1.or.nrungcm==-2.or.nrungcm==-5)

  ! Soil recycling input
  if( nrungcm<=-3 .and. nrungcm>=-5 ) then
    call histopen(ncid,surf_00,ier)
    if ( ier==0 ) then
      ! NETCDF file format
      ! clobber ifile surface data with surfin surface data
      kdate_s = kdate_sav
      ktime_s = ktime_sav
      if ( myid==0 ) then
        write(6,*) 'Replacing surface data with input from ',trim(surf_00)
      end if
      call onthefly(2,kdate,ktime,duma(:,1),duma(:,2),duma(:,3),duma(:,4),    &
                    duma(:,5),dumb(:,:,1),dumb(:,:,2),dumb(:,:,3),            &
                    dumb(:,:,4),tgg,wb,wbice,snowd,dumb(:,:,5),dumb(:,:,6),   &
                    dumb(:,:,7),dumb(:,:,8),dumb(:,:,9),tggsn,smass,ssdn,     &
                    ssdnn,snage,isflag,mlodwn,ocndwn,xtgdwn)
      call histclose
      if ( kdate/=kdate_sav .or. ktime/=ktime_sav ) then
        if ( myid==0 ) then
          write(6,*) 'WARN: Could not locate correct date/time'
          write(6,*) '      Using infile surface data instead'
        end if
        kdate = kdate_sav
        ktime = ktime_sav
      endif
    else
      ! ASCII file format
      ! for sequence of runs starting with values saved from last run
      if(ktime==1200)then
        co2in=co2_12      ! 'co2.12'
        radonin=radon_12  ! 'radon.12'
        surfin=surf_12    ! 'surf.12'
      else
        co2in=co2_00      !  'co2.00'
        radonin=radon_00  ! 'radon.00'
        surfin=surf_00    ! 'surf.00'
      endif
      if ( myid == 0 ) then
        write(6,*) 'reading previously saved wb,tgg,tss (land),snowd,sice from ',surfin
        open(87,file=surfin,form='formatted',status='old')
        read(87,'(a80)') header
        write(6,*)'header: ',header
      end if
      if(nrungcm==-5)then
        call readglobvar(87, tgg, fmt="*") ! this acts as dummy read 
      else
        call readglobvar(87, wb, fmt="*")
      endif
      call readglobvar(87, tgg, fmt="*")
      call readglobvar(87, aa, fmt="*")    ! only use land values of tss
      call readglobvar(87, snowd, fmt="*")
      if ( myid == 0 ) close(87)
      where (land(1:ifull))
        tss(1:ifull)=aa(1:ifull)
      end where
    end if  ! ier==0
  endif    !  (nrungcm<=-3.and.nrungcm>=-5)

  if(nrungcm==4)then !  wb fix for ncep input 
    ! this is related to eak's temporary fix for soil moisture
    ! - to compensate for ncep drying out, increase minimum value
    do k=1,ms
      do iq=1,ifull     
        isoil=isoilm(iq)
        wb(iq,k)=min( sfc(isoil) , max(.75*swilt(isoil)+.25*sfc(isoil),wb(iq,k)) )
        ! for antarctic snow
        if(land(iq).and.rlatt(iq)*180./pi<-60.)snowd(iq)=max(snowd(iq),400.)
      enddo   ! iq loop
    enddo    !  k loop
  endif      !  (nrungcm==4)

  if(nrungcm==5)then !  tgg, wb fix for mark 3 input
    ! unfortunately mk 3 only writes out 2 levels
    ! wb just saved as excess above wilting; top level & integrated values
    ! tgg saved for levels 2 and ms 
    do iq=1,ifull     
      isoil=isoilm(iq)
      do k=2,3
        wb(iq,k)=wb(iq,ms)
      enddo    !  k loop
      do k=1,ms
        ! wb(iq,k)=min( sfc(isoil) ,wb(iq,k)+swilt(isoil) ) ! till 22/8/02
        wb(iq,k)=wb(iq,k)+swilt(isoil) 
      enddo    !  k loop
      tgg(iq,3)=.75*tgg(iq,2)+.25*tgg(iq,6)
      tgg(iq,4)= .5*tgg(iq,2)+ .5*tgg(iq,6)
      tgg(iq,5)=.25*tgg(iq,2)+.75*tgg(iq,6)
      ! fix for antarctic snow
      if(land(iq).and.rlatt(iq)*180./pi<-60.)snowd(iq)=max(snowd(iq),400.)
    enddo   ! iq loop
    if (mydiag) then
      write(6,*) 'after nrungcm=5 fixup of mk3 soil variables:'
      write(6,*) 'tgg ',(tgg(idjd,k),k=1,ms)
      write(6,*) 'wb ',(wb(idjd,k),k=1,ms)
    end if
  endif      !  (nrungcm==5)

  if(nrungcm==1)then  ! jlm alternative wb fix for nsib runs off early mark 2 gcm
    if (mydiag ) then
      isoil = isoilm(idjd)
      write(6,"('before nrungcm=1 fix-up wb(1-ms)',9f7.3)") (wb(idjd,k),k=1,ms)
      write(6,*)'nfixwb,isoil,swilt,sfc,ssat,alb ',nfixwb,isoil,swilt(isoil),sfc(isoil),ssat(isoil),albvisnir(idjd,1)
    end if
    do ip=1,ipland  ! all land points in this nsib=1+ loop
      iq=iperm(ip)
      isoil = isoilm(iq)

      if(nfixwb==0)then
        ! very dry jlm suggestion. assume vegfrac ~.5, so try to satisfy
        ! wb0/.36=.5*(wb/sfc + (wb-swilt)/(sfc-swilt) )
        wb(iq,1)=( sfc(isoil)*(sfc(isoil)-swilt(isoil))*wb(iq,1)/.36 +.5*sfc(isoil)*swilt(isoil) )/(sfc(isoil)-.5*swilt(isoil))
        do k=2,ms
          wb(iq,k)=( sfc(isoil)*(sfc(isoil)-swilt(isoil))*wb(iq,ms)/.36+.5*sfc(isoil)*swilt(isoil) )/(sfc(isoil)-.5*swilt(isoil))
        enddo   !  k=2,ms
      endif   ! (nfixwb==0)
      if(nfixwb==1.or.nfixwb==2)then
        ! alternative simpler jlm fix-up	
        ! wb0/.36=(wb-swilt)/(sfc-swilt)
        wb(iq,1)=swilt(isoil)+(sfc(isoil)-swilt(isoil))*wb(iq,1)/.36
        do k=2,ms
          wb(iq,k)=swilt(isoil)+(sfc(isoil)-swilt(isoil))*wb(iq,ms)/.36
        enddo   !  k=2,ms
      endif   ! (nfixwb==1.or.nfixwb==2)
      if(ip==1)write(6,*)'kdate ',kdate
      if(nfixwb==2.and.kdate>3210100.and.kdate<3210200)then
        rlon_d=rlongg(iq)*180./pi
        rlat_d=rlatt(iq)*180./pi
        if(ip==1)then
          write(6,*)'kdate in nfixwb=2 ',kdate
          write(6,*)'iq,rlon_d,rlat_d ',rlon_d,rlat_d
        endif
        ! jlm fix-up for tropical oz in january 321
        if(rlon_d>130..and.rlon_d<150..and.rlat_d>-20..and.rlat_d<0.)then
          do k=1,ms
            wb(iq,k)=max(wb(iq,k),.5*(swilt(isoil)+sfc(isoil))) ! tropics
          enddo   !  k=1,ms
        endif
        ! jlm fix-up for dry interior in january 321
        if(rlon_d>117..and.rlon_d<142..and.rlat_d>-32..and.rlat_d<-22.)then
          do k=1,ms
            wb(iq,k)=swilt(isoil)  ! dry interior
          enddo   !  k=1,ms
        endif
      endif   ! (nfixwb==2)
      if(nfixwb==10)then    ! was default for nrungcm=1 till oct 2001
        ! jlm suggestion, assume vegfrac ~.5, so try to satisfy
        ! wb0/.36=.5*(wb/ssat + (wb-swilt)/(ssat-swilt) )
        wb(iq,1)=( ssat(isoil)*(ssat(isoil)-swilt(isoil))*wb(iq,1)/.36+.5*ssat(isoil)*swilt(isoil) )    &
                /(ssat(isoil)-.5*swilt(isoil))
        do k=2,ms                                                       
          wb(iq,k)=( ssat(isoil)*(ssat(isoil)-swilt(isoil))*wb(iq,ms)/.36+.5*ssat(isoil)*swilt(isoil) ) &
                  /(ssat(isoil)-.5*swilt(isoil))
        enddo   !  k=2,ms
      endif   ! (nfixwb.ne.10)

      do k=1,ms
        wb(iq,k)=max( swilt(isoil) , min(wb(iq,k),sfc(isoil)) )
      enddo     !  k=1,ms
    enddo        !  ip=1,ipland
    if (mydiag) then
      write(6,"('after nrungcm=1 fix-up wb(1-ms)',9f7.3)") (wb(idjd,k),k=1,ms)
      write(6,"('wbice(1-ms)',9f7.3)")(wbice(idjd,k),k=1,ms)
      write(6,"('wb3frac#',9f8.2)") (diagvals(wb(:,3)) - swilt(diagvals(isoilm)))     &
                                  / (sfc(diagvals(isoilm)) - swilt(diagvals(isoilm)))
    end if
  endif          !  (nrungcm==1)

  if(nrungcm==2)then  ! for nsib runs off early mark 2 gcm
    if (mydiag) then
      isoil = isoilm(idjd)
      write(6,*)'before nrungcm=2 fix-up wb(1-ms): ',wb(idjd,:)
      write(6,*)'isoil,swilt,ssat,alb ',isoil,swilt(isoil),ssat(isoil),albvisnir(idjd,1)
    end if
    do ip=1,ipland  ! all land points in this nsib=1+ loop
      iq=iperm(ip)
      isoil = isoilm(iq)
      if( albvisnir(iq,1) >= 0.25 ) then
        diffg=max(0. , wb(iq,1)-0.068)*ssat(isoil)/0.395   ! for sib3
        diffb=max(0. , wb(iq,ms)-0.068)*ssat(isoil)/0.395  ! for sib3
      else
        diffg=max(0. , wb(iq,1)-0.175)*ssat(isoil)/0.42    ! for sib3
        diffb=max(0. , wb(iq,ms)-0.175)*ssat(isoil)/0.42   ! for sib3
      endif
      wb(iq,1)=swilt(isoil)+diffg          ! for sib3
      do k=2,ms                            ! for sib3
        wb(iq,k)=swilt(isoil)+diffb         ! for sib3
      enddo     !  k=2,ms
    enddo      !  ip=1,ipland
    if(mydiag) write(6,*)'after nrungcm=2 fix-up wb(1-ms): ',wb(idjd,:)
  endif          !  (nrungcm==2)
end if ! ( .not.lrestart )


!--------------------------------------------------------------
! SPECIAL OPTIONS FOR TEMPERATURES (nspecial)

if(nspecial==34)then      ! test for Andy Pitman & Faye
  tgg(1:ifull,6)=tgg(1:ifull,6)+.1
endif
! for CAI experiment
if (nspecial==42) then
  call caispecial
endif 
if (nspecial==43) then
  call caispecial
  do iq=1,ifull
    rlongd=rlongg(iq)*180./pi
    rlatd=rlatt(iq)*180./pi
    if (rlatd.ge.-6..and.rlatd.le.6.) then
      if (rlongd.ge.180..and.rlongd.le.290.) then
        if (.not.land(iq)) then
          tgg(iq,1)=293.16
          tss(iq)=293.16
        end if
      end if
    end if
  end do
end if


!--------------------------------------------------------------
! SET-UP AMIP SSTs (namip)
if ( namip/=0 ) then
  if ( myid==0 ) then
    write(6,*) 'calling amipsst at beginning of run'
  end if
  call amipsst
endif   ! namip/=0


!--------------------------------------------------------------
! FINAL FIXES FOR SURFACE DATA

! soil moisture fixes
do iq=1,ifull
  if(.not.land(iq))then
    wb(iq,:)=0.   ! default over ocean (for plotting)
  endif    !  (.not.land(iq))
enddo     ! iq loop
      
! sea-ice fixes
do iq=1,ifull
  if(fracice(iq)<.02) then
    fracice(iq)=0.
    sicedep(iq)=0.
  end if
  if(land(iq))then
    sicedep(iq)=0.
    fracice(iq)=0.
  else
    if(fracice(iq)>0..and.sicedep(iq)==0.)then
      ! assign to 2. in NH and 1. in SH (according to spo)
      ! do this here and in nestin because of onthefly
      if(rlatt(iq)>0.)then
        sicedep(iq)=2.
      else
        sicedep(iq)=1.
      endif ! (rlatt(iq)>0.)
    elseif(fracice(iq)==0..and.sicedep(iq)>0.)then  ! e.g. from Mk3  
      fracice(iq)=1.
    endif  ! (fracice(iq)>0..and.sicedep(iq)==0.) .. elseif ..
  endif    ! (land(iq))
enddo     ! iq loop

! snow and ice fixes
snalb=.8
do iq=1,ifull
  if (nmlo==0.or.abs(nmlo)>9) then
    if(.not.land(iq))then
      ! from June '03 tgg1	holds actual sea temp, tss holds net temp 
      tgg(iq,1)=max(271.3,tss(iq)) 
      tggsn(iq,1)=tss(iq)         ! a default
    endif   ! (.not.land(iq))
    if(sicedep(iq)>0.)then
      ! at beginning of month set sice temperatures
      tggsn(iq,1)=min(271.2,tss(iq),t(iq,1)+.04*6.5) ! for 40 m level 1
      tss(iq)=tggsn(iq,1)*fracice(iq)+tgg(iq,1)*(1.-fracice(iq))
      albvisnir(iq,1)=.8*fracice(iq)+.1*(1.-fracice(iq))
      albvisnir(iq,2)=.5*fracice(iq)+.1*(1.-fracice(iq))
    endif   ! (sicedep(iq)>0.)
  endif    ! (nmlo==0.or.abs(nmlo)>9) 
  if(isoilm(iq)==9.and.(nsib==3.or.nsib==5))then
    ! also at beg. of month ensure cold deep temps over permanent ice
    do k=2,ms
      tgg(iq,k)=min(tgg(iq,k),273.1) ! was 260
      wb(iq,k)=max(wb(iq,k),sfc(9))  ! restart value may exceed sfc
      if(wbice(iq,k)<=0.)wbice(iq,k)=.8*wb(iq,k)  ! Dec 07       
    enddo
  endif   ! (isoilm(iq)==9)
enddo    ! iq loop
tpan(1:ifull)=t(1:ifull,1) ! default for land_sflux and outcdf

! albedo and roughness fixes
select case(nsib)
  case(5)
    ! MODIS input with standard surface scheme
    osnowd = snowd
    zolog=log(zmin/zolnd)   ! for land use in sflux
    sigmf=0.
    where (land)
      sigmf(:)=max(0.01,min(0.98,1.-exp(-0.4*vlai(:))))
    elsewhere
      vlai=0.
    end where
  case(3)
    ! usual input with standard surface scheme
    osnowd = snowd
    sigmf=0.
    do iq=1,ifull
      if(land(iq))then
        isoil = isoilm(iq)
        iveg  = ivegt(iq)
        sigmf(iq)=min(.8,.95*vegpsig(ivegt(iq)))  ! moved from rdnsib
        if(jlmsigmf==1)then  ! fix-up for dean's veg-fraction
          sigmf(iq)=((sfc(isoil)-wb(iq,3))*vegpmin(iveg)+(wb(iq,3)-swilt(isoil))*vegpmax(iveg))/(sfc(isoil)-swilt(isoil)) 
          sigmf(iq)=max(vegpmin(iveg),min(sigmf(iq),.8)) ! for odd wb
        endif   ! (jlmsigmf==1)
        ! following just for rsmin diags for nstn and outcdf	
        tstom=298.
        if(iveg==6+31)tstom=302.
        if(iveg>=10.and.iveg<=21.and.abs(rlatt(iq)*180./pi)<25.)tstom=302.
        tsoil=min(tstom, .5*(.3333*tgg(iq,2)+.6667*tgg(iq,3)+.95*tgg(iq,4) + .05*tgg(iq,5)))
        ftsoil=max(0.,1.-.0016*(tstom-tsoil)**2)
        ! which is same as:  ftsoil=max(0.,1.-.0016*(tstom-tsoil)**2)
        !                    if( tsoil >= tstom ) ftsoil=1.
        rlai(iq)=  max(.1,rlaim44(iveg)-slveg44(iveg)*(1.-ftsoil))
        rsmin(iq) = rsunc44(iveg)/rlai(iq)
      endif   ! (land(iq)) 
    enddo    !  iq loop
    if(jalbfix==1)then ! jlm fix for albedos, esp. for bare sandy soil
      if(mydiag)then
        isoil=isoilm(idjd)
        if (isoil.gt.0) then
          write(6,*)'before jalbfix isoil,sand,alb,rsmin ',isoil,sand(isoil),albvisnir(idjd,1),rsmin(idjd)
        else
          write(6,*)'before jalbfix isoil,sand,alb,rsmin ',isoil,0.,albvisnir(idjd,1),rsmin(idjd)
        end if
      endif
      do ip=1,ipland  
        iq=iperm(ip)
        isoil = isoilm(iq)
        albvisnir(iq,1)=max(albvisnir(iq,1),sigmf(iq)*albvisnir(iq,1)+(1.-sigmf(iq))*(sand(isoil)*.35+(1.-sand(isoil))*.06))
        albvisnir(iq,2)=albvisnir(iq,1)
      enddo                  !  ip=1,ipland
      if(mydiag)then
        write(6,*)'after jalbfix sigmf,alb ',sigmf(idjd),albvisnir(idjd,1)
      endif
    endif  ! (jalbfix==1)
    if(newrough>0)then
      call calczo
      if(mydiag)write(6,*)'after calczo zolnd ',zolnd(idjd)
      if ( mydiag ) then
        write(6,*)'after calczo with newrough = ',newrough
        write(6,"('zo#    ',9f8.2)") diagvals(zolnd)
      end if
    endif ! (newrough>0)
    zolog=log(zmin/zolnd)   ! for land use in sflux 
end select

if ( mydiag ) then
  write(6,*)'near end of indata id+-1, jd+-1'
  write(6,"(' tss#    ',9f8.2)") diagvals(tss)
  write(6,"(' tgg(1)# ',9f8.2)") diagvals(tgg(:,1))
  write(6,"(' tgg(2)# ',9f8.2)") diagvals(tgg(:,2))
  write(6,"(' tgg(3)# ',9f8.2)") diagvals(tgg(:,3))
  write(6,"(' tgg(ms)#',9f8.2)") diagvals(tgg(:,ms))
  write(6,"(' land#   ',9l8)")  diagvals(land)
  write(6,"(' sicedep#   ',9f8.2)") diagvals(sicedep)
  write(6,*)'following from rdnsib'
  write(6,"(' zo#     ',9f8.2)") diagvals(zolnd)
  write(6,"(' wb(1)#  ',9f8.3)") diagvals(wb(:,1))
  wb_tmpry(1:ms) = wb(idjd,1:ms)
  write(6,*)' wb(1-ms): ',wb_tmpry(1:ms)
  write(6,"(' wb(ms)# ',9f8.3)") diagvals(wb(:,ms))
  swilt_diag(1:9) = swilt(diagvals(isoilm))
  sfc_diag(1:9) = sfc(diagvals(isoilm))
  write(6,"(' swilt#  ',9f8.3)") swilt_diag(:)
  write(6,"(' wb3frac#',9f8.3)") (diagvals(wb(:,3)) - swilt_diag(:)) / (sfc_diag(:) - swilt_diag(:))
  write(6,"(' snowd#  ',9f8.2)") diagvals(snowd)
  write(6,"(' fracice#',9f8.3)") diagvals(fracice)
end if

! general initial checks for wb and wbice
do k=1,ms
  do iq=1,ifull
    isoil=isoilm(iq)
    wb(iq,k)=min(ssat(isoil),wb(iq,k))
    wbice(iq,k)=min(.99*wb(iq,k),wbice(iq,k)) 
    if(isoil/=9.and.wbice(iq,k)<=0.)wbice(iq,k)=min(.99,max(0.,.99*(273.1-tgg(iq,k))/5.))*wb(iq,k) ! jlm
  enddo  ! iq loop
enddo   ! ms

if ( mydiag ) then
  write(6,*)'nearer end of indata id+-1, jd+-1'
  write(6,"('tgg(2)# ',9f8.2)")  diagvals(tgg(:,2))
  write(6,"('tgg(ms)#',9f8.2)")  diagvals(tgg(:,ms))
  write(6,"('wbice(1-ms)',9f7.3)")(wbice(idjd,k),k=1,ms)
end if



!**************************************************************
!**************************************************************
! No changes to input data allowed after this point
!**************************************************************
!**************************************************************


!-----------------------------------------------------------------
! UPDATE GENERAL MODEL VARIABLES

! orography
call bounds(zs,corner=.true.)

if ( mydiag ) then
  write(6,*)'for idjd get = ',idjd,ie(idjd),iw(idjd),in(idjd),is(idjd)
  write(6,*)'with zs: ',zs(idjd),zs(ie(idjd)),zs(iw(idjd)),zs(in(idjd)),zs(is(idjd))
end if

! off-centre terms
if(abs(epsp)>1.)then   ! e.g. 20. gives epsmax=.2 for del_orog=600 m
  epsmax=abs(epsp)/100.
  do iq=1,ifull      ! sliding epsf to epsmax
    zsdiff=max(abs(zs(ie(iq))-zs(iq)),abs(zs(iw(iq))-zs(iq)),abs(zs(in(iq))-zs(iq)),abs(zs(is(iq))-zs(iq)) )
    epst(iq)=max(epsf,min(epsmax*zsdiff/(600.*grav),epsmax))
  enddo
  epsf=0.
else
  epst(1:ifull)=abs(epsp)
endif  ! (abs(epsp)>1.)
if(abs(epsp)>99.)then  ! e.g. 200. to give epsmax=.2 for orog=600 m
  epsmax=abs(epsp)/1000.
  do iq=1,ifull
    zsdiff=max(zs(iq)-zs(ie(iq)),zs(iq)-zs(iw(iq)),zs(iq)-zs(in(iq)),zs(iq)-zs(is(iq)),0. )
    epst(iq)=min(epsmax*zsdiff/(600.*grav),epsmax) ! sliding 0. to epsmax
  enddo
endif  ! (abs(epsp)>99.)
if(epsp>1..and.epsp<2.)epst(:)=epsp-1.

! saved albedo
albvissav(:)=albvisnir(:,1) ! VIS
albnirsav(:)=albvisnir(:,2) ! NIR
      
! surface pressure
ps(1:ifull)=1.e5*exp(psl(1:ifull))

!--------------------------------------------------------------
! UPDATE DAVIES NUDGING ARRAYS (nbd and nud_hrs)
! Must occur after defining initial atmosphere fields
if(nbd/=0.and.nud_hrs/=0)then
  call davset   ! as entry in subr. davies, sets psls,qgg,tt,uu,vv
  if ( myid==0 ) then
    allocate(davt_g(ifull_g))
    ! Set up the weights using global array and indexing
    ! This needs the global function indglobal for calculating the 1D index
    davt_g(:) = 0.
    if(nbd==1)then
      davt_g(:) = 1./nud_hrs !  e.g. 1/48
    endif                !  (nbd>0)
    if(nbd==-1)then    ! linearly increasing nudging, just on panel 4
      centi=.5*(il_g+1)
      do j=1,il_g
        do i=1,il_g
          dist=max(abs(i-centi),abs(j-centi)) ! dist from centre of panel
          distx=dist/(.5*il_g) ! between 0. and 1.
          davt_g(indglobal(j,i,4))=(1.-distx)/abs(nud_hrs) !  e.g. 1/24
        enddo            ! i loop
      enddo              ! j loop
    endif                !  (nbd==-1) 
    if(nbd==-2)then    ! quadr. increasing nudging, just on panel 4
      centi=.5*(il_g+1)
      do j=1,il_g  
        do i=1,il_g
          dist=max(abs(i-centi),abs(j-centi)) ! dist from centre of panel
          distx=dist/(.5*il_g) ! between 0. and 1.
          davt_g(indglobal(j,i,4))=(1.-distx**2)/abs(nud_hrs) !  e.g. 1/24
        enddo            ! i loop
      enddo              ! j loop
    endif                !  (nbd==-2) 
    if(abs(nbd)==3)then !usual far-field with no nudging on panel 1
      do n=0,5
        do j=il_g/2+1,il_g
          ! linearly between 0 (at il/2) and 1/abs(nud_hrs) (at il+1)
          rhs=(j-il_g/2)/((il_g/2+1.)*nud_hrs)
          do i=1,il_g
            if(n==0)davt_g(indglobal(i,il_g+1-j,n))=rhs
            if(n==2)davt_g(indglobal(j,i,n))=rhs
            if(n==3)davt_g(indglobal(j,i,n))=rhs
            if(n==5)davt_g(indglobal(i,il_g+1-j,n))=rhs
          enddo          ! i loop
        enddo            ! j loop
      enddo              ! n loop
      do j=1,il_g          ! full nudging on furthest panel
        do i=1,il_g
          davt_g(indglobal(j,i,4))=1./nud_hrs !  e.g. 1/48
        enddo            ! i loop
      enddo              ! j loop
    endif                !  (nbd==-3) 
    if(abs(nbd)==4)then    ! another special form with no nudging on panel 1
      do n=0,5
        do j=1,il_g
          ! linearly between 0 (at j=.5) and 1/nud_hrs (at j=il+.5)
          rhs=(j-.5)/(il_g*nud_hrs)
          do i=1,il_g
            if(n==0)davt_g(indglobal(i,il_g+1-j,n))=rhs
            if(n==2)davt_g(indglobal(j,i,n))=rhs
            if(n==3)davt_g(indglobal(j,i,n))=rhs
            if(n==5)davt_g(indglobal(i,il_g+1-j,n))=rhs
          enddo          ! i loop
        enddo            ! j loop
      enddo              ! n loop
      do j=1,il_g        ! full nudging on furthest panel
        do i=1,il_g
          davt_g(indglobal(j,i,4))=1./nud_hrs !  e.g. 1/48
        enddo            ! i loop
      enddo              ! j loop
    endif                !  (nbd==-4) 
    if(abs(nbd)==5)then    ! another special form with some nudging on panel 1
      do n=0,5
        do j=il_g/2+1,il_g
          ! linearly between 0 (at j=.5) and 1/nud_hrs (at j=il+.5)
          rhs=(.5*il_g+j-.5)/(1.5*il_g*nud_hrs)
          do i=1,il_g
            if(n==0)davt_g(indglobal(i,il_g+1-j,n))=rhs
            if(n==2)davt_g(indglobal(j,i,n))=rhs
            if(n==3)davt_g(indglobal(j,i,n))=rhs
            if(n==5)davt_g(indglobal(i,il_g+1-j,n))=rhs
          enddo          ! i loop
        enddo            ! j loop
      enddo              ! n loop
      ril2=il_g/2
      do j=1,jl
        do i=1,il
          rhs=max(abs(i-.5-ril2),abs(j-.5-ril2))/(1.5*il_g*nud_hrs)
          davt_g(indglobal(i,il_g+1-j,1))=rhs  ! panel 1
        enddo
      enddo
      do j=1,il_g        ! full nudging on furthest panel
        do i=1,il_g
          davt_g(indglobal(j,i,4))=1./nud_hrs !  e.g. 1/48
        enddo            ! i loop
      enddo              ! j loop
    endif                !  (nbd==-5)
    if(abs(nbd)==6)then ! more like 1-way nesting
      do j=1,il_g   ! full nudging on all further panels; 6 rows on 1
        do i=1,il_g
          davt_g(indglobal(j,i,0))=1./nud_hrs !  e.g. 1/48
          davt_g(indglobal(j,i,2))=1./nud_hrs !  e.g. 1/48
          davt_g(indglobal(j,i,3))=1./nud_hrs !  e.g. 1/48
          davt_g(indglobal(j,i,4))=1./nud_hrs !  e.g. 1/48
          davt_g(indglobal(j,i,5))=1./nud_hrs !  e.g. 1/48
          davt_g(indglobal(j,i,1))=0.
        enddo            ! i loop
      enddo              ! j loop
      do j=0,5
        ! linearly between 0 and 1/abs(nud_hrs) over 6 rows
        rhs=(6-j)/(6.*nud_hrs)
        do i=1+j,il_g-j
          davt_g(indglobal(i,j+1,1))=rhs
          davt_g(indglobal(i,il_g-j,1))=rhs
          davt_g(indglobal(j+1,i,1))=rhs
          davt_g(indglobal(il_g-j,i,1))=rhs
        enddo         ! i loop
      enddo           ! j loop
    endif             !  (nbd==-6)
    if(abs(nbd)==7)then    ! another special form with no nudging on panel 1
      do n=0,5
        do j=1,il_g
          ! linearly between 0 (at j=.5) and 1/nud_hrs (at j=il/2)
          rhs=min((j-.5)/(.5*il_g*nud_hrs),1./nud_hrs)
          do i=1,il_g
            if(n==0)davt_g(indglobal(i,il_g+1-j,n))=rhs
            if(n==2)davt_g(indglobal(j,i,n))=rhs
            if(n==3)davt_g(indglobal(j,i,n))=rhs
            if(n==5)davt_g(indglobal(i,il_g+1-j,n))=rhs
          enddo          ! i loop
        enddo            ! j loop
      enddo              ! n loop
      do j=1,il_g        ! full nudging on furthest panel
        do i=1,il_g
          davt_g(indglobal(j,i,4))=1./nud_hrs !  e.g. 1/48
        enddo            ! i loop
      enddo              ! j loop
    endif                !  (nbd==-7)  
    call ccmpi_distribute(davt,davt_g)
    deallocate(davt_g)
  else
    call ccmpi_distribute(davt)
  end if ! myid==0
  ! davu calc moved below next bounds call
  if(nproc==1)then
    write(6,*)'davt for i=il/2'
    write(6,'(20f6.3)') (davt(iq),iq=il/2,ifull,il)
  endif
  if(diag)call printa('davt',davt,0,0,ia,ib,ja,jb,0.,real(nud_hrs))     
endif                    ! (nbd.ne.0.and.nud_hrs.ne.0)

if(nbd>=3)then   ! separate (global) davu from (f-f) davt
  davu(:) = 1./nudu_hrs    !  e.g. 1/48
  write(6,*) 'all davu set to ',1./nudu_hrs 
elseif (nbd.ne.0) then
  davu(:) = davt(:)
endif


!--------------------------------------------------------------
! UPDATE GRAVITY WAVE DRAG DATA (lgwd)
gwdfac=.01*lgwd       ! most runs used .02 up to fri  10-10-1997
hefact=.1*abs(ngwd)   ! hal used hefact=1. (equiv to ngwd=10)
if(myid==0)write(6,*)'hefact,helim,gwdfac: ',hefact,helim,gwdfac
helo(:)=0.
if(lgwd>0)then
  do iq=1,ifull
    if(land(iq))then
      if(he(iq)==0.)write(6,*)'zero he over land for iq = ',iq
      aa(iq)=min(gwdfac*max(he(iq),.01),.8*zmin)   ! already in meters
      ! replace he by square of corresponding af value
      helo(iq)=( .4/log(zmin/aa(iq)) )**2
    endif
  enddo   ! iq loop
  if(mydiag)write(6,*)'for lgwd>0, typical zo#: ', diagvals(aa)
end if ! lgwd>0
if (ngwd/=0) then
!****    limit launching height : Palmer et al use limit on variance of
!****    (400 m)**2. we use launching height = std dev. we limit
!****    launching height to  2*400=800 m. this may be a bit severe.
!****    according to Palmer this prevents 2-grid noise at steep edge of
!****    himalayas etc.
  he(1:ifull)=min(hefact*he(1:ifull),helim)
endif     ! (ngwd.ne.0)


!--------------------------------------------------------------
! UPDATE BIOSPHERE DATA (nsib)
if (nsib==6.or.nsib==7) then
  ! Load CABLE data
  if (myid==0) write(6,*) 'Importing CABLE data'
  call loadtile
end if


!-----------------------------------------------------------------
! UPDATE MIXED LAYER OCEAN DATA (nmlo)
if (nmlo/=0.and.abs(nmlo)<=9) then
  if (any(ocndwn(:,1)>0.5)) then
    if (myid==0) write(6,*) 'Importing MLO data'
  else
    if (myid==0) write(6,*) 'Using MLO defaults'
    ocndwn(:,2)=0.
    do k=1,wlev
      call mloexpdep(0,depth,k,0)
      ! This polynomial fit is from MOM3, based on Levitus
      where (depth<2000.)
        mlodwn(:,k,1)=18.4231944        &
          -0.43030662E-1*depth(:)       &
          +0.607121504E-4*depth(:)**2   &
          -0.523806281E-7*depth(:)**3   &
          +0.272989082E-10*depth(:)**4  &
          -0.833224666E-14*depth(:)**5  &
          +0.136974583E-17*depth(:)**6  &
          -0.935923382E-22*depth(:)**7
        mlodwn(:,k,1)=mlodwn(:,k,1)*(tss-273.16)/18.4231944+273.16-wrtemp
      elsewhere
        mlodwn(:,k,1)=275.16-wrtemp
      end where
      mlodwn(:,k,2)=34.72
      mlodwn(:,k,3:4)=0.
    end do
    micdwn(:,1)=tss
    micdwn(:,2)=tss
    micdwn(:,3)=tss
    micdwn(:,4)=tss
    micdwn(:,5:6)=0.
    where (fracice>0.)
      micdwn(:,1)=min(tss,271.)
      micdwn(:,2)=min(tss,271.)
      micdwn(:,3)=min(tss,271.)
      micdwn(:,4)=min(tss,271.)
      micdwn(:,5)=1.
      micdwn(:,6)=2.
    end where
    micdwn(:,7:10)=0.
    micdwn(:,11)=10.
    where (.not.land)
      fracice=micdwn(:,5)
      sicedep=micdwn(:,6)
      snowd=micdwn(:,7)*1000.
    end where          
  end if
  !mlodwn(1:ifull,1:wlev,1)=max(mlodwn(1:ifull,1:wlev,1),271.-wrtemp)
  mlodwn(1:ifull,1:wlev,2)=max(mlodwn(1:ifull,1:wlev,2),0.)
  micdwn(1:ifull,1:4)=min(max(micdwn(1:ifull,1:4),0.),300.)
  micdwn(1:ifull,11)=max(micdwn(1:ifull,11),0.)
  call mloload(mlodwn,ocndwn(:,2),micdwn,0)
  deallocate(micdwn)
  do k=1,ms
    call mloexport(0,tgg(:,k),k,0)
    where ( tgg(:,k)<100. )
      tgg(:,k) = tgg(:,k) + wrtemp
    end where    
  end do
  do k=1,3
    call mloexpice(tggsn(:,k),k,0)
  end do 
end if


!-----------------------------------------------------------------
! UPDATE URBAN DATA (nurban)
if (nurban/=0) then
  if (myid==0) write(6,*) 'Importing ateb urban data'
  where(atebdwn(:,1)>=399.) ! must be the same as spval in onthefly.f
    atebdwn(:,1)=tss                ! roof temp 1
    atebdwn(:,2)=0.5*(tss+291.16)   ! roof temp 2
    atebdwn(:,3)=0.5*(tss+291.16)   ! roof temp 3
    atebdwn(:,4)=291.16             ! roof temp 4
    atebdwn(:,5)=tss                ! walleast temp 1
    atebdwn(:,6)=0.5*(tss+291.16)   ! walleast temp 2
    atebdwn(:,7)=0.5*(tss+291.16)   ! walleast temp 3
    atebdwn(:,8)=291.16             ! walleast temp 4
    atebdwn(:,9)=tss                ! wallwest temp 1
    atebdwn(:,10)=0.5*(tss+291.16)  ! wallwest temp 2
    atebdwn(:,11)=0.5*(tss+291.16)  ! wallwest temp 3
    atebdwn(:,12)=291.16            ! wallwest temp 4
    atebdwn(:,13)=tss               ! road temp 1
    atebdwn(:,14)=tss               ! road temp 2
    atebdwn(:,15)=tss               ! road temp 3
    atebdwn(:,16)=tss               ! road temp 4
    atebdwn(:,17)=0.5*0.26+0.5*0.18 ! Soil water road
    atebdwn(:,18)=0.18              ! Green roof water
    atebdwn(:,19)=0.   ! roof water
    atebdwn(:,20)=0.   ! road water
    atebdwn(:,21)=0.   ! canyon leaf water
    atebdwn(:,22)=0.   ! roof leaf water
    atebdwn(:,23)=0.   ! roof snow
    atebdwn(:,24)=0.   ! road snow
    atebdwn(:,25)=100. ! roof snow density
    atebdwn(:,26)=100. ! road snow density
    atebdwn(:,27)=0.85 ! roof snow albedo
    atebdwn(:,28)=0.85 ! road snow albedo
  end where
  call atebload(atebdwn,0)
  deallocate(atebdwn)
end if


!-----------------------------------------------------------------
! UPDATE AEROSOL DATA (iaero)
if ( abs(iaero)>=2 ) then
  xtg(1:ifull,:,:)=xtgdwn(:,:,:)
end if
      

!--------------------------------------------------------------     
! WRITE FORT.22 FILE FOR GRID INFO
if(nproc==1)then
  coslong=cos(rlong0*pi/180.)   
  sinlong=sin(rlong0*pi/180.)
  coslat=cos(rlat0*pi/180.)
  sinlat=sin(rlat0*pi/180.)
  polenx=-coslat
  poleny=0.
  polenz=sinlat
  write(22,920)
920     format(46x,'land            isoilm')
  write(22,921)
921     format('   iq     i    j  rlong    rlat    thet    map'          &
               '   sicedep zs(m) alb   ivegt  tss    t1    tgg2   tgg6'  &
               '   wb1   wb6   ico2  radon')
  do j=1,jl
    do i=1,il
      iq=i+(j-1)*il
      zonx=real(            -polenz*y(iq))
      zony=real(polenz*x(iq)-polenx*z(iq))
      zonz=real(polenx*y(iq)             )
      thet=atan2(-zonx*bx(iq)-zony*by(iq)-zonz*bz(iq),zonx*ax(iq)+zony*ay(iq)+zonz*az(iq))*180./pi
      if(thet<0.)thet=thet+360.
      write(22,922) iq,i,j,rlongg(iq)*180./pi,rlatt(iq)*180./pi,  &
                    thet,em(iq),land(iq),sicedep(iq),zs(iq)/grav, &
                    albvisnir(iq,1),                              &
                    isoilm(iq),ivegt(iq),                         &
                    tss(iq),t(iq,1),tgg(iq,2),tgg(iq,ms),         &
                    wb(iq,1),wb(iq,ms),                           &
                    0.
922      format(i6,2i5,3f8.3,f8.4,l2,f4.1,f7.1,f5.2,2i3,          &
                4f7.1,2f6.2,f5.2)
    enddo
  enddo
endif  ! (nproc==1)


!--------------------------------------------------------------
! INITIALISE STATION OUTPUT (nstn)
if(nstn>0)then
  if (myid==0) then
    write(6,*) 'land stations'
    write(*,"(a)") ' lu istn jstn  iq   slon   slat land rlong  rlat' &
                // ' isoil iveg zs(m) alb  wb3  wet3 vlai  zo   he'
  end if
  call ccmpi_barrier(comm_world)
  do nn=1,nstn
    call latltoij(slon(nn),slat(nn),rlong0,rlat0,schmidt,ri,rj,nface,xx4,yy4,il_g)
    ! These are global indices
    ig=nint(ri)
    jg=nint(rj) + nface*il_g
    mystn(nn) = fproc(ig,nint(rj),nface) == myid
    if ( mystn(nn) ) then
      iqg = ig + (jg-1)*il_g
      deli = nint(ri) - ri
      delj = nint(rj) - rj
      ! Local indices on this processor
      call indv_mpi(iqg,ii,jj,n)
      iq = ii + (jj-1)*ipan + (n-1)*ipan*jpan
      if(.not.land(iq))then
        ! simple search for neighbouring land point 
        ! N.B. does not search over panel/processor bndries
        ! N.B. if no land points, just returns closest point
        isav = ii
        jsav = jj
        dist=100.
        if ( isav < ipan ) then
          distnew = (deli+1)**2 + delj**2 
          if(land(iq+1).and.distnew<dist)then
            ii=isav+1
            dist=distnew
          endif
        end if
        if ( isav > 1 ) then
          distnew = (deli-1)**2 + delj**2 
          if(land(iq-1).and.distnew<dist)then
            ii=isav-1
            dist=distnew
          endif
        end if
        if ( jsav < jpan ) then
          distnew = deli**2 + (delj+1)**2 
          if(land(iq+ipan).and.distnew<dist)then
            jj=jsav+1
            dist=distnew
          endif
        end if
        if ( jsav >= 1 ) then
          distnew = deli**2 +(delj-1)**2 
          if(land(iq-ipan).and.distnew<dist)then
            jj=jsav-1
            dist=distnew
          endif
        end if
      endif              ! (.not.land(iq))
      istn(nn) = ii
      jstn(nn) = jj+(n-1)*jpan
      iq = istn(nn) + (jstn(nn)-1)*ipan
      iveg=ivegt(iq)
      isoil = isoilm(iq)
      wet3=(wb(iq,3)-swilt(isoil))/(sfc(isoil)-swilt(isoil))
      write(6,98)iunp(nn),istn(nn),jstn(nn),iq,slon(nn),slat(nn),   &
                 land(iq),rlongg(iq)*180/pi,rlatt(iq)*180/pi,       &
                 isoilm(iq),ivegt(iq),zs(iq)/grav,albvisnir(iq,1),  &
                 wb(iq,3),wet3,vlai(iq),zolnd(iq),he(iq),           &
                 myid             
    end if               ! mystn
98  format(i3,i4,i5,i6,2f7.2 ,l3,2f7.2, i3,i6,f7.1,f5.2,4f5.2,f7.1,i4)
    ! Put a barrier here to force stations to be printed in the right order
    call ccmpi_barrier(comm_world)
  enddo  ! nn=1,nstn
endif     !  (nstn>0)

    
!--------------------------------------------------------------
! OPEN MESONEST FILE
if ( mbd/=0 .or. nbd/=0 ) then
  if ( myid == 0 ) then
    write(6,*) "Opening mesonest file"
  end if
  io_in = io_nest                  ! Needs to be seen by all processors
  call histopen(ncid,mesonest,ier) ! open parallel mesonest files
  call ncmsg("mesonest",ier)       ! report error messages
  if ( myid == 0 ) then
    write(6,*) "ncid,mesonest ",ncid,trim(mesonest)
  end if
end if    ! (mbd/=0.or.nbd/=0)       


call END_LOG(indata_end)
return
end subroutine indataf


!--------------------------------------------------------------
! READ BIOSPHERIC FILES
subroutine rdnsib

use arrays_m                 ! Atmosphere dyamics prognostic arrays
use cc_mpi                   ! CC MPI routines
use infile                   ! Input file routines
use map_m                    ! Grid map arrays
use nsibd_m                  ! Land-surface arrays
use pbl_m                    ! Boundary layer arrays
use soil_m                   ! Soil and surface data
use soilsnow_m               ! Soil, snow and surface data
use tracers_m                ! Tracer data
use vegpar_m                 ! Vegetation arrays

implicit none

include 'newmpar.h'          ! Grid parameters
include 'const_phys.h'       ! Physical constants
include 'darcdf.h'           ! Netcdf data
include 'filnames.h'         ! Filenames
include 'parm.h'             ! Model configuration
include 'soilv.h'            ! Soil parameters
      
integer iq, iernc
integer ivegmin, ivegmax, ivegmax_g
integer :: idatafix = 0
integer, dimension(:,:), allocatable :: iduma
integer, dimension(ifull,2) :: idumb
integer, dimension(2) :: dumc
real sibvegver
real, dimension(:,:), allocatable :: duma
real, dimension(ifull,7) :: dumb
logical mismatch

real, parameter :: sibvegversion = 2015. ! version id for input data
real, parameter :: falbdflt = 0.
real, parameter :: frsdflt  = 990.
real, parameter :: fzodflt  = 1.
integer, parameter :: ivegdflt  = 0
integer, parameter :: isoildflt = 0

! isoilm_in holds the raw soil data.  isoilm is the data used
! by CCAM after the in-land water bodies (-1) and ocean (0)
! have been combined as water (0).  In the future, -1 will be
! used to better initialise lakes and salt emissions for
! aerosols

!------------------------------------------------------------------------
! READ BIOSPHERE FILES
! if cable, then the albedo is soil albedo only (converted to net albedo
! when cable is initialised)
if ( nsib <= 3 ) then
  if ( myid == 0 ) then
    write(6,*) "Start reading of nsib<=3 surface datafiles"
    allocate( iduma(ifull_g,2), duma(ifull_g,3) )
    write(6,*) "Reading albedo data"
    call readreal(albfile,duma(:,1),ifull_g)
    write(6,*) "Reading RSmin data"
    call readreal(rsmfile,duma(:,2),ifull_g)  ! not used these days
    write(6,*) "Reading roughness data"
    call readreal(zofile,duma(:,3),ifull_g)
    write(6,*) "Reading veg data"
    call readint(vegfile,iduma(:,1),ifull_g)
    write(6,*) "Reading soil data"
    call readint(soilfile,iduma(:,2),ifull_g)
    call ccmpi_distribute(idumb(:,1:2),iduma(:,1:2))
    call ccmpi_distribute(dumb(:,1:3),duma(:,1:3))
    deallocate( iduma, duma )
  else
    call ccmpi_distribute(idumb(:,1:2))
    call ccmpi_distribute(dumb(:,1:3))
  end if
  albvisnir(:,1) = 0.01*dumb(:,1)
  albvisnir(:,2) = 0.01*dumb(:,1) ! note VIS alb = NIR alb
  rsmin = dumb(:,2)
  zolnd = 0.01*dumb(:,3)
  ivegt = idumb(:,1)
  isoilm = idumb(:,2)
  isoilm_in = isoilm
else if ( nsib == 5 ) then
  if ( myid == 0 ) then
    write(6,*) "Start reading of nsib=5 (MODIS) surface datafiles"  
    allocate( duma(ifull_g,7) )
    if ( lncveg == 1 ) then
      call ccnf_get_attg(ncidveg,'sibvegversion',sibvegver,ierr=iernc)
      if ( iernc /= 0 ) then
        write(6,*) "Missing version of nsib=5 land-use data"
        write(6,*) "Regenerate land-use data with up-to-date version of sibveg"
        write(6,*) "or choose different land-surface option with nsib"
        call ccmpi_abort(-1)
      end if
      if ( sibvegver /= sibvegversion ) then
        write(6,*) "Wrong version of nsib=5 land-use data"
        write(6,*) "Expecting ",sibvegversion
        write(6,*) "Found     ",sibvegver
        call ccmpi_abort(-1)
      end if
      write(6,*) "Reading albedo data"
      call surfread(duma(:,1),'albvis',  netcdfid=ncidveg)
      call surfread(duma(:,2),'albnir',  netcdfid=ncidveg)
      write(6,*) "Reading RSmin data"
      call surfread(duma(:,3),'rsmin',   netcdfid=ncidveg)
      write(6,*) "Reading roughness data"
      call surfread(duma(:,4),'rough',   netcdfid=ncidveg)
      write(6,*) "Reading LAI data"
      call surfread(duma(:,5),'lai',     netcdfid=ncidveg)
      write(6,*) "Reading soil data"
      call surfread(duma(:,6),'soil',    netcdfid=ncidveg)      
      write(6,*) "Reading veg data"
      call surfread(duma(:,7),'landtype',netcdfid=ncidveg)      
    else
      write(6,*) "Cannot open vegfile as a netcdf file ",vegfile
      write(6,*) "Assuming ASCII file format"
      write(6,*) "Reading albedo data"
      call surfread(duma(:,1),'albvis',filename=albfile)
      call surfread(duma(:,2),'albnir',filename=albnirfile)
      write(6,*) "Reading RSmin data"
      call surfread(duma(:,3),'rsmin', filename=rsmfile)
      write(6,*) "Reading roughness data"
      call surfread(duma(:,4),'rough', filename=zofile)
      write(6,*) "Reading LAI data"
      call surfread(duma(:,5),'lai',   filename=laifile)
      write(6,*) "Reading soil data"
      call surfread(duma(:,6),'soilt', filename=soilfile)
      duma(:,7) = 1. ! vegt
      duma(:,1) = 0.01*duma(:,1)
      duma(:,2) = 0.01*duma(:,2)
      duma(:,4) = 0.01*duma(:,4)
      duma(:,5) = 0.01*duma(:,5)
    end if
    call ccmpi_distribute(dumb(:,1:7),duma(:,1:7))
    deallocate( duma )
  else
    call ccmpi_distribute(dumb(:,1:7))
  end if
  albvisnir(:,1) = dumb(:,1)
  albvisnir(:,2) = dumb(:,2)
  rsmin = dumb(:,3)
  zolnd = dumb(:,4)
  vlai = dumb(:,5)
  ivegt = nint(dumb(:,7))
  isoilm_in = nint(dumb(:,6))
  isoilm = max( isoilm_in, 0 )
else if ( nsib >= 6 ) then
  if ( myid == 0 ) then
    write(6,*) "Start reading of nsib>=6 (CABLE) surface datafiles"
    allocate( duma(ifull_g,3) )
    if ( lncveg == 1 ) then
      write(6,*) "Reading soil data"
      call surfread(duma(:,3),'soilt', netcdfid=ncidveg)
      write(6,*) "Reading albedo data"
      call surfread(duma(:,1),'albvis',netcdfid=ncidveg)
      call surfread(duma(:,2),'albnir',netcdfid=ncidveg)
    else
      write(6,*) "Cannot open vegfile as a netcdf file ",vegfile
      write(6,*) "Assuming ASCII file format"
      call surfread(duma(:,3),'soilt', filename=soilfile)
      call surfread(duma(:,1),'albvis',filename=albfile)
      call surfread(duma(:,2),'albnir',filename=albnirfile)
      duma(:,1:2) = 0.01*duma(:,1:2)
    end if
    call ccmpi_distribute(dumb(:,1:3),duma(:,1:3))
    deallocate( duma )
  else
    call ccmpi_distribute(dumb(:,1:3))
  end if
  ! communicate netcdf status to all processors
  albvisnir(:,1) = dumb(:,1)
  albvisnir(:,2) = dumb(:,2)
  isoilm_in = nint(dumb(:,3))
  isoilm = max( isoilm_in, 0 )
  zolnd = zobgin ! updated in cable_ccam2.f90
  ivegt = 1      ! updated in cable_ccam2.f90
end if
      
!--------------------------------------------------------------
! CHECK FOR LAND SEA MISMATCHES      
mismatch = .false.
if ( datacheck(land,albvisnir(:,1),'albv',idatafix,falbdflt) ) mismatch = .true.
if ( datacheck(land,albvisnir(:,2),'albn',idatafix,falbdflt) ) mismatch = .true.
if ( nsib < 6 ) then
  if ( datacheck(land,rsmin,'rsmin',idatafix,frsdflt) ) mismatch = .true.
end if
ivegmin = minval( ivegt, land(1:ifull) )
ivegmax = maxval( ivegt, land(1:ifull) )
if ( ivegmin<1 .or. ivegmax>44 ) then
  write(6,*) 'stopping in indata, as ivegt out of range'
  write(6,*) 'ivegmin,ivegmax ',ivegmin,ivegmax
  call ccmpi_abort(-1)
end if
if ( datacheck(land,isoilm,'isoilm',idatafix,isoildflt) ) mismatch = .true.

! --- rescale and patch up vegie data if necessary
if ( nsib/=6 .and. nsib/=7 ) then
  dumc(1)   = ivegmax
  call ccmpi_allreduce(dumc(1:1),dumc(2:2),"max",comm_world)
  ivegmax_g = dumc(2)
  if ( ivegmax_g < 14 ) then
    if ( mydiag ) write(6,*) '**** in this run veg types increased from 1-13 to 32-44'
    do iq = 1,ifull            ! add offset to sib values so 1-13 becomes 32-44
      if ( ivegt(iq) > 0 ) ivegt(iq) = ivegt(iq) + 31
    end do
  end if
end if
 
zolnd(:) = max( zolnd(:), zobgin )

return
end subroutine rdnsib


!--------------------------------------------------------------
! FUNCTIONS TO CHECK DATA
logical function rdatacheck( mask,fld,lbl,idfix,val )

implicit  none
      
include 'newmpar.h'      ! Grid parameters

integer, intent(in) :: idfix
real, intent(in) :: val
real, dimension(ifull), intent(inout) :: fld
logical, dimension(ifull), intent(in) :: mask
character(len=*), intent(in) :: lbl

rdatacheck=xrdatacheck(mask,fld,lbl,idfix,0.,val)
      
return
end function rdatacheck

    
logical function xrdatacheck(mask,fld,lbl,idfix,to,from)
      
use cc_mpi, only : myid ! CC MPI routines
      
implicit none

include 'newmpar.h'      ! Grid parameters      
      
real, intent(in) :: from,to
real, dimension(ifull), intent(inout) :: fld
integer, intent(in) :: idfix
integer iq
logical, dimension(ifull), intent(in) :: mask
logical err
character(len=*), intent(in) :: lbl

if (myid==0) write(6,*)' datacheck: verifying field ',lbl

err =.false.
do iq=1,ifull
  if( mask(iq) ) then
    if( fld(iq)==from ) then
      err = .true.
      if( idfix==1 ) then
        fld(iq) = to
        write(6,'(a,2i4,2(a,1pe12.4))') '  changing iq=',iq,' from',from,' to',to
      else
        write(6,*) '  mismatch at iq=',iq,', value',from
      end if
    end if
  end if
end do

xrdatacheck=err

return
end function xrdatacheck

    
logical function idatacheck( mask,ifld,lbl,idfix,ival )

implicit  none
      
include 'newmpar.h'      ! Grid parameters
      
integer, intent(in) :: idfix,ival
integer, dimension(ifull), intent(inout) :: ifld
logical, dimension(ifull), intent(in) :: mask
character(len=*), intent(in) :: lbl

idatacheck=xidatacheck(mask,ifld,lbl,idfix,0,ival)
      
return
end function idatacheck

    
logical function xidatacheck(mask,ifld,lbl,idfix,ifrom,ito)

use cc_mpi, only : myid ! CC MPI routines
      
implicit none

include 'newmpar.h'      ! Grid parameters      
      
integer, intent(in) :: idfix,ifrom,ito
integer iq
integer, dimension(ifull), intent(inout) :: ifld
logical, dimension(ifull), intent(in) :: mask
logical err
character(len=*), intent(in) :: lbl
      
if(myid==0) write(6,*)' datacheck: verifying field ',lbl
err =.false.
do iq=1,ifull
  if( mask(iq) ) then
    if( ifld(iq)==ifrom ) then
      err = .true.
      if( idfix==1 ) then
        ifld(iq) = ito
        write(6,'(a,2i4,2(a,i4))') '  changing iq=',iq,' from',ifrom,' to',ito
      else
        write(6,*) '  mismatch at iq=',iq,', value',ifrom
      end if
    end if
  end if
end do

xidatacheck = err

return
end function xidatacheck

    
!--------------------------------------------------------------
! READ INTEGER TEXT FILES
subroutine readint(filename,itss,ifully)
      
use cc_mpi            ! CC MPI routines
 
implicit none
      
include 'newmpar.h'   ! Grid parameters
include 'parm.h'      ! Model configuration
include 'parmgeom.h'  ! Coordinate data
            
integer ifully,ilx,jlx,ierr
integer, dimension(ifully) :: itss
integer, dimension(ifull_g) :: glob2d
real rlong0x,rlat0x,schmidtx,dsx
character(len=*) filename
character(len=47) header

if ( myid == 0 ) then
  write(6,*) 'reading data via readint from ',filename
  open(87,file=filename,status='old')
  read(87,*,iostat=ierr) ilx,jlx,rlong0x,rlat0x,schmidtx,dsx,header
  if ( ierr == 0 ) then
    write(6,*) ilx,jlx,rlong0x,rlat0x,schmidtx,dsx,header
    if(ilx/=il_g.or.jlx/=jl_g.or.abs(rlong0x-rlong0)>1.e-6.or.abs(rlat0x-rlat0)>1.E-6.or.abs(schmidtx-schmidt)>1.E-6) then
      write(6,*) 'wrong data file supplied'
      call ccmpi_abort(-1)
    end if
    read(87,*) glob2d
    close(87)
  else if ( ierr < 0 ) then ! Error, so really unformatted file
    close(87)
    write(6,*) 'now doing unformatted read'
    open(87,file=filename,status='old',form='unformatted')
    read(87) glob2d
    close(87)
  else ! ierr > 0
    write(6,*) 'End of file occurred in readint'
    call ccmpi_abort(-1)
  end if
  if (ifully==ifull) then
    call ccmpi_distribute(itss, glob2d)
  else if (ifully==ifull_g) then
    itss=glob2d
  else
    write(6,*) "ERROR: Invalid ifully for readint"
    call ccmpi_abort(-1)
  end if
  write(6,*) trim(header), glob2d(id+(jd-1)*il_g)
else
  if (ifully==ifull) then
    call ccmpi_distribute(itss)
  end if
end if
return
end subroutine readint

    
!--------------------------------------------------------------
! READ REAL TEXT FILES
subroutine readreal(filename,tss,ifully)
 
use cc_mpi            ! CC MPI routines
 
implicit none
      
include 'newmpar.h'   ! Grid parameters
include 'parm.h'      ! Model configuration
include 'parmgeom.h'  ! Coordinate data

integer ierr
integer ilx,jlx,ifully
real, dimension(ifully) :: tss
real, dimension(ifull_g) :: glob2d
real rlong0x,rlat0x,schmidtx,dsx
character(len=*) filename
character(len=47) header

if ( myid == 0 ) then
  write(6,*) 'reading data via readreal from ',trim(filename)
  open(87,file=filename,status='old')
  read(87,*,iostat=ierr) ilx,jlx,rlong0x,rlat0x,schmidtx,dsx,header
  if ( ierr == 0 ) then
    write(6,*) ilx,jlx,rlong0x,rlat0x,schmidtx,dsx,header
    if(ilx/=il_g.or.jlx/=jl_g.or.abs(rlong0x-rlong0)>1.e-6.or.abs(rlat0x-rlat0)>1.E-6.or.abs(schmidtx-schmidt)>1.E-6) then
      write(6,*) 'wrong data file supplied'
      call ccmpi_abort(-1)
    end if
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
    call ccmpi_abort(-1)
  end if
  if (ifully==ifull) then
    call ccmpi_distribute(tss, glob2d)
  else if (ifully==ifull_g) then
    tss=glob2d
  else
    write(6,*) "ERROR: Invalid ifully for readreal"
    call ccmpi_abort(-1)
  end if
  write(6,*) trim(header), glob2d(id+(jd-1)*il_g)
else
  if (ifully==ifull) then
    call ccmpi_distribute(tss)
  end if
end if
return
end subroutine readreal


!--------------------------------------------------------------
! INITALISE SOIL PARAMETERS
subroutine insoil
      
use cc_mpi, only : myid ! CC MPI routines
      
implicit none
      
include 'newmpar.h'     ! Grid parameters
include 'soilv.h'       ! Soil parameters
      
! The following common block is for soilsnow.f
real zshh, ww
common/soilzs/zshh(ms+1),ww(ms)

integer isoil, k

do isoil = 1,mxst
  cnsd(isoil)  = sand(isoil)*0.3+clay(isoil)*0.25+silt(isoil)*0.265
  hsbh(isoil)  = hyds(isoil)*abs(sucs(isoil))*bch(isoil) !difsat*etasat
  ibp2(isoil)  = nint(bch(isoil))+2
  i2bp3(isoil) = 2*nint(bch(isoil))+3
  if ( myid == 0 ) then
    write(6,"('isoil,ssat,sfc,swilt,hsbh ',i2,3f7.3,e11.4)") isoil,ssat(isoil),sfc(isoil),swilt(isoil),hsbh(isoil)
  end if
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

      
!--------------------------------------------------------------
! CALCULATE ROUGHNESS LENGTH
subroutine calczo
      
use arrays_m        ! Atmosphere dyamics prognostic arrays
use map_m           ! Grid map arrays
use nsibd_m         ! Land-surface arrays
use soil_m          ! Soil and surface data
use soilsnow_m      ! Soil, snow and surface data
      
implicit none
      
include 'newmpar.h' ! Grid parameters
include 'parm.h'    ! Model configuration
      
integer iq,iveg
real zomax,zomin,tsoil,sdep
      
real xhc(0:44)
!     vegetation height
data xhc    / 0.0,                                                   & ! 0
              30.0,28.0,25.0,17.0,12.0,10.0, 9.0, 7.0, 5.5, 3.0,     & ! 1-10
              2.5, 2.0, 1.0, 0.6, 0.5, 0.5,0.45,0.75, 0.6,0.45,      &
              0.4, 0.6, 0.6,0.24,0.25,0.35, 0.3, 2.5, 0.0, 0.0,      &
              0.0,                                                   & ! 31
              32.,20.,20.,17.,17., 1., 1., 1., 0.5, 0.6, 0., 1.,0./    !sellers 1996 j.climate

zomax=-1.e29
zomin= 1.e29
if(newrough==2)then
  do iq=1,ifull
    if(land(iq))then
      iveg=ivegt(iq)
      zolnd(iq)=max(zobgin , .1*xhc(iveg))
      zomax=max(zomax,zolnd(iq))
      zomin=min(zomin,zolnd(iq))
    endif  ! (land(iq))then
  enddo   ! iq loop
elseif(newrough==3)then
  do iq=1,ifull
    if(land(iq))then
      iveg=ivegt(iq)
      zolnd(iq)=max(zobgin , .13*xhc(iveg))  ! French factor
      zomax=max(zomax,zolnd(iq))
      zomin=min(zomin,zolnd(iq))
    endif  ! (land(iq))then
  enddo   ! iq loop
else
  do iq=1,ifull
    if(land(iq))then
      iveg=ivegt(iq)
      tsoil  = 0.5*(tgg(iq,ms)+tgg(iq,2))
      sdep=0.
      call cruf1 (iveg,tsoil,sdep,zolnd(iq),zobgin)
      zomax=max(zomax,zolnd(iq))
      zomin=min(zomin,zolnd(iq))
    endif  ! (land(iq))then
  enddo   ! iq loop
endif

write(6,*)"calczo zolnd: zomin,zomax=",zomin,zomax
return ! calczo
end subroutine calczo

subroutine cruf1(iv,tsoil,sdep,zolnd,zobgin)
      
implicit none
      
! kf, 1997
! for each vegetation type (= iv), assign veg height, total lai, albedo,
! and computed aerodynamic, radiative and interception properties.
! jmax0 assigned due to table by ray leuning and estimates  21-08-97 
! apply seasonal variations in lai and height. return via /canopy/
! nb: total lai = xrlai, veglai = xvlai, veg cover fraction = xpfc,
!     with xrlai = xvlai*xpfc
! type  0 to 31: 2d version with graetz veg types
! type 32 to 43: 2d version with gcm veg types
! type 44:       stand-alone version
!-----------------------------------------------------------------------
!   name                             symbol  code hc:cm pfc:%  veglai
!   ocean                                o!     0     0     0  0.0
!   tall dense forest                    t4     1  4200   100  4.8
!   tall mid-dense forest                t3     2  3650    85  6.3
!   dense forest                         m4     3  2500    85  5.0  (?)
!   mid-dense forest                     m3     4  1700    50  3.75
!   sparse forest (woodland)             m2     5  1200    20  2.78
!   very sparse forest (woodland)        m1     6  1000     5  2.5
!   low dense forest                     l4     7   900    85  3.9
!   low mid-dense forest                 l3     8   700    50  2.77
!   low sparse forest (woodland)         l2     9   550    20  2.04
!   tall mid-dense shrubland (scrub)     s3    10   300    50  2.6
!   tall sparse shrubland                s2    11   250    20  1.69
!   tall very sparse shrubland           s1    12   200     5  1.9
!   low mid-dense shrubland              z3    13   100    50  1.37
!   low sparse shrubland                 z2    14    60    20  1.5
!   low very sparse shrubland            z1    15    50     5  1.21
!   sparse hummock grassland             h2    16    50    20  1.58
!   very sparse hummock grassland        h1    17    45     5  1.41
!   dense tussock grassland              g4    18    75    85  2.3
!   mid-dense tussock grassland          g3    19    60    50  1.2
!   sparse tussock grassland             g2    20    45    20  1.71
!   very sparse tussock grassland        g1    21    40     5  1.21
!   dense pasture/herbfield (perennial)  f4    22    60    85  2.3
!   dense pasture/herbfield (seasonal)  f4s    23    60    85  2.3
!   mid-dense pasture/herb (perennial)   f3    24    45    50  1.2
!   mid-dense pasture/herb  (seasonal)  f3s    25    45    50  1.2
!   sparse herbfield*                    f2    26    35    20  1.87
!   very sparse herbfield                f1    27    30     5  1.0
!   littoral                             ll    28   250    50  3.0
!   permanent lake                       pl    29     0     0  0
!   ephemeral lake (salt)                sl    30     0     0  0
!   urban                                 u    31     0     0  0
!   stand alone: hc,rlai from param1      -    44     -   100  -

!   above are dean's. below are sib (added 31 to get model iveg)
!  32  1 - broadleaf evergreen trees (tropical forest)
!  33  2 - broadleaf deciduous trees
!  34  3 - broadleaf and needleaf trees
!  35  4 - needleaf evergreen trees
!  36  5 - needleaf deciduous trees 
!  37  6 - broadleaf trees with ground cover (savannah)
!  38  7 - groundcover only (perennial)
!  39  8 - broadleaf shrubs with groundcover
!  40  9 - broadleaf shrubs with bare soil
!  41 10 - dwarf trees and shrubs with groundcover
!  42 11 - bare soil
!  43 12 - agriculture or C3 grassland (newer defn)
 
!                             soil type
!       texture               
!  0   water/ocean
!  1   coarse               sand/loamy_sand
!  2   medium               clay-loam/silty-clay-loam/silt-loam
!  3   fine                 clay
!  4   coarse-medium        sandy-loam/loam
!  5   coarse-fine          sandy-clay
!  6   medium-fine          silty-clay 
!  7   coarse-medium-fine   sandy-clay-loam
!  8   organi!              peat
!  9   land ice
!-----------------------------------------------------------------------

integer iv
real ftsoil,vrlai,hc,tsoil,sdep,zolnd,zobgin
real rlai,usuh,disp,coexp
real xhc(0:44),xpfc(0:44),xvlai(0:44),xslveg(0:44)
! aerodynamic parameters, diffusivities, water density:
real vonk,a33,csw,ctl
parameter(vonk   = 0.40)     ! von karman constant
parameter(a33    = 1.25)     ! inertial sublayer sw/us
parameter(csw    = 0.50)     ! canopy sw decay (weil theory)
parameter(ctl    = 0.40)     ! wagga wheat (rdd 1992, challenges)
! vegetation height
data xhc    / 0.0,                                                   & ! 0
              30.0,28.0,25.0,17.0,12.0,10.0, 9.0, 7.0, 5.5, 3.0,     & ! 1-10
              2.5, 2.0, 1.0, 0.6, 0.5, 0.5,0.45,0.75, 0.6,0.45,      &
              0.4, 0.6, 0.6,0.24,0.25,0.35, 0.3, 2.5, 0.0, 0.0,      &
              0.0,                                                   & ! 31
              32.,20.,20.,17.,17., 1., 1., 1., 0.5, 0.6, 0., 1.,0./    !sellers 1996 j.climate

! vegetation fractional cover
data xpfc   /0.00,                                                   &
             1.00,0.85,0.85,0.50,0.20,0.05,0.85,0.50,0.20,0.50,      &
             0.20,0.05,0.50,0.20,0.05,0.20,0.05,0.85,0.50,0.20,      &
             0.05,0.85,0.85,0.50,0.50,0.20,0.05,0.50,0.00,0.00,      &
             0.00,                                                   &
             .98,.75,.75,.75,.50,.86,.65,.79,.30,.42,.02,.54,  1.0/

! veg lai from graetz table of 283 veg types (iv=0 to 31), and maximum 
! veg lai for gcm veg types (iv=32 to 43)  stand-alone: 44
data xvlai  / 0.0,                                                   &
              4.80,6.30,5.00,3.75,2.78,2.50,3.90,2.77,2.04,2.60,     &
              1.69,1.90,1.37,1.50,1.21,1.58,1.41,2.30,1.20,1.71,     &
              1.21,2.30,2.30,1.20,1.20,1.87,1.00,3.00,0.00,0.00,     &
              0.00,                                                  &
              6.0,5.0,4.0,4.0,4.0,3.0,3.0,3.0,1.0,4.0,0.5,3.0,  0.0/     ! 32-44

! for seasonally varying lai, amplitude of veg lai seasonal change
data xslveg  /0.00,                                                  &
              0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,     & 
              0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,     &
              0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,     &
              0.00,                                                  &
              2.0,2.0,2.0,2.0,2.0,1.5,1.5,1.5,1.0,0.5,0.5,0.5,  0.0/
!-----------------------------------------------------------------------
! assign aerodynamic, radiative, stomatal, interception properties
! assign total lai (xrlai) from veg lai and pfc, and assign seasonal 
!   variation in lai and veg height where necessary. this is controlled
!   by the factor season (0 =< season =< 1).
ftsoil=max(0.,1.-.0016*(298.-tsoil)**2)
if( tsoil >= 298. ) ftsoil=1.
vrlai = max(0.0,(xvlai(iv)-xslveg(iv)*(1.-ftsoil))*xpfc(iv))
hc    = max(0.0,xhc(iv) - sdep)
rlai  = vrlai*hc/max(0.01,xhc(iv))
!   find roughness length zolnd from hc and rlai:
call cruf2(hc,rlai,usuh,zolnd,disp,coexp)
!   set aerodynamic variables for bare soil and vegetated cases:
zolnd=max(zolnd, zobgin)
if (rlai<0.001 .or. hc<.05) then
  zolnd    = zobgin      ! bare soil surface
  hc     = 0.0  
  rlai   = 0.0
endif
return
end subroutine cruf1
!=======================================================================
subroutine cruf2(h,rlai,usuh,z0,d,coexp)
      
implicit none
      
!-----------------------------------------------------------------------
! m.r. raupach, 24-oct-92
! see: raupach, 1992, blm 60 375-395
!      mrr notes "simplified wind model for canopy", 23-oct-92
!      mrr draft paper "simplified expressions...", dec-92
!-----------------------------------------------------------------------
! inputs:
!   h     = roughness height
!   rlai  = leaf area index (assume rl = frontal area index = rlai/2)
! output:
!   usuh  = us/uh (us=friction velocity, uh = mean velocity at z=h)
!   z0    = roughness length
!   d     = zero-plane displacement
!   coexp = coefficient in exponential in-canopy wind profile
!           u(z) = u(h)*exp(coexp*(z/h-1)), found by gradient-matching
!           canopy and roughness-sublayer u(z) at z=h
!-----------------------------------------------------------------------
real h,rlai,usuh,z0,d,coexp
real psih,rl,usuhl,xx,dh,z0h
! preset parameters:
real cr,cs,beta,ccd,ccw,usuhm,vonk
parameter (cr    = 0.3)          ! element drag coefficient
parameter (cs    = 0.003)        ! substrate drag coefficient
parameter (beta  = cr/cs)        ! ratio cr/cs
parameter (ccd   = 15.0)         ! constant in d/h equation
parameter (ccw   = 2.0)          ! ccw=(zw-d)/(h-d)
parameter (usuhm = 0.3)          ! (max of us/uh)
parameter (vonk  = 0.4)          ! von karman constant
psih=alog(ccw)-1.0+1.0/ccw  ! i.e. .19315
rl = rlai*0.5
! find uh/us
usuhl  = sqrt(cs+cr*rl)            ! sqrt(.003 + .15*rlai)
usuh   = min(usuhl,usuhm)
! find d/h and d 
xx     = sqrt(ccd*max(rl,0.0005))  ! sqrt(7.5*rlai)
dh     = 1.0 - (1.0 - exp(-xx))/xx ! .5*xx -.166*xx*xx + .
!     dh close to 1 for large rlai (~.833 for rlai=6)    
!        equals .03 for rlai=0 
d      = dh*h                      ! not used
! find z0h and z0:
z0h    = (1.0 - dh) * exp(psih - vonk/usuh)
z0     = z0h*h
!     for rlai=   0,   .2,   .4,   .5,   .6,    1,    2,    4,    6
!     get z0h= .008, .077, .117, .128, .133, .109, .084, .058, .048      
! find coexp: see notes "simplified wind model ..." eq 34a
coexp  = usuh / (vonk*ccw*(1.0 - dh))
return ! ruff
end subroutine cruf2
!=======================================================================


!--------------------------------------------------------------
! SPECIAL FUNCTION FOR SSTs
subroutine caispecial
      
use cc_mpi
use infile
use latlong_m
use pbl_m
use soil_m
use soilsnow_m
      
implicit none
      
include 'newmpar.h'
include 'const_phys.h'  
      
integer iq,ix
integer ncid,ncs,varid
integer, dimension(3) :: spos,npos
real x,r
real, dimension(300) :: sdata,ldata
logical tst
      
if (myid==0) then
  write(6,*) "Reading nspecial=42 SSTs"
  spos=1
  call ccnf_open('sst_djf.cdf',ncid,ncs)
  if (ncs/=0) then
    write(6,*) "ERROR: Cannot open sst_djf.cdf"
    call ccmpi_abort(-1)
  end if
  npos=1
  npos(1)=300
  call ccnf_inq_varid(ncid,'SST_DJF',varid,tst)
  if (tst) then
    write(6,*) "ERROR: Cannot read SST_DJF"
    call ccmpi_abort(-1)
  end if
  call ccnf_get_vara(ncid,varid,spos(1:1),npos(1:1),sdata)
  npos=1
  npos(1)=300
  call ccnf_inq_varid(ncid,'YT_OCEAN',varid,tst)
  if (tst) then
    write(6,*) "ERROR: Cannot read SST_DJF"
    call ccmpi_abort(-1)
  end if
  call ccnf_get_vara(ncid,varid,spos(1:1),npos(1:1),ldata)
  call ccnf_close(ncid)
  sdata=sdata+273.16
end if
call ccmpi_bcast(sdata,0,comm_world)
call ccmpi_bcast(ldata,0,comm_world)
      
do iq=1,ifull
  if (.not.land(iq)) then
    r=rlatt(iq)*180./pi
    if (r.lt.ldata(2)) then
      tss(iq)=sdata(2)
    elseif (r.gt.ldata(300)) then
      tss(iq)=sdata(300)
    else
      do ix=2,300
        if (ldata(ix).gt.r) exit
      end do
      x=(r-ldata(ix))/(ldata(ix+1)-ldata(ix))
      tss(iq)=(1.-x)*sdata(ix)+x*sdata(ix+1)
    end if
    tgg(iq,1)=tss(iq)
  end if
end do
      
return
end subroutine caispecial

end module indata