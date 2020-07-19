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
    
module tracermodule

implicit none
      
private
public sitefile,shipfile
public co2em,init_tracer
public tracer_mass,interp_tracerflux,tracerlist
public writetrpm
public readtracerflux,tracunit,tracdaytime
public oh,strloss,mcfdep,jmcf
public trden, trreff, trdep, traclevel
public numtracer

integer, dimension(:), save, allocatable :: tracinterp
integer, dimension(:), save, allocatable :: nghr, igashr
integer, save :: numtracer, nhr
integer, dimension(:), save, allocatable :: traclevel
real, dimension(:,:,:), save, allocatable :: co2emhr, co2em123
real, dimension(:,:), save, allocatable :: co2hr, co2em
real, dimension(:,:), save, allocatable :: tracdaytime
real, dimension(:), save, allocatable :: tracival, trden, trreff, trdep
character(len=13), dimension(:), save, allocatable :: tracunit
character(len=50), dimension(:), save, allocatable :: tracfile
character(len=80), save :: tracerlist=' '
character(len=80), save :: sitefile=' '
character(len=80), save :: shipfile=' '
logical, save :: writetrpm=.false.
     
! rml 16/02/10 additions for Transcom methane
real, dimension(:,:,:), save, allocatable :: oh123, strloss123
real, dimension(:,:), save, allocatable :: oh, strloss
! rml 30/04/10 additions for Transcom MCF
real, dimension(:,:,:), save, allocatable :: jmcf123
real, dimension(:,:), save, allocatable :: mcfdep123, jmcf
real, dimension(:), save, allocatable :: mcfdep
logical methane, mcf

contains

! ***************************************************************************
subroutine init_tracer
use cc_mpi, only : myid, ccmpi_abort
use tracers_m
use newmpar_m
implicit none
integer nt
real tracmin,tracmax
character(len=80) :: header
logical dpoflag

! check for valid tracer data
if ( tracerlist==' ' ) then
  allocate( tr(ifull+iextra,kl,ntrac) ) ! for pgfortran  
  return
end if    


!     first read in a list of tracers, then determine what emission data is required
!     Each processor read this, perhaps not necessary?

!     rml 16/2/10 addition for Transcom methane
methane = .false.
mcf     = .false.
dpoflag = .false.

open(unit=130,file=tracerlist,form='formatted')
read(130,*) header
read(130,*) numtracer
if (numtracer<1.or.numtracer>999) then
  write(6,*) "ERROR: Invalid number of tracers"
  write(6,*) "numtracer should be between 1-999"
  write(6,*) "numtracer = ",numtracer
  call ccmpi_abort(-1)
end if
ngas=numtracer
allocate(tracname(numtracer),tractype(numtracer))
allocate(tracinterp(numtracer),tracunit(numtracer))
allocate(tracfile(numtracer),igashr(numtracer))
allocate(tracival(numtracer))
allocate(trden(numtracer),trreff(numtracer),trdep(numtracer))
allocate(traclevel(numtracer))
nhr = 0
do nt=1,numtracer
  !           name         inital value                 Emiss type   Emiss file   Density   Radius     Dry Dep?  Emiss levels
  read(130,*) tracname(nt),tracival(nt),tracmin,tracmax,tractype(nt),tracfile(nt),trden(nt),trreff(nt),trdep(nt),traclevel(nt)
  select case(tractype(nt))
   case ('monrep','month')
    tracinterp(nt)=1
   case ('hourly','3hourly','daily')
    tracinterp(nt)=2
    nhr = nhr+1
    igashr(nt) = nhr
   case ('online')
    tracunit(nt) = 'gC/m2/s'
    tracinterp(nt) = -1
   case ('monpwcb')
    tracinterp(nt)=3
   case ('mon1st')
    tracinterp(nt)=4
   case ('daypulseon')
    tracinterp(nt)=0
    dpoflag=.true.
   case default
    tracinterp(nt)=0
  end select
   if ( myid == 0 ) then
    write(6,'(a7,i5,1x,a13,a13,a50,i3)') 'Tracer ',nt,tracname(nt),tractype(nt),tracfile(nt),tracinterp(nt)
  end if
! rml 16/2/10 addition for TC methane
  if (tracname(nt)(1:7)=='methane') methane = .true.
  if (tracname(nt)(1:3)=='mcf') mcf = .true.

enddo

if (nhr>0) then
  allocate(nghr(nhr))
else
  deallocate(igashr)
end if

if (dpoflag) then
  allocate(tracdaytime(numtracer,2))    
end if

! MJT - moved here to define ntrac, ilt, jlt and klt
call tracers_init(il,jl,kl,iextra)

if (methane) then
  allocate(acloss_g(ntrac))
!       initialise accumulated loss (for methane cases)
  acloss_g = 0.
end if
!     if writing afternoon averages, initialise here
if (writetrpm) then
  allocate(trpm(il*jl,kl,ntrac),npm(il*jl))
  trpm = 0.
  npm = 0
end if

close(130)

return
end subroutine init_tracer

! *********************************************************************
subroutine readtracerflux
!     needs to happen further down the code than reading the tracer
!     list file

use cc_mpi, only : myid
use dates_m, only : kdate
use newmpar_m
use parm_m
use tracers_m

implicit none

integer nt,jyear,jmonth
real, dimension(3) :: ajunk
character(len=50) :: filename
character(len=13) :: varname

jyear=kdate/10000
jmonth=(kdate-jyear*10000)/100

do nt=1,numtracer
  select case(tracinterp(nt))
    case (0)
!           constant data
      if (.not.allocated(co2em123)) then
        allocate(co2em123(il*jl,2:2,numtracer))
      end if
      if (.not.allocated(co2em)) then
        allocate(co2em(il*jl,numtracer))
      end if
      call readrco2(nt,jyear,jmonth,1,co2em123(:,2:2,nt),ajunk(2:2))
    case (1)
!           monthly data
      if (.not.allocated(co2em123)) then
        allocate(co2em123(il*jl,3,numtracer))
      end if
      if (.not.allocated(co2em)) then
        allocate(co2em(il*jl,numtracer))
      end if
      call readrco2(nt,jyear,jmonth,3,co2em123(:,:,nt),ajunk)
    case (2)
!           daily, 3 hourly, hourly
      if (.not.allocated(co2emhr)) then
        allocate(co2emhr(il*jl,31*24+2,nhr))
        allocate(co2hr(31*24+2,nhr))
      end if
      if (.not.allocated(co2em)) then
        allocate(co2em(il*jl,numtracer))
      end if
      call readrco2(nt,jyear,jmonth,31*24+2,co2emhr(:,:,igashr(nt)),co2hr(:,igashr(nt)))
  end select
enddo
!
!     just read OH, strat loss once regardless of how many methane tracers
!     filename set here at compile time
if (methane) then
  allocate(oh123(il*jl,kl,3))
  allocate(strloss123(il*jl,kl,3))
  allocate(oh(il*jl,kl))
  allocate(strloss(il*jl,kl))
  filename = '/short/r39/TCinput/oh_c48.nc'
  varname = 'oh'
  call readoh(jmonth,3,filename,varname,oh123)
  filename = '/short/r39/TCinput/strloss_c48.nc'
  varname = 'strloss'
  call readoh(jmonth,3,filename,varname,strloss123)
endif
if (mcf) then
  allocate(mcfdep123(il*jl,3))
  allocate(mcfdep(il*jl))
  allocate(jmcf123(il*jl,kl,3))
  allocate(jmcf(il*jl,kl))
!       use standard tracer flux file read call but pass through 
!       with tracer 'ngas+1' - this will trigger MCF_loss as filename
   call readrco2(ngas+1,jyear,jmonth,3,mcfdep123,ajunk)
   filename='/short/r39/TCinput/Jmcf_cc48.nc'
   varname = 'jmcf'
   call readoh(jmonth,3,filename,varname,jmcf123)
 endif
 if (myid==0) then
   write(6,*) 'Read input fluxes/rates etc OK'
 end if

 return
 end subroutine

! *************************************************************************
 subroutine readrco2(igas,iyr,imon,nflux,fluxin,co2time)
!     rml 23/09/03 largely rewritten to use netcdf files
 use cc_mpi
 use infile
 use newmpar_m
 use parm_m
 use tracers_m
 implicit none
 
!     nflux =3 for month interp case - last month, this month, next month
!     nflux=31*24+2 for daily, hourly, 3 hourly case
 
 integer ncidfl,fluxid
 integer nregion,dayid
 integer nflux,iyr,imon,igas,ntime,lc,regnum,kount
 integer nprev,nnext,ncur,n1,n2,ntot,n,timx,gridid,ierr
 integer, dimension(3) :: start,ncount
 integer, dimension(:), allocatable :: fluxyr,fluxmon
 real timeinc, hr
 real, dimension(ifull,nflux) :: fluxin
 real, dimension(ifull_g,nflux) :: fluxin_g
 real, dimension(nflux) :: co2time
 real, dimension(2) :: dum
 real, dimension(:), allocatable :: fluxhr
 character(len=50) filename
 character(len=13) fluxtype,fluxname,fluxunit
logical gridpts,tst

n1=0
n2=0

fluxtype=tractype(igas)
fluxname=tracname(igas)
filename=tracfile(igas)

if ( trim(fluxtype)=='pulseoff' .or. trim(fluxtype)=='daypulseoff' ) then
!       no surface fluxes to read
  fluxin = 0.
  tracunit(igas)=' '    
  return
end if

if ( myid == 0 ) then ! Read on this processor and then distribute
  write(6,*)'reading ',trim(fluxname), ' with type ',trim(fluxtype),' for ',iyr,imon,' from ',filename
  call ccnf_open(filename,ncidfl,ierr)
  call ncmsg("readrco2",ierr)
  ! check for 1D or 2D formatted input
  call ccnf_inq_dimid(ncidfl,'gridpts',gridid,gridpts)
  gridpts=.not.gridpts
  call ccnf_inq_dimlen(ncidfl,'time',ntime)
  allocate(fluxyr(ntime),fluxmon(ntime))
  call ccnf_get_vara(ncidfl,'year',fluxyr)
  call ccnf_get_vara(ncidfl,'month',fluxmon)
  
!       read hours variable for daily, hourly, 3 hourly data
  if (nflux==(31*24+2)) then
    allocate(fluxhr(ntime))
    call ccnf_get_vara(ncidfl,'hour',fluxhr)
  endif

!       check fluxname first otherwise default to 'flux'
  call ccnf_inq_varid(ncidfl,fluxname,fluxid,tst)
  if (tst) then
    call ccnf_inq_varid(ncidfl,'flux',fluxid,tst)
  endif
  if (tst) then
    write(6,*) 'ERROR: flux variable not found'
    call ccmpi_abort(-1)
  end if  
!       rml 25/08/04 read flux units attribute
  call ccnf_get_att(ncidfl,fluxid,'units',fluxunit)
! rml 08/11/04 added radon units
! rml 30/4/10 exclude mcf deposition case
  if ( igas<=ngas ) then
    tracunit(igas)=fluxunit
    if ( fluxunit(1:7)/='gC/m2/s' .and. fluxunit(1:7)/='Bq/m2/s' .and. fluxunit(1:8)/='mol/m2/s' ) then
      write(6,*) 'Units for ',trim(fluxname),' are ',trim(fluxunit)
      write(6,*) 'Code not set up for units other than gC/m2/s or mol/m2/s or Bq/m2/s'
      write(6,*) 'fix flux units'
      call ccmpi_abort(-1)
    endif
  endif

  if ( trim(fluxtype)=='daypulseon' ) then
!         need to read sunset/sunrise times
    call ccnf_inq_dimlen(ncidfl,'nregion',nregion)
    call ccnf_inq_varid(ncidfl,'daylight',dayid,tst)
  endif
!
!       find required records
  if ( nflux==1 ) then
!         constant case !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!           check ntime
    if (ntime/=1) then
      write(6,*) 'flux file wrong ntime'
      call ccmpi_abort(-1)
    end if
    write(6,*)'reading ',1,fluxyr(1),fluxmon(1)

    fluxin_g=0.
    co2time=0.
    if (gridpts) then
      start(1:2)=1
      ncount(1)=ifull_g; ncount(2)=1
    else
      start(1:3)=1
      ncount(1)=il_g; ncount(2)=jl_g; ncount(3)=1
    end if
!         read current month/year
    call ccnf_get_vara(ncidfl,fluxid,start,ncount,fluxin_g(:,1))
    
  else if ( nflux==3 ) then
!         monthly/annual cases !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    nprev=0
    nnext=0
    ncur=0
    if (trim(fluxtype)=='annual') then
      do n=1,ntime
        if (fluxyr(n)==iyr) ncur=n
      enddo
    elseif (trim(fluxtype)=='monrep') then
      if (ntime/=12) then
        write(6,*) 'flux file wrong ntime'
        call ccmpi_abort(-1)
      end if
      do n=1,ntime
        if (fluxmon(n)==imon) then
          ncur=n
          nprev=n-1
          if (nprev==0) nprev=12
          nnext=n+1
          if (nnext==13) nnext=1
        endif
      enddo
    else
!           monthly case
      do n=1,ntime
        if (fluxyr(n)==iyr.and.fluxmon(n)==imon) then
          ncur=n
          nprev=n-1
!               keep flux constant at ends of data
          if (nprev==0) nprev=n
          nnext=n+1
          if (nnext==ntime+1) nnext=n
        endif
      enddo
    endif
!
    if (ncur==0) then
      write(6,*) 'ERROR: current year/month not in flux file'
      call ccmpi_abort(-1)
    end if
    write(6,*)'reading ',ncur,fluxyr(ncur),fluxmon(ncur)
!    
    fluxin_g=0.
    co2time=0.
    if (gridpts) then
      start(1)=1
      ncount(1)=ifull_g; ncount(2)=1
      timx=2
    else
      start(1:2)=1
      ncount(1)=il_g; ncount(2)=jl_g; ncount(3)=1
      timx=3
    end if
!         read preceeding month if needed
    if (nprev/=0) then
      start(timx)=nprev
      call ccnf_get_vara(ncidfl,fluxid,start,ncount,fluxin_g(:,1))
    endif
!         read current month/year
    start(timx)=ncur
    call ccnf_get_vara(ncidfl,fluxid,start,ncount,fluxin_g(:,2))
!         read next month
    if (nnext/=0) then
      start(timx)=nnext
      call ccnf_get_vara(ncidfl,fluxid,start,ncount,fluxin_g(:,3))
    endif
  else
!         daily, hourly, 3 hourly case !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!         find first time in month
!         rml 06/01/06 extend case to cope with annually repeating or
!         real year fluxes
    do n=1,ntime
      if ( ( (fluxyr(n)==0).or.(fluxyr(n)==iyr) ) .and. (fluxmon(n)==imon) ) then
        n1=n
        exit
      endif
    enddo
!         find last time in month
    do n=ntime,1,-1
      if ( ( (fluxyr(n)==0).or.(fluxyr(n)==iyr) ) .and. (fluxmon(n)==imon) ) then
        n2=n
        exit
      endif
    enddo
!         read fluxes
    ntot=n2-n1+1
    if (gridpts) then
      start(1)=1; ncount(1)=ifull_g
      start(2)=n1; ncount(2)=ntot
      timx=2
    else
      start(1)=1; ncount(1)=il_g
      start(2)=1; ncount(2)=jl_g
      start(3)=1; ncount(3)=ntot
      timx=3
    end if
    call ccnf_get_vara(ncidfl,fluxid,start,ncount,fluxin_g(:,2:ntot+1))
!         read in last time of prev month and first time of next month
    if ((n1==1).and.(fluxyr(n1)==0)) then
!           loop
      nprev=ntime
    elseif ((n1==1).and.(fluxyr(n1)/=0)) then
!           keep constant
      nprev=n1
    else
      nprev=n1-1
    endif
    start(timx)=nprev; ncount(timx)=1
    call ccnf_get_vara(ncidfl,fluxid,start,ncount,fluxin_g(:,1))
    if ((n2==ntime).and.(fluxyr(n2)==0)) then
      nnext=1
    elseif((n2==ntime).and.(fluxyr(n2)/=0)) then
      nnext=ntime
    else
      nnext=n2+1
    endif
    start(timx)=nnext; ncount(timx)=1
    call ccnf_get_vara(ncidfl,fluxid,start,ncount,fluxin_g(:,ntot+2))
!
!         need to make an array with the hour data in
    co2time(2:ntot+1)=fluxhr(n1:n2)
    timeinc = fluxhr(n1+1)-fluxhr(n1)
    co2time(1)=co2time(2)-timeinc
    co2time(ntot+2)=co2time(ntot+1)+timeinc
    deallocate(fluxhr)
  endif
!
  if (trim(fluxtype)=='daypulseon') then
!         read sunrise/sunset times for this month, region from file
    lc=len_trim(fluxname)
    read(fluxname(lc-2:lc),'(i3)') regnum
    if (regnum>nregion) then
      write(6,*) 'region number > nregion'
      call ccmpi_abort(-1)
    end if
    start(1)=regnum; start(2)=imon; start(3)=1
    ncount(1)=1; ncount(2)=1; ncount(3)=2
    call ccnf_get_vara(ncidfl,dayid,start,ncount,tracdaytime(igas,:))
  end if

  call ccnf_close(ncidfl)

  ! Should be more careful here, probably don't need the full range of
  ! the second dimension all the time.
  call ccmpi_distribute(fluxin,fluxin_g)

else ! myid /= 0
  call ccmpi_distribute(fluxin)
end if !myid == 0

!     Simple broadcast for co2 time
call ccmpi_bcast(co2time(1:nflux),0,comm_world)
!     Also need to share tracunit, if changed by input file.
call ccmpi_bcast(tracunit(igas),0,comm_world)
      
if (trim(fluxtype)=='daypulseon') then

  dum(1:2)=tracdaytime(igas,:)
  call ccmpi_bcast(dum(1:2),0,comm_world)
  tracdaytime(igas,:)=dum(1:2)

!       count number of timesteps that source emitting for
  kount=0
  do n=1,nperday
    hr = 24.*float(n)/float(nperday)
    if (tracdaytime(igas,1)<tracdaytime(igas,2) .and. tracdaytime(igas,1)<=hr .and. tracdaytime(igas,2)>=hr) kount=kount+1 
    if (tracdaytime(igas,1)>tracdaytime(igas,2) .and. (tracdaytime(igas,1)<=hr .or. tracdaytime(igas,2)>=hr)) kount=kount+1 
  enddo
!       scale flux to allow for emission over fraction of day
!       just set flux to zero if no daylight
  if (kount/=0) then
    fluxin = fluxin*float(nperday)/float(kount)
  else
    fluxin = 0.
  endif
endif

return
end subroutine readrco2

! *************************************************************************
subroutine readoh(imon,nfield,ohfile,varname,ohin)
! rml 16/2/10 New subroutine to read oh and strat loss for Transcom methane
use cc_mpi
use infile
use newmpar_m
use parm_m
implicit none
character(len=50) ohfile
character(len=13) varname
integer nfield,imon,ntime,ierr
integer nprev,nnext,ncur,n,timx,gridid
!     nflux =3 for month interp case - last month, this month, next month
real ohin(il*jl,kl,nfield)
integer ncidfl,fluxid
integer, dimension(:), allocatable :: ohmon
integer, dimension(4) :: start,ncount
real ohin_g(ifull_g,kl,nfield)
logical gridpts,tst

if ( myid == 0 ) then ! Read on this processor and then distribute
  write(6,*)'reading for ',imon,' from ',ohfile
  call ccnf_open(ohfile,ncidfl,ierr)
  call ncmsg("readoh",ierr)
  ! check for 1D or 2D formatted input
  call ccnf_inq_dimid(ncidfl,'gridpts',gridid,gridpts)
  call ccnf_inq_dimlen(ncidfl,'time',ntime)
  allocate(ohmon(ntime))
  call ccnf_get_vara(ncidfl,'month',ohmon)
!       check fluxname 
  call ccnf_inq_varid(ncidfl,varname,fluxid,tst)
!
!       find required records
  if (nfield==3) then
!         monthly case
    nprev=0
    nnext=0
    ncur=0
    if (ntime/=12) then
      write(6,*) 'oh/loss file wrong ntime'
      call ccmpi_abort(-1)
    end if
    do n=1,ntime
      if (ohmon(n)==imon) then
        ncur=n
        nprev=n-1
        if (nprev==0) nprev=12
        nnext=n+1
        if (nnext==13) nnext=1
      endif
    enddo

    if ( myid == 0 ) then
     write(6,*)'reading ',ncur,ohmon(ncur)
    end if
    if (ncur==0) then
      write(6,*) 'current month not in flux file'
      call ccmpi_abort(-1)
    end if

    ohin_g=0.
    if (gridpts) then
      start(1)=1 ; start(2)=1 ; start(3)=1
      ncount(1)=ifull_g; ncount(2)=kl ; ncount(3)=1 ! CHECK
      timx=3
    else
      start(1)=1 ; start(2)=1 ; start(3)=1 ; start(4)=1
      ncount(1)=il_g; ncount(2)=jl_g; ncount(3)=kl; ncount(4)=1
      timx=4
    end if
!         read preceeding month if needed
    if (nprev/=0) then
      start(timx)=nprev
      call ccnf_get_vara(ncidfl,fluxid,start,ncount,ohin_g(:,:,1))
    endif
!         read current month/year
    start(timx)=ncur
    call ccnf_get_vara(ncidfl,fluxid,start,ncount,ohin_g(:,:,2))
!         read next month
    if (nnext/=0) then
      start(timx)=nnext
      call ccnf_get_vara(ncidfl,fluxid,start,ncount,ohin_g(:,:,3))
    endif
  endif

  call ccnf_close(ncidfl)

  ! Should be more careful here, probably don't need the full range of
  ! the second dimension all the time.
  do n=1,nfield
    call ccmpi_distribute(ohin(:,:,n),ohin_g(:,:,n))
  enddo

else ! myid /= 0
  do n=1,nfield
    call ccmpi_distribute(ohin(:,:,n))
  enddo
end if !myid == 0

end subroutine readoh


   
! *************************************************************************
subroutine interp_tracerflux
!     interpolates tracer flux to current timestep if required
!     tracinterp 0 for no interpolation, 1 for monthly, 2 for daily/hourly
!     co2em123(:,1,:) contains prev month, co2em123(:,2,:) current month/year 
!     co2em123(:,3,:) next month
use cc_mpi, only : myid, ccmpi_abort
use dates_m, only : kdate
use newmpar_m
use parm_m

implicit none

integer, dimension(0:13) :: mdays
integer iyr,month,m1,m2,igas,igh
real a1,a2,a3,ratlm,ratlm2,hrmodel
real, dimension(ifull) :: c2,c3,c4
logical found
   
ratlm=0.
m1=0
m2=0

! this could go in the case section but then we'd do it for
! every monthly tracer.  This way we do it once but may not
! need it.
iyr=kdate/10000
month=(kdate-10000*iyr)/100
mdays=(/31,31,28,31,30,31,30,31,31,30,31,30,31,31/)
if (leap>=1) then
  if (mod(iyr,4)==0) mdays(2)=29
  if (mod(iyr,100)==0) mdays(2)=28
  if (mod(iyr,400)==0) mdays(2)=29
end if

do igas=1,numtracer
  select case(tracinterp(igas))

   ! constant
   case(0)
     co2em(:,igas)=co2em123(:,2,igas)

   ! monthly linear interpolation
   case(1)
     a1 = float(nperday*mdays(month-1))/2.
     a2 = float(nperday*mdays(month))/2.
     a3 = float(nperday*mdays(month+1))/2.
     if (ktau<a2) then
!            first half of month
       ratlm = (a1+float(ktau))/(a1+a2)
       m1 = 1; m2 = 2
     else
!            second half of month
       ratlm = (float(ktau)-a2)/(a2+a3)
       m1 = 2; m2 = 3
     endif
     co2em(:,igas) = (1.0-ratlm)*co2em123(:,m1,igas) + ratlm*co2em123(:,m2,igas)

   ! daily or hourly linear interpolation
   case(2)
       hrmodel = real(ktau-1)*dt/3600.
       igh=igashr(igas)
       if (ktau==1) nghr(igh)=1
       found=.false.
       if ( myid==0 ) then
         write(6,*) 'igas: ',igas,'igashr: ',igh
         write(6,*) 'hrmodel: ',hrmodel
       end if
       do while (.not.found)
         if ( myid==0 ) then
           write(6,*) nghr(igh),co2hr(nghr(igh),igh),co2hr(nghr(igh)+1,igh)
         end if
         if ((hrmodel.gt.co2hr(nghr(igh),igh)).and.(hrmodel.le.co2hr(nghr(igh)+1,igh))) then
           found = .true.
           ratlm2 = (hrmodel-co2hr(nghr(igh),igh))/(co2hr(nghr(igh)+1,igh)-co2hr(nghr(igh),igh))
           co2em(:,igas)= (1.-ratlm2)*co2emhr(:,nghr(igh),igh) +  ratlm2*co2emhr(:,nghr(igh)+1,igh)
         else
           nghr(igh) = nghr(igh)+1
!                this error check only useful for hourly resolution
           if (nghr(igh).eq.(31*24+2)) then
             write(6,*) 'hr flux error'
             call ccmpi_abort(-1)
           end if
         endif
       enddo

    ! monthly PWCB interpolation
    case(3)
     a2 = float(nperday*mdays(month))
     ratlm=a2/float(ktau)
     c2=co2em123(:,1,igas)
     c3=c2+co2em123(:,2,igas)
     c4=c3+co2em123(:,3,igas)
     co2em(:,igas)=.5*c3+(4.*c3-5.*c2-c4)*ratlm+1.5*(c4+3.*c2-3.*c3)*ratlm*ratlm

    ! same as case(1), but interpolate from start of month
    case(4)
     a2 = float(nperday*mdays(month))
     ratlm = float(ktau)/a2
     co2em(:,igas) = (1.-ratlm)*co2em123(:,2,igas) + ratlm*co2em123(:,3,igas)
  end select
enddo

! rml 16/2/10 interpolate oh and strat loss for methane cases
if (methane) then
!       monthly interpolation
  oh(:,:) = (1.0-ratlm)*oh123(:,:,m1) + ratlm*oh123(:,:,m2)
  strloss(:,:) = (1.0-ratlm)*strloss123(:,:,m1) + ratlm*strloss123(:,:,m2)
endif
if (mcf) then
  jmcf(:,:) = (1.0-ratlm)*jmcf123(:,:,m1) + ratlm*jmcf123(:,:,m2)
  mcfdep(:) = (1.0-ratlm)*mcfdep123(:,m1) + ratlm*mcfdep123(:,m2)
endif

return
end subroutine interp_tracerflux

! ***************************************************************************
subroutine tracer_mass
!     rml 16/10/03 check tracer mass - just write out for <= 6 tracers
use arrays_m   ! ps
use cc_mpi
use const_phys
use dates_m
use latlong_m
use newmpar_m
use sigs_m     ! dsig
use sumdd_m
use tracers_m  ! tr
use xyzinfo_m
implicit none
integer iq
real ltime

!     rml 14/5/10 code to create daily averages of afternoon concentrations
if (writetrpm) then
  do iq=1,il*jl
    ltime = timeg + rlongg(iq)*12./pi
    if (ltime>24) ltime = ltime - 24.
    if (ltime>=12.and.ltime<=15) then
      trpm(iq,1:kl,:) = trpm(iq,1:kl,:) + tr(iq,1:kl,:)
      npm(iq) = npm(iq) + 1
    endif
  enddo
endif

end subroutine tracer_mass

!subroutine tracini
!!     initial value now read from tracerlist 
!use infile
!use tracers_m
!implicit none
!integer i
!real in(il*jl,kl)
!
!do i=1,ngas
!! rml 15/11/06 facility to introduce new tracers to simulation
!! i.e. some read from restart, others initialised here
!  if (tracival(i).ne.-999) then
!! rml 16/2/10 addition for TC methane to get 3d initial condition
!      if (tracname(i)(1:7)=='methane') then
!!             read  initial condition
!        call ccnf_read('ch4in_cc48.nc','ch4in',in)
!        tr(1:il*jl,1:kl,i)=in
!      elseif (tracname(i)(1:3)=='mcf') then
!!             read  mcf initial condition
!        call ccnf_read('cmcfin_cc48.nc','mcfin',in)
!        tr(1:il*jl,1:kl,i)=in
!      else
!        tr(:,:,i)=tracival(i)
!      endif
!  endif
!enddo
!
!return
!end subroutine tracini

end module tracermodule
