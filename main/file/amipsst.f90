! Conformal Cubic Atmospheric Model
    
! Copyright 2015-2025 Commonwealth Scientific Industrial Research Organisation (CSIRO)
    
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

! Reads AMIP-style, monthly Sea Surface Temperatures, sea-ice and (optionally)
! near surface salinity.  Use onthefly for sub-monthly input data.
!
! namip=0   No input data
! namip=-1  Persisted SST anomalies
! namip=1   Use PWCB intepolation for SSTs, diagnose sea-ice
! namip=2   Use linear interpolation for SSTs and sea-ice (assumes pre-processing of monthly SSTs)
! namip=3   Use PWCB interpolation for SSTs and sea-ice equals supplied monthly value
! namip=4   Use PWCB interpolation for SSTs and sea-ice
! namip=5   Use PWCB interpolation for SSTs, sea-ice and salinity
! namip=11  Use JMc interpolation for SSTs, diagnose sea-ice
! namip=13  Use JMc interpolation for SSTs and sea-ice equals supplied monthly 
! namip=14  Use JMc interpolation for SSTs and sea-ice
! namip=15  Use JMc interpolation for SSTs, sea-ice and salinity
! namip=21  Use approx linear AMIP interpolation for SSTs and diagnose sea-ice
! namip=24  Use approx linear AMIP interpolation for SSTs and sea-ice
! namip=25  Use approx linear AMIP interpolation for SSTs, sea-ice and salinity

    
!     iday is a variable which is equal to the number of
!     day of the month.(Thus iday = 1 at the
!     start of a month, is updated to 2 at end of first day etc.)
!     imo is the number of the month we are in (1<=imo<=12)

module amipsst_m

implicit none

private
public amipsst, doints1, fill_cc1, define_grid, erase_grid, time_of_month
public time_interpolate

integer, save :: amip_mode = -1 ! 0=month (ASCII), 1=month (NetCDF), 2=day (NetCDF)
integer, save :: ncidx          ! Netcdf
integer, save :: ik
integer, dimension(:,:), allocatable :: nface4
real, dimension(:,:), allocatable :: xg4, yg4

contains
    
subroutine amipsst

use cc_mpi                ! CC MPI routines
use filnames_m            ! Filenames
use infile                ! Input file routines

implicit none

integer iernc
real timer_a, timer_b

if ( amip_mode==-1 ) then
  if ( myid==0 ) then
    call ccnf_open(sstfile,ncidx,iernc)
    if ( iernc==0 ) then
      call ccnf_get_vara(ncidx,'time',1,timer_a)
      call ccnf_get_vara(ncidx,'time',2,timer_b)
      if ( timer_b-timer_a <= 1440.01 ) then
        write(6,*) "Found daily data in sstfile for AMIPSST"  
        amip_mode = 2
      else
        write(6,*) "Found monthly data in sstfile for AMIPSST"  
        amip_mode = 1
      end if
    else
      write(6,*) "Assuming ASCII format monthly data for AMIPSST."
      amip_mode = 0
    end if    
  end if
  call ccmpi_bcast(amip_mode,0,comm_world)  
end if

select case(amip_mode)
  case(0,1)
    call amipsst_month
  case(2)
    call amipsst_day
  case default
    write(6,*) "ERROR: Cannot assign amip_mode to amipsst file"
    write(6,*) "amip_mode = ",amip_mode
    call ccmpi_abort(-1)
end select

return
end subroutine amipsst
    
subroutine amipsst_month

use arrays_m                                      ! Atmosphere dyamics prognostic arrays
use cc_mpi                                        ! CC MPI routines
use dates_m                                       ! Date data
use latlong_m                                     ! Lat/lon coordinates
use mlo_ctrl                                      ! Ocean physics control layer
use newmpar_m                                     ! Grid parameters
use nharrs_m, only : lrestart                     ! Non-hydrostatic atmosphere arrays
use parm_m                                        ! Model configuration
use pbl_m                                         ! Boundary layer arrays
use permsurf_m                                    ! Fixed surface arrays
use soil_m                                        ! Soil and surface data
use soilsnow_m                                    ! Soil, snow and surface data

implicit none

real, allocatable, save, dimension(:,:) :: ssta
real, allocatable, save, dimension(:,:) :: aice
real, allocatable, save, dimension(:,:) :: asal
real, allocatable, save, dimension(:) :: res
real, dimension(ifull) :: sssb, timelt, fraciceb
real, dimension(ifull) :: old, new, delta
real x, c2, c3, c4, rat1, rat2
real wgt
integer mdays
integer idjd_g, iq, k
integer, save :: iyr, imo, iday
integer, parameter :: curr_month = 3 ! centre of 5-points

if ( .not.allocated(ssta) ) then
  if ( myid==0 ) then
    write(6,*) 'Calling amipsst at beginning of run'
  end if
  allocate( ssta(ifull,5) )
  allocate( aice(ifull,5) )
  allocate( asal(ifull,5) )
  ssta(:,:) = 300.
  aice(:,:) = 0.
  asal(:,:) = 0.
else if ( mod(ktau,nperday)==0 ) then
  if ( myid==0 ) then
    write(6,*) "amipsst called at end of day for ktau,mtimer,namip ",ktau,mtimer,namip
  end if  
else
  ! call once per day for prescribed SSTs  
  if ( nmlo==0 ) return
end if

idjd_g = id + (jd-1)*il_g
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

call time_of_month(iyr,imo,iday,mdays,kdate,mtimer,leap)
if ( namip==-1 ) then
  iyr = 0
end if
x = real(iday-1)/real(mdays)  ! simplest at end of day

! load data for 5-point stencil.  Dummy data for unused months.
fraciceb = 0.  
if ( ktau==0 ) then
    
  if ( myid==0 ) then 
    call amiprd_month(ssta,aice,asal,iyr,imo,idjd_g)
  else
    call ccmpi_distribute(ssta(:,1:5))
    if ( namip==2 .or. (namip>=3.and.namip<=5) .or. (namip>=13.and.namip<=15) .or. &
         (namip>=24.and.namip<=25) ) then
      call ccmpi_distribute(aice(:,1:5))
    end if
    if ( namip==5 .or. namip==15 .or. namip==25 ) then
      call ccmpi_distribute(asal(:,1:5))
    end if
  end if ! myid==0
  
  ! checks
  if ( any(.not.land(1:ifull)) ) then
    do k = 1,5
      if ( minval(ssta(:,k),mask=.not.land(1:ifull))<0. ) then
        write(6,*) "ERROR: Invalid SSTs detected in AMIP file"
        write(6,*) "minval(sst) ",minval(ssta(:,k),mask=.not.land(1:ifull))
        call ccmpi_abort(-1)
      end if
      if ( maxval(ssta(:,k),mask=.not.land(1:ifull))>400. ) then
        write(6,*) "ERROR: Invalid SSTs detected in AMIP file"
        write(6,*) "maxval(sst) ",maxval(ssta(:,k),mask=.not.land(1:ifull))
        call ccmpi_abort(-1)
      end if
    end do
  end if
  
end if ! ktau==0



! namip=-1  Persisted SST anomalies
if ( namip==-1 ) then
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
  ! Each day interpolate-in-time non-land sst's    
  if ( ktau==0 ) then
    allocate(res(ifull))
    do iq = 1,ifull  
      if ( .not.land(iq) ) then
        c2 = ssta(iq,curr_month-1)
        c3 = c2 + ssta(iq,curr_month)
        c4 = c3 + ssta(iq,curr_month+1)
        res(iq) = tgg(iq,1) - ( .5*c3+(4.*c3-5.*c2-c4)*x+1.5*(c4+3.*c2-3.*c3)*x*x )
      end if      ! (.not.land(iq))
    end do
    if ( myid==0 ) then
      write(6,*)'some res values',(res(iq),iq=1,ifull,100)
    end if
  end if    ! (ktau==0)
  write(6,*) 'later_a ktau,res,tss ',ktau,res(idjd),tss(idjd)
  ! c1=0.
  do iq = 1,ifull  
    if ( .not.land(iq) ) then
      c2 = ssta(iq,curr_month-1)
      c3 = c2 + ssta(iq,curr_month)
      c4 = c3 + ssta(iq,curr_month+1)     
      tgg(iq,1) = res(iq) + ( .5*c3+(4.*c3-5.*c2-c4)*x+1.5*(c4+3.*c2-3.*c3)*x*x )
      tss(iq) = tgg(iq,1)
    end if      ! (.not.land(iq))
  end do
  write(6,*) 'later ktau,res,tss ',ktau,res(idjd),tss(iq)
  return
end if  ! (namip==-1)


! namip=1   Use PWCB intepolation for SSTs, diagnose sea-ice
if ( namip==1 ) then
  ! c1=0.
  fraciceb = 0.
  do iq = 1,ifull  
    if ( .not.land(iq) ) then
      c2 = ssta(iq,curr_month-1)
      c3 = c2 + ssta(iq,curr_month)
      c4 = c3 + ssta(iq,curr_month+1)          
      tgg(iq,1) = .5*c3+(4.*c3-5.*c2-c4)*x+1.5*(c4+3.*c2-3.*c3)*x*x
      if ( tgg(iq,1)<272. ) fraciceb(iq) = 1.
    end if      ! (.not.land(iq))
  end do
end if  ! (namip==1)


! namip=2   Use linear interpolation for SSTs and sea-ice (assumes pre-processing of monthly SSTs)
if ( namip==2 ) then
  write(6,*) "ERROR: Not currently supported"
  stop
!  !--------------------------------------------------------------------------------------------------
!  ! Linear interpolation (possibly pre-processed by AMIP tridiagonal matrix method)
!  if ( iday<mdays(imo)/2 ) then  ! 1st half of month
!    rat1 = (mdays(imo)-2.*iday)/(mdays(imo)+mdays(imo-1))
!    rat2 = (2.*iday+mdays(imo-1))/(mdays(imo)+mdays(imo-1))
!    if ( mydiag ) write(6,*)'rat1,rat2,land: ',rat1,rat2,land(idjd)
!    do iq = 1,ifull  
!      if ( .not.land(iq) ) then
!        tgg(iq,1) = rat1*ssta(iq,curr_month-1) + rat2*ssta(iq,curr_month)  ! sea water temperature
!      end if      ! (.not.land(iq))
!    end do
!    fraciceb(:) = min( .01*(rat1*aice(:,curr_month-1)+rat2*aice(:,curr_month)), 1. ) ! convert from %
!  else                             ! 2nd half of month
!    rat1 = (mdays(imo+1)+2.*mdays(imo)-2.*iday)/(mdays(imo+1)+mdays(imo))
!    rat2 = (2.*iday-mdays(imo))/(mdays(imo+1)+mdays(imo))
!    do iq = 1,ifull  
!      if ( .not.land(iq) ) then
!        tgg(iq,1) = rat1*ssta(iq,curr_month) + rat2*ssta(iq,curr_month+1)  ! sea water temperature
!      end if      ! (.not.land(iq))
!    end do
!    fraciceb(:) = min( .01*(rat1*aice(:,curr_month)+rat2*aice(:,curr_month+1)), 1. ) ! convert from %
!  end if
end if  ! (namip==2)


! namip=3   Use PWCB interpolation for SSTs and sea-ice equals supplied monthly value
! namip=4   Use PWCB interpolation for SSTs and sea-ice
! namip=5   Use PWCB interpolation for SSTs, sea-ice and salinity
if ( namip==3 .or. namip==4 .or. namip==5 ) then
  !--------------------------------------------------------------------------------------------------
  ! Piece-wise cubic bessel interpolation
  call time_interpolate(tgg(:,1),ssta(:,curr_month-2),ssta(:,curr_month-1),           &
                        ssta(:,curr_month),ssta(:,curr_month+1),ssta(:,curr_month+2), &
                        x,namip,land=land)
endif  ! (namip==3 .or. namip==4 .or. namip==5 )
if ( namip==3 ) then
  do iq=1,ifull  
    fraciceb(iq)=min(.01*aice(iq,curr_month),1.)
  enddo
else if ( namip==4 .or. namip==5 ) then
  call time_interpolate(fraciceb(:),aice(:,curr_month-2),aice(:,curr_month-1),         &
                        aice(:,curr_month),aice(:,curr_month+1),aice(:,curr_month+2),  &
                        x,namip,land=land)
  fraciceb = min( fraciceb, 1. )
endif  ! (namip==4 .or namip==5 )
if ( namip==5 ) then
  call time_interpolate(sssb(:),asal(:,curr_month-2),asal(:,curr_month-1),            &
                        asal(:,curr_month),asal(:,curr_month+1),asal(:,curr_month+2), &
                        x,namip,land=land)
  sssb = max( sssb, 0. )
end if ! namip==5


!--------------------------------------------------------------------------------------------------
! ERROR checking for invalid options
if ( (namip>5.and.namip<11) .or. namip==12 ) then
  write(6,*) "ERROR: invalid namip option ",namip
  call ccmpi_abort(-1)
end if


! namip=11  Use JMc interpolation for SSTs, diagnose sea-ice
! namip=13  Use JMc interpolation for SSTs and sea-ice equals supplied monthly 
! namip=14  Use JMc interpolation for SSTs and sea-ice
! namip=15  Use JMc interpolation for SSTs, sea-ice and salinity
if ( namip==11 .or. namip==13 .or. namip==14 .or. namip==15 ) then
  !--------------------------------------------------------------------------------------------------
  ! John McGregor 5-pt, piece-wise, cubic interpolation
  call time_interpolate(tgg(:,1),ssta(:,curr_month-2),ssta(:,curr_month-1),           &
                        ssta(:,curr_month),ssta(:,curr_month+1),ssta(:,curr_month+2), &
                        x,namip,land=land)
  if ( namip==11 ) then
    where ( tgg(1:ifull,1)<271.2 )
      fraciceb(1:ifull) = 1.
    elsewhere
      fraciceb(1:ifull) = 0.
    end where
  else if ( namip==13 ) then
    fraciceb(1:ifull)=min(.01*aice(1:ifull,curr_month),1.)
  else if ( namip==14 .or. namip==15 ) then
    call time_interpolate(fraciceb(:),aice(:,curr_month-2),aice(:,curr_month-1),         &
                          aice(:,curr_month),aice(:,curr_month+1),aice(:,curr_month+2),  &
                          x,namip,land=land)
    fraciceb = min( fraciceb, 1. )
  end if
  if ( namip==15 ) then
    call time_interpolate(sssb(:),asal(:,curr_month-2),asal(:,curr_month-1),            &
                          asal(:,curr_month),asal(:,curr_month+1),asal(:,curr_month+2), &
                          x,namip,land=land)
    sssb = max( sssb, 0. )
  end if  
end if


!--------------------------------------------------------------------------------------------------
! ERROR checking for invalid options
if ( (namip>15.and.namip<21).or.(namip>21.and.namip<24) ) then
  write(6,*) "ERROR: invalid namip option ",namip
  call ccmpi_abort(-1)
end if


! namip=21  Use approx linear AMIP interpolation for SSTs and diagnose sea-ice
! namip=24  Use approx linear AMIP interpolation for SSTs and sea-ice
! namip=25  Use approx linear AMIP interpolation for SSTs, sea-ice and salinity
if ( namip==21 .or. namip==24 .or. namip==25 ) then
  !--------------------------------------------------------------------------------------------------
  ! Approximation of piece-wise, linear AMIP interpolation
  call time_interpolate(tgg(:,1),ssta(:,curr_month-2),ssta(:,curr_month-1),           &
                        ssta(:,curr_month),ssta(:,curr_month+1),ssta(:,curr_month+2), &
                        x,namip,land=land)
  if ( namip==21 ) then
    where ( tgg(1:ifull,1)<271.2 )
      fraciceb(1:ifull) = 1.
    elsewhere
      fraciceb(1:ifull) = 0.
    end where
  else if ( namip==24 .or. namip==25 ) then
    call time_interpolate(fraciceb(:),aice(:,curr_month-2),aice(:,curr_month-1),         &
                          aice(:,curr_month),aice(:,curr_month+1),aice(:,curr_month+2),  &
                          x,namip,land=land)
    fraciceb = min( fraciceb, 1. )
  end if
  if ( namip==25 ) then
    call time_interpolate(sssb(:),asal(:,curr_month-2),asal(:,curr_month-1),            &
                          asal(:,curr_month),asal(:,curr_month+1),asal(:,curr_month+2), &
                          x,namip,land=land)
    sssb = max( sssb, 0. )
  end if      
end if


!--------------------------------------------------------------------------------------------------
! Remove small sea-ice fractions
where ( fraciceb(1:ifull)<=0.02 )
  fraciceb(1:ifull) = 0.
end where


!--------------------------------------------------------------------------------------------------
! Sea-ice and Sea-Surface-Temperature
if ( nmlo==0 ) then
  sicedep(:) = 0. 
  if ( ktau==0 .and. .not.lrestart ) then  ! will set sicedep in indata
    ! This case is for poor initial conditions ifile  
    fracice(:) = fraciceb(:)
    where ( .not.land(1:ifull) )
      tss(:) = tgg(:,1)
    end where
  else
    do iq = 1,ifull
      if ( .not.land(iq) ) then
        if ( fraciceb(iq)>0. ) then
          if ( fracice(iq)<1.e-20 ) then
            ! create values for tice, and set averaged tss
            ! N.B. if already a sice point, keep present tice
            !tggsn(iq,1)=min(271.2,tss(iq),t(iq,1)+.04*6.5) ! for 40 m lev1
            tggsn(iq,1) = 271.2
          end if  ! (fracice(iq)==0.)
          if ( rlatt(iq)>0. ) then
            sicedep(iq) = 2.
          else
            sicedep(iq) = 1.
          endif ! (rlatt(iq)>0.)
        endif   ! (fraciceb(iq)>0.)
        fracice(iq) = fraciceb(iq)
        tss(iq) = tggsn(iq,1)*fracice(iq) + tgg(iq,1)*(1.-fracice(iq))
      endif      ! (.not.land(iq))
    enddo
  end if  
elseif ( ktau>0 ) then
  if ( nud_hrs<=0 ) then
    write(6,*) "ERROR: namip/=0 has been selected with in-line ocean model (nmlo/=0)"  
    write(6,*) "nud_hrs>0 must be specified for relaxiation of SSTs"
    call ccmpi_abort(-1)
  end if
  old = tgg(:,1)
  call mloexpmelt(timelt)
  timelt = min(timelt+0.01,tgg(:,1))
  new = tgg(:,1)*(1.-fraciceb(:)) + timelt(:)*fraciceb(:) - wrtemp
  wgt = dt/real(nud_hrs*3600)
  call mloexport("temp",old,1,0)
  delta = new - old
  do k = 1,kbotmlo
    old = 273.16 - wrtemp  
    call mloexport("temp",old,k,0)
    old = old + max( wgt*delta, 271.16-wrtemp-old )
    call mloimport("temp",old,k,0)
  end do  
  do k = 1,ms
    call mloexport("temp",tgg(:,k),k,0)
    where ( tgg(:,k)<100. )
      tgg(:,k) = tgg(:,k) + wrtemp
    end where    
  end do
end if ! if (nmlo==0) ..else..


return
end subroutine amipsst_month
      
    
subroutine amiprd_month(ssta,aice,asal,iyr,imo,idjd_g)
      
use cc_mpi                ! CC MPI routines
use const_phys            ! Physical constants
use filnames_m            ! Filenames
use infile                ! Input file routines
use latlong_m             ! Lat/lon coordinates
use newmpar_m             ! Grid parameters
use parm_m                ! Model configuration
use parmgeom_m            ! Coordinate data

implicit none
      
integer, parameter :: nihead = 54
integer, parameter :: nrhead = 14
      
integer, intent(in) :: iyr, imo, idjd_g
integer imonth, iyear, il_in, jl_in, iyr_m, imo_m, ierr, leap_file
integer varid, iarchx, maxarchi, lsmid, varid_time
integer kdate_r, ktime_r
integer kdate_rsav, ktime_rsav
integer iq, mm, k
integer, dimension(3) :: spos, npos
#ifdef i8r8
integer, dimension(nihead) :: nahead
#else
integer(kind=4), dimension(nihead) :: nahead
#endif
integer(kind=8) mtimer_r
real, dimension(ifull,5), intent(out) :: ssta
real, dimension(ifull,5), intent(out) :: aice
real, dimension(ifull,5), intent(out) :: asal
real, dimension(:,:), allocatable :: ssta_g, ssta_l
real, dimension(:), allocatable :: lsmr_g
real, dimension(nrhead) :: ahead
real rlon_in, rlat_in, schmidt_in
real of, sc
real timer_r
logical ltest, tst, interpolate
logical, dimension(:), allocatable :: lsma_g
character(len=22) header
character(len=10) unitstr
character(len=10) calendarstring
character(len=80) datestring

ssta = 300.
aice = 0.
asal = 0.

! check for netcdf file format
if ( amip_mode==1 ) then

  iyr_m = iyr
  imo_m = imo
    
  ! NETCDF
  write(6,*) "Reading AMIP file in netcdf format"
  ! check grid definition
  call ccnf_get_attg(ncidx,'int_header',nahead)
  call ccnf_get_attg(ncidx,'real_header',ahead)
  il_in      = nahead(1)
  jl_in      = nahead(2)
  rlon_in    = ahead(5)
  rlat_in    = ahead(6)
  schmidt_in = ahead(7)
  if ( schmidt_in<=0. .or. schmidt_in>1. ) then
    rlon_in    = ahead(6)
    rlat_in    = ahead(7)
    schmidt_in = ahead(8)
  endif  ! (schmidtx<=0..or.schmidtx>1.)  
  call ccnf_get_attg(ncidx,'leap',leap_file,tst) ! old method
  if ( tst ) leap_file = 0                       ! old method
  call ccnf_get_att(ncidx,'time','calendar',calendarstring,ierr)                      ! new method
  if ( ierr==0 ) then                                                                 ! new method
    select case( trim(calendarstring) )                                               ! new method
      case("365_day","noleap")                                                        ! new method
        leap_file = 0                                                                 ! new method
      case("standard")                                                                ! new method
        leap_file = 1                                                                 ! new method
      case("360_day")                                                                 ! new method
        leap_file = 2                                                                 ! new method
      case default                                                                    ! new method
        write(6,*) "ERROR: Unknown calendar in AMIPSST file"                          ! new method
        write(6,*) "Found calendar = ",trim(calendarstring)                           ! new method
        write(6,*) "Valid options are 365_day, noleap, 360_day or standard"           ! new method
    end select                                                                        ! new method
  else                                                                                ! new method
    write(6,*) "ERROR: calendar is not defined in AMIPSST file"                       ! new method
    write(6,*) "Please define calendar with standard, 365_day, noleap or 360_day"     ! new method 
    call ccmpi_abort(-1)                                                              ! new method 
  end if                                                                              ! new method
      
  call ccnf_inq_dimlen(ncidx,'time',maxarchi)
         
  interpolate = ( il_g/=il_in .or. jl_g/=jl_in .or. abs(rlong0-rlon_in)>1.e-6 .or. abs(rlat0-rlat_in)>1.e-6 .or. &
       abs(schmidt-schmidt_in)>0.0002 )
  ik = il_in
  
  if ( interpolate ) then
    write(6,*) "Interpolation is required for AMIPSST"
    write(6,*) "leap_file = ",leap_file
    
    call define_grid(ik,rlon_in,rlat_in,schmidt_in)
    
    allocate( ssta_g(ik*ik*6,5), ssta_l(ifull_g,5), lsma_g(ik*ik*6) )
    lsma_g = .false.
    ssta_g = 0.
    ssta_l = 0.
    npos(1) = ik
    npos(2) = 6*ik
    spos(1:2) = 1
    allocate( lsmr_g(ik*ik*6) )  
    call ccnf_get_vara(ncidx,'lsm',spos(1:2),npos(1:2),lsmr_g(:))
    lsma_g = lsmr_g >= 0.5
    deallocate( lsmr_g )
  else
    write(6,*) "No interpolation is required for AMIPSST"
    write(6,*) "leap_file = ",leap_file
    allocate( ssta_g(ifull_g,5), ssta_l(0,0), lsma_g(0) )
    ssta_g = 0.
  end if
  
  ! search for required month
  iarchx = 0
  iyear  = -999
  imonth = -999
  ltest = .true.
  call ccnf_get_att(ncidx,'time','units',datestring)
  call processdatestring(datestring,kdate_rsav,ktime_rsav)
  do while ( ltest .and. iarchx<maxarchi )
    iarchx = iarchx + 1
    kdate_r = kdate_rsav
    ktime_r = ktime_rsav
    call ccnf_get_vara(ncidx,'time',iarchx,timer_r)
    mtimer_r = nint(timer_r,8)
    call datefix(kdate_r,ktime_r,mtimer_r,allleap=leap_file)
    iyear  = kdate_r/10000
    imonth = (kdate_r-iyear*10000)/100
    ltest  = iyr_m/=iyear .or. imo_m/=imonth
  end do    
  if ( ltest ) then
    write(6,*) "ERROR: Cannot locate year ",iyr_m
    write(6,*) "       and month ",imo_m
    write(6,*) "       in file ",trim(sstfile)
    call ccmpi_abort(-1)
  end if
  unitstr = ''
  call ccnf_get_att(ncidx,'tos','units',unitstr)
  write(6,*) "Reading SST data from amipsst file"
  if ( trim(unitstr)=='C' ) then
    write(6,*) "Converting SSTs from degC to degK"  
  end if
  npos(1) = ik
  npos(2) = 6*ik
  npos(3) = 1    
  spos(1:2) = 1
  spos(3) = max( iarchx - 2, 1 )   
  if ( namip>=11 .or. namip<=15 ) then
    if ( spos(3)>iarchx-2 .and. myid==0 ) then
      write(6,*) "Warning: Using current SSTs for previous month(s)"
    end if
    call ccnf_get_vara(ncidx,'tos',spos,npos,ssta_g(:,1))
    call ccnf_get_att(ncidx,'tos','add_offset',of,ierr=ierr)
    if ( ierr /= 0 ) of=0.
    call ccnf_get_att(ncidx,'tos','scale_factor',sc,ierr=ierr)
    if ( ierr /= 0 ) sc=1.
    ssta_g(:,1)=sc*ssta_g(:,1)+of
    if ( trim(unitstr) == 'C' ) ssta_g(:,1) = ssta_g(:,1) + 273.16
  else
    ssta_g(:,1)=0. ! dummy.  Should not be used.
  end if
  spos(3) = max( iarchx - 1, 1 )
  if ( namip==1 .or. namip==2 .or. (namip>=3.and.namip<=5) .or. (namip>=11.and.namip<=15) .or. &
       (namip>=21.and.namip<=25) ) then
    call ccnf_get_vara(ncidx,'tos',spos,npos,ssta_g(:,2))
    ssta_g(:,2)=sc*ssta_g(:,2)+of
    if ( trim(unitstr) == 'C' ) ssta_g(:,2) = ssta_g(:,2) + 273.16
  else
    ssta_g(:,2)=0. ! dummy.  Should not be used.       
  end if    
  spos(3) = iarchx
  call ccnf_get_vara(ncidx,'tos',spos,npos,ssta_g(:,3))
  ssta_g(:,3)=sc*ssta_g(:,3)+of  
  if ( trim(unitstr) == 'C' ) ssta_g(:,3) = ssta_g(:,3) + 273.16
  spos(3) = min( iarchx + 1, maxarchi )
  call ccnf_get_vara(ncidx,'tos',spos,npos,ssta_g(:,4))
  ssta_g(:,4)=sc*ssta_g(:,4)+of  
  if ( trim(unitstr) == 'C' ) ssta_g(:,4) = ssta_g(:,4) + 273.16
  spos(3) = min( iarchx + 2, maxarchi )
  if ( (namip>=11.and.namip<=15) .or. (namip>=21.and.namip<=25) ) then
    if ( spos(3)<iarchx+2 .and. myid==0 ) then
      write(6,*) "Warning: Using current SSTs for next month(s)"
    end if
    call ccnf_get_vara(ncidx,'tos',spos,npos,ssta_g(:,5))
    ssta_g(:,5)=sc*ssta_g(:,5)+of 
    if ( trim(unitstr) == 'C' ) ssta_g(:,5) = ssta_g(:,5) + 273.16
  else
    ssta_g(:,5)=0. ! dummy.  Should not be used.       
  end if 
  
  if ( interpolate ) then
    do k = 1,5  
      call fill_cc1(ssta_g(:,k),lsma_g)
      call doints1(ssta_g(:,k),ssta_l(:,k))
    end do  
    call ccmpi_distribute(ssta(:,1:5), ssta_l(:,1:5))
  else
    call ccmpi_distribute(ssta(:,1:5), ssta_g(:,1:5))
  end if  

  if ( namip==2 .or. namip==3 .or. namip==4 .or. namip==5 .or. namip==13 .or. &
       namip==14 .or. namip==15 .or. namip==24 .or. namip==25 ) then   ! sice also read at middle of month

    ! NETCDF
    write(6,*) "Reading Sea Ice data from amipsst file"
    spos(3)=max( iarchx - 2, 1 )
    if ( namip>=11 .or. namip<=15 ) then
      if ( spos(3)>iarchx-2 .and. myid==0 ) then
        write(6,*) "Warning: Using current sea-ice for previous month(s)" 
      end if
      call ccnf_get_vara(ncidx,'sic',spos,npos,ssta_g(:,1))
      call ccnf_get_att(ncidx,'sic','add_offset',of,ierr=ierr)
      if (ierr/=0) of=0.
      call ccnf_get_att(ncidx,'sic','scale_factor',sc,ierr=ierr)
      if (ierr/=0) sc=1.
      ssta_g(:,1)=sc*ssta_g(:,1)+of
      ssta_g(:,1)=100.*ssta_g(:,1)  
    else
      ssta_g(:,1)=0. ! dummy.  Should not be used.
    end if
    spos(3)=max( iarchx - 1, 1 )
    if ( namip==1 .or. namip==2 .or. (namip>=3.and.namip<=5) .or. (namip>=11.and.namip<=15) .or. &
       (namip>=21.and.namip<=25) ) then
      call ccnf_get_vara(ncidx,'sic',spos,npos,ssta_g(:,2))
      ssta_g(:,2)=sc*ssta_g(:,2)+of
      ssta_g(:,2)=100.*ssta_g(:,2)       
    else
      ssta_g(:,2)=0. ! dummy.  Should not be used.      
    end if
    spos(3)=iarchx
    call ccnf_get_vara(ncidx,'sic',spos,npos,ssta_g(:,3))
    ssta_g(:,3)=sc*ssta_g(:,3)+of
    ssta_g(:,3)=100.*ssta_g(:,3)       
    spos(3)=min( iarchx + 1, maxarchi )
    call ccnf_get_vara(ncidx,'sic',spos,npos,ssta_g(:,4))
    ssta_g(:,4)=sc*ssta_g(:,4)+of
    ssta_g(:,4)=100.*ssta_g(:,4)       
    spos(3)=min( iarchx + 2, maxarchi )
    if ( (namip>=11.and.namip<=15) .or. (namip>=21.and.namip<=25) ) then
      if ( spos(3)<iarchx+2 .and. myid==0 ) then
        write(6,*) "Warning: Using current sea-ice for next month(s)"
      end if
      call ccnf_get_vara(ncidx,'sic',spos,npos,ssta_g(:,5))
      ssta_g(:,5)=sc*ssta_g(:,5)+of
      ssta_g(:,5)=100.*ssta_g(:,5)       
    else
      ssta_g(:,5)=0. ! dummy.  Should not be used.    
    end if

    if ( interpolate ) then
      do k = 1,5  
        call fill_cc1(ssta_g(:,k),lsma_g)
        call doints1(ssta_g(:,k),ssta_l(:,k))
      end do  
      call ccmpi_distribute(aice(:,1:5), ssta_l(:,1:5))
    else
      call ccmpi_distribute(aice(:,1:5), ssta_g(:,1:5))
    end if  
    
  end if  

  if ( namip==5 .or. namip==15 .or. namip==25 ) then ! salinity also read
      
    ! NETCDF
    write(6,*) "Reading Salinity data from amipsst file"
    spos(3)=max( iarchx - 2, 1 )
    if ( namip>=11 .or. namip<=15 ) then
      if ( spos(3)>iarchx-2 .and. myid==0 ) then
        write(6,*) "Warning: Using current salinity for previous month(s)"
      end if
      call ccnf_get_vara(ncidx,'sss',spos,npos,ssta_g(:,1))
      call ccnf_get_att(ncidx,'sss','add_offset',of,ierr=ierr)
      if (ierr/=0) of=0.
      call ccnf_get_att(ncidx,'sss','scale_factor',sc,ierr=ierr)
      if (ierr/=0) sc=1.  
      ssta_g(:,1)=sc*ssta_g(:,1)+of
    else
      ssta_g(:,1)=0. ! dummy.  Should not be used.
    end if
    spos(3)=max( iarchx - 1, 1 )
    if ( namip==1 .or. namip==2 .or. (namip>=3.and.namip<=5) .or. (namip>=11.and.namip<=15) .or. &
       (namip>=21.and.namip<=25) ) then
      call ccnf_get_vara(ncidx,'sss',spos,npos,ssta_g(:,2))
      ssta_g(:,2)=sc*ssta_g(:,2)+of
    else
      ssta_g(:,2)=0. ! dummy.  Should not be used.      
    end if
    spos(3)=iarchx
    call ccnf_get_vara(ncidx,'sss',spos,npos,ssta_g(:,3))
    ssta_g(:,3)=sc*ssta_g(:,3)+of
    spos(3)=min( iarchx + 1, maxarchi )
    call ccnf_get_vara(ncidx,'sss',spos,npos,ssta_g(:,4))
    ssta_g(:,4)=sc*ssta_g(:,4)+of
    spos(3)=min( iarchx + 2, maxarchi )
    if ( (namip>=11.and.namip<=15) .or. (namip>=21.and.namip<=25) ) then
      if ( spos(3)<iarchx+2 .and. myid==0 ) then
        write(6,*) "Warning: Using current salinity for next month"
      end if
      call ccnf_get_vara(ncidx,'sss',spos,npos,ssta_g(:,5))
      ssta_g(:,5)=sc*ssta_g(:,5)+of
    else
      ssta_g(:,5)=0. ! dummy.  Should not be used.    
    end if
    
    if ( interpolate ) then
      do k = 1,5  
        call fill_cc1(ssta_g(:,k),lsma_g)
        call doints1(ssta_g(:,k),ssta_l(:,k))
      end do  
      call ccmpi_distribute(asal(:,1:5), ssta_l(:,1:5))
    else
      call ccmpi_distribute(asal(:,1:5), ssta_g(:,1:5))
    end if  

  end if

  call ccnf_close(ncidx)
  deallocate( ssta_g, ssta_l, lsma_g )
  if ( interpolate ) then
    call erase_grid
  end if  
  
else
   
  allocate( ssta_g(ifull_g,5) )
  ssta_g = 0.
 
  ! ASCII  
  if ( (namip>=11.and.namip<=15) .or. (namip>=21.and.namip<=25) ) then  
    iyr_m = iyr
    imo_m = imo - 2
  else if ( namip==1 .or. namip==2 .or. (namip>=3.and.namip<=5) ) then  
    iyr_m = iyr
    imo_m = imo - 1
  else if ( namip==-1 ) then
    iyr_m = iyr
    imo_m = imo
  else
    write(6,*) "ERROR: Invalid namip option ",namip
    call ccmpi_abort(-1)
  end if
    
  if ( imo_m<1 ) then
    imo_m = imo_m + 12
    iyr_m = iyr - 1
  endif
  
  if ( namip==-1 ) iyr_m = 0

  open(unit=75,file=sstfile,status='old',form='formatted',iostat=ierr)
  if (ierr/=0) then
    write(6,*) "ERROR: Cannot read AMIP sstfile ",trim(sstfile)
    call ccmpi_abort(-1)
  end if
  write(6,*) "Reading AMIP file in ASCII format"
  iyear=-999
  imonth=-999
  do while(iyr_m/=iyear.or.imo_m/=imonth)
    write(6,*) 'about to read amipsst file'
    read(75,*) imonth,iyear,il_in,jl_in,rlon_in,rlat_in,schmidt_in,header
    write(6,'("reading sst ",i2,i5,2i4,2f6.1,f7.4,1x,a22)') imonth,iyear,il_in,jl_in,rlon_in,rlat_in,schmidt_in,header
    if(il_g/=il_in.or.jl_g/=jl_in.or.abs(rlong0-rlon_in)>1.e-20.or.abs(rlat0-rlat_in)>1.e-20.or. &
        abs(schmidt-schmidt_in)>1.e-20)then
      write(6,*) 'il_g,il_in,jl_g,jl_in,rlong0,rlon_in',il_g,il_in,jl_g,jl_in,rlong0,rlon_in
      write(6,*) 'rlat0,rlat_in,schmidt,schmidt_in',rlat0,rlat_in,schmidt,schmidt_in
      write(6,*) 'wrong amipsst file'
      call ccmpi_abort(-1)
    endif
    read(75,*) ssta_g(:,1)
    ssta_g(:,1)=ssta_g(:,1)*.01 - 50. + 273.16
    write(6,*) 'want imo_m,iyr_m; ssta ',imo_m,iyr_m,ssta_g(idjd_g,1)
  end do
  if ( (namip>=11.and.namip<=15) .or. (namip>=21.and.namip<=25) ) then
    read(75,'(i2,i5,a22)') imonth,iyear,header
    write(6,*) 'reading sstb data:',imonth,iyear,header
    read(75,*) ssta_g(:,2)
    ssta_g(:,2)=ssta_g(:,2)*.01 -50. +273.16
  else
    ssta_g(:,2) = ssta_g(:,1)
  end if
  write(6,*) 'sstb(idjd) ',ssta_g(idjd_g,2)  
  if ( namip==1 .or. namip==2 .or. (namip>=3.and.namip<=5) .or. (namip>=11.and.namip<=15) .or. &
       (namip>=21.and.namip<=25) ) then
    read(75,'(i2,i5,a22)') imonth,iyear,header
    write(6,*) 'reading sstc data:',imonth,iyear,header
    write(6,*) 'should agree with imo,iyr ',imo,iyr
    if(iyr/=iyear.or.imo/=imonth)then
      call ccmpi_abort(-1)
    end if
    read(75,*) ssta_g(:,3)
    ssta_g(:,3)=ssta_g(:,3)*.01 -50. +273.16
  end if
  write(6,*) 'sstc(idjd) ',ssta_g(idjd_g,3)       
  read(75,'(i2,i5,a22)') imonth,iyear,header
  write(6,*) 'reading sstd data:',imonth,iyear,header
  read(75,*) ssta_g(:,4)
  ssta_g(:,4)=ssta_g(:,4)*.01 -50. +273.16
  write(6,*) 'sstd(idjd) ',ssta_g(idjd_g,4)
  if ( (namip>=11.and.namip<=15) .or. (namip>=21.and.namip<=25) ) then
    read(75,'(i2,i5,a22)') imonth,iyear,header
    write(6,*) 'reading sste data:',imonth,iyear,header
    read(75,*) ssta_g(:,5)
    ssta_g(:,5)=ssta_g(:,5)*.01 -50. +273.16
  else
    ssta_g(:,5) = ssta_g(:,4)
  end if
  write(6,*) 'sste(idjd) ',ssta_g(idjd_g,5)
  close(75)
  call ccmpi_distribute(ssta(:,1:5), ssta_g(:,1:5))

  if ( namip==2 .or. namip==3 .or. namip==4 .or. namip==5 .or. namip==13 .or. &
       namip==14 .or. namip==15 .or. namip==24 .or. namip==25 ) then   ! sice also read at middle of month

    ! ASCII
    open(unit=76,file=icefile,status='old',form='formatted',iostat=ierr)
    if (ierr/=0) then
      write(6,*) "ERROR: Cannot read AMIP icefile ",trim(icefile)
      call ccmpi_abort(-1)
    end if
    iyear=-999
    imonth=-999
    do while(iyr_m/=iyear.or.imo_m/=imonth)
      read(76,*) imonth,iyear,il_in,jl_in,rlon_in,rlat_in,schmidt_in,header
      write(6,'("reading ice ",i2,i5,2i4,2f6.1,f6.3,a22)') imonth,iyear,il_in,jl_in,rlon_in,rlat_in,schmidt_in,header
      if(il_g/=il_in.or.jl_g/=jl_in.or.abs(rlong0-rlon_in)>1.e-20.or.abs(rlat0-rlat_in)>1.e-20.or. &
          abs(schmidt-schmidt_in)>1.e-20)then
        write(6,*) 'il_g,il_in,jl_g,jl_in,rlong0,rlon_in',il_g,il_in,jl_g,jl_in,rlong0,rlon_in
        write(6,*) 'rlat0,rlat_in,schmidt,schmidt_in',rlat0,rlat_in,schmidt,schmidt_in
        write(6,*) 'wrong amipice file'
        call ccmpi_abort(-1)
      endif
      read(76,*) ssta_g(:,1)
      write(6,*) 'want imo_m,iyr_m; aice ',imo_m,iyr_m,ssta_g(idjd_g,1)
    end do
    if ( (namip>=11.and.namip<=15) .or. (namip>=21.and.namip<=25) ) then
      read(76,'(i2,i5,a22)') imonth,iyear,header
      write(6,*) 'reading b_sice data:',imonth,iyear,header
      read(76,*) ssta_g(:,2)
    else
      ssta_g(:,2) = ssta_g(:,1)
    end if
    write(6,*) 'bice(idjd) ',ssta_g(idjd_g,2)    
    if ( namip==1 .or. namip==2 .or. (namip>=3.and.namip<=5) .or. (namip>=11.and.namip<=15) .or. &
       (namip>=21.and.namip<=25) ) then
      read(76,'(i2,i5,a22)') imonth,iyear,header
      write(6,*) 'reading c_sice data:',imonth,iyear,header
      write(6,*) 'should agree with imo,iyr ',imo,iyr
      if(iyr/=iyear.or.imo/=imonth) then
        call ccmpi_abort(-1)
      end if
      read(76,*) ssta_g(:,3)
    end if
    write(6,*) 'cice(idjd) ',ssta_g(idjd_g,3)       
    read(76,'(i2,i5,a22)') imonth,iyear,header
    write(6,*) 'reading d_sice data:',imonth,iyear,header
    read(76,*) ssta_g(:,4)
    write(6,*) 'dice(idjd) ',ssta_g(idjd_g,4)
    if ( (namip>=11.and.namip<=15) .or. (namip>=21.and.namip<=25) ) then
      read(76,'(i2,i5,a22)') imonth,iyear,header
      write(6,*) 'reading e_sice data:',imonth,iyear,header
      read(76,*) ssta_g(:,5)
    else
      ssta_g(:,5) = ssta_g(:,4)
    end if
    write(6,*) 'eice(idjd) ',ssta_g(idjd_g,5)    
    close(76)
    call ccmpi_distribute(aice(:,1:5), ssta_g(:,1:5))
 
  end if

  if ( namip==5 .or. namip==15 .or. namip==25 ) then ! salinity also read
     
    ! ASCII
    open(unit=77,file=salfile,status='old',form='formatted',iostat=ierr)
    if (ierr/=0) then
      write(6,*) "ERROR: Cannot read AMIP salfile ",trim(salfile)
      call ccmpi_abort(-1)
    end if
    iyear=-999
    imonth=-999
    do while (iyr_m/=iyear.or.imo_m/=imonth)
      read(77,*) imonth,iyear,il_in,jl_in,rlon_in,rlat_in,schmidt_in,header
      write(6,'("reading sal ",i2,i5,2i4,2f6.1,f6.3,a22)') imonth,iyear,il_in,jl_in,rlon_in,rlat_in,schmidt_in,header
      if(il_g/=il_in.or.jl_g/=jl_in.or.abs(rlong0-rlon_in)>1.e-20.or.abs(rlat0-rlat_in)>1.e-20.or. &
          abs(schmidt-schmidt_in)>1.e-20)then
        write(6,*) 'il_g,il_in,jl_g,jl_in,rlong0,rlon_in',il_g,il_in,jl_g,jl_in,rlong0,rlon_in
        write(6,*) 'rlat0,rlat_in,schmidt,schmidt_in',rlat0,rlat_in,schmidt,schmidt_in
        write(6,*) 'wrong sal file'
        call ccmpi_abort(-1)
      endif
      read(77,*) ssta_g(:,1)
      write(6,*) 'want imo_m,iyr_m; asal ',imo_m,iyr_m,ssta_g(idjd_g,1)
    end do
    if ( (namip>=11.and.namip<=15) .or. (namip>=21.and.namip<=25) ) then
      read(77,'(i2,i5,a22)') imonth,iyear,header
      write(6,*) 'reading b_sal data:',imonth,iyear,header
      read(77,*) ssta_g(:,2)
    else
      ssta_g(:,2) = ssta_g(:,1)
    end if
    write(6,*) 'bsal(idjd) ',ssta_g(idjd_g,2)    
    if ( namip==1 .or. namip==2 .or. (namip>=3.and.namip<=5) .or. (namip>=11.and.namip<=15) .or. &
       (namip>=21.and.namip<=25) ) then
      read(77,'(i2,i5,a22)') imonth,iyear,header
      write(6,*) 'should agree with imo,iyr ',imo,iyr
      if(iyr/=iyear.or.imo/=imonth) then
        call ccmpi_abort(-1)
      end if
      write(6,*) 'reading c_sal data:',imonth,iyear,header
      read(77,*) ssta_g(:,3)
    end if
    write(6,*) 'csal(idjd) ',ssta_g(idjd_g,3)       
    read(77,'(i2,i5,a22)') imonth,iyear,header
    write(6,*) 'reading d_sal data:',imonth,iyear,header
    read(77,*) ssta_g(:,4)
    write(6,*) 'dsal(idjd) ',ssta_g(idjd_g,4)
    if ( (namip>=11.and.namip<=15) .or. (namip>=21.and.namip<=25) ) then
      read(77,'(i2,i5,a22)') imonth,iyear,header
      write(6,*) 'reading e_sal data:',imonth,iyear,header
      read(77,*) ssta_g(:,5)
    else
      ssta_g(:,5) = ssta_g(:,4)
    end if
    write(6,*) 'esal(idjd) ',ssta_g(idjd_g,5)    
    close(77)
    call ccmpi_distribute(asal(:,1:5), ssta_g(:,1:5))
      
  end if

  deallocate( ssta_g )
      
end if ! (iernc==0) .. else ..

return
end subroutine amiprd_month

! Read daily SST data
subroutine amipsst_day

use cc_mpi                                   ! CC MPI routines
use dates_m                                  ! Date data
use filnames_m                               ! Filenames
use infile                                   ! Input file routines
use latlong_m                                ! Lat/lon coordinates
use mlo_ctrl                                 ! Ocean physics control layer
use newmpar_m                                ! Grid parameters
use nharrs_m, only : lrestart                ! Non-hydrostatic atmosphere arrays
use parm_m                                   ! Model configuration
use parmgeom_m                               ! Coordinate data
use pbl_m                                    ! Boundary layer arrays
use soil_m                                   ! Soil and surface data
use soilsnow_m                               ! Soil, snow and surface data

implicit none

integer, parameter :: nihead = 54
integer, parameter :: nrhead = 14

integer il_in, jl_in
integer iq, k, mm, varid, varid_time, lsmid, ierr
integer iyear, imonth, iday, iyr_m, imo_m, idy_m
integer kdate_r, ktime_r
integer kdate_rsav, ktime_rsav
integer, save :: iarchx = 0
integer, save :: maxarchi = -1
integer, save :: leap_file = 0
integer, dimension(3) :: spos, npos
#ifdef i8r8
integer, dimension(nihead) :: nahead
#else
integer(kind=4), dimension(nihead) :: nahead
#endif
integer(kind=8) mtimer_r
real rlon_in, rlat_in, schmidt_in
real of, sc
real timer_r, wgt
real, dimension(nrhead) :: ahead
real, dimension(ifull) :: oldsst, newsst, deltasst, timelt
real, dimension(:), allocatable :: axs_a, ays_a, azs_a
real, dimension(:), allocatable :: bxs_a, bys_a, bzs_a
real, dimension(:), allocatable :: wts_a
real, dimension(:), allocatable :: lsmr_g
real, dimension(:,:), allocatable :: ssta_g, ssta_l
real, dimension(:), allocatable, save :: fraciceb, ssta
real(kind=8), dimension(:,:), pointer :: xx4, yy4
real(kind=8), dimension(:), pointer :: z_a, x_a, y_a
logical tst, ltest
logical, save :: firstcall = .true.
logical, save :: interpolate = .false.
logical, dimension(:), allocatable, save :: lsma_g
character(len=10), save :: unitstr
character(len=80) datestring
character(len=10) calendarstring

if ( .not.allocated(ssta) ) then
  allocate(ssta(ifull),fraciceb(ifull))
end if

if ( mod(ktau,nperday)==0 ) then
  if ( myid==0 ) then
    write(6,*) "amipsst called at end of day for ktau,mtimer,namip ",ktau,mtimer,namip

    kdate_r = kdate
    ktime_r = ktime
    mtimer_r = int(mtimer,8)
    call datefix(kdate_r,ktime_r,mtimer_r)
    
    iyr_m = kdate_r/10000
    imo_m = (kdate_r-10000*iyr_m)/100
    idy_m = kdate_r - 10000*iyr_m - 100*imo_m
    write(6,*) "Searching for iyr_m,imo_m,idy_m ",iyr_m,imo_m,idy_m
  
    if ( firstcall ) then
      ! check grid definition
      call ccnf_get_attg(ncidx,'int_header',nahead)
      call ccnf_get_attg(ncidx,'real_header',ahead)
      il_in      = nahead(1)
      jl_in      = nahead(2)
      rlon_in    = ahead(5)
      rlat_in    = ahead(6)
      schmidt_in = ahead(7)
      if ( schmidt_in<=0. .or. schmidt_in>1. ) then
        rlon_in    = ahead(6)
        rlat_in    = ahead(7)
        schmidt_in = ahead(8)
      end if  ! (schmidtx<=0..or.schmidtx>1.)  
      call ccnf_get_attg(ncidx,'leap',leap_file,tst) ! old method
      if ( tst ) leap_file = 0                       ! old method
      call ccnf_get_att(ncidx,'time','calendar',calendarstring,ierr)      ! new method
      if ( ierr==0 ) then                                                 ! new method
        if ( trim(calendarstring)=="standard" ) leap_file = 1             ! new method
        if ( trim(calendarstring)=="360_day" ) leap_file = 2              ! new method
      end if                                                              ! new method
      call ccnf_inq_dimlen(ncidx,'time',maxarchi)
      
      if ( leap==0 .or. leap==1 ) then
        if ( leap_file==2 ) then
          write(6,*) "ERROR: mismatch between amipsst file and CCAM calendar"
          write(6,*) "amipsst (leap_file=2) implies 360_day calendar"
          if ( leap==0 ) then
            write(6,*) "CCAM (leap==0) implies 365 day calendar"
          else if ( leap==1 ) then
            write(6,*) "CCAM (leap=1) implies standard calendar"
          end if
          call ccmpi_abort(-1)
        end if
      else if ( leap==2 ) then
        if ( leap_file==0 .or. leap_file==1 ) then
          write(6,*) "ERROR: mismatch between amipsst file and CCAM calendar"
          if ( leap_file==0 ) then
            write(6,*) "amipsst (leap_file==0) implies 365 day calendar"
          else if ( leap_file==1 ) then
            write(6,*) "amipsst (leap_file==1) implies standard calendar"
          end if
          write(6,*) "CCAM (leap=2) implies 360_day calendar"
          call ccmpi_abort(-1)
        end if    
      end if  
         
      interpolate = ( il_g/=il_in .or. jl_g/=jl_in .or. abs(rlong0-rlon_in)>1.e-6 .or. abs(rlat0-rlat_in)>1.e-6 .or. &
           abs(schmidt-schmidt_in)>0.0002 )
      ik = il_in
  
      if ( interpolate ) then
        write(6,*) "Interpolation is required for AMIPSST"
    
        call define_grid(ik,rlon_in,rlat_in,schmidt_in)
        
        allocate( lsma_g(ik*ik*6) )
        lsma_g = .false.
        npos(1) = ik
        npos(2) = 6*ik
        spos(1:2) = 1
        allocate( lsmr_g(ik*ik*6) )  
        call ccnf_get_vara(ncidx,'lsm',spos(1:2),npos(1:2),lsmr_g(:))
        lsma_g = lsmr_g >= 0.5
        deallocate( lsmr_g )
      else
        write(6,*) "No interpolation is required for AMIPSST"
      end if
    
      unitstr = ''
      call ccnf_get_att(ncidx,'tos','units',unitstr)
    
      firstcall = .false.
    end if
    
    ! temporary arrays
    if ( interpolate ) then
      allocate( ssta_g(ik*ik*6,1), ssta_l(ifull_g,1) )
      ssta_g = 0.
      ssta_l = 0.
    else
      allocate( ssta_g(ifull_g,1), ssta_l(0,0) )
      ssta_g = 0.
    end if
    
    ! search for required day
    iyear  = -999
    imonth = -999
    iday   = -999
    ltest = .true.
    call ccnf_get_att(ncidx,'time','units',datestring)
    call processdatestring(datestring,kdate_rsav,ktime_rsav)
    do while ( ltest .and. iarchx<maxarchi )
      iarchx = iarchx + 1
      kdate_r = kdate_rsav
      ktime_r = ktime_rsav
      call ccnf_get_vara(ncidx,'time',iarchx,timer_r)
      mtimer_r = nint(timer_r,8)
      call datefix(kdate_r,ktime_r,mtimer_r,allleap=leap_file)
      iyear  = kdate_r/10000
      imonth = (kdate_r-iyear*10000)/100
      iday   = kdate_r - 10000*iyear - 100*imonth
      ltest  = iyr_m/=iyear .or. imo_m/=imonth .or. idy_m/=iday
    end do    
    if ( ltest ) then
      write(6,*) "ERROR: Cannot locate year,month,day ",iyr_m,imo_m,idy_m
      write(6,*) "       in file ",trim(sstfile)
      call ccmpi_abort(-1)
    end if
    write(6,*) "Reading SST data from amipsst file at iarchx=",iarchx
    npos(1) = ik
    npos(2) = 6*ik
    npos(3) = 1    
    spos(1:2) = 1
    spos(3) = iarchx
    call ccnf_get_vara(ncidx,'tos',spos,npos,ssta_g(:,1))
    call ccnf_get_att(ncidx,'tos','add_offset',of,ierr=ierr)
    if ( ierr /= 0 ) of = 0.
    call ccnf_get_att(ncidx,'tos','scale_factor',sc,ierr=ierr)
    if ( ierr /= 0 ) sc = 1.
    ssta_g(:,1) = sc*ssta_g(:,1) + of
    if ( trim(unitstr) == 'C' ) then
      write(6,*) "Converting SSTs from degC to degK"    
      ssta_g(:,1) = ssta_g(:,1) + 273.16
    end if  
    if ( interpolate ) then
      call fill_cc1(ssta_g(:,1),lsma_g)
      call doints1(ssta_g(:,1),ssta_l(:,1))
      call ccmpi_distribute(ssta, ssta_l(:,1))
    else
      call ccmpi_distribute(ssta, ssta_g(:,1))
    end if  

    if ( namip==2 .or. namip==3 .or. namip==4 .or. namip==5 .or. namip==13 .or. &
         namip==14 .or. namip==15 .or. namip==24 .or. namip==25 ) then   ! sice also read at middle of month

      ! NETCDF
      write(6,*) "Reading Sea Ice data from amipsst file"
      spos(3) = iarchx
      call ccnf_get_vara(ncidx,'sic',spos,npos,ssta_g(:,1))
      call ccnf_get_att(ncidx,'sic','add_offset',of,ierr=ierr)
      if (ierr/=0) of=0.
      call ccnf_get_att(ncidx,'sic','scale_factor',sc,ierr=ierr)
      if (ierr/=0) sc=1.
      ssta_g(:,1) = sc*ssta_g(:,1) + of
      !ssta_g(:,1) = 100.*ssta_g(:,1)  
      if ( interpolate ) then
        call fill_cc1(ssta_g(:,1),lsma_g)
        call doints1(ssta_g(:,1),ssta_l(:,1))
        call ccmpi_distribute(fraciceb, ssta_l(:,1))
      else
        call ccmpi_distribute(fraciceb, ssta_g(:,1))
      end if  
      
    else   
      where( ssta<271.2 )
        fraciceb = 1.
      elsewhere
        fraciceb = 0.
      end where      
    end if       
    
    deallocate( ssta_g, ssta_l )         
         
  else      
    
    call ccmpi_distribute(ssta)  
    if ( namip==2 .or. namip==3 .or. namip==4 .or. namip==5 .or. namip==13 .or. &
         namip==14 .or. namip==15 .or. namip==24 .or. namip==25 ) then     
      call ccmpi_distribute(fraciceb)
    else
      where( ssta<271.2 )
        fraciceb = 1.
      elsewhere
        fraciceb = 0.
      end where
    end if  
      
  end if ! myid==0 ..else..        
       
else
  if ( nmlo==0 ) return  
end if    

!--------------------------------------------------------------------------------------------------
! Sea-ice and Sea-Surface-Temperature
if ( nmlo==0 ) then
  sicedep(:) = 0. 
  if ( ktau==0 .and. .not.lrestart ) then  ! will set sicedep in indata
    ! This case is for poor initial conditions ifile  
    fracice(:) = fraciceb(:)
    where ( .not.land(1:ifull) )
      tss(:)   = ssta
      tgg(:,1) = ssta
    end where
  else
    do iq = 1,ifull
      if ( .not.land(iq) ) then
        if ( fraciceb(iq)>0. ) then
          if ( fracice(iq)<1.e-20 ) then
            ! create values for tice, and set averaged tss
            ! N.B. if already a sice point, keep present tice
            tggsn(iq,1) = 271.2
          end if  ! (fracice(iq)==0.)
          if ( rlatt(iq)>0. ) then
            sicedep(iq) = 2.
          else
            sicedep(iq) = 1.
          endif ! (rlatt(iq)>0.)
        endif   ! (fraciceb(iq)>0.)
        fracice(iq) = fraciceb(iq)
        tss(iq) = tggsn(iq,1)*fracice(iq) + ssta(iq)*(1.-fracice(iq))
      endif      ! (.not.land(iq))
    enddo
  end if  
elseif ( ktau>0 ) then
  if ( nud_hrs<=0 ) then
    write(6,*) "ERROR: namip/=0 has been selected with in-line ocean model (nmlo/=0)"  
    write(6,*) "nud_hrs>0 must be specified for relaxiation of SSTs"
    call ccmpi_abort(-1)
  end if
  oldsst = ssta
  call mloexpmelt(timelt)
  timelt = min(timelt+0.01,ssta)
  newsst = ssta*(1.-fraciceb(:)) + timelt(:)*fraciceb(:) - wrtemp
  wgt = dt/real(nud_hrs*3600)
  call mloexport("temp",oldsst,1,0)
  deltasst = newsst - oldsst
  do k = 1,kbotmlo
    oldsst = 273.16 - wrtemp  
    call mloexport("temp",oldsst,k,0)
    newsst = oldsst + max( wgt*deltasst, 271.16-wrtemp-oldsst )
    call mloimport("temp",newsst,k,0)
  end do  
  do k = 1,ms
    call mloexport("temp",tgg(:,k),k,0)
    where ( tgg(:,k)<100. )
      tgg(:,k) = tgg(:,k) + wrtemp
    end where    
  end do
end if ! if (nmlo==0) ..else..

return
end subroutine amipsst_day

subroutine doints1(s,sout)
      
use newmpar_m              ! Grid parameters
use parm_m                 ! Model configuration

implicit none
      
integer mm
real, dimension(:), intent(in) :: s
real, dimension(:), intent(inout) :: sout
real, dimension(ifull_g) :: wrk
real, dimension(-1:ik+2,-1:ik+2,0:npanels) :: sx

sx(-1:ik+2,-1:ik+2,0:npanels) = 0.

sout(1:ifull_g) = 0.
sx(1:ik,1:ik,0:npanels) = reshape( s(1:6*ik*ik), (/ ik, ik, npanels+1 /) )
call sxpanelbounds(sx)
do mm = 1,m_fly
  call intsb(sx,wrk,nface4(:,mm),xg4(:,mm),yg4(:,mm))
  sout(1:ifull_g) = sout(1:ifull_g) + wrk/real(m_fly)
end do

return
end subroutine doints1

subroutine sxpanelbounds(sx_l)

use newmpar_m

implicit none

integer i, n, n_w, n_e, n_n, n_s
real, dimension(-1:ik+2,-1:ik+2,0:npanels), intent(inout) :: sx_l

do n = 0,npanels
  if ( mod(n,2)==0 ) then
    n_w = mod(n+5, 6)
    n_e = mod(n+2, 6)
    n_n = mod(n+1, 6)
    n_s = mod(n+4, 6)
    do i = 1,ik
      sx_l(-1,i,n)   = sx_l(ik-1,i,n_w)
      sx_l(0,i,n)    = sx_l(ik,i,n_w)
      sx_l(ik+1,i,n) = sx_l(ik+1-i,1,n_e)
      sx_l(ik+2,i,n) = sx_l(ik+1-i,2,n_e)
      sx_l(i,-1,n)   = sx_l(ik-1,ik+1-i,n_s)
      sx_l(i,0,n)    = sx_l(ik,ik+1-i,n_s)
      sx_l(i,ik+1,n) = sx_l(i,1,n_n)
      sx_l(i,ik+2,n) = sx_l(i,2,n_n)
    end do ! i
    sx_l(0,0,n)       = sx_l(ik,1,n_w)        ! ws
    sx_l(-1,0,n)      = sx_l(ik,2,n_w)        ! wws
    sx_l(0,-1,n)      = sx_l(ik,ik-1,n_s)     ! wss
    sx_l(ik+1,0,n)    = sx_l(ik,1,n_e)        ! es  
    sx_l(ik+2,0,n)    = sx_l(ik-1,1,n_e)      ! ees 
    sx_l(ik+1,-1,n)   = sx_l(ik,2,n_e)        ! ess        
    sx_l(0,ik+1,n)    = sx_l(ik,ik,n_w)       ! wn  
    sx_l(-1,ik+1,n)   = sx_l(ik,ik-1,n_w)     ! wwn
    sx_l(0,ik+2,n)    = sx_l(ik-1,ik,n_w)     ! wnn
    sx_l(ik+1,ik+1,n) = sx_l(1,1,n_e)         ! en  
    sx_l(ik+2,ik+1,n) = sx_l(2,1,n_e)         ! een  
    sx_l(ik+1,ik+2,n) = sx_l(1,2,n_e)         ! enn  
  else
    n_w = mod(n+4, 6)
    n_e = mod(n+1, 6)
    n_n = mod(n+2, 6)
    n_s = mod(n+5, 6)
    do i = 1,ik
      sx_l(-1,i,n)   = sx_l(ik+1-i,ik-1,n_w)  
      sx_l(0,i,n)    = sx_l(ik+1-i,ik,n_w)
      sx_l(ik+1,i,n) = sx_l(1,i,n_e)
      sx_l(ik+2,i,n) = sx_l(2,i,n_e)
      sx_l(i,-1,n)   = sx_l(i,ik-1,n_s)
      sx_l(i,0,n)    = sx_l(i,ik,n_s)
      sx_l(i,ik+1,n) = sx_l(1,ik+1-i,n_n)
      sx_l(i,ik+2,n) = sx_l(2,ik+1-i,n_n)
    end do ! i
    sx_l(0,0,n)       = sx_l(ik,ik,n_w)      ! ws
    sx_l(-1,0,n)      = sx_l(ik-1,ik,n_w)    ! wws
    sx_l(0,-1,n)      = sx_l(2,ik,n_s)       ! wss
    sx_l(ik+1,0,n)    = sx_l(1,1,n_e)        ! es
    sx_l(ik+2,0,n)    = sx_l(1,2,n_e)        ! ees
    sx_l(ik+1,-1,n)   = sx_l(2,1,n_e)        ! ess
    sx_l(0,ik+1,n)    = sx_l(1,ik,n_w)       ! wn       
    sx_l(-1,ik+1,n)   = sx_l(2,ik,n_w)       ! wwn   
    sx_l(0,ik+2,n)    = sx_l(1,ik-1,n_w)     ! wnn
    sx_l(ik+1,ik+1,n) = sx_l(1,ik,n_e)       ! en  
    sx_l(ik+2,ik+1,n) = sx_l(1,ik-1,n_e)     ! een  
    sx_l(ik+1,ik+2,n) = sx_l(2,ik,n_e)       ! enn  
  end if   ! mod(n,2)==0 ..else..
end do       ! n loop

return
end subroutine sxpanelbounds

subroutine intsb(sx_l,sout,nface_l,xg_l,yg_l)
      
!     same as subr ints, but with sout passed back and no B-S      
!     s is input; sout is output array
!     later may wish to save idel etc between array calls
!     this one does linear interp in x on outer y sides
!     doing x-interpolation before y-interpolation
!     This is a global routine 

use newmpar_m              ! Grid parameters
use parm_m                 ! Model configuration

implicit none

integer, dimension(ifull_g), intent(in) :: nface_l
integer :: idel, jdel, n, iq
real, dimension(ifull_g), intent(inout) :: sout
real, intent(in), dimension(ifull_g) :: xg_l, yg_l
real, dimension(-1:ik+2,-1:ik+2,0:npanels), intent(in) :: sx_l
real xxg, yyg, cmin, cmax
real dmul_2, dmul_3, cmul_1, cmul_2, cmul_3, cmul_4
real emul_1, emul_2, emul_3, emul_4, rmul_1, rmul_2, rmul_3, rmul_4

do iq = 1,ifull_g   ! runs through list of target points
  n = nface_l(iq)
  idel = int(xg_l(iq))
  xxg = xg_l(iq) - real(idel)
  jdel = int(yg_l(iq))
  yyg = yg_l(iq) - real(jdel)
  ! bi-cubic
  cmul_1 = (1.-xxg)*(2.-xxg)*(-xxg)/6.
  cmul_2 = (1.-xxg)*(2.-xxg)*(1.+xxg)/2.
  cmul_3 = xxg*(1.+xxg)*(2.-xxg)/2.
  cmul_4 = (1.-xxg)*(-xxg)*(1.+xxg)/6.
  dmul_2 = (1.-xxg)
  dmul_3 = xxg
  emul_1 = (1.-yyg)*(2.-yyg)*(-yyg)/6.
  emul_2 = (1.-yyg)*(2.-yyg)*(1.+yyg)/2.
  emul_3 = yyg*(1.+yyg)*(2.-yyg)/2.
  emul_4 = (1.-yyg)*(-yyg)*(1.+yyg)/6.
  cmin = min(sx_l(idel,  jdel,n),sx_l(idel+1,jdel,  n), &
             sx_l(idel,jdel+1,n),sx_l(idel+1,jdel+1,n))
  cmax = max(sx_l(idel,  jdel,n),sx_l(idel+1,jdel,  n), &
             sx_l(idel,jdel+1,n),sx_l(idel+1,jdel+1,n))
  rmul_1 = sx_l(idel,  jdel-1,n)*dmul_2 + sx_l(idel+1,jdel-1,n)*dmul_3
  rmul_2 = sx_l(idel-1,jdel,  n)*cmul_1 + sx_l(idel,  jdel,  n)*cmul_2 + &
           sx_l(idel+1,jdel,  n)*cmul_3 + sx_l(idel+2,jdel,  n)*cmul_4
  rmul_3 = sx_l(idel-1,jdel+1,n)*cmul_1 + sx_l(idel,  jdel+1,n)*cmul_2 + &
           sx_l(idel+1,jdel+1,n)*cmul_3 + sx_l(idel+2,jdel+1,n)*cmul_4
  rmul_4 = sx_l(idel,  jdel+2,n)*dmul_2 + sx_l(idel+1,jdel+2,n)*dmul_3
  sout(iq) = min( max( cmin, rmul_1*emul_1 + rmul_2*emul_2 + rmul_3*emul_3 + rmul_4*emul_4 ), cmax ) ! Bermejo & Staniforth
end do    ! iq loop

return
end subroutine intsb

! *****************************************************************************
! FILL ROUTINES

subroutine fill_cc1(a_io,land)
      
! routine fills in interior of an array which has undefined points
! this version is distributed over processes with input files

implicit none

real, dimension(:), intent(inout) :: a_io
integer nrem, j, n, cc, i
integer ccount
real, parameter :: value=999.  ! missing value flag
real, dimension(-1:ik+2,-1:ik+2,6) :: c_io
real, dimension(4) :: s_test
logical, dimension(:), intent(in) :: land
logical, dimension(4) :: l_test


where ( land(1:6*ik*ik) )
  a_io(1:6*ik*ik) = value
end where

nrem = count( abs(a_io(1:6*ik*ik)-value)<1.E-6 )
c_io = value

do while ( nrem>0 )
  c_io(1:ik,1:ik,1:6) = reshape( a_io(1:6*ik*ik), (/ ik, ik, 6 /) )  
  call sxpanelbounds(c_io(:,:,:))
  do n = 1,6
    do j = 1,ik
      cc = (j-1)*ik + (n-1)*ik**2
      do i = 1,ik
        if ( abs(a_io(cc+i)-value)<1.e-6 ) then
          s_test(1) = c_io(i,j+1,n)
          s_test(2) = c_io(i,j-1,n)
          s_test(3) = c_io(i+1,j,n)
          s_test(4) = c_io(i-1,j,n)
          l_test(:) = abs(s_test(:)-value)>=1.e-6
          ccount = count(l_test)
          if ( ccount>0 ) then
            a_io(cc+i) = sum(s_test,l_test)/real(ccount)
          end if
        end if  
      end do   
    end do
  end do
  nrem = count( abs(a_io(1:6*ik*ik)-value)<1.E-6 )
  ! test for convergence
  if ( nrem==6*ik*ik ) then
    write(6,*) "Cannot perform fill as all points are trivial"    
    a_io = 0.
    nrem = 0
  end if
end do

return
end subroutine fill_cc1

subroutine define_grid(ik_in,rlon_in,rlat_in,schmidt_in)

use latltoij_m            ! Lat/Lon to cubic ij conversion
use newmpar_m             ! Grid parameters
use parm_m                ! Model configuration
use setxyz_m              ! Define CCAM grid
use workglob_m            ! Additional grid interpolation

implicit none

integer, intent(in) :: ik_in
integer :: iq, mm
real, intent(in) :: rlon_in, rlat_in, schmidt_in
real, dimension(:), allocatable :: axs_a, ays_a, azs_a
real, dimension(:), allocatable :: bxs_a, bys_a, bzs_a
real, dimension(:), allocatable :: wts_a
real(kind=8), dimension(:,:), pointer :: xx4, yy4
real(kind=8), dimension(:), pointer :: z_a, x_a, y_a

ik = ik_in

allocate( nface4(ifull_g,4), xg4(ifull_g,4), yg4(ifull_g,4) )
allocate( xx4(1+4*ik,1+4*ik), yy4(1+4*ik,1+4*ik) )
allocate( axs_a(ik*ik*6), ays_a(ik*ik*6), azs_a(ik*ik*6) )
allocate( bxs_a(ik*ik*6), bys_a(ik*ik*6), bzs_a(ik*ik*6) )
allocate( x_a(ik*ik*6), y_a(ik*ik*6), z_a(ik*ik*6) )
allocate( wts_a(ik*ik*6) )
!   following setxyz call is for source data geom    ****   
do iq = 1,ik*ik*6
  axs_a(iq) = real(iq)
  ays_a(iq) = real(iq)
  azs_a(iq) = real(iq)
end do 
call setxyz(ik,rlon_in,rlat_in,-schmidt_in,x_a,y_a,z_a,wts_a,axs_a,ays_a,azs_a,bxs_a,bys_a,bzs_a,xx4,yy4, &
            id,jd,ktau,ds)
deallocate( x_a, y_a, z_a ) 
deallocate( axs_a, ays_a, azs_a, bxs_a, bys_a, bzs_a )
deallocate( wts_a ) 
nullify( x_a, y_a, z_a )    
! setup interpolation arrays
do mm = 1,m_fly  !  was 4, now may be set to 1 in namelist
  call latltoij(rlong4(:,mm),rlat4(:,mm),          & !input
                rlon_in,rlat_in,schmidt_in,        & !input
                xg4(:,mm),yg4(:,mm),nface4(:,mm),  & !output (source)
                xx4,yy4,ik)
end do
deallocate( xx4, yy4 )    
nullify( xx4, yy4 )  

return
end subroutine define_grid

subroutine erase_grid

implicit none

deallocate( nface4, xg4, yg4 )

end subroutine erase_grid

subroutine time_of_month(iyr,imo,iday,mdays,kdate,mtimer,leap)

implicit none

integer, intent(in) :: kdate, leap, mtimer
integer, intent(out) :: iyr, imo, iday, mdays
integer, dimension(0:13) :: mdays_a

iyr = kdate/10000
imo = (kdate-10000*iyr)/100
iday = kdate - 10000*iyr - 100*imo + mtimer/(60*24)
if ( leap==0 ) then ! 365 day calendar
  mdays_a = (/ 31, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31, 31 /)
else if ( leap==1 ) then ! 365/366 day calendar
  mdays_a = (/ 31, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31, 31 /)
  if ( mod(iyr,4)==0 ) mdays_a(2) = 29
  if ( mod(iyr,100)==0 ) mdays_a(2) = 28
  if ( mod(iyr,400)==0 ) mdays_a(2) = 29
else if ( leap==2 ) then ! 360 day calendar
  mdays_a = (/ 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30 /)  
else
  write(6,*) "ERROR: Invalid leap = ",leap
  stop
end if
do while ( iday>mdays_a(imo) )
  iday = iday - mdays_a(imo)
  imo = imo + 1
  if ( imo>12 ) then
    imo = 1
    iyr = iyr + 1
    if ( leap==1 ) then
      if ( mod(iyr,4)==0 ) mdays_a(2) = 29
      if ( mod(iyr,100)==0 ) mdays_a(2) = 28
      if ( mod(iyr,400)==0 ) mdays_a(2) = 29
    end if
  end if
end do
mdays = mdays_a(imo)

return
end subroutine time_of_month

subroutine time_interpolate(datout,dat_pm,dat_p,dat_c,dat_n,dat_np,x,method,land)

integer, intent(in) :: method
integer iq
real, intent(in) :: x
real a0, a1, a2, aa, bb, cc, mp1, mp2, c2, c3, c4
real ssta2, ssta3, ssta4
real, dimension(:), intent(out) :: datout
real, dimension(:), intent(in) :: dat_pm, dat_p, dat_c, dat_n, dat_np
logical, dimension(:), intent(in), optional :: land
logical, dimension(size(datout)) :: land_l

! dat_pm - previous month -1 (current month - 2)
! dat_p  - previous month (current month - 1)
! dat_c  - current month
! dat_n  - next month ( current month + 1)
! dat_np - next month +1 ( current month + 2 )

land_l = .false.
if ( present(land) ) then
  land_l = land
end if

! method=3,4,5  Use PWCB interpolation
if ( method==3 .or. method==4 .or. method==5 ) then
  !--------------------------------------------------------------------------------------------------
  ! Piece-wise cubic bessel interpolation
  do iq=1,size(datout)
    if( .not.land_l(iq) )then
      c2=dat_p(iq)
      c3=dat_p(iq)+dat_c(iq)
      c4=c3+dat_n(iq)          
      datout(iq)=.5*c3+(4.*c3-5.*c2-c4)*x+1.5*(c4+3.*c2-3.*c3)*x*x
    endif      ! (.not.land(iq))
  enddo
endif

! method=11,13,14,15  Use JMc interpolation
if ( method==11 .or. method==13 .or. method==14 .or. method==15 ) then
  !--------------------------------------------------------------------------------------------------
  ! John McGregor 5-pt, piece-wise, cubic interpolation
  do iq=1,size(datout)
    if ( .not.land_l(iq) ) then
      ssta2=(24.*dat_p(iq)-dat_pm(iq)-dat_c(iq))/22. 
      ssta3=(24.*dat_c(iq)-dat_p(iq)-dat_n(iq))/22. 
      ssta4=(24.*dat_n(iq)-dat_c(iq)-dat_np(iq))/22. 
      c2=(-dat_pm(iq)+9*ssta2+9*ssta3-dat_n(iq))/16.
      c3=(-dat_p(iq)+9*ssta3+9*ssta4-dat_np(iq))/16.
      datout(iq)=c2+(6.*dat_c(iq)-4.*c2-2.*c3)*x    &
                 +(3.*c2+3.*c3-6.*dat_c(iq))*x*x 
    endif      ! (.not.land(iq))
  enddo
end if

! method=21,24,25  Use approx linear AMIP interpolation
if ( method==21 .or. method==24 .or. method==25 ) then
  !--------------------------------------------------------------------------------------------------
  ! Approximation of piece-wise, linear AMIP interpolation
  if ( x<0.5 ) then
    do iq = 1,size(datout)
      if ( .not.land_l(iq) ) then
        a0 = 0.5*dat_p(iq)
        a1 = -dat_c(iq)
        a2 = 0.5*dat_n(iq)
        aa = a0 + a1 + a2
        bb = -3.*a0 - 2.*a1 - a2
        cc = 2.*a0
        mp1 = 0.25*aa + 0.5*bb + cc ! start of month value
        a0 = 0.5*dat_c(iq)
        a1 = -dat_n(iq)
        a2 = 0.5*dat_np(iq)
        aa = a0 + a1 + a2
        bb = -3.*a0 - 2.*a1 - a2
        cc = 2.*a0
        mp2 = 0.25*aa + 0.5*bb + cc ! end of month value
        c4 = 2.*dat_c(iq) - 0.5*mp1 - 0.5*mp2 ! mid-point value
        c2 = mp1                              ! intercept
        c3 = 2.*(c4-c2)                       ! gradient
        datout(iq) = c3*x + c2
      end if
    end do
  else
    do iq = 1,size(datout)
      if ( .not.land_l(iq) ) then
        a0 = 0.5*dat_p(iq)
        a1 = -dat_c(iq)
        a2 = 0.5*dat_n(iq)
        aa = a0 + a1 + a2
        bb = -3.*a0 - 2.*a1 - a2
        cc = 2.*a0
        mp1 = 0.25*aa + 0.5*bb + cc ! start of month value
        a0 = 0.5*dat_c(iq)
        a1 = -dat_n(iq)
        a2 = 0.5*dat_np(iq)
        aa = a0 + a1 + a2
        bb = -3.*a0 - 2.*a1 - a2
        cc = 2.*a0
        mp2 = 0.25*aa + 0.5*bb + cc ! end of month value
        c4 = 2.*dat_c(iq) - 0.5*mp1 - 0.5*mp2 ! mid-point value
        c3 = 2.*(mp2 - c4)                    ! gradient
        c2 = 2.*c4 - mp2                      ! intercept
        datout(iq) = c3*x + c2
      end if
    end do
  end if
end if

return
end subroutine time_interpolate

end module amipsst_m