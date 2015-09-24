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

! This module reads 3D input fields (e.g., ozone) and interpolates
! them onto the model grid

module ozoneread

implicit none

private      
public o3_read, o3set, fieldinterpolate, o3regrid
public o3read_amip, o3set_amip

! CMIP3 / CMIP5 ozone fields
integer, save :: ii, jj, kk
real, dimension(:,:), allocatable, save :: o3pre, o3mth, o3nxt
real, dimension(:), allocatable, save :: o3pres
real, dimension(:,:), allocatable, save :: dduo3n, ddo3n2
real, dimension(:,:), allocatable, save :: ddo3n3, ddo3n4

! o3amip parameters
integer, parameter :: mo = 12    ! Months
integer, parameter :: jg = 64    ! Latitudes
integer, parameter :: kg = 59    ! Levels
integer, parameter :: lg = kg+1  ! Layer Interfaces
real, dimension(:), allocatable, save :: glat, dp, gpri
real, dimension(:,:,:), allocatable, save :: gdat

contains

!--------------------------------------------------------------
! This subroutine reads ozone fields

! This version supports CMIP3 (ASCII) and CMIP5 (netcdf) file formats            
subroutine o3_read(sigma,jyear,jmonth)

use cc_mpi
use infile

implicit none

include 'newmpar.h'
include 'filnames.h'
      
integer, intent(in) :: jyear, jmonth
integer nlev, i, l, k
integer ncstatus, ncid, tt
integer valident, yy, mm, nn
integer, dimension(1) :: iti
integer, dimension(4) :: spos, npos
integer, dimension(3) :: dum
real, dimension(:,:,:,:), allocatable, save :: o3dum
real, dimension(:), allocatable, save :: o3lon, o3lat
real, dimension(:), allocatable, save :: o3pack
real, dimension(kl) :: sigma, sigin
character(len=32) cdate
logical tst

real, parameter :: sigtol = 1.e-3
      
!--------------------------------------------------------------
! Read montly ozone data for interpolation
if ( myid == 0 ) then
  write(6,*) "Reading ",trim(o3file)
  call ccnf_open(o3file,ncid,ncstatus)
  if ( ncstatus == 0 ) then
    write(6,*) "Ozone in NetCDF format (CMIP5)"
    call ccnf_inq_dimlen(ncid,'lon',ii)
    allocate( o3lon(ii) )
    call ccnf_inq_varid(ncid,'lon',valident,tst)
    if ( tst ) then
      write(6,*) "lon variable not found"
      call ccmpi_abort(-1)
    end if
    spos(1) = 1
    npos(1) = ii
    call ccnf_get_vara(ncid,valident,spos(1:1),npos(1:1),o3lon)
    call ccnf_inq_dimlen(ncid,'lat',jj)
    allocate( o3lat(jj) )
    call ccnf_inq_varid(ncid,'lat',valident,tst)
    if ( tst ) then
      write(6,*) "lat variable not found"
      call ccmpi_abort(-1)
    end if
    npos(1) = jj
    call ccnf_get_vara(ncid,valident,spos(1:1),npos(1:1),o3lat)
    call ccnf_inq_dimlen(ncid,'plev',kk)
    allocate( o3pres(kk) )
    call ccnf_inq_varid(ncid,'plev',valident,tst)
    if ( tst ) then
      write(6,*) "plev variable not found"
      call ccmpi_abort(-1)
    end if
    npos(1) = kk
    call ccnf_get_vara(ncid,valident,spos(1:1),npos(1:1),o3pres)
    call ccnf_inq_dimlen(ncid,'time',tt)
    call ccnf_inq_varid(ncid,'time',valident,tst)
    if ( tst ) then
      write(6,*) "time variable not found"
      call ccmpi_abort(-1)
    end if
    call ccnf_get_att(ncid,valident,'units',cdate)
    npos(1) = 1
    call ccnf_get_vara(ncid,valident,spos(1:1),npos(1:1),iti)
    write(6,*) "Found ozone dimensions ",ii,jj,kk,tt
    allocate( o3dum(ii,jj,kk,3) )
    read(cdate(14:17),*) yy
    read(cdate(19:20),*) mm
    yy = yy + iti(1)/12
    mm = mm + mod(iti(1),12)
    write(6,*) "Requested date ",jyear,jmonth
    write(6,*) "Initial ozone date ",yy,mm
    nn = (jyear-yy)*12 + (jmonth-mm) + 1
    if ( nn < 1 .or. nn > tt ) then
      write(6,*) "ERROR: Cannot find date in ozone data"
      call ccmpi_abort(-1)
    end if
    write(6,*) "Found ozone data at index ",nn
    spos = 1
    npos(1) = ii
    npos(2) = jj
    npos(3) = kk
    npos(4) = 1
    write(6,*) "Reading O3"
    call ccnf_inq_varid(ncid,'O3',valident,tst)
    if ( tst ) then
      write(6,*) "O3 variable not found"
      call ccmpi_abort(-1)
    end if
    spos(4) = max(nn-1,1)
    call ccnf_get_vara(ncid,valident,spos,npos,o3dum(:,:,:,1))
    spos(4) = nn
    call ccnf_get_vara(ncid,valident,spos,npos,o3dum(:,:,:,2))
    spos(4) = min(nn+1,tt)
    call ccnf_get_vara(ncid,valident,spos,npos,o3dum(:,:,:,3))
    call ccnf_close(ncid)
    ! Here we fix missing values by filling down
    ! If we try to neglect these values in the
    ! vertical column integration, then block
    ! artifacts are apparent
     write(6,*) "Fix missing values in ozone"
     do l = 1,3
       do k = kk-1,1,-1
         where ( o3dum(:,:,k,l) > 1.E34 )
           o3dum(:,:,k,l) = o3dum(:,:,k+1,l)
         end where
       end do
     end do
  else
    write(6,*) "Ozone in ASCII format (CMIP3)"
    ii = 0
    jj = 0
    kk = 0
    open(16,file=o3file,form='formatted',status='old')
    read(16,*) nlev
    if ( nlev /= kl ) then
      write(6,*) ' ERROR - Number of levels wrong in o3_data file'
      call ccmpi_abort(-1)
    end if
    ! Check that the sigma levels are the same
    ! Note that the radiation data has the levels in the reverse order
    read(16,*) (sigin(i),i=kl,1,-1)
    do k = 1,kl
      if ( abs(sigma(k)-sigin(k)) > sigtol ) then
        write(6,*) ' ERROR - sigma level wrong in o3_data file'
        write(6,*) k, sigma(k), sigin(k)
        call ccmpi_abort(-1)
      end if
    end do
          
    allocate( dduo3n(37,kl), ddo3n2(37,kl) )
    allocate( ddo3n3(37,kl), ddo3n4(37,kl) )
    ! Note that the data is written as MAM, SON, DJF, JJA. The arrays in
    ! o3dat are in the order DJF, MAM, JJA, SON
    read(16,"(9f8.5)") ddo3n2
    read(16,"(9f8.5)") ddo3n4
    read(16,"(9f8.5)") dduo3n
    read(16,"(9f8.5)") ddo3n3
    close(16)
  end if
  write(6,*) "Finished reading ozone data"
  dum(1) = ii
  dum(2) = jj
  dum(3) = kk
end if
call ccmpi_bcast(dum(1:3),0,comm_world)
ii = dum(1)
jj = dum(2)
kk = dum(3)
if ( ii > 0 ) then
  allocate(o3pre(ifull,kk),o3mth(ifull,kk),o3nxt(ifull,kk))
  if ( myid == 0 ) then
    allocate( o3pack(ii+jj+kk) )
    o3pack(1:ii) = o3lon(1:ii)
    o3pack(ii+1:ii+jj) = o3lat(1:jj)
    o3pack(ii+jj+1:ii+jj+kk) = o3pres(1:kk)
  else
    allocate( o3lon(ii), o3lat(jj), o3pres(kk) )
    allocate( o3dum(ii,jj,kk,3) )
    allocate( o3pack(ii+jj+kk) )
  end if
  call ccmpi_bcast(o3dum,0,comm_world)
  call ccmpi_bcast(o3pack,0,comm_world)
  o3lon(1:ii) = o3pack(1:ii)
  o3lat(1:jj) = o3pack(ii+1:ii+jj)
  o3pres(1:kk) = o3pack(ii+jj+1:ii+jj+kk)
  deallocate ( o3pack )
  if ( myid == 0 ) write(6,*) "Interpolate ozone data to CC grid"
  call o3regrid(o3pre,o3mth,o3nxt,o3dum,o3lon,o3lat)
  deallocate( o3dum, o3lat, o3lon )
else
  if ( myid /= 0 ) then
    allocate(dduo3n(37,kl),ddo3n2(37,kl))
    allocate(ddo3n3(37,kl),ddo3n4(37,kl))
  end if
  call ccmpi_bcast(ddo3n2,0,comm_world)
  call ccmpi_bcast(ddo3n4,0,comm_world)
  call ccmpi_bcast(dduo3n,0,comm_world)
  call ccmpi_bcast(ddo3n3,0,comm_world)
  call resetd(dduo3n,ddo3n2,ddo3n3,ddo3n4,37*kl)
end if
      
if ( myid == 0 ) write(6,*) "Finished processing ozone data"
      
return
end subroutine o3_read

!--------------------------------------------------------------------
! Interpolate ozone to model grid
subroutine o3set(npts,istart,mins,duo3n,sig,ps)

use latlong_m

implicit none

include 'newmpar.h'
include 'const_phys.h'
      
integer, intent(in) :: npts,mins,istart
integer j,ilat,m,iend,ilatp
real date,rang,rsin1,rcos1,rcos2,theta,angle,than
real do3,do3p,t5
real, dimension(npts), intent(in) :: ps
real, dimension(kl), intent(in) :: sig
real, dimension(npts,kl), intent(out) :: duo3n
real, dimension(npts,kk) :: duma,dumb,dumc

real, parameter :: rlag=14.8125
real, parameter :: year=365.
real, parameter :: amd=28.9644 ! molecular mass of air g/mol
real, parameter :: amo=47.9982 ! molecular mass of o3 g/mol
real, parameter :: dobson=6.022e3/2.69/48.e-3
      
if (allocated(o3mth)) then ! CMIP5 ozone
      
  iend=istart+npts-1
  duma=o3pre(istart:iend,:)
  dumb=o3mth(istart:iend,:)
  dumc=o3nxt(istart:iend,:)
  ! note this inverts levels
  call fieldinterpolate(duo3n,duma,dumb,dumc,o3pres,npts,kl,kk,mins,sig,ps,interpmeth=1)

  ! convert units from mol/mol to g/g
  where (duo3n<1.)
    duo3n=duo3n*amo/amd
  end where
         
else ! CMIP3 ozone

  ! Convert time to day number
  ! date = amod( float(mins)/1440., year)
  ! Use year of exactly 365 days
  date = float(mod(mins,525600))/1440.
  rang = tpi*(date-rlag)/year
  rsin1 = sin(rang)
  rcos1 = cos(rang)
  rcos2 = cos(2.*rang)

  do j=1,npts
    theta = 90. - rlatt(istart+j-1)*180./pi
    t5 = theta/5.
    ilat = int(t5) + 1
    angle = real(ilat)
    than = (t5 - angle)
    ilatp = min(ilat + 1,37)
    ! levels are already in reverse order
    do m = 1,kl
      do3  = dduo3n(ilat,m) + rsin1*ddo3n2(ilat,m) + rcos1*ddo3n3(ilat,m) + rcos2*ddo3n4(ilat,m)
      do3p = dduo3n(ilatp,m) + rsin1*ddo3n2(ilatp,m) + rcos1*ddo3n3(ilatp,m) + rcos2*ddo3n4(ilatp,m)
      duo3n(j,m) = do3 + than*(do3p - do3)
      ! convert from cm stp to gm/gm
      duo3n(j,m) = duo3n(j,m)*1.01325e2*0.1/ps(j)
    end do
  end do      
end if
duo3n=max(1.e-10,duo3n)
      
return
end subroutine o3set

!--------------------------------------------------------------------
! This subroutine interoplates monthly intput fields to the current
! time step.  It also interpolates from pressure levels to CCAM
! vertical levels
subroutine fieldinterpolate(out,fpre,fmth,fnxt,fpres,ipts,ilev,nlev,mins,sig,ps,interpmeth)
      
implicit none

include 'dates.h'

integer leap
common/leap_yr/leap  ! 1 to allow leap years
      
integer, intent(in) :: ipts,ilev,nlev,mins
integer, intent(in), optional :: interpmeth
integer ozoneintp ! ozone interpolation (0=simple, 1=integrate column)
integer date,iq,m,k1,jyear,jmonth
real, dimension(ipts,ilev), intent(out) :: out
real, dimension(ipts,nlev), intent(in) :: fpre,fmth,fnxt
real, dimension(nlev), intent(in) :: fpres
real, dimension(ipts), intent(in) :: ps
real, dimension(ilev), intent(in) :: sig
real, dimension(ilev) :: prf,o3new
real, dimension(nlev) :: o3inp,o3sum,b,c,d
real, dimension(nlev,3) :: o3tmp
real, dimension(12) :: monlen
real, dimension(12), parameter :: oldlen= (/ 0.,31.,59.,90.,120.,151.,181.,212.,243.,273.,304.,334. /)
real rang,fp,mino3

if (present(interpmeth)) then
  ozoneintp=interpmeth
else
  ozoneintp=1
end if

monlen=oldlen
jyear =kdate/10000
jmonth=(kdate-jyear*10000)/100
if (leap==1) then
  if (mod(jyear,4)  ==0) monlen(3:12)=oldlen(3:12)+1
  if (mod(jyear,100)==0) monlen(3:12)=oldlen(3:12)
  if (mod(jyear,400)==0) monlen(3:12)=oldlen(3:12)+1
end if

!     Time interpolation factors
date = nint(mins/1440.)
if (jmonth==12) then
  rang=(date-monlen(12))/31.
else
  rang=(date-monlen(jmonth))/(monlen(jmonth+1)-monlen(jmonth))
end if
if (rang>1.4) then ! use 1.4 to give 40% tolerance
  write(6,*) "WARN: fieldinterpolation is outside input range"
end if
      
do iq=1,ipts

  if (rang<=1.) then
    ! temporal interpolation (PWCB)
    o3tmp(:,1)=fpre(iq,:)
    o3tmp(:,2)=fmth(iq,:)+o3tmp(:,1)
    o3tmp(:,3)=fnxt(iq,:)+o3tmp(:,2)
    b=0.5*o3tmp(:,2)
    c=4.*o3tmp(:,2)-5.*o3tmp(:,1)-o3tmp(:,3)
    d=1.5*o3tmp(:,3)+4.5*o3tmp(:,1)-4.5*o3tmp(:,2)
    o3inp=b+c*rang+d*rang*rang
  else
    ! linear interpolation when rang is out-of-range
    o3inp=max(3.-2.*rang,0.)*0.5*(fmth(iq,:)+fnxt(iq,:))+min(2.*rang-2.,1.)*fnxt(iq,:)
  end if
         
  !-----------------------------------------------------------
  ! Simple interpolation on pressure levels
  ! vertical interpolation (from LDR - Mk3.6)
  ! Note inverted levels
  if (ozoneintp==0) then
    do m=nlev-1,1,-1
      if (o3inp(m)>1.E34) o3inp(m)=o3inp(m+1)
    end do
    prf=0.01*ps(iq)*sig
    do m=1,ilev
      if (prf(m)>fpres(1)) then
        out(iq,ilev-m+1)=o3inp(1)
      elseif (prf(m)<fpres(nlev)) then
        out(iq,ilev-m+1)=o3inp(nlev)
      else
        do k1=2,nlev-1
          if (prf(m)>fpres(k1)) exit
        end do
        fp=(prf(m)-fpres(k1))/(fpres(k1-1)-fpres(k1))
        out(iq,ilev-m+1)=(1.-fp)*o3inp(k1)+fp*o3inp(k1-1)
      end if
    end do

  else
    !-----------------------------------------------------------
    ! Approximate integral of ozone column

    ! pressure levels on CCAM grid
    prf=0.01*ps(iq)*sig
         
    mino3=minval(o3inp)
         
    ! calculate total column of ozone
    o3sum=0.
    o3sum(nlev)=o3inp(nlev)*0.5*sum(fpres(nlev-1:nlev))
    do m=nlev-1,2,-1
      if (o3inp(m)>1.E34) then
        o3sum(m)=o3sum(m+1)
      else
        o3sum(m)=o3sum(m+1)+o3inp(m)*0.5*(fpres(m-1)-fpres(m+1))
      end if
    end do
    if (o3inp(1)>1.E34) then
      o3sum(1)=o3sum(2)
    else
      o3sum(1)=o3sum(2)+o3inp(1)*(fpres(1)+0.01*ps(iq)-prf(1)-0.5*sum(fpres(1:2)))
    end if
        
    ! vertical interpolation
    o3new=0.
    do m=1,ilev
      if (prf(m)>fpres(1)) then
        o3new(m)=o3sum(1)
      elseif (prf(m)<fpres(nlev)) then
        o3new(m)=o3sum(nlev)
      else
        do k1=2,nlev-1
          if (prf(m)>fpres(k1)) exit
        end do
        fp=(prf(m)-fpres(k1))/(fpres(k1-1)-fpres(k1))
        o3new(m)=(1.-fp)*o3sum(k1)+fp*o3sum(k1-1)
      end if
    end do        
         
    ! output ozone (invert levels)
    out(iq,ilev)=(o3sum(1)-o3new(2))/(0.01*ps(iq)-0.5*sum(prf(1:2)))
    do m=2,ilev-1
      out(iq,ilev-m+1)=2.*(o3new(m)-o3new(m+1))/(prf(m-1)-prf(m+1))
    end do
    out(iq,1)=2.*o3new(ilev)/sum(prf(ilev-1:ilev))
    out(iq,:)=max(out(iq,:),mino3)
  end if
end do
      
return
end subroutine fieldinterpolate

!--------------------------------------------------------------------
! This subroutine interpolates input fields to the CCAM model grid
subroutine o3regrid(o3pre,o3mth,o3nxt,o3dum,o3lon,o3lat)
      
use latlong_m

implicit none

include 'newmpar.h'
include 'const_phys.h'

real, dimension(:,:), intent(out) :: o3pre, o3mth, o3nxt
real, dimension(:,:,:,:), intent(in) :: o3dum
real, dimension(:), intent(in) :: o3lon
real, dimension(:), intent(in) :: o3lat
real, dimension(size(o3dum,3)) :: o3tmp, b, c, d
real, dimension(ifull) :: blon, blat
real alonx, lonadj, serlon, serlat
integer nlon, nlat, nlev
integer iq, l, ilon, ilat, ip

nlon = size(o3dum,1)
nlat = size(o3dum,2)
nlev = size(o3dum,3)

blon(1:ifull) = rlongg(1:ifull)*180./pi
where ( blon(1:ifull) < 0. )
  blon(1:ifull) = blon(1:ifull) + 360.
end where
blat(1:ifull) = rlatt(1:ifull)*180./pi

do iq = 1,ifull
        
  alonx = blon(iq)
  if ( alonx < o3lon(1) ) then
    alonx = alonx + 360.
    ilon = nlon
  else
    do ilon = 1,nlon-1
      if ( o3lon(ilon+1) > alonx ) exit
    end do
  end if
  ip = ilon + 1
  lonadj = 0.
  if ( ip > nlon ) then
    ip = 1
    lonadj = 360.
  end if
  serlon = (alonx-o3lon(ilon))/(o3lon(ip)+lonadj-o3lon(ilon))

  if ( blat(iq) < o3lat(1) ) then
    ilat = 1
    serlat = 0.
  else if ( blat(iq) > o3lat(nlat) ) then
    ilat = nlat - 1
    serlat = 1.
  else
    do ilat = 1,nlat-1
      if ( o3lat(ilat+1) > blat(iq) ) exit
    end do
    serlat = (blat(iq)-o3lat(ilat))/(o3lat(ilat+1)-o3lat(ilat))  
  end if

  ! spatial interpolation
  do l = 1,3
    d(1:nlev) = o3dum(ip,ilat+1,1:nlev,l) - o3dum(ilon,ilat+1,1:nlev,l) &
              - o3dum(ip,ilat,1:nlev,l) + o3dum(ilon,ilat,1:nlev,l)
    b(1:nlev) = o3dum(ip,ilat,1:nlev,l)-o3dum(ilon,ilat,1:nlev,l)
    c(1:nlev) = o3dum(ilon,ilat+1,1:nlev,l)-o3dum(ilon,ilat,1:nlev,l)
    o3tmp(1:nlev) = b(1:nlev)*serlon + c(1:nlev)*serlat + d(1:nlev)*serlon*serlat &
                  + o3dum(ilon,ilat,1:nlev,l)
    where ( o3dum(ilon,ilat,1:nlev,l) > 1.E34 .or. o3dum(ip,ilat,1:nlev,l) > 1.E34 .or.   &
            o3dum(ilon,ilat+1,1:nlev,l) > 1.E34 .or. o3dum(ip,ilat+1,1:nlev,l) > 1.E34 )
      o3tmp(1:nlev) = 1.E35
    end where
             
    select case(l)
      case(1)
        o3pre(iq,1:nlev) = o3tmp(1:nlev)
      case(2)
        o3mth(iq,1:nlev) = o3tmp(1:nlev)
      case(3)
        o3nxt(iq,1:nlev) = o3tmp(1:nlev)
    end select
        
  end do

  ! avoid interpolating missing values
  where ( o3pre(iq,1:nlev) > 1.E34 .or. o3mth(iq,1:nlev) > 1.E34 .or. o3nxt(iq,1:nlev) > 1.E34 )
    o3pre(iq,1:nlev) = 1.E35
    o3mth(iq,1:nlev) = 1.E35
    o3nxt(iq,1:nlev) = 1.E35
  end where

end do

return
end subroutine o3regrid

subroutine resetd(x1,x2,x3,x4,n)

implicit none

integer n, i
real x1(n), x2(n), x3(n), x4(n)
real avg, a1, b1, b2

do i = 1,n
  avg = 0.25*(x1(i)+x2(i)+x3(i)+x4(i))
  a1 = 0.5*(x2(i)-x4(i))
  b1 = 0.5*(x1(i)-x3(i))
  b2 = 0.25*((x1(i)+x3(i))-(x2(i)+x4(i)))
  x1(i) = avg
  x2(i) = a1
  x3(i) = b1
  x4(i) = b2
end do

return
end subroutine resetd

subroutine o3read_amip

!  Read the AMIP2 lat-height ozone data (volume mixing ratio)
!  Based on the supplied example program to read the data.
!  This routine is fixed format f90

use cc_mpi

implicit none
      
include 'filnames.h'

character(len=120) :: label
real, parameter :: amd = 28.9644, amo = 48.0000
real, parameter :: massratio = amo / amd
integer :: lato3d  ! = jg  number of data latitudinal grid
integer :: layo3d  ! = kg  number of data vertical layers
integer :: lvlo3d  ! = lg  number of data vertical layer interfaces
integer :: j, k, month
integer, parameter :: un = 16  ! Unit number for file
real, dimension(kg) :: galt, gprs
real, dimension(lg) :: gali
real :: o3vubc  ! =     upper boundary o3 vmr (ppmv)
real :: prsubc  ! =     upper boundary pressure (mb)
real :: altubc  ! =     upper boundary altitude (km)

allocate(glat(jg),dp(kg),gpri(lg),gdat(jg,kg,mo))

if (myid==0) then

  open(un,file=o3file,form='formatted',status='old')
  read (un,"(i3)") lato3d
  read (un,"(i3)") layo3d
  lvlo3d  =    1 + layo3d
  if ( lato3d /= jg ) then
    write(6,*) "Error in horizontal resolution of ozone data"
    call ccmpi_abort(-1)
  end if
  if ( layo3d /= kg ) then
    write(6,*) "Error in vertical resolution of ozone data"
    call ccmpi_abort(-1)
  end if

  read (un,"(3(1pe12.5))") o3vubc, prsubc, altubc
  read (un,"(a)") label   ! grid latitudes (deg) --------->
  read (un,"(10(1pe12.5))") (glat(j),j=1,lato3d)
  read (un,"(a)") label   ! layer pressure (mb) ---------->
  read (un,"(10(1pe12.5))") (gprs(k),k=1,layo3d)
  read (un,"(a)") label   ! layer altitude (km) ---------->
  read (un,"(10(1pe12.5))") (galt(k),k=1,layo3d)
  read (un,"(a)") label   ! interface pressure (mb) ------>
  read (un,"(10(1pe12.5))") (gpri(k),k=1,lvlo3d)
  read (un,"(a)") label   ! interface altitude (km) ------>
  read (un,"(10(1pe12.5))") (gali(k),k=1,lvlo3d)
  do month = 1,mo
    read (un,"(a)") label   ! monthly o3 vmr (ppmv) lat-alt
    read (un,"(10(1pe12.5))") ((gdat(j,k,month),j=1,lato3d),k=1,layo3d)
  enddo
  close(un)
end if

call ccmpi_bcast(glat,0,comm_world)
call ccmpi_bcast(gpri,0,comm_world)
call ccmpi_bcast(gdat,0,comm_world)

dp = gpri(2:lg) - gpri(1:lg-1)

! Convert from ppmv to mass mixing ratio
gdat = gdat * 1.e-6 * massratio

return
end subroutine o3read_amip

subroutine o3set_amip ( alat, npts, mins, sigh, ps, qo3 )

! Interpolate the AMIP2 ozone data in time, latitude and pressure to 
! produce ozone mixing ratio for model radiation code. 

! In CCAM version latitude may be different for every point.

implicit none

include 'newmpar.h'
include 'const_phys.h'
include 'dates.h'

integer, intent(in) :: npts
real, dimension(npts), intent(in) :: alat     ! Latitude
integer, intent(in) :: mins  ! Time
real, intent(in),  dimension(kl+1) :: sigh ! Half level sigma
real, intent(in),  dimension(npts) :: ps  ! Surface pressure
real, intent(out), dimension(npts,kl) :: qo3  ! Ozone mixing ratio
! Day no of middle of month
real, dimension(12) :: monmid
real, parameter, dimension(12) :: oldmid = (/ 15.5, 45.0, 74.5, 105.0, 135.5, 166.0, 196.5, 227.5, 258.0, 288.5, 319.0, 349.5 /)
real, dimension(kg) :: oz    ! Column for this date and lat.
real, dimension(lg) :: ozcol ! Integrated column
real, dimension(kl+1) :: qo3p
integer :: k1, jyear
integer :: j, m1, m2, j1, j2, k
real, dimension(kl+1) :: prh ! Half level pressure
real :: fac1, fac2, tfac1, tfac2, date, theta

integer leap
common/leap_yr/leap  ! 1 to allow leap years

monmid=oldmid
if ( leap==1 ) then
  jyear =kdate/10000
  if ( mod(jyear,4)==0 ) then
    monmid(2)=oldmid(2)+0.5
    monmid(3:12)=oldmid(3:12)+1.
  end if
  if ( mod(jyear,100)==0 ) then
    monmid(2:12)=oldmid(2:12)
  end if
  if ( mod(jyear,400)==0 ) then
    monmid(2)=oldmid(2)+0.5
    monmid(3:12)=oldmid(3:12)+1.
  end if
end if

! Time interpolation factors
date = real(modulo(mins,nint(monmid(12)+15.5)*1440)) / 1440.0
if ( date <= monmid(1) ) then ! First half of Jan
  m1 = 12
  m2 = 1
  tfac1 = (monmid(1) - date)/31.0
  tfac2 = 1.0 - tfac1
else if ( date >= monmid(12) ) then
  m1 = 12
  m2 = 1
  tfac2 = (date - monmid(12))/31.0
  tfac1 = 1 - tfac2
else
! Search for bracketing dates, i such that monmid(i) >= date
  do m2=2,12
    if ( monmid(m2) >= date ) exit
  end do
  m1 = m2 - 1
  tfac1 = ( monmid(m2) - date ) / ( monmid(m2) - monmid(m1) )
  tfac2 = 1.0 - tfac1
end if

do j=1,npts
! Factors for interpolation in latitude, 
! Note that ozone grid latitudes run S to N.
  theta = alat(j)*180.0/pi
  if ( theta <= glat(1) ) then
    j1 = 1
    j2 = 1
    fac1 = 1.0
    fac2 = 0.0
  else if ( theta >= glat(jg) ) then
    j1 = jg
    j2 = jg
    fac1 = 1.0
    fac2 = 0.0
  else
!   Input data isn't exactly equally spaced, so search to find the spanning
!   latitudes. Want to find first j such that glat(j) >= theta
    j2 = 1 + count ( glat < theta )
    j1 = j2 - 1
    fac1 = ( glat(j2) - theta ) / ( glat(j2) - glat(j1) )
    fac2 = 1.0 - fac1
  end if

! Now interpolate in latitude and time.
  oz = tfac1 * ( fac1*gdat(j1,:,m1) + fac2*gdat(j2,:,m1 ) ) + tfac2 * ( fac1*gdat(j1,:,m2) + fac2*gdat(j2,:,m2 ) ) 

! Each longitude has same vertical profile as a function of pressure,
! but model levels depend on surface pressure so each must be 
! interpolated separately. To ensure the column amount is maintained
! correctly, first calculate ozone column between every half level and
! TOA.
  ozcol(1) = 0.0
  do k=1,kg
    ozcol(k+1) = ozcol(k) + oz(k)*dp(k)
  end do

! To calculate model amounts, find the data half level immediately above
! the given model half level. Note that model levels start at surface while
! data levels start at the top.
!
  prh = sigh * ps(j)*0.01 ! Input PS in Pa

  qo3p = 0.
  qo3p(kl+1) = 0.0      ! TOA
  do k=kl,1,-1          ! Start at the top model level

!   Find largest k1 such that gpri(k1) <= prh(k)
!   Restrict to range 1:kg so as not to go out of bounds of oz
    k1 = max(1, count( gpri(1:kg) <= prh(k) ) )
    qo3p(k) = ozcol(k1) + (prh(k)-gpri(k1))*oz(k1)
  end do

! Finally get the model ozone mixing ratios by differencing the half level
! path amounts. The vertical order of the qo3 array is reversed to 
! match the radiation code.
  do k=1,kl
    qo3(j,kl+1-k) = ( qo3p(k+1) - qo3p(k) ) / ( prh(k+1) - prh(k) )
  end do
end do

return
end subroutine o3set_amip

end module ozoneread
