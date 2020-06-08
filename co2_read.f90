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
    
! Read CO2 and other GHG concentrations.
    
! These routines are designed for processing CMIP3, CMIP5 and CMIP6
! emission scenarios.
    
module co2_read_m

implicit none

private
public co2_read
    
contains

subroutine co2_read(sigma,jyear,csolar)

!  This routine reads the CO2 transmission coefficients from the
!  co2_datafile (filename set in namelist)

use cc_mpi
use co2dta_m
use filnames_m
use infile
use newmpar_m
use parm_m
use radisw_m
      
implicit none
      
include 'rdparm.h'
      
real, intent(inout) :: csolar
real, dimension(kl), intent(in) :: sigma
real, parameter :: sigtol=1e-3
real, dimension(kl) :: sigin
real, dimension(35) :: rcn
real, dimension(8) :: rdum
integer, intent(in) :: jyear
integer k, i, ierr, nlev, iyr
integer ncidrad, ncidsolar, ncidch4, ncidn2o
integer ncidcfc11, ncidcfc12, ncidcfc113, ncidhcfc22
integer, parameter :: lu=15
      
if (nrad==5) then
    
  ! SEA-ESF radiation code
  if (myid==0) then
    if ( radfile==' ' ) then ! default
      write(6,*) "Using default CO2 for radiative data"  
      rrvco2  = 0.33e-3
      rrvch4  = 0.
      rrvn2o  = 0.
      rrvf11  = 0.
      rrvf12  = 0.
      rrvf113 = 0.
      rrvf22  = 0.
    else  

      ! Attempt to open netcdf file
      call ccnf_open(radfile,ncidrad,ierr)
      if ( ierr==0 ) then
          
        ! CMIP6 format  
        write(6,*) 'CO2 data read from file ',trim(radfile)
        call readrad2d(ncidrad,'mole_fraction_of_carbon_dioxide_in_air',rrvco2)
        rrvco2 = rrvco2*1.e-6
        call ccnf_close(ncidrad)
        
        write(6,*) 'CH4 data read from file ',trim(ch4file)
        call ccnf_open(ch4file,ncidch4,ierr)
        if ( ierr/=0 ) then
          write(6,*) "ERROR: Cannot read CH4 data from file ",trim(ch4file)
          call ccmpi_abort(-1)
        end if
        call readrad2d(ncidch4,'mole_fraction_of_methane_in_air',rrvch4)
        rrvch4 = rrvch4*1.e-9
        call ccnf_close(ncidch4)
        
        write(6,*) 'N2O data read from file ',trim(n2ofile)
        call ccnf_open(n2ofile,ncidn2o,ierr)
        if ( ierr/=0 ) then
          write(6,*) "ERROR: Cannot read N2O data from file ",trim(n2ofile)
          call ccmpi_abort(-1)
        end if
        call readrad2d(ncidn2o,'mole_fraction_of_nitrous_oxide_in_air',rrvn2o)
        rrvn2o = rrvn2o*1.e-9
        call ccnf_close(ncidn2o)
        
        write(6,*) 'CFC11 data read from file ',trim(cfc11file)
        call ccnf_open(cfc11file,ncidcfc11,ierr)
        if ( ierr/=0 ) then
          write(6,*) "ERROR: Cannot read CFC11 data from file ",trim(cfc11file)
          call ccmpi_abort(-1)
        end if
        call readrad2d(ncidcfc11,'mole_fraction_of_cfc11_in_air',rrvf11)
        rrvf11 = rrvf11*1.e-12
        call ccnf_close(ncidcfc11)
                
        write(6,*) 'CFC12 data read from file ',trim(cfc12file)
        call ccnf_open(cfc12file,ncidcfc12,ierr)
        if ( ierr/=0 ) then
          write(6,*) "ERROR: Cannot read CFC12 data from file ",trim(cfc12file)
          call ccmpi_abort(-1)
        end if
        call readrad2d(ncidcfc12,'mole_fraction_of_cfc12_in_air',rrvf12)
        rrvf12 = rrvf12*1.e-12
        call ccnf_close(ncidcfc12)
        
        write(6,*) 'CFC113 data read from file ',trim(cfc113file)
        call ccnf_open(cfc113file,ncidcfc113,ierr)
        if ( ierr/=0 ) then
          write(6,*) "ERROR: Cannot read CFC113 data from file ",trim(cfc113file)
          call ccmpi_abort(-1)
        end if
        call readrad2d(ncidcfc113,'mole_fraction_of_cfc113_in_air',rrvf113)
        rrvf113 = rrvf113*1.e-12
        call ccnf_close(ncidcfc113)
        
        write(6,*) 'HCFC22 data read from file ',trim(hcfc22file)
        call ccnf_open(hcfc22file,ncidhcfc22,ierr)
        if ( ierr/=0 ) then
          write(6,*) "ERROR: Cannot read HCFC22 data from file ",trim(hcfc22file)
          call ccmpi_abort(-1)
        end if
        call readrad2d(ncidhcfc22,'mole_fraction_of_hcfc22_in_air',rrvf22)
        rrvf22 = rrvf22*1.e-12
        call ccnf_close(ncidhcfc22)
        
        write(6,*) 'Solar data read from file ',trim(solarfile)
        call ccnf_open(solarfile,ncidsolar,ierr)
        if ( ierr/=0 ) then
          write(6,*) "ERROR: Cannot read solar data from file ",trim(solarfile)
          call ccmpi_abort(-1)
        end if
        call readrad1d(ncidsolar,'tsi',csolar)
        call ccnf_close(ncidsolar)

      else    
        ! check for ASCII file  
        write(6,*) 'Radiative data read from file ',trim(radfile)
        open(lu,file=radfile,form='formatted',status='old')
        
        nlev=0
        read(lu,*,iostat=ierr) nlev
        if (nlev>0) then ! old format
        
          ! CMIP3 emission data
          read(lu,*) (sigin(1),i=nlev,1,-1)
          read(lu,*) rrvco2
          rrvch4=0.
          rrvn2o=0.
          rrvf11=0.
          rrvf12=0.
          rrvf113=0.
          rrvf22=0.
        else
        
          ! CMIP5 emission data
          iyr=-9999
          do while (iyr<jyear)
            read(lu,*,iostat=ierr) iyr,rcn(1:35)
            if (ierr<0) then
              write(6,*) "ERROR: Cannot find concentration data"
              call ccmpi_abort(-1)
            end if
          end do
          rrvco2=rcn(3)*1.E-6
          rrvch4=rcn(4)*1.E-9
          rrvn2o=rcn(5)*1.E-9
          rrvf11=rcn(20)*1.E-12
          rrvf12=rcn(21)*1.E-12
          rrvf113=rcn(22)*1.E-12
          rrvf22=rcn(27)*1.E-12            
        end if ! nlev>0 ..else..
        close(lu)
      end if ! ierr==0 ..else..  
    end if ! radfile==' ' ..else..  
    write(6,*) ' CO2  mixing ratio is ',rrvco2*1e6,' ppmv'
    write(6,*) ' CH4  mixing ratio is ',rrvch4*1e9,' ppbv'
    write(6,*) ' N2O  mixing ratio is ',rrvn2o*1e9,' ppbv'
    write(6,*) ' F11  mixing ratio is ',rrvf11*1e12,' pptv'
    write(6,*) ' F12  mixing ratio is ',rrvf12*1e12,' pptv'
    write(6,*) ' F113 mixing ratio is ',rrvf113*1e12,' pptv'
    write(6,*) ' F22  mixing ratio is ',rrvf22*1e12,' pptv'
    write(6,*) ' Solar constant is    ',csolar,' W/m2'
    rdum(1)=rrvco2
    rdum(2)=rrvch4
    rdum(3)=rrvn2o
    rdum(4)=rrvf11
    rdum(5)=rrvf12
    rdum(6)=rrvf113
    rdum(7)=rrvf22
    rdum(8)=csolar
  end if ! myid==0
  call ccmpi_bcast(rdum(1:8),0,comm_world)
  rrvco2=rdum(1)
  rrvch4=rdum(2)
  rrvn2o=rdum(3)
  rrvf11=rdum(4)
  rrvf12=rdum(5)
  rrvf113=rdum(6)
  rrvf22=rdum(7)
  csolar=rdum(8)

else
      
  ! Older GFDL radiation code (LH/SF)
  if (myid==0) then 
    if ( radfile==' ' ) then
      write(6,*) "ERROR: Must specify radfile for nrad=4"
      call ccmpi_abort(-1)
    end if    
    
    write(6,*) 'Radiative data read from file ',trim(radfile)
    open(lu,file=radfile,form='formatted',status='old')
    
    read(lu,*) nlev
    write(6,*)'co2_read nlev=',nlev
!   Check that the number of levels is the same
    if ( nlev/=kl ) then
      write(6,*) ' ERROR - Number of levels wrong in co2_data file'
      call ccmpi_abort(-1)
    end if
!   Check that the sigma levels are the same
!   Note that the radiation data has the levels in the reverse order
    read(lu,*) (sigin(i),i=kl,1,-1)
    write(6,*)'co2_read sigin=',sigin
    do k=1,kl
      if ( abs(sigma(k)-sigin(k)) > sigtol ) then
        write(6,*) ' ERROR - sigma level wrong in co2_data file'
        write(6,*) k, sigma(k), sigin(k)
        call ccmpi_abort(-1)
      end if
    end do
    read(lu,*) rrvco2
    write(6,*) ' CO2 mixing ratio is ', rrvco2*1e6,' ppmv'
    read(lu,*) stemp
    read(lu,*) gtemp
    read(lu,*) cdt51
    read(lu,*) co251
    read(lu,*) c2d51
    read(lu,*) cdt58
    read(lu,*) co258
    read(lu,*) c2d58
    read(lu,*) cdtm51
    read(lu,*) co2m51
    read(lu,*) c2dm51
    read(lu,*) cdtm58
    read(lu,*) co2m58
    read(lu,*) c2dm58
    read(lu,*) cdt31
    read(lu,*) co231
    read(lu,*) c2d31
    read(lu,*) cdt38
    read(lu,*) co238
    read(lu,*) c2d38
    read(lu,*) cdt71
    read(lu,*) co271
    read(lu,*) c2d71
    read(lu,*) cdt78
    read(lu,*) co278
    read(lu,*) c2d78
    read(lu,*) co211
    read(lu,*) co218
    close(lu)
  end if
  rdum(1)=rrvco2
  call ccmpi_bcast(rdum(1:1),0,comm_world)
  rrvco2=rdum(1)
  call ccmpi_bcast(stemp,0,comm_world)
  call ccmpi_bcast(gtemp,0,comm_world)
  call ccmpi_bcast(cdt51,0,comm_world)
  call ccmpi_bcast(co251,0,comm_world)
  call ccmpi_bcast(c2d51,0,comm_world)
  call ccmpi_bcast(cdt58,0,comm_world)
  call ccmpi_bcast(co258,0,comm_world)
  call ccmpi_bcast(c2d58,0,comm_world)
  call ccmpi_bcast(cdtm51,0,comm_world)
  call ccmpi_bcast(co2m51,0,comm_world)
  call ccmpi_bcast(c2dm51,0,comm_world)
  call ccmpi_bcast(cdtm58,0,comm_world)
  call ccmpi_bcast(co2m58,0,comm_world)
  call ccmpi_bcast(c2dm58,0,comm_world)
  call ccmpi_bcast(cdt31,0,comm_world)
  call ccmpi_bcast(co231,0,comm_world)
  call ccmpi_bcast(c2d31,0,comm_world)
  call ccmpi_bcast(cdt38,0,comm_world)
  call ccmpi_bcast(co238,0,comm_world)
  call ccmpi_bcast(c2d38,0,comm_world)
  call ccmpi_bcast(cdt71,0,comm_world)
  call ccmpi_bcast(co271,0,comm_world)
  call ccmpi_bcast(c2d71,0,comm_world)
  call ccmpi_bcast(cdt78,0,comm_world)
  call ccmpi_bcast(co278,0,comm_world)
  call ccmpi_bcast(c2d78,0,comm_world)
  call ccmpi_bcast(co211,0,comm_world)
  call ccmpi_bcast(co218,0,comm_world)
end if      
      
return
end subroutine co2_read

subroutine readrad2d(ncid,vname,rrvc)

use cc_mpi
use infile
use stime_m

implicit none

integer, intent(in) :: ncid
integer kdate_rsav, ktime_rsav
integer idvtime, iarchi, kdate_r, ktime_r
integer maxarchi
integer(kind=8) :: mtimer
integer, dimension(2) :: nstart, ncount
real, dimension(1) :: rdat
real, intent(out) :: rrvc
real timer
logical ltest
character(len=*), intent(in) :: vname
character(len=80) datestring

call ccnf_inq_dimlen(ncid,'time',maxarchi)
call ccnf_inq_varid(ncid,'time',idvtime)
call ccnf_get_att(ncid,idvtime,'units',datestring)
call processdatestring(datestring,kdate_rsav,ktime_rsav)
ltest = .true.
iarchi = 0
do while ( ltest .and. iarchi<maxarchi )
  iarchi = iarchi + 1  
  kdate_r = kdate_rsav
  ktime_r = ktime_rsav
  call ccnf_get_vara(ncid,idvtime,iarchi,timer)
  mtimer = nint(timer,8)*1440_8 ! units=days
  call datefix(kdate_r,ktime_r,mtimer,allleap=0,silent=.true.)
  ltest = (kdate_r/100-kdate_s/100)<0
end do
if ( ltest ) then
  write(6,*) "ERROR: Search failed with ltest,iarchi = ",ltest,iarchi
  write(6,*) "kdate_r = ",kdate_r
  call ccmpi_abort(-1)
end if
nstart(1:2) = (/ 1, iarchi  /)
ncount(1:2) = (/ 1, 1 /)
call ccnf_get_vara(ncid,vname,nstart,ncount,rdat)
rrvc = rdat(1)

end subroutine readrad2d
    
subroutine readrad1d(ncid,vname,csolar)

use cc_mpi
use infile
use stime_m

implicit none

integer, intent(in) :: ncid
integer kdate_rsav, ktime_rsav
integer idvtime, iarchi, kdate_r, ktime_r
integer maxarchi
integer(kind=8) mtimer
integer, dimension(1) :: nstart, ncount
real, dimension(1) :: rdat
real, intent(out) :: csolar
real timer
logical ltest
character(len=*), intent(in) :: vname
character(len=80) datestring

call ccnf_inq_dimlen(ncid,'time',maxarchi)
call ccnf_inq_varid(ncid,'time',idvtime)
call ccnf_get_att(ncid,idvtime,'units',datestring)
call processdatestring(datestring,kdate_rsav,ktime_rsav)
ltest = .true.
iarchi = 0
do while ( ltest .and. iarchi<maxarchi )
  iarchi = iarchi + 1  
  kdate_r = kdate_rsav
  ktime_r = ktime_rsav
  call ccnf_get_vara(ncid,idvtime,iarchi,timer)
  mtimer = nint(timer,8)*1440_8 ! units=days
  call datefix(kdate_r,ktime_r,mtimer,allleap=0,silent=.true.)
  ltest = (kdate_r/100-kdate_s/100)<0
end do
if ( ltest ) then
  write(6,*) "ERROR: Search failed with ltest,iarchi = ",ltest,iarchi
  write(6,*) "kdate_r = ",kdate_r
  call ccmpi_abort(-1)
end if
nstart(1) = iarchi
ncount(1) = 1
call ccnf_get_vara(ncid,vname,nstart,ncount,rdat)
csolar = rdat(1)

end subroutine readrad1d

end module co2_read_m