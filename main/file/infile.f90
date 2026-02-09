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
    
module infile
      
! This module contains routines for reading netcdf files, vertical interpolation and some calendar functions.
! This module also contains the interface to the netcdf library.
      
! This version of infile.f90 supports parallel input files.  Multiple processors read as many input files as
! supplied in parallel and then makes this data avaliable to all processors for interpolation. The code can also identify
! restart files, in which case no additional message passing is required.

use netcdf_m

implicit none
            
private
public datefix, datefix_month, getzinp, ncmsg, processdatestring
public ptest, pfall, ncidold, resprocformat, pncid
public histopen, histclose, histrd, surfread
public attrib, histwrt, lsmask
public ccnf_open, ccnf_create, ccnf_close, ccnf_sync, ccnf_enddef
public ccnf_redef, ccnf_nofill, ccnf_read
public ccnf_inq_varid, ccnf_inq_dimid, ccnf_inq_dimlen, ccnf_inq_varndims
public ccnf_inq_exist
public ccnf_def_dim, ccnf_def_dimu, ccnf_def_var
public ccnf_get_vara, ccnf_get_att, ccnf_get_attg
public ccnf_put_vara, ccnf_put_att, ccnf_put_attg
public comm_ip, pil_single

public driving_model_id, driving_model_ensemble_number, driving_experiment_name
public driving_institution_id

public any_m, daily_m, sixhr_m, point_m, tmean_m, max_m, min_m, fixed_m, sum_m
public short_m, double_m, float_m, anotdef_m, amean_m, land_m, sea_m, seaice_m

public time_of_month, time_interpolate

integer(kind=4), dimension(:), allocatable, save :: pncid
integer, dimension(:), allocatable, save :: pprid
integer, dimension(:), allocatable, save :: ppanid, ppiid, ppjid
integer, save :: ncidold = -1
integer, save :: comm_ip
integer, save :: pil_single = 48 ! grid size for single file decomposition
logical, dimension(:), allocatable, save :: pfown
logical, save :: ptest, pfall, resprocformat
integer(kind=2), parameter :: minv = -32500
integer(kind=2), parameter :: maxv =  32500
integer(kind=2), parameter :: missval = -32501

! host metadata for CORDEX
character(len=256), save :: driving_model_id = ' '
character(len=256), save :: driving_institution_id = ' '
character(len=256), save :: driving_model_ensemble_number = ' '
character(len=256), save :: driving_experiment_name = ' '

! define configuration for attrib
integer, parameter :: any_m = 0      ! output is valid at any time-step (usually hourly)
integer, parameter :: daily_m = 1    ! output is daily only
integer, parameter :: sixhr_m = 2    ! output is 6-hourly only
integer, parameter :: point_m = 0    ! output is for a point in time
integer, parameter :: tmean_m = 1    ! output is averaged in time
integer, parameter :: max_m = 2      ! output is maximum over a time interval
integer, parameter :: min_m = 3      ! output is minimum over a time interval
integer, parameter :: fixed_m = 4    ! output is fixed (independent of time)
integer, parameter :: sum_m = 5      ! output is a sum over time interval
integer, parameter :: short_m = 1    ! output is compressed as a short
integer, parameter :: double_m = 2   ! output is double precision
integer, parameter :: float_m = -1   ! output is single precision (or default model precision)
integer, parameter :: anotdef_m = -1 ! output has no area mean (e.g., index)
integer, parameter :: amean_m = 0    ! output is mean over area
integer, parameter :: land_m = 1     ! output is mean over land
integer, parameter :: sea_m = 2      ! output is mean over sea
integer, parameter :: seaice_m = 3   ! output is mean over sea_ice

interface histrd
  module procedure histrd3r4, histrd4r4, histrd5r4
#ifndef i8r8  
  module procedure histrd3r8, histrd4r8, histrd5r8
#endif
end interface histrd

interface histwrt
  module procedure histwrt3r4, histwrt4r4, histwrt5r4
#ifndef i8r8
  module procedure histwrt3r8, histwrt4r8, histwrt5r8
#endif
end interface histwrt

interface ccnf_get_att
  module procedure ccnf_get_att_text, ccnf_get_att_real
end interface ccnf_get_att

interface ccnf_get_attg
  module procedure ccnf_get_att_intg1i, ccnf_get_att_intg2i
  module procedure ccnf_get_att_realg1r, ccnf_get_att_realg2r
  module procedure ccnf_get_att_textg
end interface ccnf_get_attg

interface ccnf_get_vara
  module procedure ccnf_get_vara_text1t_t
  module procedure ccnf_get_var_real2r_t
  module procedure ccnf_get_vara_real2r_0
  module procedure ccnf_get_vara_real1r_s, ccnf_get_vara_real2r_s, ccnf_get_vara_real3r_s, ccnf_get_vara_real4r_s
  module procedure ccnf_get_vara_real1r_t, ccnf_get_vara_real2r_t, ccnf_get_vara_real3r_t
  module procedure ccnf_get_var_int2i_t
  module procedure ccnf_get_vara_int1i_s, ccnf_get_vara_int2i_s
  module procedure ccnf_get_vara_int2i_t, ccnf_get_vara_int3i_t
#ifndef i8r8
  module procedure ccnf_get_vara_double1r_s, ccnf_get_vara_double4r_s 
  module procedure ccnf_get_vara_double2r_t
#endif
end interface ccnf_get_vara

interface ccnf_def_var
  module procedure ccnf_def_var_s, ccnf_def_var_v
end interface ccnf_def_var

interface ccnf_put_att
  module procedure ccnf_put_att_text, ccnf_put_att_rs
end interface ccnf_put_att

interface ccnf_put_attg
  module procedure ccnf_put_att_realg1, ccnf_put_att_realg2
  module procedure ccnf_put_att_intg1, ccnf_put_att_intg2
  module procedure ccnf_put_att_textg
end interface ccnf_put_attg

interface ccnf_put_vara
  module procedure ccnf_put_var_text2t_0
  module procedure ccnf_put_var_int2i_0, ccnf_put_var_int3i_0
  module procedure ccnf_put_vara_int1i_s, ccnf_put_vara_int2i_s
  module procedure ccnf_put_vara_int1i_t
  module procedure ccnf_put_vara_real2r_0
  module procedure ccnf_put_vara_real1r_s, ccnf_put_vara_real2r_s, ccnf_put_vara_real3r_s
  module procedure ccnf_put_vara_real1r_t, ccnf_put_vara_real2r_t, ccnf_put_vara_real3r_t
#ifndef i8r8
  module procedure ccnf_put_vara_double1r_s, ccnf_put_vara_double2r_s
  module procedure ccnf_put_vara_double1r_t
#endif
end interface ccnf_put_vara

contains

!--------------------------------------------------------------
! Interface for reading 2D+time fields
subroutine histrd3r4(iarchi,ier,name,var,ifull)
      
use cc_mpi, only : myid, ccmpi_reduce, histrd_begin, histrd_end, fnresid, &
                   start_log, end_log, ccmpi_distribute, pil_g
use parm_m

integer, intent(in) :: iarchi,ifull
integer, intent(out) :: ier
real, dimension(:), intent(inout) :: var ! may be dummy argument from myid/=0
real, dimension(:,:,:), allocatable :: lvar
real, dimension(:,:,:), allocatable :: globvar
character(len=*), intent(in) :: name

call START_LOG(histrd_begin)

allocate( lvar(size(var),1,1) )
lvar(:,1,1) = var(:)  

if ( ifull==6*pil_g**2 .or. ptest ) then
  ! read local arrays
  call hr5p(iarchi,ier,name,1,1,.true.,lvar)
else
  ! gather and distribute (i.e., change in number of processors) 
  if ( myid==0 ) then
     allocate( globvar(6*pil_g**2,1,1) )
     globvar(:,:,:) = 0.
     call hr5p(iarchi,ier,name,1,1,.false.,globvar)
     call ccmpi_distribute(lvar,globvar)
     deallocate( globvar )
  else
    call hr5p(iarchi,ier,name,1,1,.false.)
    call ccmpi_distribute(lvar)
  end if
end if ! ifull==6*pil_g**2 .or. ptest ..else..

var(:) = lvar(:,1,1)

if ( ier==0 .and. myid==0 ) then
  write(6,'(" done histrd3 ",a8,i4,i3)') trim(name),ier,iarchi
else if ( myid==0 ) then
  write(6,'("***absent field for name,ier: ",a8,i4)') trim(name),ier
end if

deallocate( lvar )

call END_LOG(histrd_end)

return
end subroutine histrd3r4   

#ifndef i8r8
!--------------------------------------------------------------
! Interface for reading 2D+time fields (double precision version)
subroutine histrd3r8(iarchi,ier,name,var,ifull)
      
use cc_mpi, only : myid, ccmpi_reducer8, histrd_begin, histrd_end, fnresid, &
                   start_log, end_log, ccmpi_distributer8, pil_g
use parm_m

integer, intent(in) :: iarchi, ifull
integer, intent(out) :: ier
real(kind=8), dimension(:), intent(inout) :: var ! may be dummy argument from myid/=0
real(kind=8), dimension(:,:,:), allocatable :: lvar
real(kind=8), dimension(:,:,:), allocatable :: globvar
character(len=*), intent(in) :: name

call START_LOG(histrd_begin)

allocate( lvar(size(var),1,1) )
lvar(:,1,1) = var(:)

if ( ifull==6*pil_g**2 .or. ptest ) then
  ! read local arrays without gather and distribute (e.g., restart file)
  call hr5pr8(iarchi,ier,name,1,1,.true.,lvar)
else
  ! read local arrays with gather and distribute (i.e., change in number of processors) 
  if ( myid==0 ) then
    allocate( globvar(6*pil_g**2,1,1) )
    globvar(:,:,:) = 0.
    call hr5pr8(iarchi,ier,name,1,1,.false.,globvar)
    call ccmpi_distributer8(lvar,globvar)
    deallocate( globvar )
  else
    call hr5pr8(iarchi,ier,name,1,1,.false.)
    call ccmpi_distributer8(lvar)
  end if
end if ! ifull==6*pil_g**2 .or. ptest ..else..

var(:) = lvar(:,1,1)

if ( ier==0 .and. myid==0 ) then
  write(6,'(" done histrd3r8 ",a8,i4,i3)') trim(name),ier,iarchi
else if ( myid==0 ) then
  write(6,'("***absent field for name,ier: ",a8,i4)') trim(name),ier
end if

deallocate( lvar )

call END_LOG(histrd_end)

return
end subroutine histrd3r8
#endif

!--------------------------------------------------------------   
! Interface for reading 3D+time fields
subroutine histrd4r4(iarchi,ier,name,var,ifull)
      
use cc_mpi, only : myid, ccmpi_reduce, histrd_begin, histrd_end, fnresid, &
                   start_log, end_log, ccmpi_distribute, ccmpi_abort, pil_g
use parm_m
      
integer, intent(in) :: iarchi, ifull
integer, intent(out) :: ier
integer kk
real, dimension(:,:), intent(inout) :: var ! may be dummy argument from myid/=0
real, dimensioN(:,:,:), allocatable :: lvar
real, dimension(:,:,:), allocatable :: globvar
character(len=*), intent(in) :: name

call START_LOG(histrd_begin)

kk = size(var,2)

allocate( lvar(size(var,1),kk,1) )
lvar(:,:,1) = var(:,:)

if ( ifull==6*pil_g**2 .or. ptest ) then
  ! read local arrays without gather and distribute
  call hr5p(iarchi,ier,name,kk,1,.true.,lvar)
else 
  ! read local arrays with gather and distribute
  if ( myid==0 ) then
    allocate( globvar(6*pil_g**2,kk,1) )
    globvar(:,:,:) = 0.
    call hr5p(iarchi,ier,name,kk,1,.false.,globvar)     
    call ccmpi_distribute(lvar,globvar)
    deallocate( globvar )
  else 
    call hr5p(iarchi,ier,name,kk,1,.false.)
    call ccmpi_distribute(lvar)
  end if
end if

var(:,:) = lvar(:,:,1)

if ( ier==0 .and. myid==0 ) then
  write(6,'(" done histrd4 ",a8,i4,i3)') trim(name),ier,iarchi
else if ( myid==0 ) then
  write(6,'("***absent field for name,ier: ",a8,i4)') trim(name),ier
end if

deallocate( lvar )

call END_LOG(histrd_end)

return
end subroutine histrd4r4

#ifndef i8r8
!--------------------------------------------------------------   
! Interface for reading 3D+time fields (double precision version)
subroutine histrd4r8(iarchi,ier,name,var,ifull)
      
use cc_mpi, only : myid, ccmpi_reducer8, histrd_begin, histrd_end, fnresid, &
                   start_log, end_log, ccmpi_distributer8, ccmpi_abort, pil_g
use parm_m
      
integer, intent(in) :: iarchi, ifull
integer, intent(out) :: ier
integer kk
real(kind=8), dimension(:,:), intent(inout) :: var ! may be dummy argument from myid/=0
real(kind=8), dimensioN(:,:,:), allocatable :: lvar
real(kind=8), dimension(:,:,:), allocatable :: globvar
character(len=*), intent(in) :: name

call START_LOG(histrd_begin)

kk = size(var,2)

allocate( lvar(size(var,1),kk,1) )
lvar(:,:,1) = var(:,:)

if ( ifull==6*pil_g**2 .or. ptest ) then
  ! read local arrays without gather and distribute
  call hr5pr8(iarchi,ier,name,kk,1,.true.,lvar)
else
  ! gather and distribute
  if ( myid==0 ) then
    allocate( globvar(6*pil_g**2,kk,1) )
    globvar(:,:,:) = 0._8
    call hr5pr8(iarchi,ier,name,kk,1,.false.,globvar)
    call ccmpi_distributer8(lvar,globvar)
    deallocate( globvar )
  else 
    call hr5pr8(iarchi,ier,name,kk,1,.false.)
    call ccmpi_distributer8(lvar)
  end if
end if

var(:,:) = lvar(:,:,1)

if ( ier==0 .and. myid==0 ) then
  write(6,'(" done histrd4r8 ",a8,i4,i3)') trim(name),ier,iarchi
else if ( myid==0 ) then
  write(6,'("***absent field for name,ier: ",a8,i4)') trim(name),ier
end if

deallocate( lvar )

call END_LOG(histrd_end)

return
end subroutine histrd4r8
#endif

!--------------------------------------------------------------   
! Interface for reading 4D+time fields
subroutine histrd5r4(iarchi,ier,name,var,ifull)
      
use cc_mpi, only : myid, ccmpi_reduce, histrd_begin, histrd_end, start_log, end_log, &
                   ccmpi_distribute, ccmpi_abort, pil_g
use parm_m
      
integer, intent(in) :: iarchi, ifull
integer, intent(out) :: ier
integer kk, ll
real, dimension(:,:,:), intent(inout) :: var ! may be dummy argument from myid/=0
real, dimension(:,:,:), allocatable :: globvar
character(len=*), intent(in) :: name

call START_LOG(histrd_begin)

kk = size(var,2)
ll = size(var,3)

if ( ifull==6*pil_g**2 .or. ptest ) then
  ! read local arrays without gather and distribute
  call hr5p(iarchi,ier,name,kk,ll,.true.,var)
else    
  ! read local arrays with gather and distribute
  if ( myid==0 ) then
    allocate( globvar(6*pil_g**2,kk,ll) )
    globvar(:,:,:) = 0.
    call hr5p(iarchi,ier,name,kk,ll,.false.,globvar)
    call ccmpi_distribute(var,globvar)
    deallocate( globvar )
  else 
    call hr5p(iarchi,ier,name,kk,ll,.false.)
    call ccmpi_distribute(var)
  end if
end if

if ( ier==0 .and. myid==0 ) then
  write(6,'(" done histrd5 ",a8,i4,i3)') trim(name),ier,iarchi
else if ( myid==0 ) then
  write(6,'("***absent field for name,ier: ",a8,i4)') trim(name),ier
end if

call END_LOG(histrd_end)

return
end subroutine histrd5r4

#ifndef i8r8
!--------------------------------------------------------------   
! Interface for reading 4D+time fields (double precision version)
subroutine histrd5r8(iarchi,ier,name,var,ifull)
      
use cc_mpi, only : myid, ccmpi_reducer8, histrd_begin, histrd_end, start_log, end_log, &
                   ccmpi_distributer8, ccmpi_abort, pil_g
use parm_m
      
integer, intent(in) :: iarchi, ifull
integer, intent(out) :: ier
integer kk, ll
real(kind=8), dimension(:,:,:), intent(inout) :: var
real(kind=8), dimension(:,:,:), allocatable :: globvar
character(len=*), intent(in) :: name

call START_LOG(histrd_begin)

kk = size(var,2)
ll = size(var,3)

if ( ifull==6*pil_g**2 .or. ptest ) then
  ! read local arrays without gather and distribute
  call hr5pr8(iarchi,ier,name,kk,ll,.true.,var)
else    
  ! read local arrays with gather and distribute
  if ( myid==0 ) then
    allocate( globvar(6*pil_g**2,kk,ll) )
    globvar(:,:,:) = 0._8
    call hr5pr8(iarchi,ier,name,kk,ll,.false.,globvar)
    call ccmpi_distributer8(var,globvar)
    deallocate( globvar )
  else 
    call hr5pr8(iarchi,ier,name,kk,ll,.false.)
    call ccmpi_distributer8(var)
  end if
end if

if ( ier==0 .and. myid==0 ) then
  write(6,'(" done histrd5r8 ",a8,i4,i3)') trim(name),ier,iarchi
else if ( myid==0 ) then
  write(6,'("***absent field for name,ier: ",a8,i4)') trim(name),ier
end if

call END_LOG(histrd_end)

return
end subroutine histrd5r8
#endif

!--------------------------------------------------------------
! This subroutine reads 2D/3D/4D+time input files
      
! when qtest=.true. the input grid decomposition should
! match the current processor decomposition.  We can then
! skip the MPI gather and distribute steps.
subroutine hr5p(iarchi,ier,name,kk,ll,qtest,var)

use cc_mpi
      
integer, intent(in) :: iarchi, kk, ll
integer, intent(out) :: ier
integer :: ipf, ca, jpf, ip, n, no, cc, j, k, l
integer(kind=4), dimension(6) :: start, ncount
integer(kind=4) :: idv, ndims
real, dimension(:,:,:), intent(inout), optional :: var
real, dimension(:,:,:), allocatable :: rvar
real, dimension(:,:,:,:), allocatable :: gvar
real(kind=4) :: laddoff, lsf
logical, intent(in) :: qtest
character(len=*), intent(in) :: name
character(len=80) :: newname

ier = 0

if ( mynproc>0 ) then

  allocate( rvar(pipan*pjpan*pnpan,kk,ll) )
    
  do ipf = 0,mynproc-1
      
    rvar(:,:,:) = 0. ! default for missing field  
    
    ! get variable idv
    ier = nf90_inq_varid(pncid(ipf),name,idv)
    if ( ier == nf90_noerr ) then
      ier = nf90_inquire_variable(pncid(ipf),idv,ndims=ndims)
      call ncmsg(name,ier)        
      if ( resprocformat ) then
        select case(ndims)
          case(3)
            start(1:3)  = (/ 1, 1, pprid(ipf) /)
            ncount(1:3) = (/ pipan, pjpan*pnpan, 1 /)
          case(4)
            start(1:4)  = (/ 1, 1, pprid(ipf), iarchi /)
            ncount(1:4) = (/ pipan, pjpan*pnpan, 1, 1 /)
          case(5)
            start(1:5)  = (/ 1, 1, 1, pprid(ipf), iarchi /)
            ncount(1:5) = (/ pipan, pjpan*pnpan, kk, 1, 1 /)
          case(6)
            start(1:6)  = (/ 1, 1, 1, 1, pprid(ipf), iarchi /)
            ncount(1:6) = (/ pipan, pjpan*pnpan, kk, ll, 1, 1 /)
          case default
            write(6,*) "ERROR: Unexpected number of dimensions reading ",trim(name)
            write(6,*) "ndims = ",ndims
            call ccmpi_abort(-1)
        end select    
      else
        select case(ndims)
          case(2)
            start(1:2)  = (/ ppiid(ipf), ppjid(ipf)+ppanid(ipf)*pil_g /)
            ncount(1:2) = (/ pipan, pjpan*pnpan /)   
          case(3)
            start(1:3)  = (/ ppiid(ipf), ppjid(ipf)+ppanid(ipf)*pil_g, iarchi /)
            ncount(1:3) = (/ pipan, pjpan*pnpan, 1 /)   
          case(4)
            start(1:4)  = (/ ppiid(ipf), ppjid(ipf)+ppanid(ipf)*pil_g, 1, iarchi /)
            ncount(1:4) = (/ pipan, pjpan*pnpan, kk, 1 /)   
          case(5)
            start(1:5)  = (/ ppiid(ipf), ppjid(ipf)+ppanid(ipf)*pil_g, 1, 1, iarchi /)
            ncount(1:5) = (/ pipan, pjpan*pnpan, kk, ll, 1 /)   
          case default
            write(6,*) "ERROR: Unexpected number of dimensions reading ",trim(name)
            write(6,*) "ndims = ",ndims
            call ccmpi_abort(-1)
        end select            
      end if    
      ! obtain scaling factors and offsets from attributes
      ier = nf90_get_att(pncid(ipf),idv,'add_offset',laddoff)
      if ( ier/=nf90_noerr ) laddoff = 0.
      ier = nf90_get_att(pncid(ipf),idv,'scale_factor',lsf)
      if ( ier/=nf90_noerr ) lsf = 1.
      ier = nf90_get_var(pncid(ipf),idv,rvar,start=start(1:ndims),count=ncount(1:ndims))
      call ncmsg(name,ier)
      ! unpack data
      rvar(:,:,:) = rvar(:,:,:)*real(lsf) + real(laddoff)
    else ! ier==nf_noerr ..else..
      if ( resprocformat ) then  
        start(1:4) = (/ 1, 1, pprid(ipf), iarchi /)
        ncount(1:4) = (/ pipan, pjpan*pnpan, 1, 1 /)
        ndims = 4
      else
        start(1:3) = (/ ppiid(ipf), ppjid(ipf)+ppanid(ipf)*pil_g, iarchi /)
        ncount(1:3) = (/ pipan, pjpan*pnpan, 1 /)
        ndims = 3
      end if    
      do k = 1,kk        
        write(newname,'("'//trim(name)//'",I3.3)') k
        ier = nf90_inq_varid(pncid(ipf),newname,idv)
        if ( ier/=nf90_noerr .and. k<100 ) then
          write(newname,'("'//trim(name)//'",I2.2)') k
          ier = nf90_inq_varid(pncid(ipf),newname,idv)          
        end if
        if ( ier/=nf90_noerr .and. k<10 ) then
          write(newname,'("'//trim(name)//'",I1.1)') k
          ier = nf90_inq_varid(pncid(ipf),newname,idv)          
        end if
        if ( ier/=nf90_noerr ) then
          exit
        end if
        ! obtain scaling factors and offsets from attributes
        ier = nf90_get_att(pncid(ipf),idv,'add_offset',laddoff)
        if ( ier/=nf90_noerr ) laddoff = 0.
        ier = nf90_get_att(pncid(ipf),idv,'scale_factor',lsf)
        if ( ier/=nf90_noerr ) lsf = 1.
        ier = nf90_get_var(pncid(ipf),idv,rvar(:,k,1),start=start(1:ndims),count=ncount(1:ndims))
        call ncmsg(name,ier)
        ! unpack data
        rvar(:,k,1) = rvar(:,k,1)*real(lsf) + real(laddoff)      
      end do
    end if ! ier==nf_noerr ..else. 

    if ( qtest ) then
      ! e.g., restart file or nogather=.true.
      ca = pipan*pjpan*pnpan*ipf
      var(1+ca:pipan*pjpan*pnpan+ca,1:kk,1:ll) = rvar(:,:,:)
    else
      ! e.g., mesonest file
      if ( myid==0 ) then
        allocate( gvar(pipan*pjpan*pnpan,kk,ll,fnresid) )
        call ccmpi_gatherx(gvar,rvar,0,comm_ip)
        do jpf = 1,fnresid
          ip = ipf*fnresid + jpf - 1   ! local file number
          do l = 1,ll
            do k = 1,kk
              do n = 0,pnpan-1
                no = n - pnoff(ip) + 1   ! global panel number of local file
                ca = pioff(ip,no) + pjoff(ip,no)*pil_g + no*pil_g**2 - pil_g
                cc = n*pipan*pjpan - pipan
                do j = 1,pjpan
                  var(1+j*pil_g+ca:pipan+j*pil_g+ca,k,l) = gvar(1+j*pipan+cc:pipan+j*pipan+cc,k,l,jpf)
                end do
              end do
            end do
          end do
        end do
        deallocate( gvar )
      else
        allocate( gvar(1,1,1,1) )  
        call ccmpi_gatherx(gvar,rvar,0,comm_ip)
        deallocate( gvar )
      end if
    end if ! qtest

  end do ! ipf

  deallocate(rvar)
  
end if ! mynproc>0

return
end subroutine hr5p

#ifndef i8r8
subroutine hr5pr8(iarchi,ier,name,kk,ll,qtest,var)

use cc_mpi
      
integer, intent(in) :: iarchi, kk, ll
integer, intent(out) :: ier
integer :: ipf, ca,jpf, ip, n, no, cc, j, k, l
integer(kind=4), dimension(6) :: start, ncount
integer(kind=4) idv, ndims, dimlen
real(kind=8), dimension(:,:,:), intent(inout), optional :: var
real(kind=8), dimension(:,:,:), allocatable :: rvar
real(kind=8), dimension(:,:,:,:), allocatable :: gvar
real(kind=4) laddoff, lsf
logical, intent(in) :: qtest
character(len=*), intent(in) :: name
character(len=80) :: newname

ier = 0

if ( mynproc>0 ) then
      
  allocate( rvar(pipan*pjpan*pnpan,kk,ll) )  
    
  do ipf = 0,mynproc-1
    
    rvar(:,:,:) = 0._8  
    
    ! get variable idv
    ier = nf90_inq_varid(pncid(ipf),name,idv)
    if ( ier==nf90_noerr ) then
      ier = nf90_inquire_variable(pncid(ipf),idv,ndims=ndims)
      call ncmsg(name,ier)        
      if ( resprocformat ) then
        select case(ndims)  
          case(3)
            start(1:3)  = (/ 1, 1, pprid(ipf) /)
            ncount(1:3) = (/ pipan, pjpan*pnpan, 1 /)
          case(4)
            start(1:4)  = (/ 1, 1, pprid(ipf), iarchi /)
            ncount(1:4) = (/ pipan, pjpan*pnpan, 1, 1 /)
          case(5)  
            start(1:5)  = (/ 1, 1, 1, pprid(ipf), iarchi /)
            ncount(1:5) = (/ pipan, pjpan*pnpan, kk, 1, 1 /)
          case(6)  
            start(1:6)  = (/ 1, 1, 1, 1, pprid(ipf), iarchi /)
            ncount(1:6) = (/ pipan, pjpan*pnpan, kk, ll, 1, 1 /)
          case default
            write(6,*) "ERROR: Unexpected number of dimensions reading ",trim(name)
            write(6,*) "ndims = ",ndims
            call ccmpi_abort(-1)
        end select            
      else
        select case(ndims)
          case(2)
            start(1:2)  = (/ ppiid(ipf), ppjid(ipf)+ppanid(ipf)*pil_g /)
            ncount(1:2) = (/ pipan, pjpan*pnpan /)
          case(3)  
            start(1:3)  = (/ ppiid(ipf), ppjid(ipf)+ppanid(ipf)*pil_g, iarchi /)
            ncount(1:3) = (/ pipan, pjpan*pnpan, 1 /)   
          case(4)  
            start(1:4)  = (/ ppiid(ipf), ppjid(ipf)+ppanid(ipf)*pil_g, 1, iarchi /)
            ncount(1:4) = (/ pipan, pjpan*pnpan, kk, 1 /)   
          case(5)  
            start(1:5)  = (/ ppiid(ipf), ppjid(ipf)+ppanid(ipf)*pil_g, 1, 1, iarchi /)
            ncount(1:5) = (/ pipan, pjpan*pnpan, kk, ll, 1 /)   
          case default
            write(6,*) "ERROR: Unexpected number of dimensions reading ",trim(name)
            write(6,*) "ndims = ",ndims
            call ccmpi_abort(-1)
        end select            
      end if    
      ! obtain scaling factors and offsets from attributes
      ier = nf90_get_att(pncid(ipf),idv,'add_offset',laddoff)
      if ( ier/=nf90_noerr ) laddoff = 0.
      ier = nf90_get_att(pncid(ipf),idv,'scale_factor',lsf)
      if ( ier/=nf90_noerr ) lsf = 1.
      ier = nf90_get_var(pncid(ipf),idv,rvar,start=start(1:ndims),count=ncount(1:ndims))
      call ncmsg(name,ier)
      ! unpack data
      rvar(:,:,:) = rvar(:,:,:)*real(lsf,8) + real(laddoff,8)
    else
      if ( resprocformat ) then  
        start(1:4) = (/ 1, 1, pprid(ipf), iarchi /)
        ncount(1:4) = (/ pipan, pjpan*pnpan, 1, 1 /)
        ndims = 4
      else
        start(1:3) = (/ ppiid(ipf), ppjid(ipf)+ppanid(ipf)*pil_g, iarchi /)
        ncount(1:3) = (/ pipan, pjpan*pnpan, 1 /)
        ndims = 3
      end if    
      do k = 1,kk        
        write(newname,'("'//trim(name)//'",I3.3)') k
        ier = nf90_inq_varid(pncid(ipf),newname,idv)
        if ( ier/=nf90_noerr .and. k<100 ) then
          write(newname,'("'//trim(name)//'",I2.2)') k
          ier = nf90_inq_varid(pncid(ipf),newname,idv)          
        end if
        if ( ier/=nf90_noerr .and. k<10 ) then
          write(newname,'("'//trim(name)//'",I1.1)') k
          ier = nf90_inq_varid(pncid(ipf),newname,idv)          
        end if
        if ( ier/=nf90_noerr ) then
          exit
        end if
        ! obtain scaling factors and offsets from attributes
        ier = nf90_get_att(pncid(ipf),idv,'add_offset',laddoff)
        if ( ier/=nf90_noerr ) laddoff = 0.
        ier = nf90_get_att(pncid(ipf),idv,'scale_factor',lsf)
        if ( ier/=nf90_noerr ) lsf = 1.
        ier = nf90_get_var(pncid(ipf),idv,rvar(:,k,1),start=start(1:ndims),count=ncount(1:ndims))
        call ncmsg(name,ier)
        ! unpack data
        rvar(:,k,1) = rvar(:,k,1)*real(lsf,8) + real(laddoff,8)
      end do
    end if ! ier==nf_noerr

    if ( qtest ) then
      ! e.g., restart file or nogather=.true.
      ca = pipan*pjpan*pnpan*ipf
      var(1+ca:pipan*pjpan*pnpan+ca,1:kk,1:ll) = rvar(:,:,:)
    else
      ! e.g., mesonest file
      if ( myid==0 ) then
        allocate( gvar(pipan*pjpan*pnpan,kk,ll,fnresid) )
        call ccmpi_gatherxr8(gvar,rvar,0,comm_ip)
        do jpf = 1,fnresid
          ip = ipf*fnresid + jpf - 1   ! local file number
          do l = 1,ll
            do k = 1,kk
              do n = 0,pnpan-1
                no = n - pnoff(ip) + 1   ! global panel number of local file
                ca = pioff(ip,no) + pjoff(ip,no)*pil_g + no*pil_g**2 - pil_g
                cc = n*pipan*pjpan - pipan
                do j = 1,pjpan
                  var(1+j*pil_g+ca:pipan+j*pil_g+ca,k,l) = gvar(1+j*pipan+cc:pipan+j*pipan+cc,k,l,jpf)
                end do
              end do  
            end do
          end do
        end do
        deallocate( gvar )
      else
        allocate( gvar(1,1,1,1) ) 
        call ccmpi_gatherxr8(gvar,rvar,0,comm_ip)
        deallocate( gvar )
      end if
    end if ! qtest

  end do ! ipf
  
  deallocate( rvar )
  
end if ! mynproc>0

return
end subroutine hr5pr8
#endif

!--------------------------------------------------------------
! This subroutine opens parallel input files
subroutine histopen(ncid,ifile,ier,fileerror)
      
use cc_mpi
use newmpar_m
use parm_m
      
integer, parameter :: nihead = 54
      
integer, dimension(0:5) :: duma, dumb
integer, dimension(12) :: idum
integer, intent(out) :: ncid, ier
integer is, ipf, dmode
integer ipin, ipin_f, ipin_new, nxpr, nypr
integer ltst, der, myrank, resprocmode
integer nxp_test
integer, dimension(:,:), allocatable, save :: dum_off
integer, dimension(:,:), allocatable, save :: resprocdata_inv
integer, dimension(:), allocatable, save :: resprocmap_inv
integer, dimension(:), allocatable, save :: procfileowner
integer(kind=4), dimension(1) :: start, ncount
integer(kind=4), dimension(nihead) :: lahead
integer(kind=4) lncid, lidum, ldid, lvid, llen
logical, intent(in), optional :: fileerror
logical ferror
character(len=*), intent(in) :: ifile
character(len=170) pfile
character(len=8) fdecomp

resprocmode = 0
ferror = .false.
if ( present(fileerror) ) then
  ferror = fileerror
end if

if ( myid==0 ) then
  ! attempt to open single file with myid==0
  ier = nf90_open(ifile,nf90_nowrite,lncid)
  ncid = lncid
  fnproc = 1              ! number of files to be read over all processors
  dmode = 0               ! Single file (dmode=0), Face decomposition (dmode=1),
                          ! Depreciated (dmode=2) or Uniform decomposition (dmode=3)
  pipan = 0               ! Number of X grid points within a file panel
  pjpan = 0               ! Number of Y grid points within a file panel
  pnpan = 0               ! Number of panels in file
  ptest = .false.         ! Files match current processor (e.g., Restart file), allowing MPI gather/scatter to be avoided
  pfall = .false.         ! Every processor has been assigned at least one file, no need to Bcast metadata data
  resprocformat = .false. ! procformat file format with multiple processes per file
      
  ! attempt to open parallel files
  if ( ier/=nf90_noerr ) then
    write(pfile,"(a,'.',i6.6)") trim(ifile), 0
    ier = nf90_open(pfile,nf90_nowrite,lncid)
    ncid = lncid
    if ( ier/=nf90_noerr ) then
      if ( ferror ) then
        write(6,*) "ERROR: Cannot open ",trim(pfile)," or ",trim(ifile)
        call ccmpi_abort(-1)
      else  
        write(6,*) "WARN: Cannot open ",trim(pfile)," or ",trim(ifile)
      end if  
    else
      ! found parallel input file  
      der = nf90_get_att(lncid,nf90_global,"nproc",lidum)
      fnproc = lidum
      call ncmsg("nproc",der)
      der = nf90_get_att(lncid,nf90_global,"procmode",lidum)
      if ( der==nf90_noerr ) then
        write(6,*) "Found procformat input file ",trim(ifile)
        resprocformat = .true.
        resprocmode = lidum
      else
        write(6,*) "Found parallel input file ",trim(ifile)
        resprocformat = .false.
        resprocmode = 0
      end if
      fdecomp = ''      
      der = nf90_get_att(lncid,nf90_global,"decomp",fdecomp)
      call ncmsg("decomp",der)
      select case(fdecomp)
        case('face')
          dmode = 1
        case('uniform')  ! old uniform (Depreciated)
          dmode = 2
        case('uniform1') ! new uniform (Dix style)
          dmode = 3
        case default
          write(6,*) "ERROR: Unknown decomposition ",trim(fdecomp)
          call ccmpi_abort(-1)
      end select
    end if
  else
    ! nproc should only exist in multi-file input
    der = nf90_get_att(lncid,nf90_global,"nproc",lidum)
    if ( der==nf90_noerr ) then
      write(6,*) "ERROR: Incorrect base filename"
      write(6,*) "Possibly incorrectly included 000000 in filename"
      call ccmpi_abort(-1)
    end if
    write(6,*) "Found single input file ",trim(ifile)
  end if

  ! Read grid metadata
  if ( ier==nf90_noerr ) then
    ! Newer global attributes method
    der = nf90_get_att(lncid,nf90_global,"il_g",lidum)
    if ( der==nf90_noerr ) then
      pil_g = lidum
      der = nf90_get_att(lncid,nf90_global,"jl_g",lidum)
      pjl_g = lidum
    else
      ! Older int_header method
      der = nf90_get_att(lncid,nf90_global,"int_header",lahead)
      pil_g = lahead(1)
      pjl_g = lahead(2)
      call ncmsg("int_header",der)
    end if  
    ! Atmosphere vertical levels
    der = nf90_inq_dimid(lncid,"lev",ldid)
    if ( der==nf90_noerr ) then
      der = nf90_inquire_dimension(lncid,ldid,len=llen)
      pka_g = llen
      call ncmsg("lev",der)
    else
      pka_g = 0  
    end if
    ! Ocean vertical levels if present
    der = nf90_inq_dimid(lncid,"olev",ldid)
    if ( der==nf90_noerr ) then
      der = nf90_inquire_dimension(lncid,ldid,len=llen)
      pko_g = llen
      call ncmsg("olev",der)
    else
      pko_g = 0
    end if
        
    ! special case for single file input
    if ( dmode==0 ) then
      nxp_test = max( min( pil_g/pil_single, int(sqrt(real(nproc)/6.)) ), 1 )
      do while ( mod(pil_g,nxp_test)/=0 .and. nxp_test>1 )
        nxp_test = nxp_test - 1
      end do  
      fnproc = nxp_test**2*6
      if ( myid==0 ) then
        write(6,*) "-> Decompose single file into sections = ",fnproc
      end if      
    end if  
    
    if ( allocated(pioff) ) then
      write(6,*) "ERROR: Cannot open new input file until old file is closed"
      call ccmpi_abort(-1)
    end if
    allocate( pioff(0:fnproc-1,0:5), pjoff(0:fnproc-1,0:5) )
    allocate( pnoff(0:fnproc-1) )
        
    select case(dmode)
      case(0) ! single file - decompose into multiple sections
        pnpan = max(1,6/fnproc)
        do ipin = 0,fnproc-1
          call face_set(pipan,pjpan,pnoff(ipin),duma,dumb,pnpan,pil_g,ipin,fnproc,nxpr,nypr)
          pioff(ipin,:) = duma(:)
          pjoff(ipin,:) = dumb(:)
        end do
      case(1) ! face decomposition
        pnpan = max(1,6/fnproc)
        do ipin = 0,fnproc-1
          call face_set(pipan,pjpan,pnoff(ipin),duma,dumb,pnpan,pil_g,ipin,fnproc,nxpr,nypr)
          pioff(ipin,:) = duma(:)
          pjoff(ipin,:) = dumb(:)
        end do
      case(2) ! old uniform decomposition - depreciated
        pnpan = 6
        do ipin = 0,fnproc-1
          call uniform_set(pipan,pjpan,pnoff(ipin),duma,dumb,pnpan,pil_g,ipin,fnproc,nxpr,nypr)
          pioff(ipin,:) = duma(:)
          pjoff(ipin,:) = dumb(:)
        end do
      case(3) ! new uniform decomposition - depreciated
        pnpan = 6
        do ipin = 0,fnproc-1
          call dix_set(pipan,pjpan,pnoff(ipin),duma,dumb,pnpan,pil_g,ipin,fnproc,nxpr,nypr)
          pioff(ipin,:) = duma(:)
          pjoff(ipin,:) = dumb(:)
        end do
    end select

    ptest = .false.
    if ( nproc==fnproc ) then
      if ( pil_g==il_g .and. pjl_g==jl_g ) then
        if ( dmode==1 ) then
          ptest = .true.
        end if
      end if
    end if

    write(6,*) "-> dmode,ptest,resprocformat ",dmode,ptest,resprocformat
    write(6,*) "-> fnproc,pil_g,pjl_g        ",fnproc,pil_g,pjl_g
    write(6,*) "-> pipan,pjpan,pnpan         ",pipan,pjpan,pnpan
    
  end if

  idum(1) = fnproc
  if ( resprocformat ) then
    idum(2) = 1
  else
    idum(2) = 0
  end if
  idum(3) = pipan
  idum(4) = pjpan
  idum(5) = pnpan
  if (ptest) then
    idum(6) = 1
  else
    idum(6) = 0      
  end if
  idum(7) = ier
  idum(8) = pka_g
  idum(9) = pko_g
  idum(10) = pil_g
  idum(11) = pjl_g
  idum(12) = dmode
  
  if ( resprocformat ) then
    start(1) = 1
    ncount(1) = fnproc
    allocate( resprocdata_inv(0:fnproc-1,2) )
    der = nf90_inq_varid(lncid,'gprocnode',lvid)  
    if ( der==nf90_noerr ) then
      ! procformat v2 format  
      write(6,*) "-> Found procformat v2"  
      der = nf90_get_var(lncid,lvid,resprocdata_inv(:,1),start,ncount)
      der = nf90_inq_varid(lncid,'gprocoffset',lvid)    
      der = nf90_get_var(lncid,lvid,resprocdata_inv(:,2),start,ncount)
    else
      ! procformat v1 format  
      write(6,*) "-> Found procformat v1"  
      allocate( resprocmap_inv(0:fnproc-1) )  
      der = nf90_inq_varid(lncid,'gprocessor',lvid)
      if ( der/=nf90_noerr ) then
        write(6,*) "ERROR: Corrupted procformat file"
        call ccmpi_abort(-1)
      end if    
      der = nf90_get_var(lncid,lvid,resprocmap_inv,start,ncount)
      do ipin = 0,fnproc-1
        ipin_new = resprocmap_inv(ipin)  
        resprocdata_inv(ipin,1) = ipin_new/resprocmode   
        resprocdata_inv(ipin,2) = mod(ipin_new, resprocmode)  
      end do    
      deallocate( resprocmap_inv )
    end if
  end if
  
  write(6,*) "-> Broadcasting file metadata"
end if

! Broadcast file metadata
call ccmpi_bcast(idum(1:12),0,comm_world)
fnproc        = idum(1)      ! number of files to be read
resprocformat = idum(2)==1   ! test for procformat file format
pipan         = idum(3)      ! width of panel in each file
pjpan         = idum(4)      ! length of panel in each file
pnpan         = idum(5)      ! number of panels in each file
ptest         = idum(6)==1   ! test for match between files and processes
ier           = idum(7)      ! file error flag
pka_g         = idum(8)      ! number of atmosphere levels
pko_g         = idum(9)      ! number of ocean levels
pil_g         = idum(10)     ! global grid size
pjl_g         = idum(11)     ! global grid size
dmode         = idum(12)     ! file decomposition

if ( ier/=nf90_noerr ) return

if ( myid==0 ) then
  write(6,*) "-> Opening data files"
end if

! calculate number of files to be read on this processor
fnresid = min( fnproc, nproc ) 
do while ( mod(fnproc,fnresid)/=0 )
  fnresid = fnresid - 1     ! limit on processor ranks that will read files    
end do
fncount = fnproc/fnresid
if ( myid<fnresid ) then
  mynproc = fncount  ! calculate the number of files to be read per process
else
  mynproc = 0
end if

! allocate array of file handles.  Unallocated implies no files to be read on this processor
if ( allocated(pncid) ) then
  write(6,*) "ERROR: Cannot open new input file until old file is closed"
  call ccmpi_abort(-1)
end if
if ( mynproc>0 ) then
  allocate( pncid(0:mynproc-1), pfown(0:mynproc-1) )
  allocate( ppanid(0:mynproc-1), ppiid(0:mynproc-1), ppjid(0:mynproc-1) )
  pfown(:) = .false.
  ppanid(:) = 0
  ppiid(:) = 1
  ppjid(:) = 1
end if
! Rank 0 can start with the second file, because the first file has already been opened
if ( myid==0 ) then 
  is = 1
  pncid(0) = ncid
else
  is = 0
end if


! distribute comms
if ( myid==0 ) then
  write(6,*) "-> Splitting comms for distributing file data with fnresid ",fnresid
  write(6,*) "-> Number of files to be read with mynproc ",mynproc
end if

! define comm group to read the residual files
if ( myid<fnresid ) then
  ltst = 0
  myrank = myid
else
  ltst = -1 ! undefined
  myrank = myid - fnresid
end if
call ccmpi_commsplit(comm_ip,comm_world,ltst,myrank)


! Open files
if ( mynproc>0 ) then
  if ( resprocformat ) then

    ! procformat ----------------------------------------------------
      
    allocate( pprid(0:mynproc-1) )
      
    ! copy process map
    allocate( procfileowner(0:fnproc-1) )
    procfileowner(:) = -1
    if ( .not.allocated(resprocdata_inv) ) then
      allocate( resprocdata_inv(0:fnproc-1,2) )
    end if
    call ccmpi_bcast(resprocdata_inv,0,comm_ip)
 
    ! update required process and load files
    do ipf = 0,mynproc-1
      ipin = ipf*fnresid + myid                ! parallel file number
      pprid(ipf) = resprocdata_inv(ipin,2) + 1 ! procformat process
    end do
    
    if ( myid==0 ) then
      procfileowner(0) = 0
      pfown(0) = .true.
    end if
    
    do ipf = is,mynproc-1
      ipin = ipf*fnresid + myid               ! parallel file number
      ipin_f = resprocdata_inv(ipin,1)        ! procformat file
      if ( procfileowner(ipin_f)==-1 ) then
        procfileowner(ipin_f) = ipf ! which ipf is responsible for opening this file
        pfown(ipf) = .true.         ! which ipf is responsible for closing this file
        write(pfile,"(a,'.',i6.6)") trim(ifile), ipin_f
        der = nf90_open(pfile,nf90_nowrite,pncid(ipf))
        if ( der/=nf90_noerr ) then
          write(6,*) "ERROR: Cannot open ",pfile
          call ncmsg("open",der)
        end if
      else
        pncid(ipf) = pncid(procfileowner(ipin_f)) ! file is already open
      end if
    end do
    deallocate( procfileowner )
    deallocate( resprocdata_inv )    
   
  else if ( dmode/=0 ) then
    
    ! original parallel file without resprocformat
    ! loop through files to be opened by this processor
    do ipf = is,mynproc-1
      ipin = ipf*fnresid + myid
      pfown(ipf) = .true.
      write(pfile,"(a,'.',i6.6)") trim(ifile), ipin
      der = nf90_open(pfile,nf90_nowrite,pncid(ipf))
      if ( der/=nf90_noerr ) then
        write(6,*) "ERROR: Cannot open ",pfile
        call ncmsg("open",der)
      end if
    end do
   
  else  
    
    ! single input file decomposed into multiple sections
    do ipf = is,mynproc-1
      pfown(ipf) = .true.
      der = nf90_open(ifile,nf90_nowrite,pncid(ipf))
      if ( der/=nf90_noerr ) then
        write(6,*) "ERROR: Cannot open ",ifile
        call ncmsg("open",der)
      end if
    end do 
      
  end if ! resprocformat ..else..
end if   ! mynproc>0



pfall = fnresid==nproc  ! are all processes associated with a file?
                        ! this means we do not need to Bcst file metadata
if ( mynproc>0 ) then
  ncid = pncid(0)       ! set ncid to the first file handle as onthefly
                        ! assumes changes in ncid reflect a new file
                        ! and hence updates the metadata
end if


allocate( dum_off(0:fnproc-1,1:13) )
if ( myid==0 ) then
  write(6,*) "-> Broadcast file coordinate data"
  dum_off(0:fnproc-1,1:6)  = pioff(0:fnproc-1,0:5)
  dum_off(0:fnproc-1,7:12) = pjoff(0:fnproc-1,0:5)
  dum_off(0:fnproc-1,13)   = pnoff(0:fnproc-1)
else
  allocate( pioff(0:fnproc-1,0:5), pjoff(0:fnproc-1,0:5) )
  allocate( pnoff(0:fnproc-1) )
end if
call ccmpi_bcast(dum_off,0,comm_world)
pioff(0:fnproc-1,0:5) = dum_off(0:fnproc-1,1:6)
pjoff(0:fnproc-1,0:5) = dum_off(0:fnproc-1,7:12)
pnoff(0:fnproc-1)     = dum_off(0:fnproc-1,13)
deallocate( dum_off )


! single input file decomposed into multiple sections
if ( dmode==0 ) then
  if ( mynproc>0 ) then
    if ( pnpan/=1 ) then
      write(6,*) "ERROR: Decomposing single input file requires pnpan=1"
      call ccmpi_abort(-1)
    end if
    do ipf = 0,mynproc-1
      ipin = ipf*fnresid + myid 
      ppanid(ipf) = 1 - pnoff(ipin)
      ppiid(ipf) = pioff(ipin,0) + 1 ! assumes face decomposition
      ppjid(ipf) = pjoff(ipin,0) + 1 ! assumes face decomposition
    end do 
  end if
end if

! final checks
if ( .not.resprocformat ) then
  do ipf = 0,mynproc-1  
    if ( pnpan>1 .and. ppjid(ipf)>1 ) then
      write(6,*) "ERROR: segment of input file requested over multiple panels"
      call ccmpi_abort(-1)
    end if    
  end do    
end if

if ( myid==0 ) then
  write(6,*) "-> Ready to read data from input file"
end if

return
end subroutine histopen
      
!--------------------------------------------------------------
! This subroutine closes parallel input files
subroutine histclose

use cc_mpi
      
implicit none

integer ipf, plen
integer(kind=4) ierr

if ( myid==0 ) then
  write(6,*) 'Closing input file'
end if

if ( allocated(pncid) ) then
  plen = size(pncid)
  do ipf = 0,plen-1
    if ( pfown(ipf) ) then
      ierr = nf90_close(pncid(ipf))
    end if
  end do
  deallocate( pfown, pncid )
  deallocate( ppanid, ppiid, ppjid )
end if

if ( allocated(pprid) ) then
  deallocate(pprid)
end if
if (allocated(pioff)) then
  deallocate(pioff,pjoff,pnoff)
end if
if ( allocated(filemap_req) ) then
  deallocate( filemap_req, filemap_qmod )
end if
if ( allocated(filemap_recv) ) then
  deallocate( filemap_recv, filemap_rmod )
end if
if ( allocated(filemap_send) ) then
  deallocate( filemap_send, filemap_smod )
end if
if ( allocated(filemap_indx) ) then
  deallocate( filemap_indx )
  filemap_indxlen = 0
end if

call ccmpi_filewinfinalize

ncidold = -1 ! flag onthefly to load metadata

return
end subroutine histclose

!--------------------------------------------------------------------
! This subroutine advances input date by the amount of time defined by mtimer_r
subroutine datefix(kdate_r,ktime_r,mtimer_r,allleap,silent)

use cc_mpi
use dates_m
use parm_m

implicit none

integer, intent(inout) :: kdate_r,ktime_r
integer(kind=8), intent(inout) :: mtimer_r
integer, intent(in), optional :: allleap
integer(kind=8), dimension(12) :: mdays
integer, dimension(12) :: mdays4
integer leap_l
integer(kind=8) iyr,imo,iday,ihr,imins
integer(kind=8) mtimerh,mtimerm
integer(kind=8) mdays_save
integer(kind=8), parameter :: minsday = 1440
logical, intent(in), optional :: silent
logical quiet

if ( present(allleap) ) then
  leap_l = allleap
else
  leap_l = leap
end if

if ( present(silent) ) then
  quiet = silent
else
  quiet = .false.
end if

iyr=int(kdate_r,8)/10000_8
imo=(int(kdate_r,8)-10000_8*iyr)/100_8
iday=int(kdate_r,8)-10000_8*iyr-100_8*imo
ihr=int(ktime_r,8)/100_8
imins=int(ktime_r,8)-100_8*ihr
if ( myid==0 .and. .not.quiet ) then
  write(6,*) 'entering datefix'
  write(6,'(A,I6,I4,I4)') ' -> iyr,imo,iday:       ',iyr,imo,iday
  write(6,'(A,I4,I4,I8)') ' -> ihr,imins,mtimer_r: ',ihr,imins,int(mtimer_r)
end if

call calendar_function(mdays4,kdate_r,leap_l)
mdays(:) = int( mdays4(:) )
do while ( mtimer_r>minsday*mdays(imo) )
  mtimer_r=mtimer_r-minsday*mdays(imo)
  imo=imo+1_8
  if ( imo>12_8 ) then
    imo=1_8
    iyr=iyr+1_8
    if ( leap_l==1 ) then
      mdays(2)=28_8      
      if ( mod(iyr,4_8)==0   ) mdays(2)=29_8
      if ( mod(iyr,100_8)==0 ) mdays(2)=28_8
      if ( mod(iyr,400_8)==0 ) mdays(2)=29_8
    end if
  end if
end do
if ( diag .and. .not.quiet ) then
  write(6,*)'b datefix iyr,imo,iday,ihr,imins,mtimer_r: ', &
                       iyr,imo,iday,ihr,imins,int(mtimer_r)
end if
  
iday=iday+mtimer_r/minsday
mtimer_r=mod(mtimer_r,minsday)
if ( diag .and. .not.quiet ) then
  write(6,*)'c datefix iyr,imo,iday,ihr,imins,mtimer_r: ', &
                       iyr,imo,iday,ihr,imins,int(mtimer_r)
end if
  
! at this point mtimer_r has been reduced to fraction of a day
mtimerh=mtimer_r/60_8
mtimerm=mtimer_r-mtimerh*60_8  ! minutes left over
ihr=ihr+mtimerh
imins=imins+mtimerm

ihr=ihr+imins/60_8
imins=mod(imins,60_8)
if ( diag .and. .not.quiet ) then
  write(6,*)'d datefix iyr,imo,iday,ihr,imins,mtimer_r: ', &
                       iyr,imo,iday,ihr,imins,int(mtimer_r)
end if
  
iday=iday+ihr/24_8
ihr=mod(ihr,24_8)
if ( diag .and. .not.quiet ) then
  write(6,*)'e datefix iyr,imo,iday,ihr,imins,mtimer_r: ', &
                       iyr,imo,iday,ihr,imins,int(mtimer_r)
end if
  
mdays_save=mdays(imo)
imo=imo+(iday-1_8)/mdays(imo)
iday=mod(iday-1_8,mdays_save)+1_8

iyr=iyr+(imo-1_8)/12_8
imo=mod(imo-1_8,12_8)+1_8

kdate_r=int(iday+100_8*(imo+100_8*iyr),4)
ktime_r=int(ihr*100_8+imins,4)
mtimer_r = 0.
if ( diag .and. .not.quiet ) then
  write(6,*)'end datefix iyr,imo,iday,ihr,imins,mtimer_r: ', &
                         iyr,imo,iday,ihr,imins,int(mtimer_r)
end if
  
if ( myid==0 .and. .not.quiet ) then
  write(6,*)'leaving datefix kdate_r,ktime_r: ',kdate_r,ktime_r
end if
  
return
end subroutine datefix

!--------------------------------------------------------------------
! Same as datefix, but for time with units of months
subroutine datefix_month(kdate_r,mtimer_r)

implicit none

integer, intent(inout) :: kdate_r
integer(kind=8), intent(inout) :: mtimer_r
integer(kind=8) iyr,imo,iday

iyr=int(kdate_r,8)/10000_8
imo=(int(kdate_r,8)-10000_8*iyr)/100_8
iday=15_8

do while ( mtimer_r>0_8 )
  mtimer_r = mtimer_r-1_8
  imo = imo + 1_8
  if ( imo>12_8 ) then
    imo = 1_8
    iyr = iyr + 1_8
  end if
end do
  
kdate_r = int(iday + 100_8*imo + 10000_8*iyr)
mtimer_r = 0_8
  
return
end subroutine datefix_month

!--------------------------------------------------------------------
! This subroutine converts time units into a date
subroutine processdatestring(datestring,kdate_rsav,ktime_rsav)

use cc_mpi            ! CC MPI routines

implicit none

integer, intent(out) :: kdate_rsav, ktime_rsav
integer iposa, iposb, ierx
integer yyyy, mm, dd, hh, mt
character(len=*), intent(in) :: datestring

! process year
iposa = index(trim(datestring),'since')
iposa = iposa + 5 ! skip 'since'
iposb = index(trim(datestring(iposa:)),'-')
iposb = iposa + iposb - 2 ! remove '-'
read(datestring(iposa:iposb),FMT=*,iostat=ierx) yyyy
if ( ierx/=0 ) then
  write(6,*) "ERROR reading time units.  Expecting year but found ",datestring(iposa:iposb)
  call ccmpi_abort(-1)
end if

! process month
iposa = iposb + 2 ! skip '-'
iposb = index(trim(datestring(iposa:)),'-')
iposb = iposa + iposb - 2 ! remove '-'
read(datestring(iposa:iposb),FMT=*,iostat=ierx) mm
if ( ierx/=0 ) then
  write(6,*) "ERROR reading time units.  Expecting month but found ",datestring(iposa:iposb)
  call ccmpi_abort(-1)
end if

! process day
iposa = iposb + 2 ! skip '-'
iposb = index(trim(datestring(iposa:)),' ')
iposb = iposa + iposb - 2 ! remove ' '
if ( iposb<iposa ) then
  read(datestring(iposa:),FMT=*,iostat=ierx) dd
  if ( ierx/=0 ) then
    write(6,*) "ERROR reading time units.  Expecting day but found ",datestring(iposa:)
    call ccmpi_abort(-1)
  end if
  kdate_rsav = yyyy*10000 + mm*100 + dd
  ktime_rsav = 0
  return
end if
read(datestring(iposa:iposb),FMT=*,iostat=ierx) dd
if ( ierx/=0 ) then
  write(6,*) "ERROR reading time units.  Expecting day but found ",datestring(iposa:iposb)
  call ccmpi_abort(-1)
end if

! process hour
iposa = iposb + 2 ! skip ' '
iposb = index(trim(datestring(iposa:)),':')
iposb = iposa + iposb - 2 ! remove ':'
read(datestring(iposa:iposb),FMT=*,iostat=ierx) hh
if ( ierx/=0 ) then
  write(6,*) "ERROR reading time units.  Expecting hour but found ",datestring(iposa:iposb)
  call ccmpi_abort(-1)
end if

! process mins
iposa = iposb + 2 ! skip ':'
iposb = index(trim(datestring(iposa:)),':')
iposb = iposa + iposb - 2 ! remove ':'
if ( iposb<iposa ) then
  read(datestring(iposa:),FMT=*,iostat=ierx) mt
  if ( ierx/=0 ) then
    write(6,*) "ERROR reading time units.  Expecting minutes but found ",datestring(iposa:)
    call ccmpi_abort(-1)
  end if
  kdate_rsav = yyyy*10000 + mm*100 + dd
  ktime_rsav = hh*100 + mt
  return
end if
read(datestring(iposa:iposb),FMT=*,iostat=ierx) mt
if ( ierx/=0 ) then
  write(6,*) "ERROR reading time units.  Expecting minutes but found ",datestring(iposa:iposb)
  call ccmpi_abort(-1)
end if

! final date and time
kdate_rsav = yyyy*10000 + mm*100 + dd
ktime_rsav = hh*100 + mt

return
end subroutine processdatestring

!--------------------------------------------------------------------
! Set up number of minutes from beginning of year
subroutine getzinp(jyear,jmonth,jday,jhour,jmin,mins,allleap)

use cc_mpi
use dates_m
use parm_m

implicit none

integer, intent(out) :: jyear,jmonth,jday,jhour,jmin ! start date of run
integer, intent(out) :: mins                         ! elapsed time from start of year
integer mstart, n1
integer, dimension(12) :: ndoy
logical, intent(in), optional :: allleap ! force use of leap days even if leap=0
logical lleap

if ( present(allleap) ) then
  lleap = allleap
else
  lleap = .false.    
end if

jyear  = kdate/10000
jmonth = (kdate-jyear*10000)/100
jday   = kdate - jyear*10000 - jmonth*100
jhour  = ktime/100
jmin   = ktime - jhour*100
      
if ( jmonth<1 .or. jmonth>12 ) then
  write(6,*) "ERROR: Invalid month ",jmonth," for kdate ",kdate
  call ccmpi_abort(-1)
end if 

if ( lleap .or. leap==1 ) then ! 365/366 day calendar
  ndoy(:) = (/ 0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334 /)
  n1 = 0
  if ( mod(jyear,4)  ==0 ) n1 = 1
  if ( mod(jyear,100)==0 ) n1 = 0
  if ( mod(jyear,400)==0 ) n1 = 1
  ndoy(3:12)=ndoy(3:12)+n1
else if ( leap==0 ) then ! 365 day calendar
  ndoy(:) = (/ 0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334 /)  
else if ( leap==2 ) then ! 360 day calendar  
  ndoy(:) = (/ 0, 30, 60, 90, 120, 150, 180, 210, 240, 270, 300, 330 /)
else
  write(6,*) "ERROR: Unknown option for leap = ",leap
  call ccmpi_abort(-1)
end if  

mstart = 1440*(ndoy(jmonth)+jday-1) + 60*jhour + jmin ! mins from start of year

! mtimer contains number of minutes since the start of the run.
mins = mtimer + mstart

return
end subroutine getzinp

!--------------------------------------------------------------------
! DEFINE ATTRIBUTES
subroutine attrib(cdfid,dim,ndim,name,lname,units,xmin,xmax,time_freq,time_method,area_method,itype)

use cc_mpi
use newmpar_m
use parm_m

implicit none

integer, intent(in) :: cdfid, itype, ndim
integer, intent(in) :: time_freq, time_method, area_method
integer, dimension(ndim), intent(in) :: dim
integer ier
integer(kind=4) vtype, idv, lcdfid, lsize, lcompression
integer(kind=4) ldlen, ltlen, tlen, lil, ljl, lvnode
integer(kind=4), dimension(ndim) :: ldim
integer(kind=4), dimension(ndim) :: chunks
real, intent(in) :: xmin, xmax
real(kind=4) lscalef, laddoff
character(len=*), intent(in) :: name
character(len=*), intent(in) :: lname
character(len=*), intent(in) :: units
character(len=80) :: cell_methods, cell_methods_area, cell_methods_time

if ( itype==short_m ) then
  vtype = nf90_short
else if ( itype==double_m ) then
  vtype = nf90_double  
else
#ifdef i8r8
  vtype = nf90_double
#else
  vtype = nf90_float
#endif
end if
lcdfid = cdfid
ldim   = dim
if ( localhist .and. ndim>3 ) then
  ! MJT notes - PR identified (/il, jl, kl, vnode_nproc, min(10, tlen)/) as optimal.
  ! However, here we simplify the code and PR reports that the performance is
  ! similar
  lil = int(il,kind=4)
  ljl = int(jl,kind=4)
  lvnode = int( min(vnode_nproc,node_nx) ,kind=4)
  select case(ndim)
    case(6)
      ier = nf90_inquire_dimension(lcdfid,ldim(3),len=ldlen)
      chunks = (/ lil, ljl, ldlen, 1_4, lvnode, 1_4 /)
    case(5)
      ier = nf90_inquire_dimension(lcdfid,ldim(3),len=ldlen)
      chunks = (/ lil, ljl, ldlen, lvnode, 1_4 /)
    case(4)
      ier = nf90_inquire_dimension(lcdfid,ldim(4),len=tlen)
      ltlen = int( chunk_time, kind=4)
      do while ( mod(tlen,ltlen)/=0 .and. ltlen>1 )
        ltlen = ltlen - 1
      end do
      chunks = (/ lil, ljl, lvnode, ltlen /)
    case default
      write(6,*) "ERROR: Invalid ndim in attrib ",ndim
      call ccmpi_abort(-1)
  end select
  ! MJT notes - compression can degrade model runtime when using a subset of processes
  ! for file IO.  However, file sizes can increase significantly without compression.
  lcompression = compression   
  ier = nf90_def_var(lcdfid, name, vtype, ldim, idv)
  call ncmsg("def_var - "//trim(name),ier)
  ier = nf90_def_var_deflate(lcdfid, idv, 1_4, 1_4, lcompression)
  ier = nf90_def_var_chunking(lcdfid, idv, NF90_CHUNKED, chunks )
else
  lcompression = compression
  ier = nf90_def_var(lcdfid, name, vtype, ldim, idv)
  call ncmsg("def_var - "//trim(name),ier)
  ier = nf90_def_var_deflate(lcdfid, idv, 1_4, 1_4, lcompression)
end if
lsize = len_trim(lname)
ier = nf90_put_att(lcdfid,idv,'long_name',lname)
call ncmsg("long_name",ier)
lsize = len_trim(units)
if ( lsize>0 ) then
  ier = nf90_put_att(lcdfid,idv,'units',units)
  call ncmsg("units",ier)
endif
if ( vtype==nf90_short ) then
  ier = nf90_put_att(lcdfid,idv,'valid_min',minv)
  call ncmsg("valid_min",ier)
  ier = nf90_put_att(lcdfid,idv,'valid_max',maxv)
  call ncmsg("valid_max",ier)
  ier = nf90_put_att(lcdfid,idv,'missing_value',missval)
  call ncmsg("missing_value",ier)
  lscalef = (xmax-xmin)/(real(maxv)-real(minv))
  laddoff = xmin-real(lscalef)*real(minv)
  ier = nf90_put_att(lcdfid,idv,'add_offset',laddoff)
  call ncmsg("add_offset",ier)
  ier = nf90_put_att(lcdfid,idv,'scale_factor',lscalef)
  call ncmsg("scale_factor",ier)
else
  ier = nf90_put_att(lcdfid,idv,'missing_value',nf90_fill_float)
  call ncmsg("missing_value",ier)
endif
ier = nf90_put_att(lcdfid,idv,'FORTRAN_format','G11.4')
call ncmsg("FORTRAN_format",ier)
if ( time_freq==daily_m ) then
  ier = nf90_put_att(lcdfid,idv,'valid_time','daily')
  call ncmsg("valid_time",ier)
else if ( time_freq==sixhr_m ) then
  ier = nf90_put_att(lcdfid,idv,'valid_time','6hr')
  call ncmsg("valid_time",ier)
endif

cell_methods_area = ''
select case(area_method)
  case(amean_m)
    cell_methods_area = 'area: mean'
  case(land_m)
    cell_methods_area = 'area: mean where land'
  case(sea_m)
    cell_methods_area = 'area: mean where sea'
  case(seaice_m)
    cell_methods_area = 'area: mean where sea_ice'
 end select

cell_methods_time = ''
select case(time_method)
  case(point_m)
    cell_methods_time = 'time: point'
  case(tmean_m)  
    cell_methods_time = 'time: mean'
  case(max_m)    
    cell_methods_time = 'time: maximum'
  case(min_m)  
    cell_methods_time = 'time: minimum'
  case(fixed_m)    
    cell_methods_time = 'time: fixed'
  case(sum_m)  
    cell_methods_time = 'time: sum'
end select    

cell_methods = ''
if ( len_trim(cell_methods_area)>0 .and. len_trim(cell_methods_time)>0 ) then
  cell_methods = trim(cell_methods_area) // ' ' // trim(cell_methods_time)
else if ( len_trim(cell_methods_area)>0 ) then
  cell_methods = trim(cell_methods_area)
else if ( len_trim(cell_methods_time)>0 ) then
  cell_methods = trim(cell_methods_time)
end if
if ( len_trim(cell_methods)>0 ) then
  ier = nf90_put_att(lcdfid,idv,'cell_methods',cell_methods)
  call ncmsg("cell_methods",ier) 
end if    
      
return
end subroutine attrib

!--------------------------------------------------------------------
! 3D NETCDF WRITE ARRAY ROUTINES
subroutine histwrt3r4(var,sname,idnc,iarch,local,lwrite)

use cc_mpi              ! CC MPI routines
use newmpar_m           ! Grid parameters
use parm_m              ! Model configuration

implicit none

integer, intent(in) :: idnc, iarch
real, dimension(ifull), intent(in) :: var
real, dimension(ifull,1,1) :: wvar
real, dimension(1,1,1,1) :: var_t
character(len=*), intent(in) :: sname
logical, intent(in) :: local, lwrite

if (.not.lwrite) then
  wvar(1:ifull,1,1) = nf90_fill_float
else
  wvar(1:ifull,1,1) = var(1:ifull)
end if

if ( local ) then
  call hw5lp(wvar,sname,idnc,iarch)  
else if ( localhist ) then
  call ccmpi_gatherx(var_t,wvar,0,comm_vnode)
else if ( myid==0 ) then
  call hw5a(wvar,sname,idnc,iarch)  
else
  call ccmpi_gather(wvar)
end if

if ( myid==0 .and. nmaxpr==1 ) then
  if ( any(abs(wvar-real(nf90_fill_float))<1.e-20) ) then
    write(6,'(" histwrt3 ",a20,i8,a7)') sname,iarch,"missing"
  else
    write(6,'(" histwrt3 ",a20,i8)') sname,iarch
  end if
end if

return
end subroutine histwrt3r4

#ifndef i8r8
!--------------------------------------------------------------------
! 3D NETCDF WRITE ARRAY ROUTINES (double precision version)
subroutine histwrt3r8(var,sname,idnc,iarch,local,lwrite)

use cc_mpi              ! CC MPI routines
use newmpar_m           ! Grid parameters
use parm_m              ! Model configuration

implicit none

integer, intent(in) :: idnc, iarch
real(kind=8), dimension(ifull), intent(in) :: var
real(kind=8), dimension(ifull,1,1) :: wvar
real(kind=8), dimension(1,1,1,1) :: var_t
character(len=*), intent(in) :: sname
logical, intent(in) :: local, lwrite

if (.not.lwrite) then
  wvar(1:ifull,1,1) = real(nf90_fill_float,8)
else
  wvar(1:ifull,1,1) = var(1:ifull)
end if

if ( local ) then
  call hw5lpr8(wvar,sname,idnc,iarch)  
else if ( localhist ) then
  call ccmpi_gatherxr8(var_t,wvar,0,comm_vnode)
else if ( myid==0 ) then
  call hw5ar8(wvar,sname,idnc,iarch)
else
  call ccmpi_gatherr8(wvar)
end if

if ( myid==0 .and. nmaxpr==1 ) then
  if ( any(abs(wvar-real(nf90_fill_float,8))<1.e-20_8) ) then
    write(6,'(" histwrt3r8 ",a20,i8,a7)') sname,iarch,"missing"
  else
    write(6,'(" histwrt3r8 ",a20,i8)') sname,iarch
  end if
end if

return
end subroutine histwrt3r8
#endif

!--------------------------------------------------------------------
! 4D NETCDF WRITE ARRAY ROUTINES
subroutine histwrt4r4(var,sname,idnc,iarch,local,lwrite)

use cc_mpi              ! CC MPI routines
use newmpar_m           ! Grid parameters
use parm_m              ! Model configuration

implicit none

integer, intent(in) :: idnc, iarch
integer kk
real, dimension(:,:), intent(in) :: var
real, dimension(ifull,size(var,2),1) :: wvar
real, dimension(1,1,1,1) :: var_g
character(len=*), intent(in) :: sname
logical, intent(in) :: local, lwrite

kk = size(var,2)

if ( .not.lwrite ) then
  wvar(1:ifull,1:kk,1) = real(nf90_fill_float)
else
  wvar(1:ifull,1:kk,1) = var(1:ifull,1:kk)
endif

if ( local ) then
  call hw5lp(wvar,sname,idnc,iarch)  
else if ( localhist ) then
  call ccmpi_gatherx(var_g,wvar,0,comm_vnode)
else if ( myid==0  ) then
  call hw5a(wvar,sname,idnc,iarch)  
else
  call ccmpi_gather(wvar)  
endif

if ( myid==0 .and. nmaxpr==1 ) then
  if ( any(abs(wvar-real(nf90_fill_float))<1.e-20) ) then
    write(6,'(" histwrt4 ",a20,i4,a7)') sname,iarch,"missing"
  else
    write(6,'(" histwrt4 ",a20,i4)') sname,iarch
  end if
end if

return
end subroutine histwrt4r4

#ifndef i8r8
subroutine histwrt4r8(var,sname,idnc,iarch,local,lwrite)

use cc_mpi              ! CC MPI routines
use newmpar_m           ! Grid parameters
use parm_m              ! Model configuration

implicit none

integer, intent(in) :: idnc, iarch
integer kk
real(kind=8), dimension(:,:), intent(in) :: var
real(kind=8), dimension(ifull,size(var,2),1) :: wvar
real(kind=8), dimension(1,1,1,1) :: var_g
character(len=*), intent(in) :: sname
logical, intent(in) :: local, lwrite

kk = size(var,2)

if ( .not.lwrite ) then
  wvar(1:ifull,1:kk,1)=real(nf90_fill_float)
else
  wvar(1:ifull,1:kk,1)=var(1:ifull,1:kk)
endif

if ( local ) then
  call hw5lpr8(wvar,sname,idnc,iarch)  
else if ( localhist ) then
  call ccmpi_gatherxr8(var_g,wvar,0,comm_vnode)
else if ( myid==0 ) then
  call hw5ar8(wvar,sname,idnc,iarch)
else
  call ccmpi_gatherr8(wvar)
endif

if ( myid==0 .and. nmaxpr==1 ) then
  if ( any(abs(wvar-real(nf90_fill_float,8))<1.e-20_8) ) then
    write(6,'(" histwrt4r8 ",a20,i4,a7)') sname,iarch,"missing"
  else
    write(6,'(" histwrt4r8 ",a20,i4)') sname,iarch
  end if
end if

return
end subroutine histwrt4r8
#endif

!--------------------------------------------------------------------
! 5D NETCDF WRITE ARRAY ROUTINES
subroutine histwrt5r4(var,sname,idnc,iarch,local,lwrite)

use cc_mpi              ! CC MPI routines
use newmpar_m           ! Grid parameters
use parm_m              ! Model configuration

implicit none

integer, intent(in) :: idnc, iarch
integer ll, kk
real, dimension(:,:,:), intent(in) :: var
real, dimension(ifull,size(var,2),size(var,3)) :: wvar
real, dimension(1,1,1,1) :: var_g
character(len=*), intent(in) :: sname
logical, intent(in) :: local, lwrite

kk = size(var,2)
ll = size(var,3)

if ( .not.lwrite ) then
  wvar(1:ifull,1:kk,1:ll) = real(nf90_fill_float)
else
  wvar(1:ifull,1:kk,1:ll) = var(1:ifull,1:kk,1:ll)
endif

if ( local ) then
  call hw5lp(wvar,sname,idnc,iarch)  
else if ( localhist ) then
  call ccmpi_gatherx(var_g,wvar,0,comm_vnode)
else if ( myid==0 ) then
  call hw5a(wvar,sname,idnc,iarch)
else
  call ccmpi_gather(wvar)
endif

if ( myid==0 .and. nmaxpr==1 ) then
  if ( any(abs(wvar-real(nf90_fill_float))<1.e-20) ) then
    write(6,'(" histwrt5 ",a20,i4,a7)') sname,iarch,"missing"
  else
    write(6,'(" histwrt5 ",a20,i4)') sname,iarch
  end if
end if

return
end subroutine histwrt5r4

#ifndef i8r8
subroutine histwrt5r8(var,sname,idnc,iarch,local,lwrite)

use cc_mpi              ! CC MPI routines
use newmpar_m           ! Grid parameters
use parm_m              ! Model configuration

implicit none

integer, intent(in) :: idnc, iarch
integer kk, ll
real(kind=8), dimension(:,:,:), intent(in) :: var
real(kind=8), dimension(ifull,size(var,2),size(var,3)) :: wvar
real(kind=8), dimension(1,1,1,1) :: var_g
character(len=*), intent(in) :: sname
logical, intent(in) :: local, lwrite

kk = size(var,2)
ll = size(var,3)

if ( .not.lwrite ) then
  wvar(1:ifull,1:kk,1:ll) = real(nf90_fill_float)
else
  wvar(1:ifull,1:kk,1:ll) = var(1:ifull,1:kk,1:ll)
endif

if ( local ) then
  call hw5lpr8(wvar,sname,idnc,iarch)  
else if ( localhist ) then
  call ccmpi_gatherxr8(var_g,wvar,0,comm_vnode)
else if ( myid==0 ) then
  call hw5ar8(wvar,sname,idnc,iarch)
else
  call ccmpi_gatherr8(wvar)
endif

if ( myid==0 .and. nmaxpr==1 ) then
  if ( any(abs(wvar-real(nf90_fill_float,8))<1.e-20_8) ) then
    write(6,'(" histwrt5r8 ",a20,i4,a7)') sname,iarch,"missing"
  else
    write(6,'(" histwrt5r8 ",a20,i4)') sname,iarch
  end if
end if

return
end subroutine histwrt5r8
#endif

! procformat and local(write)
subroutine hw5lp(var,sname,idnc,iarch)

use cc_mpi               ! CC MPI routines
use newmpar_m            ! Grid parameters
use parm_m               ! Model configuration

implicit none

integer, intent(in) :: idnc, iarch
integer iq, k, ier, v, kk, l, ll
integer(kind=4) mid, vtype, lidnc, ndims
integer(kind=4), dimension(6) :: start, ncount
integer(kind=2), dimension(:,:,:,:), allocatable :: ipack_g
real, dimension(:,:,:), intent(in) :: var
real, dimension(:,:,:,:), allocatable :: var_g
real(kind=4) laddoff, lscale_f
character(len=*), intent(in) :: sname
#ifndef i8r8  
real(kind=8), dimension(:,:,:,:), allocatable :: dpack_g
#endif


kk = size(var,2)
ll = size(var,3)

allocate( var_g(ifull,kk,ll,vnode_nproc) )
call ccmpi_gatherx(var_g,var,0,comm_vnode)

lidnc = idnc
ier = nf90_inq_varid(lidnc,sname,mid)
call ncmsg(sname,ier)
ier = nf90_inquire_variable(lidnc,mid,xtype=vtype,ndims=ndims)
select case(ndims)
  case(3)
    start(1:ndims) = (/ 1, 1, 1 /)
    ncount(1:ndims) = (/ il, jl, vnode_nproc /)
  case(4)
    start(1:ndims) = (/ 1, 1, 1, iarch /)
    ncount(1:ndims) = (/ il, jl, vnode_nproc, 1 /)
  case(5)
    start(1:ndims) = (/ 1, 1, 1, 1, iarch /)
    ncount(1:ndims) = (/ il, jl, kk, vnode_nproc, 1 /)
  case(6)
    start(1:ndims) = (/ 1, 1, 1, 1, 1, iarch /)
    ncount(1:ndims) = (/ il, jl, kk, ll, vnode_nproc, 1 /)
  case default
    write(6,*) "ERROR: Variable ",trim(sname)," expected to have 6 or less dimensions,"
    write(6,*) "but was created with ndims = ",ndims
    call ccmpi_abort(-1)
end select

if ( vtype==nf90_short ) then
  allocate( ipack_g(ifull,kk,ll,vnode_nproc) )  
  if ( all(var_g>9.8e36) ) then
    ipack_g(:,:,:,:) = missval
  else
    ier = nf90_get_att(lidnc,mid,'add_offset',laddoff)
    ier = nf90_get_att(lidnc,mid,'scale_factor',lscale_f)
    do v = 1,vnode_nproc
      do l = 1,ll  
        do k = 1,kk
          do iq = 1,ifull
            ipack_g(iq,k,l,v) = nint(max(min((var_g(iq,k,l,v)-real(laddoff))/real(lscale_f),real(maxv)),real(minv)),2)
          end do  
        end do
      end do
    end do
  end if
  ier = nf90_put_var(lidnc,mid,ipack_g,start=start(1:ndims),count=ncount(1:ndims))
  deallocate( ipack_g )
#ifndef i8r8  
else if ( vtype==nf90_double ) then
  allocate( dpack_g(ifull,kk,ll,vnode_nproc) )
  dpack_g = real(var_g,8)
  ier = nf90_put_var(lidnc,mid,dpack_g,start=start(1:ndims),count=ncount(1:ndims))
  deallocate( dpack_g )
#endif  
else
  ier = nf90_put_var(lidnc,mid,var_g,start=start(1:ndims),count=ncount(1:ndims))
end if
call ncmsg(sname,ier)

deallocate( var_g )

return
end subroutine hw5lp           

subroutine hw5a(var,sname,idnc,iarch)

use cc_mpi               ! CC MPI routines
use newmpar_m            ! Grid parameters
use parm_m               ! Model configuration

implicit none

integer, intent(in) :: idnc, iarch
integer iq, k, ier, kk, l, ll
integer(kind=4) mid, vtype, lidnc, ndims
integer(kind=4), dimension(5) :: start, ncount
integer(kind=2), dimension(:,:,:), allocatable :: ipack_g
real, dimension(:,:,:), intent(in) :: var
real, dimension(:,:,:), allocatable :: var_g
real(kind=4) laddoff, lscale_f
character(len=*), intent(in) :: sname
#ifndef i8r8
real(kind=8), dimension(:,:,:), allocatable :: dpack_g
#endif

kk = size(var,2)
ll = size(var,3)

allocate( var_g(ifull_g,kk,ll) )
call ccmpi_gather(var,var_g)

lidnc = idnc
ier = nf90_inq_varid(lidnc,sname,mid)
call ncmsg(sname,ier)
ier = nf90_inquire_variable(lidnc,mid,xtype=vtype,ndims=ndims)
select case(ndims)
  case(2)
    start(1:ndims) = (/ 1, 1 /)
    ncount(1:ndims) = (/ il_g, jl_g /)
  case(3)
    start(1:ndims) = (/ 1, 1, iarch /)
    ncount(1:ndims) = (/ il_g, jl_g, 1 /)
  case(4)
    start(1:ndims) = (/ 1, 1, 1, iarch /)
    ncount(1:ndims) = (/ il_g, jl_g, kk, 1 /)
  case(5)
    start(1:ndims) = (/ 1, 1, 1, 1, iarch /)
    ncount(1:ndims) = (/ il_g, jl_g, kk, ll, 1 /)
  case default
    write(6,*) "ERROR: Variable ",trim(sname)," expected to have 5 or less dimensions,"
    write(6,*) "but was created with ndims = ",ndims
    call ccmpi_abort(-1)
end select

if ( vtype==nf90_short ) then
  allocate( ipack_g(ifull_g,kk,ll) )
  if ( all(var_g>9.8e36) ) then
    ipack_g(:,:,:) = missval
  else
    ier = nf90_get_att(lidnc,mid,'add_offset',laddoff)
    ier = nf90_get_att(lidnc,mid,'scale_factor',lscale_f)
    do l = 1,ll  
      do k = 1,kk
        do iq = 1,ifull_g
          ipack_g(iq,k,l) = nint(max(min((var_g(iq,k,l)-real(laddoff))/real(lscale_f),real(maxv)),real(minv)),2)
        end do
      end do
    end do
  end if
  ier = nf90_put_var(lidnc,mid,ipack_g,start=start(1:ndims),count=ncount(1:ndims))
  deallocate( ipack_g )
#ifndef i8r8  
else if ( vtype==nf90_double ) then
  allocate( dpack_g(ifull,kk,ll) )
  dpack_g = real(var_g,8)
  ier = nf90_put_var(lidnc,mid,dpack_g,start=start(1:ndims),count=ncount(1:ndims))
  deallocate( dpack_g )
#endif  
else
  ier = nf90_put_var(lidnc,mid,var_g,start=start(1:ndims),count=ncount(1:ndims))
end if
call ncmsg(sname,ier)

deallocate( var_g )

return
end subroutine hw5a

#ifndef i8r8
! procformat and local(write)
subroutine hw5lpr8(var,sname,idnc,iarch)

use cc_mpi               ! CC MPI routines
use newmpar_m            ! Grid parameters
use parm_m               ! Model configuration

implicit none

integer, intent(in) :: idnc, iarch
integer :: iq, k, ier, v, kk, l, ll
integer(kind=4) mid, vtype, lidnc, ndims
integer(kind=4), dimension(6) :: start, ncount
integer(kind=2), dimension(:,:,:,:), allocatable :: ipack_g
real(kind=8), dimension(:,:,:), intent(in) :: var
real(kind=8), dimension(:,:,:,:), allocatable :: var_g
real(kind=4) :: laddoff, lscale_f
real(kind=4), dimension(:,:,:,:), allocatable :: fpack_g
character(len=*), intent(in) :: sname


kk = size(var,2)
ll = size(var,3)

allocate( var_g(ifull,kk,ll,vnode_nproc) )
call ccmpi_gatherxr8(var_g,var,0,comm_vnode)

lidnc = idnc
ier = nf90_inq_varid(lidnc,sname,mid)
call ncmsg(sname,ier)
ier = nf90_inquire_variable(lidnc,mid,xtype=vtype,ndims=ndims)
select case(ndims)
  case(3)
    start(1:ndims) = (/ 1, 1, 1 /)
    ncount(1:ndims) = (/ il, jl, vnode_nproc /)
  case(4)
    start(1:ndims) = (/ 1, 1, 1, iarch /)
    ncount(1:ndims) = (/ il, jl, vnode_nproc, 1 /)
  case(5)
    start(1:ndims) = (/ 1, 1, 1, 1, iarch /)
    ncount(1:ndims) = (/ il, jl, kk, vnode_nproc, 1 /)
  case(6)
    start(1:ndims) = (/ 1, 1, 1, 1, 1, iarch /)
    ncount(1:ndims) = (/ il, jl, kk, ll, vnode_nproc, 1 /)
  case default
    write(6,*) "ERROR: Variable ",trim(sname)," expected to have 6 or less dimensions,"
    write(6,*) "but was created with ndims = ",ndims
    call ccmpi_abort(-1)
end select

if ( vtype==nf90_short ) then
  allocate( ipack_g(ifull,kk,ll,vnode_nproc) )  
  if ( all(var>9.8e36_8) ) then
    ipack_g(:,:,:,:) = missval
  else
    ier = nf90_get_att(lidnc,mid,'add_offset',laddoff)
    ier = nf90_get_att(lidnc,mid,'scale_factor',lscale_f)
    do v = 1,vnode_nproc
      do l = 1,ll  
        do k = 1,kk
          do iq = 1,ifull
            ipack_g(iq,k,l,v) = nint(max(min((var_g(iq,k,l,v)-real(laddoff,8))/real(lscale_f,8),real(maxv,8)), &
                                     real(minv,8)),2)
          end do  
        end do
      end do
    end do
  end if
  ier = nf90_put_var(lidnc,mid,ipack_g,start=start(1:ndims),count=ncount(1:ndims))
  deallocate( ipack_g )
else if ( vtype==nf90_float ) then
  allocate( fpack_g(ifull,kk,ll,vnode_nproc) )
  fpack_g = real(var_g,4)
  ier = nf90_put_var(lidnc,mid,fpack_g,start=start(1:ndims),count=ncount(1:ndims))
  deallocate( fpack_g )
else
  ier = nf90_put_var(lidnc,mid,var_g,start=start(1:ndims),count=ncount(1:ndims))
end if
call ncmsg(sname,ier)

deallocate( var_g )

return
end subroutine hw5lpr8    

subroutine hw5ar8(var,sname,idnc,iarch)

use cc_mpi               ! CC MPI routines
use newmpar_m            ! Grid parameters
use parm_m               ! Model configuration

implicit none

integer, intent(in) :: idnc, iarch
integer :: iq, k, ier, kk, l, ll
integer(kind=4) mid, vtype, lidnc, ndims
integer(kind=4), dimension(5) :: start, ncount
integer(kind=2), dimension(:,:,:), allocatable :: ipack_g
real(kind=8), dimension(:,:,:), intent(in) :: var
real(kind=8), dimension(:,:,:), allocatable :: var_g
real(kind=4) :: laddoff, lscale_f
real(kind=4), dimension(:,:,:), allocatable :: fpack_g
character(len=*), intent(in) :: sname

kk = size(var,2)
ll = size(var,3)

allocate( var_g(ifull_g,kk,ll) )
call ccmpi_gatherr8(var,var_g)

lidnc = idnc
ier = nf90_inq_varid(lidnc,sname,mid)
call ncmsg(sname,ier)
ier = nf90_inquire_variable(lidnc,mid,xtype=vtype,ndims=ndims)
select case(ndims)
  case(2)
    start(1:ndims) = (/ 1, 1 /)
    ncount(1:ndims) = (/ il_g, jl_g /)
  case(3)
    start(1:ndims) = (/ 1, 1, iarch /)
    ncount(1:ndims) = (/ il_g, jl_g, 1 /)
  case(4)
    start(1:ndims) = (/ 1, 1, 1, iarch /)
    ncount(1:ndims) = (/ il_g, jl_g, kk, 1 /)
  case(5)
    start(1:ndims) = (/ 1, 1, 1, 1, iarch /)
    ncount(1:ndims) = (/ il_g, jl_g, kk, ll, 1 /)
  case default
    write(6,*) "ERROR: Variable ",trim(sname)," expected to have 5 or less dimensions,"
    write(6,*) "but was created with ndims = ",ndims
    call ccmpi_abort(-1)
end select

if ( vtype==nf90_short ) then
  allocate( ipack_g(ifull_g,kk,ll) )
  if ( all(var>9.8e36_8) ) then
    ipack_g(:,:,:) = missval
  else
    ier = nf90_get_att(lidnc,mid,'add_offset',laddoff)
    ier = nf90_get_att(lidnc,mid,'scale_factor',lscale_f)
    do l = 1,ll  
      do k = 1,kk
        do iq = 1,ifull_g
          ipack_g(iq,k,l) = nint(max(min((var_g(iq,k,l)-real(laddoff,8))/real(lscale_f,8),real(maxv,8)), &
                                   real(minv,8)),2)
        end do
      end do
    end do
  end if
  ier = nf90_put_var(lidnc,mid,ipack_g,start=start(1:ndims),count=ncount(1:ndims))
  deallocate( ipack_g )
else if ( vtype==nf90_float ) then
  allocate( fpack_g(ifull,kk,ll) )
  fpack_g = real(var_g,4)
  ier = nf90_put_var(lidnc,mid,fpack_g,start=start(1:ndims),count=ncount(1:ndims))
  deallocate( fpack_g )
else
  ier = nf90_put_var(lidnc,mid,var_g,start=start(1:ndims),count=ncount(1:ndims))
end if
call ncmsg(sname,ier)

deallocate( var_g )

return
end subroutine hw5ar8
#endif

!--------------------------------------------------------------------

subroutine ccnf_open(fname,ncid,status)

integer, intent(out) :: ncid
integer, intent(out), optional :: status
integer ncstatus
integer(kind=4) lncid
character(len=*), intent(in) :: fname

lncid = 0
ncstatus = nf90_open(fname,nf90_nowrite,lncid)
ncid = lncid

if ( present(status) ) then
  status = ncstatus
else
  call ncmsg("open",ncstatus)
end if

return
end subroutine ccnf_open

subroutine ccnf_create(fname,ncid)

integer, intent(out) :: ncid
integer ncstatus
integer(kind=4) :: lncid
character(len=*), intent(in) :: fname

ncstatus = nf90_create(fname,nf90_netcdf4,lncid)
ncid = lncid
if ( ncstatus/=nf90_noerr ) then
  write(6,*) "ERROR: Cannot create fname = ",trim(fname)
end if
call ncmsg("create",ncstatus)

return
end subroutine ccnf_create

subroutine ccnf_close(ncid)

integer, intent(in) :: ncid
integer ncstatus
integer(kind=4) lncid

lncid = ncid
ncstatus = nf90_close(lncid)
call ncmsg("close",ncstatus)

return
end subroutine ccnf_close

subroutine ccnf_nofill(ncid)

integer, intent(in) :: ncid
integer ncstatus
integer(kind=4) lncid, lomode

lncid = ncid
ncstatus = nf90_set_fill(lncid,nf90_nofill,lomode)
call ncmsg("nofill",ncstatus)

return
end subroutine ccnf_nofill

subroutine ccnf_inq_dimid(ncid,dname,did,tst)

integer, intent(in) :: ncid
integer, intent(out) :: did
integer ncstatus
integer(kind=4) lncid, ldid
logical ltst
logical, intent(out), optional :: tst
character(len=*), intent(in) :: dname

lncid = ncid
ncstatus = nf90_inq_dimid(lncid,dname,ldid)
ltst = (ncstatus/=nf90_noerr)
did = ldid

if ( present(tst) ) then
  tst = ltst
else
  call ncmsg("dimid",ncstatus)
end if

return
end subroutine ccnf_inq_dimid

subroutine ccnf_inq_dimlen(ncid,dname,dlen,failok)

use cc_mpi

integer, intent(in) :: ncid
integer, intent(inout) :: dlen
integer(kind=4) lncid, lncstatus, ldid, ldlen
character(len=*), intent(in) :: dname
logical, intent(in), optional :: failok
logical ftest

if ( present(failok) ) then
  ftest = failok
else
  ftest = .false.
end if

lncid = ncid
ldlen = dlen
lncstatus = nf90_inq_dimid(lncid,dname,ldid)
if ( lncstatus/=nf90_noerr ) then
  if ( ftest ) return
  write(6,*) nf90_strerror(lncstatus)
  call ccmpi_abort(-1)
end if
lncstatus = nf90_inquire_dimension(lncid,ldid,len=ldlen)
if ( lncstatus/=nf90_noerr ) then
  if ( ftest ) return
  write(6,*) nf90_strerror(lncstatus)
  call ccmpi_abort(-1)
end if
dlen = ldlen

return
end subroutine ccnf_inq_dimlen

subroutine ccnf_inq_varid(ncid,vname,vid,tst)

integer, intent(in) :: ncid
integer, intent(out) :: vid
integer ncstatus
integer(kind=4) lncid, lvid
character(len=*), intent(in) :: vname
logical, intent(out), optional :: tst
logical ltst

lncid = ncid
ncstatus = nf90_inq_varid(lncid,vname,lvid)
ltst = (ncstatus/=nf90_noerr)
vid = lvid

if (present(tst)) then
  tst = ltst
else
  call ncmsg(vname,ncstatus)
end if

return
end subroutine ccnf_inq_varid

subroutine ccnf_inq_exist(ncid,vname,tst)

integer, intent(in) :: ncid
integer ncstatus
integer(kind=4) lncid, lvid
character(len=*), intent(in) :: vname
logical, intent(out) :: tst

lncid = ncid
ncstatus = nf90_inq_varid(lncid,vname,lvid)
! true if varid is found
tst = (ncstatus==nf90_noerr)

return
end subroutine ccnf_inq_exist

subroutine ccnf_inq_varndims(ncid,vid,ndims)

integer, intent(in) :: ncid, vid
integer, intent(out) :: ndims
integer ncstatus
integer(kind=4) lncid, lvid, lndims

lncid = ncid
lvid = vid
ncstatus = nf90_inquire_variable(lncid,lvid,ndims=lndims)
call ncmsg('ccnf_inq_varndims',ncstatus)
ndims = lndims

return
end subroutine ccnf_inq_varndims

subroutine ccnf_def_dim(ncid,dname,nsize,did)

integer, intent(in) :: ncid, nsize
integer, intent(out) :: did
integer ncstatus
integer(kind=4) lncid, lnsize, ldid
character(len=*), intent(in) :: dname

lncid = ncid
lnsize = nsize
ncstatus = nf90_def_dim(lncid,dname,lnsize,ldid)
did = ldid
call ncmsg("def_dim",ncstatus)

return
end subroutine ccnf_def_dim

subroutine ccnf_def_dimu(ncid,dname,did)

integer, intent(in) :: ncid
integer, intent(out) :: did
integer ncstatus
integer(kind=4) lncid, ldid
character(len=*), intent(in) :: dname

lncid = ncid
ncstatus = nf90_def_dim(lncid,dname,nf90_unlimited,ldid)
did = ldid
call ncmsg("def_dimu",ncstatus)

return
end subroutine ccnf_def_dimu

subroutine ccnf_def_var_v(ncid,vname,vtype,vndim,dims,vid)

use cc_mpi

integer, intent(in) :: ncid, vndim
integer, intent(out) :: vid
integer, dimension(vndim), intent(in) :: dims
integer ncstatus
integer(kind=4) lncid, ltype, lvid
integer(kind=4), dimension(vndim) :: ldims
character(len=*), intent(in) :: vname
character(len=*), intent(in) :: vtype

select case(vtype)
  case('int')
    ltype = nf90_int
  case('float')
    ltype = nf90_float
  case('double')
    ltype = nf90_double
  case('char')
    ltype = nf90_char
  case default
    write(6,*) "ERROR: Unknown option for ccnf_def_var ",vtype
    call ccmpi_abort(-1)
end select

lncid = ncid
ldims = dims
ncstatus = nf90_def_var(lncid,vname,ltype,ldims,lvid,deflate_level=1_4)
vid = lvid
call ncmsg("def_var - "//trim(vname),ncstatus)

return
end subroutine ccnf_def_var_v

subroutine ccnf_def_var_s(ncid,vname,vtype,vid)

use cc_mpi

integer, intent(in) :: ncid
integer, intent(out) :: vid
integer ncstatus
integer(kind=4) lncid, lvid, ltype
character(len=*), intent(in) :: vname
character(len=*), intent(in) :: vtype

select case(vtype)
  case('int')
    ltype = nf90_int
  case('float')
    ltype = nf90_float
  case('double')
    ltype = nf90_double
  case('char')
    ltype = nf90_char
  case default
    write(6,*) "ERROR: Unknown option for ccnf_def_var ",vtype
    call ccmpi_abort(-1)
end select

lncid = ncid
ncstatus = nf90_def_var(lncid,vname,ltype,lvid)
vid = lvid
call ncmsg("def_var - "//trim(vname),ncstatus)

return
end subroutine ccnf_def_var_s

subroutine ccnf_enddef(ncid)

integer, intent(in) :: ncid
integer ncstatus
integer(kind=4) lncid

lncid = ncid
ncstatus = nf90_enddef(lncid)
call ncmsg("enddef",ncstatus)

return
end subroutine ccnf_enddef

subroutine ccnf_redef(ncid)

integer, intent(in) :: ncid
integer ncstatus
integer(kind=4) lncid

lncid = ncid
ncstatus = nf90_redef(lncid)
call ncmsg("redef",ncstatus)

return
end subroutine ccnf_redef

subroutine ccnf_sync(ncid)

integer, intent(in) :: ncid
integer ncstatus
integer(kind=4) lncid

lncid = ncid
ncstatus = nf90_sync(lncid)
call ncmsg("sync",ncstatus)

return
end subroutine ccnf_sync

subroutine ccnf_get_var_real2r_t(ncid,vname,vdat)

integer, intent(in) :: ncid
integer ncstatus
integer(kind=4) :: lvid, lncid
real, dimension(:), intent(out) :: vdat
character(len=*), intent(in) :: vname

lncid = ncid
ncstatus = nf90_inq_varid(lncid,vname,lvid)
call ncmsg("get_var_varid",ncstatus)
ncstatus = nf90_get_var(lncid,lvid,vdat)
call ncmsg("get_var",ncstatus)

return
end subroutine ccnf_get_var_real2r_t

subroutine ccnf_get_vara_text1t_t(ncid,name,start,ncount,vdat)

integer, intent(in) :: ncid
integer, dimension(:), intent(in) :: start, ncount
integer ncstatus
integer(kind=4) lncid, lvid
integer(kind=4), dimension(size(start)) :: lstart
integer(kind=4), dimension(size(ncount)) :: lncount
character(len=*), intent(out) :: vdat
character(len=*), intent(in) :: name

lncid = ncid
ncstatus = nf90_inq_varid(lncid,name,lvid)
call ncmsg(name,ncstatus)
lstart(:) = start(:)
lncount(:) = ncount(:)
ncstatus = nf90_get_var(lncid,lvid,vdat,start=lstart,count=lncount)
call ncmsg("get_vara_text_t",ncstatus)

return
end subroutine ccnf_get_vara_text1t_t 

subroutine ccnf_get_var_int2i_t(ncid,vname,vdat)

integer, intent(in) :: ncid
integer ncstatus
integer, dimension(:), intent(out) :: vdat
integer(kind=4) :: lncid, lvid
character(len=*), intent(in) :: vname

lncid = ncid
ncstatus = nf90_inq_varid(lncid,vname,lvid)
call ncmsg("get_var_varid",ncstatus)
ncstatus = nf90_get_var(lncid,lvid,vdat)
call ncmsg("get_var",ncstatus)

return
end subroutine ccnf_get_var_int2i_t

subroutine ccnf_get_vara_real1r_s(ncid,vid,start,vdat)

integer, intent(in) :: ncid, vid, start
integer ncstatus
integer(kind=4) lncid, lvid
integer(kind=4), dimension(1) :: lstart
integer(kind=4), dimension(1) :: lncount
real, intent(out) :: vdat
real, dimension(1) :: lvdat

lncid = ncid
lvid = vid
lstart = start
lncount = 1
ncstatus = nf90_get_var(lncid,lvid,lvdat,start=lstart,count=lncount)
call ncmsg("get_vara_real1r",ncstatus)
vdat = lvdat(1)

return
end subroutine ccnf_get_vara_real1r_s

subroutine ccnf_get_vara_real1r_t(ncid,name,start,vdat)

integer, intent(in) :: ncid, start
integer ncstatus
integer(kind=4) lncid, lvid
integer(kind=4), dimension(1) :: lstart
integer(kind=4), dimension(1) :: lncount
real, intent(out) :: vdat
real, dimension(1) :: lvdat
character(len=*), intent(in) :: name

lncid = ncid
ncstatus = nf90_inq_varid(lncid,name,lvid)
call ncmsg(name,ncstatus)
lstart = start
lncount = 1
ncstatus = nf90_get_var(lncid,lvid,lvdat,start=lstart,count=lncount)
call ncmsg("get_vara_real1r",ncstatus)
vdat = lvdat(1)

return
end subroutine ccnf_get_vara_real1r_t

subroutine ccnf_get_vara_real2r_0(ncid,vid,start,ncount,vdat)

integer, intent(in) :: ncid, vid
integer, intent(in) :: start, ncount
integer ncstatus
integer(kind=4) lncid, lvid
integer(kind=4), dimension(1) :: lstart
integer(kind=4), dimension(1) :: lncount
real, dimension(:), intent(out) :: vdat

lncid = ncid
lvid = vid
lstart = start
lncount = ncount
ncstatus = nf90_get_var(lncid,lvid,vdat,start=lstart,count=lncount)
call ncmsg("get_vara_real2r",ncstatus)

return
end subroutine ccnf_get_vara_real2r_0 

subroutine ccnf_get_vara_real2r_t(ncid,name,start,ncount,vdat)

integer, intent(in) :: ncid
integer, dimension(:), intent(in) :: start, ncount
integer ncstatus
integer(kind=4) lncid, lvid
integer(kind=4), dimension(size(start)) :: lstart
integer(kind=4), dimension(size(ncount)) :: lncount
real, dimension(:), intent(out) :: vdat
character(len=*), intent(in) :: name

lncid = ncid
ncstatus = nf90_inq_varid(lncid,name,lvid)
call ncmsg(name,ncstatus)
lstart(:) = start(:)
lncount(:) = ncount(:)
ncstatus = nf90_get_var(lncid,lvid,vdat,start=lstart,count=lncount)
call ncmsg("get_vara_real2r_t",ncstatus)

return
end subroutine ccnf_get_vara_real2r_t 

subroutine ccnf_get_vara_real2r_s(ncid,vid,start,ncount,vdat)

integer, intent(in) :: ncid, vid
integer, dimension(:), intent(in) :: start, ncount
integer ncstatus
integer(kind=4) lncid, lvid
integer(kind=4), dimension(size(start)) :: lstart
integer(kind=4), dimension(size(ncount)) :: lncount
real, dimension(:), intent(out) :: vdat

lncid = ncid
lvid = vid
lstart = start
lncount = ncount
ncstatus = nf90_get_var(lncid,lvid,vdat,start=lstart,count=lncount)
call ncmsg("get_vara_real2r",ncstatus)

return
end subroutine ccnf_get_vara_real2r_s

subroutine ccnf_get_vara_real3r_s(ncid,vid,start,ncount,vdat)

integer, intent(in) :: ncid, vid
integer, dimension(:) :: start, ncount
integer ncstatus
integer(kind=4) lncid, lvid
integer(kind=4), dimension(size(start)) :: lstart
integer(kind=4), dimension(size(ncount)) :: lncount
real, dimension(:,:), intent(out) :: vdat

lncid = ncid
lvid = vid
lstart = start
lncount = ncount
ncstatus = nf90_get_var(lncid,lvid,vdat,start=lstart,count=lncount)
call ncmsg("get_vara_real3r",ncstatus)

return
end subroutine ccnf_get_vara_real3r_s

subroutine ccnf_get_vara_real3r_t(ncid,name,start,ncount,vdat)

integer, intent(in) :: ncid
integer, dimension(:), intent(in) :: start, ncount
integer ncstatus
integer(kind=4) lncid, lvid
integer(kind=4), dimension(size(start)) :: lstart
integer(kind=4), dimension(size(ncount)) :: lncount
real, dimension(:,:), intent(out) :: vdat
character(len=*), intent(in) :: name

lncid = ncid
ncstatus = nf90_inq_varid(lncid,name,lvid)
call ncmsg(name,ncstatus)
lstart(:) = start(:)
lncount(:) = ncount(:)
ncstatus = nf90_get_var(lncid,lvid,vdat,start=lstart,count=lncount)
call ncmsg("get_vara_real3r_t",ncstatus)

return
end subroutine ccnf_get_vara_real3r_t 

subroutine ccnf_get_vara_real4r_s(ncid,vid,start,ncount,vdat)

integer, intent(in) :: ncid, vid
integer, dimension(:) :: start, ncount
integer ncstatus
integer(kind=4) :: lncid, lvid
integer(kind=4), dimension(size(start)) :: lstart
integer(kind=4), dimension(size(ncount)) :: lncount
real, dimension(:,:,:), intent(out) :: vdat

lncid = ncid
lvid = vid
lstart = start
lncount = ncount
ncstatus = nf90_get_var(lncid,lvid,vdat,start=lstart,count=lncount)
call ncmsg("get_vara_real4r",ncstatus)

return
end subroutine ccnf_get_vara_real4r_s

subroutine ccnf_get_vara_int1i_s(ncid,vid,start,vdat)

integer, intent(in) :: ncid, vid, start
integer ncstatus
integer, intent(out) :: vdat
integer(kind=4) :: lncid, lvid
integer(kind=4), dimension(1) :: lstart
integer(kind=4), dimension(1) :: lncount
integer(kind=4), dimension(1) :: lvdat

lncid = ncid
lvid = vid
lstart = start
lncount = 1
ncstatus = nf90_get_var(lncid,lvid,lvdat,start=lstart,count=lncount)
call ncmsg("get_vara_int1i",ncstatus)
vdat = lvdat(1)

return
end subroutine ccnf_get_vara_int1i_s

subroutine ccnf_get_vara_int2i_t(ncid,name,start,ncount,vdat)

integer, intent(in) :: ncid
integer, dimension(:) :: start, ncount
integer ncstatus
integer, dimension(:), intent(out) :: vdat
integer(kind=4) :: lncid, lvid
integer(kind=4), dimension(size(start)) :: lstart
integer(kind=4), dimension(size(ncount)) :: lncount
character(len=*), intent(in) :: name

lncid = ncid
ncstatus = nf90_inq_varid(lncid,name,lvid)
lstart(:) = start(:)
lncount(:) = ncount(:)
ncstatus = nf90_get_var(lncid,lvid,vdat,start=lstart,count=lncount)
call ncmsg("get_vara_int2i_t",ncstatus)

return
end subroutine ccnf_get_vara_int2i_t

subroutine ccnf_get_vara_int2i_s(ncid,vid,start,ncount,vdat)

integer, intent(in) :: ncid, vid
integer, dimension(:) :: start, ncount
integer ncstatus
integer, dimension(:), intent(out) :: vdat
integer(kind=4) :: lncid, lvid
integer(kind=4), dimension(size(start)) :: lstart
integer(kind=4), dimension(size(ncount)) :: lncount

lncid = ncid
lvid = vid
lstart = start
lncount = ncount
ncstatus = nf90_get_var(lncid,lvid,vdat,start=lstart,count=lncount)
call ncmsg("get_vara_int2i",ncstatus)

return
end subroutine ccnf_get_vara_int2i_s

subroutine ccnf_get_vara_int3i_t(ncid,name,start,ncount,vdat)

integer, intent(in) :: ncid
integer, dimension(:) :: start, ncount
integer ncstatus
integer, dimension(:,:), intent(out) :: vdat
integer(kind=4) :: lncid, lvid
integer(kind=4), dimension(size(start)) :: lstart
integer(kind=4), dimension(size(ncount)) :: lncount
character(len=*), intent(in) :: name

lncid = ncid
ncstatus = nf90_inq_varid(lncid,name,lvid)
lstart(:) = start(:)
lncount(:) = ncount(:)
ncstatus = nf90_get_var(lncid,lvid,vdat,start=lstart,count=lncount)
call ncmsg("get_vara_int3i_t",ncstatus)

return
end subroutine ccnf_get_vara_int3i_t

#ifndef i8r8
subroutine ccnf_get_vara_double1r_s(ncid,vid,start,vdat)

integer, intent(in) :: ncid, vid, start
integer ncstatus
integer(kind=4) lncid, lvid
integer(kind=4), dimension(1) :: lstart
integer(kind=4), dimension(1) :: lncount
real(kind=8), intent(out) :: vdat
real(kind=8), dimension(1) :: lvdat

lncid = ncid
lvid = vid
lstart = start
lncount = 1
ncstatus = nf90_get_var(lncid,lvid,lvdat,start=lstart,count=lncount)
call ncmsg("get_vara_real1r",ncstatus)
vdat = lvdat(1)

return
end subroutine ccnf_get_vara_double1r_s

subroutine ccnf_get_vara_double2r_t(ncid,name,start,ncount,vdat)

integer, intent(in) :: ncid
integer, dimension(:), intent(in) :: start, ncount
integer ncstatus
integer(kind=4) lncid, lvid
integer(kind=4), dimension(size(start)) :: lstart
integer(kind=4), dimension(size(ncount)) :: lncount
real(kind=8), dimension(:), intent(out) :: vdat
character(len=*), intent(in) :: name

lncid = ncid
ncstatus = nf90_inq_varid(lncid,name,lvid)
call ncmsg(name,ncstatus)
lstart(:) = start(:)
lncount(:) = ncount(:)
ncstatus = nf90_get_var(lncid,lvid,vdat,start=lstart,count=lncount)
call ncmsg("get_vara_double2r_t",ncstatus)

return
end subroutine ccnf_get_vara_double2r_t 

subroutine ccnf_get_vara_double4r_s(ncid,vid,start,ncount,vdat)

integer, intent(in) :: ncid, vid
integer, dimension(:) :: start, ncount
integer ncstatus
integer(kind=4) :: lncid, lvid
integer(kind=4), dimension(size(start)) :: lstart
integer(kind=4), dimension(size(ncount)) :: lncount
real(kind=8), dimension(:,:,:), intent(out) :: vdat

lncid = ncid
lvid = vid
lstart = start
lncount = ncount
ncstatus = nf90_get_var(lncid,lvid,vdat,start=lstart,count=lncount)
call ncmsg("get_vara_double4d",ncstatus)

return
end subroutine ccnf_get_vara_double4r_s
#endif

subroutine ccnf_get_att_text(ncid,vid,aname,atext,ierr)

integer, intent(in) :: ncid, vid
integer, intent(out), optional :: ierr
integer ncstatus
integer(kind=4) lncid, lvid
character(len=*), intent(in) :: aname
character(len=*), intent(out) :: atext

atext = ''
lncid = ncid
lvid = vid
ncstatus = nf90_get_att(lncid,lvid,aname,atext)
if (present(ierr)) then
  ierr = ncstatus
else
  if ( ncstatus/=nf90_noerr ) then  
    write(6,*) "ERROR: Cannot read ",trim(aname)  
  end if     
  call ncmsg(aname,ncstatus)
end if

return
end subroutine ccnf_get_att_text

subroutine ccnf_get_att_real(ncid,vid,aname,vdat,ierr)

integer, intent(in) :: ncid, vid
integer, intent(out), optional :: ierr
integer ncstatus
integer(kind=4) lncid, lvid
character(len=*), intent(in) :: aname
real, intent(out) :: vdat

lncid = ncid
lvid = vid
ncstatus = nf90_get_att(lncid,lvid,aname,vdat)
if (present(ierr)) then
  ierr = ncstatus
else
  if ( ncstatus/=nf90_noerr ) then  
    write(6,*) "ERROR: Cannot read ",trim(aname)  
  end if     
  call ncmsg(aname,ncstatus)
end if

return
end subroutine ccnf_get_att_real

subroutine ccnf_get_att_realg1r(ncid,aname,vdat,ierr)

integer, intent(in) :: ncid
integer, intent(out), optional :: ierr
integer ncstatus
integer(kind=4) lncid
character(len=*), intent(in) :: aname
real, intent(out) :: vdat

lncid = ncid
ncstatus = nf90_get_att(lncid,nf90_global,aname,vdat)
if (present(ierr)) then
  ierr = ncstatus
else
  if ( ncstatus/=nf90_noerr ) then  
    write(6,*) "ERROR: Cannot read ",trim(aname)  
  end if  
  call ncmsg("get_attg",ncstatus)
end if

return
end subroutine ccnf_get_att_realg1r

subroutine ccnf_get_att_realg2r(ncid,aname,vdat)

integer, intent(in) :: ncid
integer ncstatus
integer(kind=4) lncid
character(len=*), intent(in) :: aname
real, dimension(:), intent(out) :: vdat

lncid = ncid
ncstatus = nf90_get_att(lncid,nf90_global,aname,vdat(:))
if ( ncstatus/=nf90_noerr ) then  
  write(6,*) "ERROR: Cannot read ",aname  
end if   
call ncmsg("get_attg",ncstatus)

return
end subroutine ccnf_get_att_realg2r

subroutine ccnf_get_att_intg1i(ncid,aname,vdat,tst)

integer, intent(in) :: ncid
integer, intent(out) :: vdat
integer ncstatus
integer(kind=4) :: lncid
integer(kind=4) :: lvdat
logical, intent(out), optional :: tst
character(len=*), intent(in) :: aname

lncid = ncid
ncstatus = nf90_get_att(lncid,nf90_global,aname,lvdat)
vdat = lvdat
if (present(tst)) then
  tst=ncstatus/=nf90_noerr
else
  if ( ncstatus/=nf90_noerr ) then  
    write(6,*) "ERROR: Cannot read ",trim(aname)  
  end if  
  call ncmsg("get_attg",ncstatus)
end if

return
end subroutine ccnf_get_att_intg1i

subroutine ccnf_get_att_intg2i(ncid,aname,vdat)

integer, intent(in) :: ncid
integer, dimension(:), intent(out) :: vdat
integer ncstatus
integer(kind=4) lncid
integer(kind=4), dimension(size(vdat)) :: lvdat
character(len=*), intent(in) :: aname

lncid = ncid
ncstatus = nf90_get_att(lncid,nf90_global,aname,lvdat)
vdat = lvdat
if ( ncstatus/=nf90_noerr ) then  
  write(6,*) "ERROR: Cannot read ",trim(aname)  
end if 
call ncmsg("get_attg",ncstatus)

return
end subroutine ccnf_get_att_intg2i

subroutine ccnf_get_att_textg(ncid,aname,atext,ierr)

integer, intent(in) :: ncid
integer, intent(out), optional :: ierr
integer ncstatus
integer(kind=4) lncid
character(len=*), intent(in) :: aname
character(len=*), intent(out) :: atext

atext = ''
lncid = ncid
ncstatus = nf90_get_att(lncid,nf90_global,aname,atext)
if (present(ierr)) then
  ierr = ncstatus
else
  if ( ncstatus/=nf90_noerr ) then  
    write(6,*) "ERROR: Cannot read ",trim(aname)  
  end if     
  call ncmsg(aname,ncstatus)
end if

return
end subroutine ccnf_get_att_textg

subroutine ccnf_read(fname,vname,vdat)

use cc_mpi
use newmpar_m

integer ncstatus
integer(kind=4) lncid, lvid
real, dimension(ifull), intent(out) :: vdat
real, dimension(:), allocatable :: vdat_g
character(len=*), intent(in) :: fname
character(len=*), intent(in) :: vname

if (myid==0) then
  allocate( vdat_g(ifull_g) )  
  ncstatus = nf90_open(fname,nf90_nowrite,lncid)
  call ncmsg(fname,ncstatus)
  ncstatus = nf90_inq_varid(lncid,vname,lvid)
  call ncmsg(fname,ncstatus)
  ncstatus = nf90_get_var(lncid,lvid,vdat_g)
  call ncmsg(fname,ncstatus)
  ncstatus = nf90_close(lncid)
  call ncmsg(fname,ncstatus)
  call ccmpi_distribute(vdat,vdat_g)
  deallocate( vdat_g )
else
  call ccmpi_distribute(vdat)
end if

return
end subroutine ccnf_read

subroutine ccnf_put_var_text2t_0(ncid,vid,vtxt)

integer, intent(in) :: ncid, vid
integer ncstatus
integer(kind=4) lncid, lvid
character(len=*), dimension(:), intent(in) :: vtxt

lncid = ncid
lvid = vid
ncstatus = nf90_put_var(lncid,lvid,vtxt)
call ncmsg("put_var",ncstatus)

return
end subroutine ccnf_put_var_text2t_0

subroutine ccnf_put_var_int2i_0(ncid,vid,vdat)

integer, intent(in) :: ncid,vid
integer ncstatus
integer, dimension(:), intent(in) :: vdat
integer(kind=4) lncid, lvid

lncid = ncid
lvid = vid
ncstatus = nf90_put_var(lncid,lvid,vdat)
call ncmsg("put_var",ncstatus)

return
end subroutine ccnf_put_var_int2i_0

subroutine ccnf_put_var_int3i_0(ncid,vid,vdat)

integer, intent(in) :: ncid, vid
integer ncstatus
integer, dimension(:,:), intent(in) :: vdat
integer(kind=4) lncid, lvid

lncid = ncid
lvid = vid
ncstatus = nf90_put_var(lncid,lvid,vdat)
call ncmsg("put_var",ncstatus)

return
end subroutine ccnf_put_var_int3i_0

subroutine ccnf_put_vara_real1r_t(ncid,name,start,vdat)

integer, intent(in) :: ncid
integer ncstatus
integer, intent(in) :: start
integer(kind=4) lncid, lvid
integer(kind=4), dimension(1) :: lstart
integer(kind=4), dimension(1) :: lncount
real, intent(in) :: vdat
real, dimension(1) :: lvdat
character(len=*), intent(in) :: name

lncid = ncid
ncstatus = nf90_inq_varid(lncid,name,lvid)
call ncmsg(name,ncstatus)
lstart = start
lncount = 1
lvdat = vdat
ncstatus = nf90_put_var(lncid,lvid,lvdat,start=lstart,count=lncount)
call ncmsg("put_vara_real1r",ncstatus)

return
end subroutine ccnf_put_vara_real1r_t

subroutine ccnf_put_vara_real1r_s(ncid,vid,start,vdat)

integer, intent(in) :: ncid, vid
integer ncstatus
integer, intent(in) :: start
integer(kind=4) lncid, lvid
integer(kind=4), dimension(1) :: lstart
integer(kind=4), dimension(1) :: lncount
real, intent(in) :: vdat
real, dimension(1) :: lvdat

lncid = ncid
lvid = vid
lstart = start
lncount = 1
lvdat = vdat
ncstatus = nf90_put_var(lncid,lvid,lvdat,start=lstart,count=lncount)
call ncmsg("put_vara_real1r",ncstatus)

return
end subroutine ccnf_put_vara_real1r_s

subroutine ccnf_put_vara_real2r_t(ncid,name,start,ncount,vdat)

integer, intent(in) :: ncid
integer ncstatus
integer, intent(in) :: start, ncount
integer(kind=4) lncid, lvid
integer(kind=4), dimension(1) :: lstart
integer(kind=4), dimension(1) :: lncount
real, dimension(:), intent(in) :: vdat
character(len=*), intent(in) :: name

lncid = ncid
ncstatus = nf90_inq_varid(lncid,name,lvid)
call ncmsg("put_vara_real2r_t"//trim(name),ncstatus)
lstart = start
lncount = ncount
ncstatus = nf90_put_var(lncid,lvid,vdat,start=lstart,count=lncount)
call ncmsg("put_vara_real2r_t "//trim(name),ncstatus)

return
end subroutine ccnf_put_vara_real2r_t

subroutine ccnf_put_vara_real2r_0(ncid,vid,start,ncount,vdat)

integer, intent(in) :: ncid, vid
integer ncstatus
integer, intent(in) :: start, ncount
integer(kind=4) lncid, lvid
integer(kind=4), dimension(1) :: lstart
integer(kind=4), dimension(1) :: lncount
real, dimension(:), intent(in) :: vdat

lncid = ncid
lvid = vid
lstart = start
lncount = ncount
ncstatus = nf90_put_var(lncid,lvid,vdat,start=lstart,count=lncount)
call ncmsg("put_vara_real2r_s",ncstatus)

return
end subroutine ccnf_put_vara_real2r_0

subroutine ccnf_put_vara_real3r_t(ncid,name,start,ncount,vdat)

integer, intent(in) :: ncid
integer ncstatus
integer, dimension(2), intent(in) :: start, ncount
integer(kind=4) lncid, lvid
integer(kind=4), dimension(2) :: lstart
integer(kind=4), dimension(2) :: lncount
real, dimension(:,:), intent(in) :: vdat
character(len=*), intent(in) :: name

lncid = ncid
ncstatus = nf90_inq_varid(lncid,name,lvid)
call ncmsg("put_vara_real3r_t",ncstatus)
lstart(:) = start(:)
lncount(:) = ncount(:)
ncstatus = nf90_put_var(lncid,lvid,vdat,start=lstart,count=lncount)
call ncmsg("put_vara_real3r_t",ncstatus)

return
end subroutine ccnf_put_vara_real3r_t

subroutine ccnf_put_vara_real2r_s(ncid,vid,start,ncount,vdat)

integer, intent(in) :: ncid, vid
integer ncstatus
integer, dimension(:), intent(in) :: start, ncount
integer(kind=4) lncid, lvid
integer(kind=4), dimension(size(start)) :: lstart
integer(kind=4), dimension(size(ncount)) :: lncount
real, dimension(:), intent(in) :: vdat

lncid = ncid
lvid = vid
lstart = start
lncount = ncount
ncstatus = nf90_put_var(lncid,lvid,vdat,start=lstart,count=lncount)
call ncmsg("put_vara_real2r",ncstatus)

return
end subroutine ccnf_put_vara_real2r_s

subroutine ccnf_put_vara_real3r_s(ncid,vid,start,ncount,vdat)

integer, intent(in) :: ncid, vid
integer ncstatus
integer, dimension(:), intent(in) :: start, ncount
integer(kind=4) :: lncid, lvid
integer(kind=4), dimension(size(start)) :: lstart
integer(kind=4), dimension(size(ncount)) :: lncount
real, dimension(:,:), intent(in) :: vdat

lncid = ncid
lvid = vid
lstart = start
lncount = ncount
ncstatus = nf90_put_var(lncid,lvid,vdat,start=lstart,count=lncount)
call ncmsg("put_vara",ncstatus)

return
end subroutine ccnf_put_vara_real3r_s

#ifndef i8r8
subroutine ccnf_put_vara_double1r_s(ncid,vid,start,vdat)

integer, intent(in) :: ncid, vid, start
integer ncstatus
integer(kind=4) lncid, lvid
integer(kind=4), dimension(1) :: lstart
integer(kind=4), dimension(1) :: lncount
real(kind=8), intent(in) :: vdat
real(kind=8), dimension(1) :: lvdat

lncid = ncid
lvid = vid
lstart = start
lncount = 1
lvdat = vdat
ncstatus = nf90_put_var(lncid,lvid,lvdat,start=lstart,count=lncount)
call ncmsg("put_vara",ncstatus)

return
end subroutine ccnf_put_vara_double1r_s

subroutine ccnf_put_vara_double1r_t(ncid,name,start,vdat)

integer, intent(in) :: ncid
integer ncstatus
integer, intent(in) :: start
integer(kind=4) lncid, lvid
integer(kind=4), dimension(1) :: lstart
integer(kind=4), dimension(1) :: lncount
real(kind=8), intent(in) :: vdat
real(kind=8), dimension(1) :: lvdat
character(len=*), intent(in) :: name

lncid = ncid
ncstatus = nf90_inq_varid(lncid,name,lvid)
call ncmsg(name,ncstatus)
lstart = start
lncount = 1
lvdat = vdat
ncstatus = nf90_put_var(lncid,lvid,lvdat,start=lstart,count=lncount)
call ncmsg("put_vara_double1r",ncstatus)

return
end subroutine ccnf_put_vara_double1r_t

subroutine ccnf_put_vara_double2r_s(ncid,vid,start,ncount,vdat)

integer, intent(in) :: ncid, vid
integer ncstatus
integer, dimension(:), intent(in) :: start, ncount
integer(kind=4) lncid, lvid
integer(kind=4), dimension(size(start)) :: lstart
integer(kind=4), dimension(size(ncount)) :: lncount
real(kind=8), dimension(:), intent(in) :: vdat

lncid = ncid
lvid = vid
lstart = start
lncount = ncount
ncstatus = nf90_put_var(lncid,lvid,vdat,start=lstart,count=lncount)
call ncmsg("put_vara",ncstatus)

return
end subroutine ccnf_put_vara_double2r_s
#endif

subroutine ccnf_put_vara_int1i_t(ncid,name,start,vdat)

integer, intent(in) :: ncid
integer ncstatus
integer, intent(in) :: start
integer(kind=4) lncid, lvid
integer(kind=4), dimension(1) :: lstart
integer(kind=4), dimension(1) :: lncount
integer, intent(in) :: vdat
integer, dimension(1) :: lvdat
character(len=*), intent(in) :: name

lncid = ncid
ncstatus = nf90_inq_varid(lncid,name,lvid)
call ncmsg("put_vara_int1i",ncstatus)
lstart = start
lncount = 1
lvdat = vdat
ncstatus = nf90_put_var(lncid,lvid,lvdat,start=lstart,count=lncount)
call ncmsg("put_vara_int1i",ncstatus)

return
end subroutine ccnf_put_vara_int1i_t

subroutine ccnf_put_vara_int1i_s(ncid,vid,start,vdat)

integer, intent(in) :: ncid, vid
integer ncstatus
integer, intent(in) :: start
integer(kind=4) lncid, lvid
integer(kind=4), dimension(1) :: lstart
integer(kind=4), dimension(1) :: lncount
integer, intent(in) :: vdat
integer, dimension(1) :: lvdat

lncid = ncid
lvid = vid
lstart = start
lncount = 1
lvdat = vdat
ncstatus = nf90_put_var(lncid,lvid,lvdat,start=lstart,count=lncount)
call ncmsg("put_vara_int1i",ncstatus)

return
end subroutine ccnf_put_vara_int1i_s

subroutine ccnf_put_vara_int2i_s(ncid,vid,start,ncount,vdat)

integer, intent(in) :: ncid, vid
integer ncstatus
integer, dimension(:), intent(in) :: start, ncount
integer, dimension(:), intent(in) :: vdat
integer(kind=4) lncid, lvid
integer(kind=4), dimension(size(start)) :: lstart
integer(kind=4), dimension(size(ncount)) :: lncount

lncid = ncid
lvid = vid
lstart = start
lncount = ncount
ncstatus = nf90_put_var(lncid,lvid,vdat,start=lstart,count=lncount)
call ncmsg("put_vara",ncstatus)

return
end subroutine ccnf_put_vara_int2i_s

subroutine ccnf_put_att_text(ncid,vid,aname,atext)

integer, intent(in) :: ncid, vid
integer ncstatus
integer(kind=4) lncid, lvid
character(len=*), intent(in) :: aname
character(len=*), intent(in) :: atext

lncid = ncid
lvid = vid
ncstatus = nf90_put_att(lncid,lvid,aname,atext)
call ncmsg("put_att",ncstatus)

return
end subroutine ccnf_put_att_text

subroutine ccnf_put_att_rs(ncid,vid,aname,rval)

integer, intent(in) :: ncid, vid
integer ncstatus
integer(kind=4) lncid, lvid
real, intent(in) :: rval
character(len=*), intent(in) :: aname

lncid = ncid
lvid = vid
ncstatus = nf90_put_att(lncid,lvid,aname,rval)
call ncmsg("put_att",ncstatus)

return
end subroutine ccnf_put_att_rs

subroutine ccnf_put_att_textg(ncid,aname,atext)

integer, intent(in) :: ncid
integer ncstatus
integer(kind=4) lncid
character(len=*), intent(in) :: aname
character(len=*), intent(in) :: atext

lncid = ncid
ncstatus = nf90_put_att(lncid,nf90_global,aname,atext)
call ncmsg("put_attg",ncstatus)

return
end subroutine ccnf_put_att_textg

subroutine ccnf_put_att_intg1(ncid,aname,vdat)

integer, intent(in) :: ncid
integer ncstatus
integer, intent(in) :: vdat
integer(kind=4) lncid
integer(kind=4), dimension(1) :: ldat
character(len=*), intent(in) :: aname

lncid = ncid
ldat(1) = vdat
ncstatus = nf90_put_att(lncid,nf90_global,aname,ldat)
call ncmsg("put_attg",ncstatus)

return
end subroutine ccnf_put_att_intg1

subroutine ccnf_put_att_intg2(ncid,aname,vdat)

integer, intent(in) :: ncid
integer ncstatus
integer, dimension(:), intent(in) :: vdat
integer(kind=4) lncid
integer(kind=4), dimension(size(vdat)) :: lvdat
character(len=*), intent(in) :: aname

lncid = ncid
lvdat = vdat
ncstatus = nf90_put_att(lncid,nf90_global,aname,lvdat)
call ncmsg("put_att",ncstatus)

return
end subroutine ccnf_put_att_intg2

subroutine ccnf_put_att_realg1(ncid,aname,vdat)

integer, intent(in) :: ncid
integer ncstatus
integer(kind=4) lncid
real, intent(in) :: vdat
real, dimension(1) :: ldat
character(len=*), intent(in) :: aname

lncid = ncid
ldat(1) = vdat
ncstatus = nf90_put_att(lncid,nf90_global,aname,ldat)
call ncmsg("put_attg",ncstatus)

return
end subroutine ccnf_put_att_realg1

subroutine ccnf_put_att_realg2(ncid,aname,vdat)

integer, intent(in) :: ncid
integer ncstatus
integer(kind=4) lncid
real, dimension(:), intent(in) :: vdat
character(len=*), intent(in) :: aname

lncid = ncid
ncstatus = nf90_put_att(lncid,nf90_global,aname,vdat)
call ncmsg("put_attg",ncstatus)

return
end subroutine ccnf_put_att_realg2

subroutine surfread(dat,varname,netcdfid,filename)

use cc_mpi
use newmpar_m

integer, intent(in), optional :: netcdfid
integer ifull_l
character(len=*), intent(in), optional :: filename
character(len=*), intent(in) :: varname
real, dimension(:), intent(out) :: dat

ifull_l = size(dat)
if ( myid==0 ) then
  if ( present(filename) ) then
    call surfreadglob(dat,varname,filename=filename)  
  else if ( present(netcdfid) ) then
    call surfreadglob(dat,varname,netcdfid=netcdfid)  
  else
    write(6,*) 'Failed to specify input file'
    call ccmpi_abort(-1)
  end if
else
  if ( ifull_l==ifull ) then
    call ccmpi_distribute(dat)
  end if
end if

return
end subroutine surfread

! Read surface data and distribute over processors
! This version suports both netcdf and text file formats
subroutine surfreadglob(dat,vname,netcdfid,filename)

use cc_mpi
use newmpar_m
use parm_m
use parmgeom_m

integer, intent(in), optional :: netcdfid
integer, dimension(3) :: spos, npos
integer ifull_l, ncidx, iernc, varid, ierr
integer ilx, jlx, ndims
character(len=*), intent(in), optional :: filename
character(len=*), intent(in) :: vname
character(len=47) header
real, dimension(:), intent(out) :: dat
real, dimension(:), allocatable :: glob2d
real rlong0x, rlat0x, schmidtx, dsx
logical tst

ifull_l = size(dat)
allocate( glob2d(ifull_g) )

iernc = 0
ncidx = 0
if (present(filename)) then
  call ccnf_open(filename,ncidx,iernc)
else if (present(netcdfid)) then
  ncidx = netcdfid
else
  write(6,*) "ERROR: call to surfreadglobal without filename or netcdfid"
  call ccmpi_abort(-1)
end if

if ( iernc==0 ) then ! Netcdf file

  call ccnf_inq_dimlen(ncidx,'longitude',ilx)
  call ccnf_inq_dimlen(ncidx,'latitude',jlx)
  call ccnf_get_attg(ncidx,'lon0',rlong0x)
  call ccnf_get_attg(ncidx,'lat0',rlat0x)
  call ccnf_get_attg(ncidx,'schmidt',schmidtx)
  if(ilx/=il_g.or.jlx/=jl_g.or.abs(rlong0x-rlong0)>=1.e-20.or.abs(rlat0x-rlat0)>=1.e-20.or. &
      abs(schmidtx-schmidt)>=1.e-20) then
    write(6,*) 'wrong data file supplied'
    if ( ilx/=il_g ) write(6,*) "ilx,il_g ",ilx,il_g
    if ( jlx/=jl_g ) write(6,*) "jlx,jl_g ",jlx,jl_g
    if ( rlong0x/=rlong0 ) write(6,*) "rlong0x,rlong0 ",rlong0x,rlong0
    if ( rlat0x/=rlat0 ) write(6,*) "rlat0x,rlat0 ",rlat0x,rlat0
    if ( schmidtx/=schmidt ) write(6,*) "schmidtx,schmidt ",schmidtx,schmidt
    call ccmpi_abort(-1)
  end if
  spos(1:3)=1
  npos(1)=il_g
  npos(2)=6*il_g
  npos(3)=1
  call ccnf_inq_varid(ncidx,vname,varid,tst)
  if ( .not.tst ) then
    call ccnf_inq_varndims(ncidx,varid,ndims)
    call ccnf_get_vara(ncidx,varid,spos(1:ndims),npos(1:ndims),glob2d)
  else
    glob2d(:)=0.
  end if
  if (present(filename)) then
    call ccnf_close(ncidx)
  end if
  
else ! ASCII file

  open(87,file=filename,status='old')
  read(87,*,iostat=ierr) ilx,jlx,rlong0x,rlat0x,schmidtx,dsx,header
  if ( ierr == 0 ) then
    write(6,*) ilx,jlx,rlong0x,rlat0x,schmidtx,dsx,header
    if(ilx/=il_g.or.jlx/=jl_g.or.abs(rlong0x-rlong0)>=1.e-20.or.abs(rlat0x-rlat0)>=1.e-20.or. &
      abs(schmidtx-schmidt)>=1.e-20) then
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
    write(6,*) "error in surfreadglob",trim(filename),ierr
    call ccmpi_abort(-1)
  end if
  
end if

! distrubte data over processors
if ( ifull_l==ifull ) then
  call ccmpi_distribute(dat, glob2d)
else if ( ifull_l==ifull_g ) then
  dat = glob2d
else
  write(6,*) "ERROR: Invalid array size in surfreadglob"
  call ccmpi_abort(-1)
end if
  
deallocate( glob2d )

return
end subroutine surfreadglob

!--------------------------------------------------------------
! Trap netcdf error messages
subroutine ncmsg(txt,ierr)

use cc_mpi

integer, intent(in) :: ierr
integer(kind=4) :: lierr
character(len=*), intent(in) :: txt

if (ierr/=nf90_noerr) then
  lierr=ierr
  write(6,*) "ERROR: Netcdf error = ",ierr," on rank = ",myid
  write(6,*) txt," ",nf90_strerror(lierr)
  call ccmpi_abort(-1)
end if

return
end subroutine ncmsg

!--------------------------------------------------------------
! define mask for output
subroutine lsmask(inp,outp,msk)

use newmpar_m

real, dimension(ifull), intent(in) :: inp
real, dimension(ifull), intent(out) :: outp
logical, dimension(ifull), intent(in) :: msk

where ( msk )
  outp = inp
elsewhere
  outp = nf90_fill_float
end where

return
end subroutine lsmask

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

end module infile