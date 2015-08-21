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
    
module infile
      
! This module contains routines for reading netcdf files,
! vertical interpolation and some calendar functions.
! This module also contains the interface to the netcdf
! library.
      
! This version of infile.f90 supports parallel localhist input
! files.  Multiple processors read as many input files as
! supplied in parallel and then send this data to all processors
! for interpolation. The code can also identify restart files,
! in which case no additional message passing is required.

#ifdef usenc3
use netcdf_m
#else
use netcdf
#endif

implicit none
            
private
public vertint, datefix, getzinp, ncmsg
public histopen, histclose, histrd1, histrd4, pfall, ncidold
public attrib, histwrt3, histwrt4, freqwrite, surfread
public ccnf_open, ccnf_create, ccnf_close, ccnf_sync, ccnf_enddef
public ccnf_redef, ccnf_nofill, ccnf_inq_varid, ccnf_inq_dimid
public ccnf_inq_dimlen, ccnf_inq_varndims, ccnf_def_dim
public ccnf_def_dimu, ccnf_def_var, ccnf_def_var0, ccnf_get_var
public ccnf_get_vara, ccnf_get_var1, ccnf_get_att, ccnf_get_attg
public ccnf_read, ccnf_put_var, ccnf_put_var1, ccnf_put_vara
public ccnf_put_att, ccnf_put_attg
public pil_g, pjl_g, pka_g, pko_g, mynproc

interface ccnf_get_att
  module procedure ccnf_get_att_text, ccnf_get_att_real
end interface ccnf_get_att
interface ccnf_get_attg
  module procedure ccnf_get_att_intg1i, ccnf_get_att_intg2i
  module procedure ccnf_get_att_realg1r, ccnf_get_att_realg2r
end interface ccnf_get_attg
interface ccnf_get_vara
  module procedure ccnf_get_vara_real2r, ccnf_get_vara_real3r, ccnf_get_vara_real4r 
  module procedure ccnf_get_vara_int2i
#ifndef i8r8
  module procedure ccnf_get_vara_double4d
#endif
end interface ccnf_get_vara
interface ccnf_get_var1
  module procedure ccnf_get_var1_real, ccnf_get_var1_int
end interface ccnf_get_var1
interface ccnf_get_var
  module procedure ccnf_get_var_real, ccnf_get_var_int
end interface ccnf_get_var
interface ccnf_put_att
  module procedure ccnf_put_att_text
end interface ccnf_put_att
interface ccnf_put_attg
  module procedure ccnf_put_att_realg1, ccnf_put_att_realg2
  module procedure ccnf_put_att_intg1, ccnf_put_att_intg2
  module procedure ccnf_put_att_textg
end interface ccnf_put_attg
interface ccnf_put_vara
  module procedure ccnf_put_vara_real2r, ccnf_put_vara_real3r
  module procedure ccnf_put_vara_int2i
#ifndef i8r8
  module procedure ccnf_put_vara_double2r
#endif
end interface ccnf_put_vara
interface ccnf_put_var
  module procedure ccnf_put_var_text2r
  module procedure ccnf_put_var_int2i, ccnf_put_var_int3i
end interface ccnf_put_var
interface ccnf_put_var1
  module procedure ccnf_put_var1_real, ccnf_put_var1_int
#ifndef i8r8
  module procedure ccnf_put_var1_double
#endif
end interface ccnf_put_var1

integer(kind=4), dimension(:), allocatable, save :: pncid
integer, save :: ncidold = -1
integer, save :: mynproc
integer, save :: pil_g, pjl_g, pka_g, pko_g
integer, save :: comm_ip
logical, save :: ptest, pfall

integer(kind=2), parameter :: minv = -32500
integer(kind=2), parameter :: maxv =  32500
integer(kind=2), parameter :: missval = -32501
      
contains

!--------------------------------------------------------------
! Interface for reading 2D+time fields
subroutine histrd1(iarchi,ier,name,ik,var,ifull,nogather)
      
use cc_mpi

implicit none

include 'parm.h'
      
integer iarchi,ier,ik,ifull
character(len=*) name
real, dimension(:), intent(inout) :: var ! may be dummy argument from myid/=0
logical, intent(in), optional :: nogather
logical ngflag

call START_LOG(histrd1_begin)

ngflag = .false.
if ( present(nogather) ) then
  ngflag = nogather
end if

if ( (ifull/=6*ik*ik.and.ptest) .or. ngflag ) then
  ! read local arrays without gather and distribute (e.g., restart file)
  call hr1p(iarchi,ier,name,.true.,var)
  if ( ier==0 .and. mod(ktau,nmaxpr)==0 .and. myid==0 ) then
    write(6,'("done histrd1 ",a8,i4,i3)') name,ier,iarchi
  end if

else if ( myid==0 ) then
  ! read global arrays for myid==0 (possibly distribute or possibly keep on myid==0)
  ! split up processors to save memory.  No need to allocate global arrays on myid/=0.
  call hr1a(iarchi,ier,name,ik,var,ifull)

else if ( ifull/=6*ik*ik ) then
  ! read local arrays with gather and distribute (i.e., change in number of processors)
  call hr1p(iarchi,ier,name,.false.)
  call ccmpi_distribute(var)

else
  ! read global arrays and gather on myid==0
  call hr1p(iarchi,ier,name,.false.)

end if

call END_LOG(histrd1_end)

return
end subroutine histrd1    

!--------------------------------------------------------------------
! Gather 2D+time fields on myid==0
subroutine hr1a(iarchi,ier,name,ik,var,ifull)
      
use cc_mpi
      
implicit none
      
include 'parm.h'
      
integer, intent(in) :: iarchi, ik, ifull
integer, intent(out) :: ier
character(len=*), intent(in) :: name
real, dimension(:), intent(inout) :: var
real, dimension(6*ik*ik) :: globvar
real vmax, vmin

call hr1p(iarchi,ier,name,.false.,globvar)

if ( ier==0 .and. mod(ktau,nmaxpr)==0 ) then
  vmax = maxval(globvar)
  vmin = minval(globvar)
  write(6,'("done histrd1 ",a8,i4,i3,3e14.6)') name,ier,iarchi,vmin,vmax,globvar(id+(jd-1)*ik)
end if

if ( ifull==6*ik*ik ) then
  ! read global arrays for myid==0
  var(1:ifull)=globvar(:) ! really ifull_g
else
  ! read local arrays with gather and distribute (i.e., change in number of processors)
  call ccmpi_distribute(var,globvar)
endif

return
end subroutine hr1a  

!--------------------------------------------------------------
! This subroutine reads 2D+time input files
      
! when qtest=.true. the input grid decomposition should
! match the current processor decomposition.  We can then
! skip the MPI gather and distribute steps.
subroutine hr1p(iarchi,ier,name,qtest,var)
      
use cc_mpi
      
implicit none

include 'newmpar.h'

integer, intent(in) :: iarchi
integer, intent(out) :: ier
integer(kind=4), dimension(3) :: start, ncount
integer ipf, ca
integer(kind=4) idv, ndims
real, dimension(:), intent(inout), optional :: var
real, dimension(pil*pjl*pnpan) :: rvar
real(kind=4) laddoff, lsf
logical, intent(in) :: qtest
character(len=*), intent(in) :: name

start  = (/ 1, 1, iarchi /)
ncount = (/ pil, pjl*pnpan, 1 /)
ier = 0
      
do ipf=0,mynproc-1

  rvar(:)=0. ! default for missing field
  
  ! get variable idv
#ifdef usenc3
  ier=nf_inq_varid(pncid(ipf),name,idv)
  if(ier/=nf_noerr)then
    if (myid==0.and.ipf==0) then
      write(6,*) '***absent field for ncid,name,ier: ',pncid(0),name,ier
    end if
  else
    ! obtain scaling factors and offsets from attributes
    ier=nf_get_att_real(pncid(ipf),idv,'add_offset',laddoff)
    if (ier/=nf_noerr) laddoff=0.
    ier=nf_get_att_real(pncid(ipf),idv,'scale_factor',lsf)
    if (ier/=nf_noerr) lsf=1.
    ier=nf_inq_varndims(pncid(ipf),idv,ndims)
    call START_LOG(ncgetv_begin)
#ifdef i8r8
    ier=nf_get_vara_double(pncid(ipf),idv,start(1:ndims),ncount(1:ndims),rvar)
#else
    ier=nf_get_vara_real(pncid(ipf),idv,start(1:ndims),ncount(1:ndims),rvar)
#endif
    call END_LOG(ncgetv_end)
    call ncmsg(name,ier)
    ! unpack compressed data
    rvar(:)=rvar(:)*real(lsf)+real(laddoff)
  end if ! ier
#else
  ier=nf90_inq_varid(pncid(ipf),name,idv)
  if(ier/=nf90_noerr)then
    if (myid==0.and.ipf==0) then
      write(6,*) '***absent field for ncid,name,ier: ',pncid(0),name,ier
    end if
  else
    ! obtain scaling factors and offsets from attributes
    ier=nf90_get_att(pncid(ipf),idv,'add_offset',laddoff)
    if (ier/=nf90_noerr) laddoff=0.
    ier=nf90_get_att(pncid(ipf),idv,'scale_factor',lsf)
    if (ier/=nf90_noerr) lsf=1.
    ier=nf90_inquire_variable(pncid(ipf),idv,ndims=ndims)
    call START_LOG(ncgetv_begin)
    ier=nf90_get_var(pncid(ipf),idv,rvar,start=start(1:ndims),count=ncount(1:ndims))
    call END_LOG(ncgetv_end)
    call ncmsg(name,ier)
    ! unpack compressed data
    rvar(:)=rvar(:)*real(lsf)+real(laddoff)
  end if ! ier
#endif
      
  if (qtest) then
    ! e.g., restart file or nogather=.true.
    ca = pil*pjl*pnpan*ipf
    var(1+ca:pil*pjl*pnpan+ca)=rvar(:)
  else
    ! e.g., mesonest file
    if ( myid==0 ) then
      call host_hr1p(ipf,rvar,var)
    else
      call proc_hr1p(rvar)
    end if
  end if ! qtest

end do ! ipf

return
end subroutine hr1p

subroutine host_hr1p(ipf,rvar,var)

use cc_mpi

implicit none

include 'newmpar.h'

integer, intent(in) :: ipf
integer jpf, ip, n, no, ca, cc, j
real, dimension(:), intent(inout) :: var
real, dimension(pil*pjl*pnpan), intent(in) :: rvar
real, dimension(pil*pjl*pnpan,fnresid) :: gvar 

call ccmpi_gatherx(gvar,rvar,0,comm_ip)
do jpf = 1,fnresid
  ip = ipf*fnresid + jpf - 1
  do n = 0,pnpan-1
    no = n - pnoff(ip) + 1
    ca = pioff(ip,no) + (pjoff(ip,no)-1)*pil_g + no*pil_g*pil_g
    cc = n*pil*pjl - pil
    do j = 1,pjl
      var(1+j*pil_g+ca:pil+j*pil_g+ca) = gvar(1+j*pil+cc:pil+j*pil+cc,jpf)
    end do
  end do
end do

return
end subroutine host_hr1p

subroutine proc_hr1p(rvar)

use cc_mpi

implicit none

real, dimension(pil*pjl*pnpan), intent(in) :: rvar
real, dimension(0) :: gvar 

call ccmpi_gatherx(gvar,rvar,0,comm_ip)

return
end subroutine proc_hr1p

!--------------------------------------------------------------   
! Interface for reading 3D+time fields
subroutine histrd4(iarchi,ier,name,ik,kk,var,ifull,nogather)
      
use cc_mpi
      
implicit none
      
include 'parm.h'
      
integer, intent(in) :: iarchi,ik,kk,ifull
integer, intent(out) :: ier
character(len=*), intent(in) :: name
real, dimension(:,:), intent(inout) :: var ! may be dummy argument from myid/=0
logical, intent(in), optional :: nogather
logical ngflag

call START_LOG(histrd4_begin)

ngflag = .false.
if ( present(nogather) ) then
  ngflag = nogather
end if

if ( (ifull/=6*ik*ik.and.ptest) .or. ngflag ) then
  ! read local arrays without gather and distribute
  call hr4p(iarchi,ier,name,kk,.true.,var)
  if ( ier==0 .and. mod(ktau,nmaxpr)==0 .and. myid==0 ) then
    write(6,'("done histrd4 ",a8,i4,i3)') name,ier,iarchi
  end if

else if ( myid==0 ) then  
  ! split up processors to save memory.  No need to allocate global
  ! arrays on myid/=0.
  call hr4sa(iarchi,ier,name,ik,kk,var,ifull)

else if ( ifull/=6*ik*ik ) then
  ! read local arrays with gather and distribute
  call hr4p(iarchi,ier,name,kk,.false.)
  call ccmpi_distribute(var)

else
  ! read global arrays and gather on myid==0
  call hr4p(iarchi,ier,name,kk,.false.)

end if

call END_LOG(histrd4_end)

return
end subroutine histrd4

!--------------------------------------------------------------------
! Gather 3D+time fields on myid==0

subroutine hr4sa(iarchi,ier,name,ik,kk,var,ifull)
      
use cc_mpi
      
implicit none
      
include 'parm.h'
      
integer, intent(in) :: iarchi, ik, kk, ifull
integer, intent(out) :: ier
character(len=*), intent(in) :: name
real, dimension(:,:) :: var
real, dimension(6*ik*ik,kk) :: globvar
real vmax, vmin
      
call hr4p(iarchi,ier,name,kk,.false.,globvar)     

if( ier==0 .and. mod(ktau,nmaxpr)==0 ) then
  vmax = maxval(globvar)
  vmin = minval(globvar)
  write(6,'("done histrd4 ",a6,i3,i4,i3,3f12.4)') name,kk,ier,iarchi,vmin,vmax,globvar(id+(jd-1)*ik,nlv)
end if

! Have to return correct value of ier on all processes because it's 
! used for initialisation in calling routine
if ( ifull==6*ik*ik ) then
  ! read global arrays for myid==0
  var(1:ifull,1:kk)=globvar(:,:)
else
  ! read local arrays with gather and distribute (i.e., change in number of processors)
  call ccmpi_distribute(var,globvar)
endif

return
end subroutine hr4sa      

!--------------------------------------------------------------
! This subroutine reads 3D+time input files

! when qtest=.true. the input grid decomposition should
! match the current processor decomposition.  We can then
! skip the MPI gather and distribute steps.
subroutine hr4p(iarchi,ier,name,kk,qtest,var)
      
use cc_mpi
      
implicit none

include 'newmpar.h'

integer, intent(in) :: iarchi, kk
integer, intent(out) :: ier
integer(kind=4), dimension(4) :: start, ncount
integer ipf, k, ca
integer(kind=4) idv, ndims
real, dimension(:,:), intent(inout), optional :: var
real, dimension(pil*pjl*pnpan,kk) :: rvar
real(kind=4) laddoff, lsf
logical, intent(in) :: qtest
character(len=*), intent(in) :: name
character(len=80) :: newname

ier = 0
      
do ipf = 0,mynproc-1

  ! get variable idv
#ifdef usenc3
  ier = nf_inq_varid(pncid(ipf),name,idv)
  if ( ier==nf_noerr ) then
    start = (/ 1, 1, 1, iarchi /)
    ncount = (/ pil, pjl*pnpan, kk, 1 /)      
    ! obtain scaling factors and offsets from attributes
    ier = nf_get_att_real(pncid(ipf),idv,'add_offset',laddoff)
    if ( ier/=nf_noerr ) laddoff=0.
    ier = nf_get_att_real(pncid(ipf),idv,'scale_factor',lsf)
    if ( ier/=nf_noerr ) lsf=1.
    ier = nf_inq_varndims(pncid(ipf),idv,ndims)
    call START_LOG(ncgetv_begin)
#ifdef i8r8
    ier = nf_get_vara_double(pncid(ipf),idv,start(1:ndims),ncount(1:ndims),rvar)
#else
    ier = nf_get_vara_real(pncid(ipf),idv,start(1:ndims),ncount(1:ndims),rvar)
#endif
    call END_LOG(ncgetv_end)
    call ncmsg(name,ier)
    ! unpack data
    rvar(:,:) = rvar(:,:)*real(lsf)+real(laddoff)    
  else
    start(1:3) = (/ 1, 1, iarchi /)
    ncount(1:3) = (/ pil, pjl*pnpan, 1 /)
    do k = 1,kk
      write(newname,'("'//trim(name)//'",I3.3)') k
      ier = nf_inq_varid(pncid(ipf),newname,idv)
      if ( ier/=nf_noerr .and. k<100 ) then
        write(newname,'("'//trim(name)//'",I2.2)') k
        ier = nf_inq_varid(pncid(ipf),newname,idv)
      end if
      if ( ier/=nf_noerr .and. k<10 ) then
        write(newname,'("'//trim(name)//'",I1.1)') k
        ier = nf_inq_varid(pncid(ipf),newname,idv)          
      end if
      if ( ier/=nf_noerr ) exit
      ! obtain scaling factors and offsets from attributes
      ier = nf_get_att_real(pncid(ipf),idv,'add_offset',laddoff)
      if ( ier/=nf_noerr ) laddoff=0.
      ier = nf_get_att_real(pncid(ipf),idv,'scale_factor',lsf)
      if ( ier/=nf_noerr ) lsf=1.
      ier = nf_inq_varndims(pncid(ipf),idv,ndims)
      call START_LOG(ncgetv_begin)
#ifdef i8r8
      ier = nf_get_vara_double(pncid(ipf),idv,start(1:ndims),ncount(1:ndims),rvar(:,k))
#else
      ier = nf_get_vara_real(pncid(ipf),idv,start(1:ndims),ncount(1:ndims),rvar(:,k))
#endif
      call END_LOG(ncgetv_end)
      call ncmsg(name,ier)
      ! unpack data
      rvar(:,k) = rvar(:,k)*real(lsf)+real(laddoff)  
    end do
    if ( ier/=nf_noerr ) then
      if ( myid==0 .and. ipf==0 ) then
        write(6,*) '***absent field for ncid,name,ier: ',pncid(0),name,ier
      end if
      rvar(:,:) = 0. ! default value for missing field
    end if
  end if ! ier
#else
  ier = nf90_inq_varid(pncid(ipf),name,idv)
  if ( ier==nf90_noerr ) then
    start = (/ 1, 1, 1, iarchi /)
    ncount = (/ pil, pjl*pnpan, kk, 1 /)   
    ! obtain scaling factors and offsets from attributes
    ier = nf90_get_att(pncid(ipf),idv,'add_offset',laddoff)
    if ( ier/=nf90_noerr ) laddoff=0.
    ier = nf90_get_att(pncid(ipf),idv,'scale_factor',lsf)
    if ( ier/=nf90_noerr ) lsf=1.
    ier = nf90_inquire_variable(pncid(ipf),idv,ndims=ndims)
    call START_LOG(ncgetv_begin)
    ier = nf90_get_var(pncid(ipf),idv,rvar,start=start(1:ndims),count=ncount(1:ndims))
    call END_LOG(ncgetv_end)
    call ncmsg(name,ier)
    ! unpack data
    rvar(:,:) = rvar(:,:)*real(lsf)+real(laddoff)
  else
    start(1:3) = (/ 1, 1, iarchi /)
    ncount(1:3) = (/ pil, pjl*pnpan, 1 /)
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
      if ( ier/=nf90_noerr ) exit
      ! obtain scaling factors and offsets from attributes
      ier = nf90_get_att(pncid(ipf),idv,'add_offset',laddoff)
      if ( ier/=nf90_noerr ) laddoff=0.
      ier = nf90_get_att(pncid(ipf),idv,'scale_factor',lsf)
      if ( ier/=nf90_noerr ) lsf=1.
      ier = nf90_inquire_variable(pncid(ipf),idv,ndims=ndims)
      call START_LOG(ncgetv_begin)
      ier = nf90_get_var(pncid(ipf),idv,rvar(:,k),start=start(1:ndims),count=ncount(1:ndims))
      call END_LOG(ncgetv_end)
      call ncmsg(name,ier)
      ! unpack data
      rvar(:,k) = rvar(:,k)*real(lsf)+real(laddoff)      
    end do
    if ( ier/=nf90_noerr ) then
      if ( myid==0 .and. ipf==0 ) then
        write(6,*) '***absent field for ncid,name,ier: ',pncid(0),name,ier
      end if
      rvar(:,:) = 0. ! default value for missing field
    end if
  end if ! ier
#endif

  if ( qtest ) then
    ! e.g., restart file or nogather=.true.
    ca = pil*pjl*pnpan*ipf
    var(1+ca:pil*pjl*pnpan+ca,1:kk) = rvar(:,:)
  else
    ! e.g., mesonest file
    if ( myid==0 ) then
      call host_hr4p(ipf,kk,rvar,var)
    else
      call proc_hr4p(kk,rvar)
    end if
  end if ! qtest

end do ! ipf

return
end subroutine hr4p

subroutine host_hr4p(ipf,kk,rvar,var)

use cc_mpi

implicit none

include 'newmpar.h'

integer, intent(in) :: ipf, kk
integer jpf, ip, n, no, ca, cc, j, k
real, dimension(:,:), intent(inout) :: var
real, dimension(pil*pjl*pnpan,kk), intent(in) :: rvar
real, dimension(pil*pjl*pnpan,kk,fnresid) :: gvar 

call ccmpi_gatherx(gvar,rvar,0,comm_ip)
do jpf = 1,fnresid
  ip = ipf*fnresid + jpf - 1   ! local file number
  do k = 1,kk
    do n = 0,pnpan-1
      no = n - pnoff(ip) + 1   ! global panel number of local file
      ca = pioff(ip,no) + pjoff(ip,no)*pil_g + no*pil_g*pil_g - pil_g
      cc = n*pil*pjl - pil
      do j = 1,pjl
        var(1+j*pil_g+ca:pil+j*pil_g+ca,k) = gvar(1+j*pil+cc:pil+j*pil+cc,k,jpf)
      end do
    end do
  end do
end do

return
end subroutine host_hr4p

subroutine proc_hr4p(kk,rvar)

use cc_mpi

implicit none

integer, intent(in) :: kk
real, dimension(pil*pjl*pnpan,kk), intent(in) :: rvar
real, dimension(0,0) :: gvar 

call ccmpi_gatherx(gvar,rvar,0,comm_ip)

return
end subroutine proc_hr4p

!--------------------------------------------------------------
! This subroutine opens parallel input files
subroutine histopen(ncid,ifile,ier)
      
use cc_mpi
      
implicit none
      
include 'newmpar.h'
      
integer, parameter :: nihead = 54
      
integer, dimension(nihead) :: ahead
integer, dimension(0:5) :: duma, dumb
integer, dimension(10) :: idum
integer, intent(out) :: ncid, ier
integer is, ipf, dmode
integer ipin, nxpr, nypr
integer ltst, der, myrank, pkl
integer(kind=4), dimension(nihead) :: lahead
integer(kind=4) lncid, lidum, ldid
character(len=*), intent(in) :: ifile
character(len=170) pfile
character(len=8) fdecomp

if (myid==0) then
  ! attempt to open single file with myid==0
#ifdef usenc3
  ier=nf_open(ifile,nf_nowrite,lncid)
#else
  ier=nf90_open(ifile,nf90_nowrite,lncid)
#endif
  ncid=lncid
  fnproc=1      ! number of files to be read over all processors
  dmode=0       ! Single file (dmode=0), Face decomposition (dmode=1), Depreciated (dmode=2) or Uniform decomposition (dmode=3)
  pil=0         ! Number of X grid points within a file panel
  pjl=0         ! Number of Y grid points within a file panel
  pnpan=0       ! Number of panels in file
  ptest=.false. ! Files match current processor (e.g., Restart file), allowing MPI gather/scatter to be avoided
  pfall=.false. ! Every processor has been assigned at least one file, no need to Bcast metadata data
      
  ! attempt to open parallel files
#ifdef usenc3
  if (ier/=nf_noerr) then
    write(pfile,"(a,'.',i6.6)") trim(ifile), 0
    ier=nf_open(pfile,nf_nowrite,lncid)
    ncid=lncid
    if (ier/=nf_noerr) then
      write(6,*) "WARN: Cannot open ",trim(pfile)
      write(6,*) "WARN: Cannot open ",trim(ifile)
    else  
      write(6,*) "Found parallel input file ",trim(ifile)
      fdecomp=''
      der=nf_get_att_int(lncid,nf_global,"nproc",lidum)
      fnproc=lidum
      call ncmsg("nproc",der)
      der=nf_get_att_text(lncid,nf_global,"decomp",fdecomp)
      call ncmsg("decomp",der)
#else
  if (ier/=nf90_noerr) then
    write(pfile,"(a,'.',i6.6)") trim(ifile), 0
    ier=nf90_open(pfile,nf90_nowrite,lncid)
    ncid=lncid
    if (ier/=nf90_noerr) then
      write(6,*) "WARN: Cannot open ",trim(pfile)
      write(6,*) "WARN: Cannot open ",trim(ifile)
    else  
      write(6,*) "Found parallel input file ",trim(ifile)
      fdecomp=''
      der=nf90_get_att(lncid,nf90_global,"nproc",lidum)
      fnproc=lidum
      call ncmsg("nproc",der)
      der=nf90_get_att(lncid,nf90_global,"decomp",fdecomp)
      call ncmsg("decomp",der)
#endif
      select case(fdecomp)
        case('face')
          dmode=1
        case('uniform')  ! old uniform
          dmode=2
        case('uniform1') ! new uniform (Dix style)
          dmode=3
        case default
          write(6,*) "ERROR: Unknown decomposition ",trim(fdecomp)
          call ccmpi_abort(-1)
      end select
    end if
  else
    ! nproc should only exist in multi-file input
#ifdef usenc3
    der=nf_get_att_int(lncid,nf_global,"nproc",lidum)
    if (der==nf_noerr) then
#else
    der=nf90_get_att(lncid,nf90_global,"nproc",lidum)
    if (der==nf90_noerr) then
#endif
      write(6,*) "ERROR: Incorrect base filename"
      call ccmpi_abort(-1)
    end if
    write(6,*) "Found single input file ",trim(ifile)
  end if

#ifdef usenc3
  if (ier==nf_noerr) then
    der=nf_get_att_int(lncid,nf_global,"int_header",lahead)
    ahead=lahead
    call ncmsg("int_header",der)
    der=nf_inq_dimid(lncid,"olev",ldid)
    if ( der==nf_noerr ) then
      der=nf_inq_dimlen(lncid,ldid,pko_g)
      call ncmsg("olev",der)
    else
      pko_g=0
    end if
#else
  if (ier==nf90_noerr) then
    der=nf90_get_att(lncid,nf90_global,"int_header",lahead)
    ahead=lahead
    call ncmsg("int_header",der)
    der=nf90_inq_dimid(lncid,"olev",ldid)
    if ( der==nf90_noerr ) then
      der=nf90_inquire_dimension(lncid,ldid,len=pko_g)
      call ncmsg("olev",der)
    else
      pko_g=0
    end if
#endif
    pil_g=ahead(1)
    pjl_g=ahead(2)
    pka_g=ahead(3)
        
    if (allocated(pioff)) then
      deallocate(pioff,pjoff,pnoff)
    end if
    allocate(pioff(0:fnproc-1,0:5),pjoff(0:fnproc-1,0:5))
    allocate(pnoff(0:fnproc-1))
        
    select case(dmode)
      case(0) ! no decomposition
        pnpan=1
        pnoff=1
        pioff=0
        pjoff=0
        pil=pil_g
        pjl=pjl_g
      case(1) ! face decomposition
        pnpan=max(1,6/fnproc)
        do ipf=0,fnproc-1
          call face_set(pil,pjl,pnoff(ipf),duma,dumb,pnpan,pil_g,ipf,fnproc,nxpr,nypr)
          pioff(ipf,:)=duma
          pjoff(ipf,:)=dumb
        end do
      case(2) ! old uniform decomposition
        pnpan=6
        do ipf=0,fnproc-1
          call uniform_set(pil,pjl,pnoff(ipf),duma,dumb,pnpan,pil_g,ipf,fnproc,nxpr,nypr)
          pioff(ipf,:)=duma
          pjoff(ipf,:)=dumb
        end do
      case(3) ! new uniform decomposition
        pnpan=6
        do ipf=0,fnproc-1
          call dix_set(pil,pjl,pnoff(ipf),duma,dumb,pnpan,pil_g,ipf,fnproc,nxpr,nypr)
          pioff(ipf,:)=duma
          pjoff(ipf,:)=dumb
        end do
    end select

    ptest=.false.
#ifdef uniform_decomp
    if (dmode==3) then
      if (nproc==fnproc) then
        if (pil_g==il_g.and.pjl_g==jl_g) then
          ptest=.true.
        end if
      end if
    end if
#else
    if (dmode==1) then
      if (nproc==fnproc) then
        if (pil_g==il_g.and.pjl_g==jl_g) then
          ptest=.true.
        end if
      end if
    end if
#endif

    write(6,*) "Found pil_g,pjl_g,fnproc,dmode,ptest ",pil_g,pjl_g,fnproc,dmode,ptest
    
  end if

  idum(1)=fnproc
  idum(2)=pil
  idum(3)=pjl
  idum(4)=pnpan
  if (ptest) then
    idum(5)=1
  else
    idum(5)=0      
  end if
  idum(6)=ier
  idum(7)=pka_g
  idum(8)=pko_g
  idum(9)=pil_g
  idum(10)=pjl_g
  
  write(6,*) "Broadcasting file metadata"
end if

! Broadcast file metadata
call ccmpi_bcast(idum(1:10),0,comm_world)
fnproc=idum(1)      ! number of files to be read
pil   =idum(2)      ! width of panel in each file
pjl   =idum(3)      ! length of panel in each file
pnpan =idum(4)      ! number of panels in each file
ptest =(idum(5)==1) ! test for match between files and processes
ier   =idum(6)      ! file error flag
pka_g =idum(7)      ! number of atmosphere levels
pko_g =idum(8)      ! number of atmosphere levels
pil_g =idum(9)      ! grid size
pjl_g =idum(10)     ! grid size

#ifdef usenc3
if (ier/=nf_noerr) return
#else
if (ier/=nf90_noerr) return
#endif

if ( myid==0 ) then
  write(6,*) "Opening data files"
end if

! calculate number of files to be read on this processor
fnresid = min( fnproc, nproc ) 
do while ( mod(fnproc,fnresid)/=0 )
  fnresid = fnresid - 1     ! limit on processor ranks that will read files    
end do
if ( myid<fnresid) then
  mynproc = fnproc/fnresid  ! calculate the number of files to be read per process
else
  mynproc = 0
end if

! allocate array of file handles.  Unallocated implies no files to be read on this processor
if ( allocated(pncid) ) then
  write(6,*) "ERROR: Cannot open new parallel output until old file is closed"
  call ccmpi_abort(-1)
end if
if ( mynproc>0 ) then
  allocate(pncid(0:mynproc-1))
end if

! Rank 0 can start with the second file, because the first file has already been opened
if ( myid==0 ) then 
  is=1
  pncid(0)=ncid
else
  is=0
end if

! loop through files to be opened by this processor
do ipf = is,mynproc-1
  ipin=ipf*fnresid+myid
  write(pfile,"(a,'.',i6.6)") trim(ifile), ipin
#ifdef usenc3
  der=nf_open(pfile,nf_nowrite,pncid(ipf))
  if ( der/=nf_noerr ) then
#else
  der=nf90_open(pfile,nf90_nowrite,pncid(ipf))
  if ( der/=nf90_noerr ) then
#endif
    write(6,*) "ERROR: Cannot open ",pfile
    call ncmsg("open",der)
  end if
end do

if ( myid==0 ) then
  write(6,*) "Splitting comms for distributing file data with fnresid ",fnresid
end if

! define comm group to read the residual files
if ( myid<fnresid ) then
  ltst=0
  myrank=myid
else
  ltst=-1 ! undefined
  myrank=myid-fnresid
end if
call ccmpi_commsplit(comm_ip,comm_world,ltst,myrank)

pfall=fnresid==nproc  ! are all processes associated with a file?
                      ! this means we do not need to Bcst file metadata
if ( mynproc>0 ) then
  ncid=pncid(0)       ! set ncid to the first file handle as onthefly
                      ! assumes changes in ncid reflect a new file
                      ! and hence updates the metadata
end if

if ( myid==0 ) then
  write(6,*) "Broadcast file coordinate data"
else
  allocate(pioff(0:fnproc-1,0:5),pjoff(0:fnproc-1,0:5))
  allocate(pnoff(0:fnproc-1))
end if

call ccmpi_bcast(pioff,0,comm_world)
call ccmpi_bcast(pjoff,0,comm_world)
call ccmpi_bcast(pnoff,0,comm_world)

if ( myid==0 ) then
  write(6,*) "Create file RMA windows"
end if
pkl=max(pka_g,pko_g)
call ccmpi_filewincreate(pkl)

if ( myid==0 ) then
  write(6,*) "Ready to read data from input file"
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
#ifdef usenc3
    ierr = nf_close(pncid(ipf))
#else
    ierr = nf90_close(pncid(ipf))
#endif
  end do
  deallocate(pncid)
end if
call ccmpi_commfree(comm_ip)
if (allocated(pioff)) then
  deallocate(pioff,pjoff,pnoff)
end if
call ccmpi_filewinfree
ncidold = -1 ! flag onthefly to load metadata

return
end subroutine histclose

!--------------------------------------------------------------
! transforms 3d array from dimension kk in vertical to kl   jlm

! This version of vertint can interpolate from host models with
! a greater number of vertical levels than the nested model.
subroutine vertint(told,t,n,kk,sigin)

use sigs_m

implicit none

include 'newmpar.h'
include 'parm.h'

integer, intent(in) :: kk, n
integer k, kin
integer, dimension(:), allocatable, save :: ka
integer, save :: kk_save = -1
integer, save :: klapse = 0
real, dimension(:,:), intent(out) :: t
real, dimension(ifull,kk), intent(in) :: told
real, dimension(kk), intent(in) :: sigin
real, dimension(:), allocatable, save :: sigin_save
real, dimension(:), allocatable, save :: wta
      
if ( kk==kl ) then
  if ( all(abs(sig-sigin)<0.0001) ) then
    t(1:ifull,1:kl)=told(:,:)
    return
  end if
end if

if ( kk_save/=kk ) then
  if ( allocated(sigin_save) ) then
    deallocate(sigin_save)
  end if
  allocate(sigin_save(kk))
  sigin_save = 0.
  kk_save    = kk
  if ( .not.allocated(wta) ) then
    allocate(wta(kl),ka(kl))
  end if
end if

if ( any(abs(sigin-sigin_save)>=0.0001) ) then
  sigin_save = sigin
  klapse = 0
  kin    = 2
  do k=1,kl
    if ( sig(k)>=sigin(1) ) then
      ka(k)=2
      wta(k)=0.
      klapse=k   ! i.e. T lapse correction for k<=klapse
    else if ( sig(k)<=sigin(kk) ) then   ! at top
      ka(k)=kk
      wta(k)=1.
    else
      do while ( sig(k)<=sigin(kin) .and. kin < kk )
        kin = kin + 1
      end do
      ka(k)=kin
      wta(k)=(sigin(kin-1)-sig(k))/(sigin(kin-1)-sigin(kin))
    endif  !  (sig(k)>=sigin(1)) ... ...
  enddo   !  k loop
end if
#ifdef debug
if ( myid==0 ) then
  write(6,*) 'in vertint kk,kl ',kk,kl
  write(6,"('sigin',10f7.4)") (sigin(k),k=1,kk)
  write(6,"('sig  ',10f7.4)") sig
  write(6,*) 'ka ',ka
  write(6,*) 'kb ',kb
  write(6,"('wta',10f7.4)") wta
  write(6,"('wtb',10f7.4)") wtb
endif   !  (myid==0)
#endif

do k=1,kl
  ! N.B. "a" denotes "above", "b" denotes "below"
  t(1:ifull,k)=wta(k)*told(:,ka(k))+(1.-wta(k))*told(:,ka(k)-1)
enddo    ! k loop

if ( n==1 .and. klapse/=0 ) then  ! for T lapse correction
  do k=1,klapse
    ! assume 6.5 deg/km, with dsig=.1 corresponding to 1 km
    t(1:ifull,k)=t(1:ifull,k)+(sig(k)-sigin(1))*6.5/.1
  enddo    ! k loop
else if ( n==2 ) then  ! for qg do a -ve fix
  t(1:ifull,:)=max(t(1:ifull,:),qgmin)
else if ( n==5 ) then  ! for qfg, qlg do a -ve fix
  t(1:ifull,:)=max(t(1:ifull,:),0.)
endif
      
return
end subroutine vertint

!--------------------------------------------------------------------
! This subroutine advances input date by the amount of time defined by mtimer_r
subroutine datefix(kdate_r,ktime_r,mtimer_r)

use cc_mpi

implicit none

include 'newmpar.h'
include 'parm.h'

integer leap
common/leap_yr/leap  ! 1 to allow leap years

integer, intent(inout) :: kdate_r,ktime_r,mtimer_r
integer, dimension(12) :: mdays = (/31,28,31,30,31,30,31,31,30,31,30,31/)
integer iyr,imo,iday,ihr,imins
integer mtimerh,mtimerm,mtimer
integer mdays_save
integer, parameter :: minsday = 1440

if ( kdate_r>=00600000 .and. kdate_r<=00991231 ) then   ! old 1960-1999
  kdate_r=kdate_r+19000000
  if ( myid==0 ) write(6,*) 'For Y2K kdate_r altered to: ',kdate_r
endif

iyr=kdate_r/10000
imo=(kdate_r-10000*iyr)/100
iday=kdate_r-10000*iyr-100*imo
ihr=ktime_r/100
imins=ktime_r-100*ihr
if ( myid==0 ) write(6,*) 'entering datefix iyr,imo,iday,ihr,imins,mtimer_r: ', &
                                            iyr,imo,iday,ihr,imins,mtimer_r

mdays(2)=28
if ( leap==1 ) then
  if ( mod(iyr,4)==0   ) mdays(2)=29
  if ( mod(iyr,100)==0 ) mdays(2)=28
  if ( mod(iyr,400)==0 ) mdays(2)=29
end if
do while ( mtimer_r>minsday*mdays(imo) )
  mtimer_r=mtimer_r-minsday*mdays(imo)
  imo=imo+1
  if ( imo>12 ) then
    imo=1
    iyr=iyr+1
    if ( leap==1 ) then
      mdays(2)=28      
      if ( mod(iyr,4)==0   ) mdays(2)=29
      if ( mod(iyr,100)==0 ) mdays(2)=28
      if ( mod(iyr,400)==0 ) mdays(2)=29
    end if
  end if
end do
if(diag)write(6,*)'b datefix iyr,imo,iday,ihr,imins,mtimer_r: ', &
                             iyr,imo,iday,ihr,imins,mtimer_r

iday=iday+mtimer_r/minsday
mtimer_r=mod(mtimer_r,minsday)
if(diag)write(6,*)'c datefix iyr,imo,iday,ihr,imins,mtimer_r: ', &
                             iyr,imo,iday,ihr,imins,mtimer_r

! at this point mtimer_r has been reduced to fraction of a day
mtimerh=mtimer_r/60
mtimerm=mtimer_r-mtimerh*60  ! minutes left over
ihr=ihr+mtimerh
imins=imins+mtimerm

if ( imins==58 .or. imins==59 ) then
  ! allow for roundoff for real timer from old runs
  if ( myid==0 ) write(6,*)'*** imins increased to 60 from imins = ',imins
  imins=60
endif

ihr=ihr+imins/60
imins=mod(imins,60)
if(diag)write(6,*)'d datefix iyr,imo,iday,ihr,imins,mtimer_r: ', &
                             iyr,imo,iday,ihr,imins,mtimer_r

iday=iday+ihr/24
ihr=mod(ihr,24)
if(diag)write(6,*)'e datefix iyr,imo,iday,ihr,imins,mtimer_r: ', &
                             iyr,imo,iday,ihr,imins,mtimer_r

mdays_save=mdays(imo)
imo=imo+(iday-1)/mdays(imo)
iday=mod(iday-1,mdays_save)+1

iyr=iyr+(imo-1)/12
imo=mod(imo-1,12)+1

kdate_r=iday+100*(imo+100*iyr)
ktime_r=ihr*100+imins
mtimer=0
if(diag)write(6,*)'end datefix iyr,imo,iday,ihr,imins,mtimer_r: ', &
                               iyr,imo,iday,ihr,imins,mtimer_r

if ( myid==0 ) write(6,*)'leaving datefix kdate_r,ktime_r: ',kdate_r,ktime_r

return
end subroutine datefix

!--------------------------------------------------------------------
! Set up number of minutes from beginning of year
subroutine getzinp(fjd,jyear,jmonth,jday,jhour,jmin,mins,allleap)

use cc_mpi

implicit none

include 'dates.h'
include 'parm.h'

integer leap
common/leap_yr/leap  ! 1 to allow leap years

integer, intent(out) :: jyear,jmonth,jday,jhour,jmin ! start date of run
integer, intent(out) :: mins                         ! elapsed time from start of year
integer mstart
integer, dimension(12) :: ndoy
! days from beginning of year (1st Jan is 0)
integer, dimension(12), parameter :: odoy=(/ 0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334 /) 
real, intent(out) :: fjd
logical, intent(in), optional :: allleap
logical lleap

if ( present(allleap) ) then
  lleap=allleap
else
  lleap=.false.    
end if

jyear =kdate/10000
jmonth=(kdate-jyear*10000)/100
jday  =kdate-jyear*10000-jmonth*100
jhour =ktime/100
jmin  =ktime-jhour*100
      
if ( jmonth<1 .or. jmonth>12 ) then
  write(6,*) "ERROR: Invalid month ",jmonth," for kdate ",kdate
  call ccmpi_abort(-1)
end if 

ndoy=odoy
if ( leap==1 .or. lleap ) then
  if ( mod(jyear,4)  ==0 ) ndoy(3:12)=odoy(3:12)+1
  if ( mod(jyear,100)==0 ) ndoy(3:12)=odoy(3:12)
  if ( mod(jyear,400)==0 ) ndoy(3:12)=odoy(3:12)+1
end if

mstart=1440*(ndoy(jmonth)+jday-1) + 60*jhour + jmin ! mins from start of year
! mtimer contains number of minutes since the start of the run.
mins = mtimer + mstart

if ( nhstest<0 ) then  ! aquaplanet test
  fjd = 79.+float(mod(mins,1440))/1440.  ! set to 21 March +frac of day
  mins=nint(fjd*1440.)
else
  fjd = float(mod(mins,(ndoy(12)+31)*1440))/1440.    ! 525600 = 1440*365
endif

return
end subroutine getzinp

!--------------------------------------------------------------------
! DEFINE ATTRIBUTES
subroutine attrib(cdfid,dim,ndim,name,lname,units,xmin,xmax,daily,itype)

use cc_mpi

implicit none

integer, intent(in) :: cdfid, itype, ndim
integer, intent(in) :: daily
integer, dimension(ndim), intent(in) :: dim
integer ier
integer(kind=4) vtype, idv, lcdfid, lndim, lsize
integer(kind=4), dimension(ndim) :: ldim
integer(kind=2), dimension(1) :: ivals2
real, intent(in) :: xmin, xmax
real(kind=4) lscalef, laddoff
real(kind=4), dimension(1) :: rvals
character(len=*), intent(in) :: name
character(len=*), intent(in) :: lname
character(len=*), intent(in) :: units

#ifdef usenc3
if ( itype==1 ) then
  vtype = nf_short
else
#ifdef i8r8
  vtype = nf_double
#else
  vtype = nf_float
#endif
end if
lcdfid = cdfid
lndim  = ndim
ldim   = dim
ier = nf_def_var(lcdfid, name, vtype, lndim, ldim, idv)
call ncmsg("def_var",ier)
lsize = len_trim(lname)
ier = nf_put_att_text(lcdfid,idv,'long_name',lsize,lname)
call ncmsg("long_name",ier)
lsize = len_trim(units)
if ( lsize>0 ) then
  ier = nf_put_att_text(lcdfid,idv,'units',lsize,units)
  call ncmsg("units",ier)
endif
if ( vtype==nf_short ) then
  ivals2 = minv
  ier = nf_put_att_int2(lcdfid,idv,'valid_min',nf_short,1_4,ivals2)
  call ncmsg("valid_min",ier)
  ivals2 = maxv
  ier = nf_put_att_int2(lcdfid,idv,'valid_max',nf_short,1_4,ivals2)
  call ncmsg("valid_max",ier)
  ivals2 = missval
  ier = nf_put_att_int2(lcdfid,idv,'missing_value',nf_short,1_4,ivals2)
  call ncmsg("missing_value",ier)
  lscalef=(xmax-xmin)/(real(maxv)-real(minv))
  laddoff=xmin-real(lscalef)*real(minv)
  rvals=laddoff
  ier = nf_put_att_real(lcdfid,idv,'add_offset',nf_float,1_4,rvals)
  call ncmsg("add_offset",ier)
  rvals=lscalef
  ier = nf_put_att_real(lcdfid,idv,'scale_factor',nf_float,1_4,rvals)
  call ncmsg("scale_factor",ier)
else
  rvals=nf_fill_float
  ier = nf_put_att_real(lcdfid,idv,'missing_value',nf_float,1_4,rvals)
  call ncmsg("missing_value",ier)
endif
ier = nf_put_att_text(lcdfid,idv,'FORTRAN_format',5_4,'G11.4')
call ncmsg("FORTRAN_format",ier)
if( daily>0 )then
  ier = nf_put_att_text(lcdfid,idv,'valid_time',5_4,'daily')
  call ncmsg("valid_time",ier)
endif
#else
if ( itype==1 ) then
  vtype = nf90_short
else
#ifdef i8r8
  vtype = nf90_double
#else
  vtype = nf90_float
#endif
end if
lcdfid = cdfid
ldim   = dim
ier = nf90_def_var(lcdfid, name, vtype, ldim, idv, deflate_level=1_4)
call ncmsg("def_var",ier)
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
if ( daily>0 ) then
  ier = nf90_put_att(lcdfid,idv,'valid_time','daily')
  call ncmsg("valid_time",ier)
endif
#endif
      
return
end subroutine attrib

!--------------------------------------------------------------------
! 3D NETCDF WRITE ARRAY ROUTINES
subroutine histwrt3(var,sname,idnc,iarch,local,lwrite)

use cc_mpi              ! CC MPI routines

implicit none

include 'newmpar.h'     ! Grid parameters

integer, intent(in) :: idnc, iarch
real, dimension(ifull), intent(in) :: var
real, dimension(ifull) :: wvar
character(len=*), intent(in) :: sname
logical, intent(in) :: local, lwrite

if (.not.lwrite) then
#ifdef usenc3
  wvar(:)=real(nf_fill_float)
#else
  wvar(:)=real(nf90_fill_float)
#endif
else
  wvar(:)=var(:)
end if

if ( local ) then
  call fw3l(wvar,sname,idnc,iarch)
else if ( myid==0 ) then
  call fw3a(wvar,sname,idnc,iarch)
else
  call ccmpi_gather(wvar)
end if

return
end subroutine histwrt3

subroutine freqwrite(fncid,cname,fiarch,local,datain)

use cc_mpi               ! CC MPI routines
      
implicit none
      
include 'newmpar.h'      ! Grid parameters
      
integer, intent(in) :: fncid, fiarch
real, dimension(ifull), intent(in) :: datain
logical, intent(in) :: local
character(len=*), intent(in) :: cname
      
if ( local ) then
  call fw3l(datain,cname,fncid,fiarch)
elseif ( myid==0 ) then
  call fw3a(datain,cname,fncid,fiarch)
else
  call ccmpi_gather(datain)
endif
     
return
end subroutine freqwrite

subroutine fw3l(var,sname,idnc,iarch)

use cc_mpi               ! CC MPI routines
      
implicit none
      
include 'newmpar.h'      ! Grid parameters
include 'parm.h'         ! Model configuration
      
integer, intent(in) :: idnc, iarch
integer ier
integer(kind=4) :: lidnc, mid, vtype, ndims
integer(kind=4), dimension(3) :: start, ncount
integer(kind=2), dimension(ifull) :: ipack
real, dimension(ifull), intent(in) :: var
real(kind=4) laddoff, lscale_f
character(len=*), intent(in) :: sname

start = (/ 1, 1, iarch /)
ncount = (/ il, jl, 1 /)

lidnc = idnc
#ifdef usenc3
ier = nf_inq_varid(lidnc,sname,mid)
call ncmsg(sname,ier)
ier = nf_inq_varndims(lidnc,mid,ndims)
ier = nf_inq_vartype(lidnc,mid,vtype)
if ( vtype==nf_short ) then
  if ( all(var>9.8E36) ) then
    ipack = missval
  else
    ier = nf_get_att_real(lidnc,mid,'add_offset',laddoff)
    ier = nf_get_att_real(lidnc,mid,'scale_factor',lscale_f)
    ipack = nint(max(min((var-real(laddoff))/real(lscale_f),real(maxv)),real(minv)),2)
  end if
  ier = nf_put_vara_int2(lidnc,mid,start(1:ndims),ncount(1:ndims),ipack)
else
#ifdef i8r8
  ier = nf_put_vara_double(lidnc,mid,start(1:ndims),ncount(1:ndims),var)
#else
  ier = nf_put_vara_real(lidnc,mid,start(1:ndims),ncount(1:ndims),var)
#endif
end if
call ncmsg(sname,ier)
if ( mod(ktau,nmaxpr)==0 .and. myid==0 ) then
  if ( any(var==real(nf_fill_float)) ) then
    write(6,'("histwrt3 ",a7,i8,a7)') sname,iarch,"missing"
  else
    write(6,'("histwrt3 ",a7,i8)') sname,iarch
  end if
end if
#else
ier = nf90_inq_varid(lidnc,sname,mid)
call ncmsg(sname,ier)
ier = nf90_inquire_variable(lidnc,mid,xtype=vtype,ndims=ndims)
if ( vtype==nf90_short ) then
  if ( all(var>9.8E36) ) then
    ipack = missval
  else
    ier = nf90_get_att(lidnc,mid,'add_offset',laddoff)
    ier = nf90_get_att(lidnc,mid,'scale_factor',lscale_f)
    ipack = nint(max(min((var-real(laddoff))/real(lscale_f),real(maxv)),real(minv)),2)
  end if
  ier = nf90_put_var(lidnc,mid,ipack,start=start(1:ndims),count=ncount(1:ndims))
else
  ier = nf90_put_var(lidnc,mid,var,start=start(1:ndims),count=ncount(1:ndims))
end if
call ncmsg(sname,ier)
if ( mod(ktau,nmaxpr)==0 .and. myid==0 ) then
  if ( any(var==real(nf90_fill_float)) ) then
    write(6,'("histwrt3 ",a7,i8,a7)') sname,iarch,"missing"
  else
    write(6,'("histwrt3 ",a7,i8)') sname,iarch
  end if
end if
#endif

return
end subroutine fw3l
      
subroutine fw3a(var,sname,idnc,iarch)

use cc_mpi               ! CC MPI routines

implicit none

include 'newmpar.h'      ! Grid parameters
include 'parm.h'         ! Model configuration

integer, intent(in) :: idnc, iarch
integer ier, imn, imx, jmn, jmx, iq
integer(kind=4) lidnc, mid, vtype, ndims
integer(kind=4), dimension(3) :: start, ncount
integer(kind=2), dimension(ifull_g) :: ipack
real, dimension(ifull), intent(in) :: var
real, dimension(ifull_g) :: globvar
real varn, varx
real(kind=4) laddoff, lscale_f
character(len=*), intent(in) :: sname
      
call ccmpi_gather(var, globvar)

start = (/ 1, 1, iarch /)
ncount = (/ il_g, jl_g, 1 /)

!     find variable index
lidnc = idnc
#ifdef usenc3
ier = nf_inq_varid(lidnc,sname,mid)
call ncmsg(sname,ier)
ier = nf_inq_varndims(lidnc,mid,ndims)
ier = nf_inq_vartype(lidnc,mid,vtype)
if ( vtype==nf_short ) then
  if ( all(globvar>9.8e36) ) then
    ipack = missval
  else
    ier = nf_get_att_real(lidnc,mid,'add_offset',laddoff)
    ier = nf_get_att_real(lidnc,mid,'scale_factor',lscale_f)
    ipack = nint(max(min((globvar-real(laddoff))/real(lscale_f),real(maxv)),real(minv)),2)
  endif
  ier = nf_put_vara_int2(lidnc,mid,start(1:ndims),ncount(1:ndims),ipack)
else
#ifdef i8r8
  ier = nf_put_vara_double(lidnc,mid,start(1:ndims),ncount(1:ndims),globvar)
#else
  ier = nf_put_vara_real(lidnc,mid,start(1:ndims),ncount(1:ndims),globvar)
#endif
endif
call ncmsg(sname,ier)
#else
ier = nf90_inq_varid(lidnc,sname,mid)
call ncmsg(sname,ier)
ier = nf90_inquire_variable(lidnc,mid,xtype=vtype,ndims=ndims)
if ( vtype==nf90_short ) then
  if ( all(globvar>9.8e36) ) then
    ipack = missval
  else
    ier = nf90_get_att(lidnc,mid,'add_offset',laddoff)
    ier = nf90_get_att(lidnc,mid,'scale_factor',lscale_f)
    ipack = nint(max(min((globvar-real(laddoff))/real(lscale_f),real(maxv)),real(minv)),2)
  endif
  ier = nf90_put_var(lidnc,mid,ipack,start=start(1:ndims),count=ncount(1:ndims))
else
  ier = nf90_put_var(lidnc,mid,globvar,start=start(1:ndims),count=ncount(1:ndims))
endif
call ncmsg(sname,ier)
#endif

if ( mod(ktau,nmaxpr)==0 ) then
#ifdef usenc3
  if ( any(globvar==real(nf_fill_float)) ) then
#else
  if ( any(globvar==real(nf90_fill_float)) ) then
#endif
    write(6,'("histwrt3 ",a7,i4,a7)') sname,iarch,"missing"
  else
    varn = minval(globvar(:))
    varx = maxval(globvar(:))
    ! This should work ???? but sum trick is more portable???
    ! iq = minloc(globvar,dim=1)
    iq = sum(minloc(globvar(:)))
    ! Convert this 1D index to 2D
    imn = 1 + modulo(iq-1,il_g)
    jmn = 1 + (iq-1)/il_g
    iq = sum(maxloc(globvar(:)))
    ! Convert this 1D index to 2D
    imx = 1 + modulo(iq-1,il_g)
    jmx = 1 + (iq-1)/il_g
    write(6,'("histwrt3 ",a7,i4,f12.4,2i4,f12.4,2i4,f12.4)') &
                   sname,iarch,varn,imn,jmn,varx,imx,jmx,    &
                   globvar(id+(jd-1)*il_g)
  end if
end if

return
end subroutine fw3a

!--------------------------------------------------------------------
! 4D NETCDF WRITE ARRAY ROUTINES
subroutine histwrt4(var,sname,idnc,iarch,local,lwrite)

use cc_mpi              ! CC MPI routines

implicit none

include 'newmpar.h'     ! Grid parameters

integer, intent(in) :: idnc, iarch
real, dimension(:,:), intent(in) :: var
real, dimension(ifull,kl) :: wvar
character(len=*), intent(in) :: sname
logical, intent(in) :: local, lwrite

if ( .not.lwrite ) then
#ifdef usenc3
  wvar=real(nf_fill_float)
#else
  wvar=real(nf90_fill_float)
#endif
else
  wvar(:,:)=var(1:ifull,1:kl)
endif

if ( local ) then
  call hw4l(wvar,sname,idnc,iarch)
elseif ( myid==0 ) then
  call hw4a(wvar,sname,idnc,iarch)
else
  call ccmpi_gather(wvar)
endif

return
end subroutine histwrt4
      
subroutine hw4l(var,sname,idnc,iarch)

use cc_mpi               ! CC MPI routines

implicit none

include 'newmpar.h'      ! Grid parameters
include 'parm.h'         ! Model configuration

integer, intent(in) :: idnc, iarch
integer iq, k, ier
integer(kind=4) mid, vtype, lidnc, ndims
integer(kind=4), dimension(4) :: start, ncount
integer(kind=2), dimension(ifull,kl) :: ipack
real, dimension(ifull,kl), intent(in) :: var
real(kind=4) laddoff, lscale_f
character(len=*), intent(in) :: sname

start = (/ 1, 1, 1, iarch /)
ncount = (/ il, jl, kl, 1 /)

lidnc = idnc
#ifdef usenc3
ier = nf_inq_varid(lidnc,sname,mid)
call ncmsg(sname,ier)
ier = nf_inq_varndims(lidnc,mid,ndims)
ier = nf_inq_vartype(lidnc,mid,vtype)
if ( vtype==nf_short ) then
  if ( all(var>9.8e36) ) then
    ipack = missval
  else
    ier = nf_get_att_real(lidnc,mid,'add_offset',laddoff)
    ier = nf_get_att_real(lidnc,mid,'scale_factor',lscale_f)
    do k = 1,kl
      do iq = 1,ifull
        ipack(iq,k) = nint(max(min((var(iq,k)-real(laddoff))/real(lscale_f),real(maxv)),real(minv)),2)
      end do
    end do
  end if
  ier = nf_put_vara_int2(lidnc,mid,start(1:ndims),ncount(1:ndims),ipack)
else
#ifdef i8r8
  ier = nf_put_vara_double(lidnc,mid,start(1:ndims),ncount(1:ndims),var)
#else
  ier = nf_put_vara_real(lidnc,mid,start(1:ndims),ncount(1:ndims),var)
#endif
end if
call ncmsg(sname,ier)
#else
ier = nf90_inq_varid(lidnc,sname,mid)
call ncmsg(sname,ier)
ier = nf90_inquire_variable(lidnc,mid,xtype=vtype,ndims=ndims)
if ( vtype==nf90_short ) then
  if ( all(var>9.8e36) ) then
    ipack = missval
  else
    ier = nf90_get_att(lidnc,mid,'add_offset',laddoff)
    ier = nf90_get_att(lidnc,mid,'scale_factor',lscale_f)
    do k = 1,kl
      do iq = 1,ifull
        ipack(iq,k) = nint(max(min((var(iq,k)-real(laddoff))/real(lscale_f),real(maxv)),real(minv)),2)
      end do
    end do
  end if
  ier = nf90_put_var(lidnc,mid,ipack,start=start(1:ndims),count=ncount(1:ndims))
else
  ier = nf90_put_var(lidnc,mid,var,start=start(1:ndims),count=ncount(1:ndims))
endif
call ncmsg(sname,ier)
#endif

if ( mod(ktau,nmaxpr)==0 .and. myid==0 ) then
#ifdef usenc3
  if ( any(var==real(nf_fill_float)) ) then
#else
  if ( any(var==real(nf90_fill_float)) ) then
#endif
    write(6,'("histwrt4 ",a7,i4,a7)') sname,iarch,"missing"
  else
    write(6,'("histwrt4 ",a7,i4)') sname,iarch
  end if
end if

return
end subroutine hw4l      

subroutine hw4a(var,sname,idnc,iarch)

use cc_mpi              ! CC MPI routines

implicit none

include 'newmpar.h'     ! Grid parameters
include 'parm.h'        ! Model configuration

integer, intent(in) :: idnc, iarch
integer ier, imx, jmx, kmx, iq, k
integer(kind=4) mid, vtype, lidnc, ndims
integer(kind=4), dimension(4) :: start, ncount
integer, dimension(2) :: max_result
integer(kind=2), dimension(ifull_g,kl) :: ipack
real varn, varx
real(kind=4) laddoff, lscale_f
real, dimension(ifull,kl), intent(in) :: var
real, dimension(ifull_g,kl) :: globvar
character(len=*), intent(in) :: sname
      
call ccmpi_gather(var, globvar)
start = (/ 1, 1, 1, iarch /)
ncount = (/ il_g, jl_g, kl, 1 /)

!     find variable index
lidnc = idnc
#ifdef usenc3
ier = nf_inq_varid(lidnc,sname,mid)
call ncmsg(sname,ier)
ier = nf_inq_varndims(lidnc, mid, ndims)
ier = nf_inq_vartype(lidnc, mid, vtype)
if ( vtype==nf_short ) then
  if ( all(globvar>9.8e36) ) then
    ipack = missval
  else
    ier = nf_get_att_real(lidnc,mid,'add_offset',laddoff)
    ier = nf_get_att_real(lidnc,mid,'scale_factor',lscale_f)
    do k = 1,kl
      do iq = 1,ifull_g
        ipack(iq,k) = nint(max(min((globvar(iq,k)-real(laddoff))/real(lscale_f),real(maxv)),real(minv)),2)
      end do
    end do
  endif
  ier = nf_put_vara_int2(lidnc,mid,start(1:ndims),ncount(1:ndims),ipack)
else
#ifdef i8r8
  ier = nf_put_vara_double(lidnc,mid,start(1:ndims),ncount(1:ndims),globvar)
#else
  ier = nf_put_vara_real(lidnc,mid,start(1:ndims),ncount(1:ndims),globvar)
#endif
endif
call ncmsg(sname,ier)
#else
ier = nf90_inq_varid(lidnc,sname,mid)
call ncmsg(sname,ier)
ier = nf90_inquire_variable(lidnc, mid, xtype=vtype, ndims=ndims)
if ( vtype==nf90_short ) then
  if ( all(globvar>9.8e36) )then
    ipack = missval
  else
    ier = nf90_get_att(lidnc,mid,'add_offset',laddoff)
    ier = nf90_get_att(lidnc,mid,'scale_factor',lscale_f)
    do k = 1,kl
      do iq = 1,ifull_g
        ipack(iq,k) = nint(max(min((globvar(iq,k)-real(laddoff))/real(lscale_f),real(maxv)),real(minv)),2)
      end do
    end do
  end if
  ier = nf90_put_var(lidnc,mid,ipack,start=start(1:ndims),count=ncount(1:ndims))
else
  ier = nf90_put_var(lidnc,mid,globvar,start=start(1:ndims),count=ncount(1:ndims))
endif
call ncmsg(sname,ier)
#endif

if ( mod(ktau,nmaxpr)==0 ) then
#ifdef usenc3
  if ( any(globvar==real(nf_fill_float)) ) then
#else
  if ( any(globvar==real(nf90_fill_float)) ) then
#endif
    write(6,'("histwrt4 ",a7,i4,a7)') sname,iarch,"missing"
  else
    varn = minval(globvar)
    varx = maxval(globvar)
    max_result = maxloc(globvar)
    kmx = max_result(2)
    iq = max_result(1)
    ! Convert this 1D index to 2D
    imx = 1 + modulo(iq-1,il_g)
    jmx = 1 + (iq-1)/il_g
    write(6,'("histwrt4 ",a7,i4,2f12.4,3i4,f12.4)') sname,iarch,varn,varx,imx,jmx,kmx,globvar(id+(jd-1)*il_g,nlv)
  end if
end if

return
end subroutine hw4a      

!--------------------------------------------------------------------

subroutine ccnf_open(fname,ncid,status)

use cc_mpi

implicit none

integer, intent(out) :: ncid
integer, intent(out), optional :: status
integer ncstatus
integer(kind=4) lncid
character(len=*), intent(in) :: fname

#ifdef usenc3
ncstatus = nf_open(fname,nf_nowrite,lncid)
#else
ncstatus = nf90_open(fname,nf90_nowrite,lncid)
#endif
ncid=lncid

if ( present(status) ) then
  status=ncstatus
else
  call ncmsg("open",ncstatus)
endif

return
end subroutine ccnf_open

subroutine ccnf_create(fname,ncid)

use cc_mpi

implicit none

integer, intent(out) :: ncid
integer ncstatus
integer(kind=4) lncid
character(len=*), intent(in) :: fname

#ifdef usenc3
  !ncstatus = nf_create(fname,nf_clobber,lncid)
  ncstatus = nf_create(fname,nf_64bit_offset,lncid)
#else
  ncstatus = nf90_create(fname,nf90_netcdf4,lncid)
#endif
ncid=lncid
call ncmsg("create",ncstatus)

return
end subroutine ccnf_create

subroutine ccnf_close(ncid)

use cc_mpi

implicit none

integer, intent(in) :: ncid
integer ncstatus
integer(kind=4) lncid

lncid=ncid
#if usenc3
ncstatus = nf_close(lncid)
#else
ncstatus = nf90_close(lncid)
#endif
call ncmsg("close",ncstatus)

return
end subroutine ccnf_close

subroutine ccnf_nofill(ncid)

use cc_mpi

implicit none

integer, intent(in) :: ncid
integer ncstatus
integer(kind=4) lncid, lomode

lncid = ncid
#ifdef usenc3
ncstatus = nf_set_fill(lncid,nf_nofill,lomode)
#else
ncstatus = nf90_set_fill(lncid,nf90_nofill,lomode)
#endif
call ncmsg("nofill",ncstatus)

return
end subroutine ccnf_nofill

subroutine ccnf_inq_dimid(ncid,dname,did,tst)

use cc_mpi

implicit none

integer, intent(in) :: ncid
integer, intent(out) :: did
integer ncstatus
integer(kind=4) lncid, ldid
logical ltst
logical, intent(out), optional :: tst
character(len=*), intent(in) :: dname

lncid=ncid
#ifdef usenc3
ncstatus=nf_inq_dimid(lncid,dname,ldid)
ltst=(ncstatus/=nf_noerr)
#else
ncstatus=nf90_inq_dimid(lncid,dname,ldid)
ltst=(ncstatus/=nf90_noerr)
#endif
did=ldid

if (present(tst)) then
  tst=ltst
else
  call ncmsg("dimid",ncstatus)
end if

return
end subroutine ccnf_inq_dimid

subroutine ccnf_inq_dimlen(ncid,dname,dlen,failok)

use cc_mpi

implicit none

integer, intent(in) :: ncid
integer, intent(inout) :: dlen
integer(kind=4) lncid, lncstatus, ldid, ldlen
character(len=*), intent(in) :: dname
logical, intent(in), optional :: failok
logical ftest

if (present(failok)) then
  ftest=failok
else
  ftest=.false.
end if

lncid=ncid
ldlen=dlen
#ifdef usenc3
lncstatus=nf_inq_dimid(lncid,dname,ldid)
if (lncstatus/=nf_noerr) then
  if (ftest) return
  write(6,*) nf_strerror(lncstatus)
  call ccmpi_abort(-1)
end if
lncstatus=nf_inq_dimlen(lncid,ldid,ldlen)
if (lncstatus/=nf_noerr) then
  if (ftest) return
  write(6,*) nf_strerror(lncstatus)
  call ccmpi_abort(-1)
end if
#else
lncstatus=nf90_inq_dimid(lncid,dname,ldid)
if (lncstatus/=nf90_noerr) then
  if (ftest) return
  write(6,*) nf90_strerror(lncstatus)
  call ccmpi_abort(-1)
end if
lncstatus=nf90_inquire_dimension(lncid,ldid,len=ldlen)
if (lncstatus/=nf90_noerr) then
  if (ftest) return
  write(6,*) nf90_strerror(lncstatus)
  call ccmpi_abort(-1)
end if
#endif
dlen=ldlen

return
end subroutine ccnf_inq_dimlen

subroutine ccnf_inq_varid(ncid,vname,vid,tst)

use cc_mpi

implicit none

integer, intent(in) :: ncid
integer, intent(out) :: vid
integer ncstatus
integer(kind=4) lncid, lvid
character(len=*), intent(in) :: vname
logical, intent(out), optional :: tst
logical ltst

lncid=ncid
#ifdef usenc3
ncstatus = nf_inq_varid(lncid,vname,lvid)
ltst=(ncstatus/=nf_noerr)
#else
ncstatus = nf90_inq_varid(lncid,vname,lvid)
ltst=(ncstatus/=nf90_noerr)
#endif
vid=lvid

if (present(tst)) then
  tst=ltst
else
  call ncmsg(vname,ncstatus)
end if

return
end subroutine ccnf_inq_varid


subroutine ccnf_inq_varndims(ncid,vid,ndims)

use cc_mpi

implicit none

integer, intent(in) :: ncid, vid
integer, intent(out) :: ndims
integer ncstatus
integer(kind=4) lncid, lvid, lndims

lncid=ncid
lvid=vid
#ifdef usenc3
ncstatus = nf_inq_varndims(lncid,lvid,lndims)
#else
ncstatus = nf90_inquire_variable(lncid,lvid,ndims=lndims)
#endif
call ncmsg('ccnf_inq_varndims',ncstatus)
ndims = lndims

return
end subroutine ccnf_inq_varndims

subroutine ccnf_def_dim(ncid,dname,nsize,did)

use cc_mpi

implicit none

integer, intent(in) :: ncid, nsize
integer, intent(out) :: did
integer ncstatus
integer(kind=4) lncid, lnsize, ldid
character(len=*), intent(in) :: dname

lncid=ncid
lnsize=nsize
#ifdef usenc3
ncstatus=nf_def_dim(lncid,dname,lnsize,ldid)
#else
ncstatus=nf90_def_dim(lncid,dname,lnsize,ldid)
#endif
did=ldid
call ncmsg("def_dim",ncstatus)

return
end subroutine ccnf_def_dim

subroutine ccnf_def_dimu(ncid,dname,did)

use cc_mpi

implicit none

integer, intent(in) :: ncid
integer, intent(out) :: did
integer ncstatus
integer(kind=4) lncid, ldid
character(len=*), intent(in) :: dname

lncid=ncid
#ifdef usenc3
ncstatus=nf_def_dim(lncid,dname,nf_unlimited,ldid)
#else
ncstatus=nf90_def_dim(lncid,dname,nf90_unlimited,ldid)
#endif
did=ldid
call ncmsg("def_dimu",ncstatus)

return
end subroutine ccnf_def_dimu

subroutine ccnf_def_var(ncid,vname,vtype,vndim,dims,vid)

use cc_mpi

implicit none

integer, intent(in) :: ncid, vndim
integer, intent(out) :: vid
integer, dimension(vndim), intent(in) :: dims
integer ncstatus
integer(kind=4) lncid, lvndim, ltype, lvid
integer(kind=4), dimension(vndim) :: ldims
character(len=*), intent(in) :: vname
character(len=*), intent(in) :: vtype

select case(vtype)
#ifdef usenc3
  case('int')
    ltype=nf_int
  case('float')
    ltype=nf_float
  case('double')
    ltype=nf_double
  case('char')
    ltype=nf_char
#else
  case('int')
    ltype=nf90_int
  case('float')
    ltype=nf90_float
  case('double')
    ltype=nf90_double
  case('char')
    ltype=nf90_char
#endif
  case default
    write(6,*) "ERROR: Unknown option for ccnf_def_var ",vtype
    call ccmpi_abort(-1)
end select

lncid=ncid
ldims=dims
#ifdef usenc3
lvndim=vndim
ncstatus = nf_def_var(lncid,vname,ltype,lvndim,ldims,lvid)
#else
ncstatus = nf90_def_var(lncid,vname,ltype,ldims,lvid,deflate_level=1_4)
#endif
vid=lvid
call ncmsg("def_var",ncstatus)

return
end subroutine ccnf_def_var

subroutine ccnf_def_var0(ncid,vname,vtype,vid)

use cc_mpi

implicit none

integer, intent(in) :: ncid
integer, intent(out) :: vid
integer ncstatus
integer(kind=4) lncid, lvid, ltype
integer(kind=4), dimension(1) :: dimids
character(len=*), intent(in) :: vname
character(len=*), intent(in) :: vtype

select case(vtype)
#ifdef usenc3
  case('int')
    ltype=nf_int
  case('float')
    ltype=nf_float
  case('double')
    ltype=nf_double
  case('char')
    ltype=nf_char
#else
  case('int')
    ltype=nf90_int
  case('float')
    ltype=nf90_float
  case('double')
    ltype=nf90_double
  case('char')
    ltype=nf90_char
#endif
  case default
    write(6,*) "ERROR: Unknown option for ccnf_def_var ",vtype
    call ccmpi_abort(-1)
end select

lncid=ncid
#ifdef usenc3
dimids=1
ncstatus = nf_def_var(lncid,vname,ltype,0_4,dimids,lvid)
#else
ncstatus = nf90_def_var(lncid,vname,ltype,lvid)
#endif
vid=lvid
call ncmsg("def_var0",ncstatus)

return
end subroutine ccnf_def_var0

subroutine ccnf_enddef(ncid)

use cc_mpi

implicit none

integer, intent(in) :: ncid
integer ncstatus
integer(kind=4) lncid

lncid=ncid
#ifdef usenc3
ncstatus=nf_enddef(lncid)
#else
ncstatus=nf90_enddef(lncid)
#endif
call ncmsg("enddef",ncstatus)

return
end subroutine ccnf_enddef

subroutine ccnf_redef(ncid)

use cc_mpi

implicit none

integer, intent(in) :: ncid
integer ncstatus
integer(kind=4) lncid

lncid=ncid
#ifdef usenc3
ncstatus=nf_redef(lncid)
#else
ncstatus=nf90_redef(lncid)
#endif
call ncmsg("redef",ncstatus)

return
end subroutine ccnf_redef

subroutine ccnf_sync(ncid)

use cc_mpi

implicit none

integer, intent(in) :: ncid
#ifdef outsync
integer ncstatus
integer(kind=4) lncid

lncid=ncid
#ifdef usenc3
ncstatus = nf_sync(lncid)
#else
ncstatus = nf90_sync(lncid)
#endif
call ncmsg("sync",ncstatus)
#endif

return
end subroutine ccnf_sync

subroutine ccnf_get_var_real(ncid,vname,vdat)

use cc_mpi

implicit none

integer, intent(in) :: ncid
integer ncstatus
integer(kind=4) :: lvid, lncid
real, dimension(:), intent(out) :: vdat
character(len=*), intent(in) :: vname

lncid=ncid
#ifdef usenc3
ncstatus = nf_inq_varid(lncid,vname,lvid)
call ncmsg("get_var_varid",ncstatus)
call START_LOG(ncgetv_begin)
#ifdef i8r8
ncstatus = nf_get_var_double(lncid,lvid,vdat)
#else
ncstatus = nf_get_var_real(lncid,lvid,vdat)
#endif
call END_LOG(ncgetv_end)
call ncmsg("get_var",ncstatus)
#else
ncstatus = nf90_inq_varid(lncid,vname,lvid)
call ncmsg("get_var_varid",ncstatus)
call START_LOG(ncgetv_begin)
ncstatus = nf90_get_var(lncid,lvid,vdat)
call END_LOG(ncgetv_end)
call ncmsg("get_var",ncstatus)
#endif

return
end subroutine ccnf_get_var_real

subroutine ccnf_get_var_int(ncid,vname,vdat)

use cc_mpi

implicit none

integer, intent(in) :: ncid
integer ncstatus
integer, dimension(:), intent(out) :: vdat
integer(kind=4) :: lncid, lvid
integer(kind=4), dimension(size(vdat)) :: lvdat
integer(kind=4), dimension(1) :: lstart, lncount
character(len=*), intent(in) :: vname

lncid=ncid
#ifdef usenc3
ncstatus = nf_inq_varid(lncid,vname,lvid)
call ncmsg("get_var_varid",ncstatus)
lstart = 1
lncount = size(vdat)
call START_LOG(ncgetv_begin)
ncstatus = nf_get_vara_int(lncid,lvid,lstart,lncount,lvdat)
vdat=lvdat
call END_LOG(ncgetv_end)
call ncmsg("get_var",ncstatus)
#else
ncstatus = nf90_inq_varid(lncid,vname,lvid)
call ncmsg("get_var_varid",ncstatus)
call START_LOG(ncgetv_begin)
ncstatus = nf90_get_var(lncid,lvid,vdat)
call END_LOG(ncgetv_end)
call ncmsg("get_var",ncstatus)
#endif

return
end subroutine ccnf_get_var_int

subroutine ccnf_get_var1_real(ncid,vid,start,vdat)

use cc_mpi

implicit none

integer, intent(in) :: ncid, vid, start
integer ncstatus
integer(kind=4) lncid, lvid
#ifdef usenc3
integer(kind=4), dimension(1) :: lstart
#else
integer(kind=4), dimension(1) :: lstart, lcount
real, dimension(1) :: ldat
#endif
real, intent(out) :: vdat
#ifdef i8r8
real(kind=8), dimension(1) :: rvals
#else
real(kind=4), dimension(1) :: rvals
#endif

lncid=ncid
lvid=vid
#ifdef usenc3
lstart = start
call START_LOG(ncgetv_begin)
#ifdef i8r8
ncstatus = nf_get_var1_double(lncid,lvid,lstart(1),rvals(1))
vdat = rvals(1)
#else
ncstatus = nf_get_var1_real(lncid,lvid,lstart(1),rvals(1))
vdat = rvals(1)
#endif
call END_LOG(ncgetv_end)
#else
lstart = start
lcount = 1
call START_LOG(ncgetv_begin)
ncstatus = nf90_get_var(lncid,lvid,ldat,start=lstart,count=lcount)
vdat=ldat(1)
call END_LOG(ncgetv_end)
#endif
call ncmsg("get_var1",ncstatus)

return
end subroutine ccnf_get_var1_real

subroutine ccnf_get_var1_int(ncid,vid,start,vdat)

use cc_mpi

implicit none

integer, intent(in) :: ncid, vid, start
integer ncstatus
integer(kind=4) :: lncid, lvid
#ifdef usenc3
integer(kind=4), dimension(1) :: lstart, ldat
#else
integer(kind=4), dimension(1) :: lstart, lcount
integer, dimension(1) :: ldat
#endif
integer, intent(out) :: vdat

lncid=ncid
lvid=vid
#ifdef usenc3
lstart=start
call START_LOG(ncgetv_begin)
#ifdef i8r8
ncstatus = nf_get_var1_int(lncid,lvid,lstart(1),ldat(1))
#else
ncstatus = nf_get_var1_int(lncid,lvid,lstart(1),ldat(1))
#endif
vdat=ldat(1)
call END_LOG(ncgetv_end)
#else
lstart=start
lcount=1
call START_LOG(ncgetv_begin)
ncstatus = nf90_get_var(lncid,lvid,ldat,start=lstart,count=lcount)
vdat=ldat(1)
call END_LOG(ncgetv_end)
#endif
call ncmsg("get_var1",ncstatus)

return
end subroutine ccnf_get_var1_int

subroutine ccnf_get_vara_real2r(ncid,vid,start,ncount,vdat)

use cc_mpi

implicit none

integer, intent(in) :: ncid, vid
integer, dimension(:), intent(in) :: start, ncount
integer ncstatus
integer(kind=4) lncid, lvid
integer(kind=4), dimension(size(start)) :: lstart
integer(kind=4), dimension(size(ncount)) :: lncount
real, dimension(:), intent(out) :: vdat

lncid=ncid
lvid=vid
lstart=start
lncount=ncount
call START_LOG(ncgetv_begin)
#ifdef usenc3
#ifdef i8r8
ncstatus=nf_get_vara_double(lncid,lvid,lstart,lncount,vdat)
#else
ncstatus=nf_get_vara_real(lncid,lvid,lstart,lncount,vdat)
#endif
#else
ncstatus=nf90_get_var(lncid,lvid,vdat,start=lstart,count=lncount)
#endif
call END_LOG(ncgetv_end)
call ncmsg("get_vara_real2r",ncstatus)

return
end subroutine ccnf_get_vara_real2r 

subroutine ccnf_get_vara_real3r(ncid,vid,start,ncount,vdat)

use cc_mpi

implicit none

integer, intent(in) :: ncid, vid
integer, dimension(:) :: start, ncount
integer ncstatus
integer(kind=4) lncid, lvid
integer(kind=4), dimension(size(start)) :: lstart
integer(kind=4), dimension(size(ncount)) :: lncount
real, dimension(:,:), intent(out) :: vdat

lncid=ncid
lvid=vid
lstart=start
lncount=ncount
call START_LOG(ncgetv_begin)
#ifdef usenc3
#ifdef i8r8
ncstatus=nf_get_vara_double(lncid,lvid,lstart,lncount,vdat)
#else
ncstatus=nf_get_vara_real(lncid,lvid,lstart,lncount,vdat)
#endif
#else
ncstatus=nf90_get_var(lncid,lvid,vdat,start=lstart,count=lncount)
#endif
call END_LOG(ncgetv_end)
call ncmsg("get_vara_real3r",ncstatus)

return
end subroutine ccnf_get_vara_real3r

subroutine ccnf_get_vara_real4r(ncid,vid,start,ncount,vdat)

use cc_mpi

implicit none

integer, intent(in) :: ncid, vid
integer, dimension(:) :: start, ncount
integer ncstatus
integer(kind=4) :: lncid, lvid
integer(kind=4), dimension(size(start)) :: lstart
integer(kind=4), dimension(size(ncount)) :: lncount
real, dimension(:,:,:), intent(out) :: vdat

lncid=ncid
lvid=vid
lstart=start
lncount=ncount
call START_LOG(ncgetv_begin)
#ifdef usenc3
#ifdef i8r8
ncstatus=nf_get_vara_double(lncid,lvid,lstart,lncount,vdat)
#else
ncstatus=nf_get_vara_real(lncid,lvid,lstart,lncount,vdat)
#endif
#else
ncstatus=nf90_get_var(lncid,lvid,vdat,start=lstart,count=lncount)
#endif
call END_LOG(ncgetv_end)
call ncmsg("get_vara_real4r",ncstatus)

return
end subroutine ccnf_get_vara_real4r

subroutine ccnf_get_vara_int2i(ncid,vid,start,ncount,vdat)

use cc_mpi

implicit none

integer, intent(in) :: ncid, vid
integer, dimension(:) :: start, ncount
integer ncstatus
integer, dimension(:), intent(out) :: vdat
integer(kind=4) :: lncid, lvid
integer(kind=4), dimension(size(start)) :: lstart
integer(kind=4), dimension(size(ncount)) :: lncount
integer(kind=4), dimension(size(vdat)) :: lvdat

lncid=ncid
lvid=vid
lstart=start
lncount=ncount
call START_LOG(ncgetv_begin)
#ifdef usenc3
ncstatus=nf_get_vara_int(lncid,lvid,lstart,lncount,lvdat)
vdat=lvdat
#else
ncstatus=nf90_get_var(lncid,lvid,vdat,start=lstart,count=lncount)
#endif
call END_LOG(ncgetv_end)
call ncmsg("get_vara_int2i",ncstatus)

return
end subroutine ccnf_get_vara_int2i

#ifndef i8r8
subroutine ccnf_get_vara_double4d(ncid,vid,start,ncount,vdat)

use cc_mpi

implicit none

integer, intent(in) :: ncid, vid
integer, dimension(:) :: start, ncount
integer ncstatus
integer(kind=4) :: lncid, lvid
integer(kind=4), dimension(size(start)) :: lstart
integer(kind=4), dimension(size(ncount)) :: lncount
real(kind=8), dimension(:,:,:), intent(out) :: vdat

lncid=ncid
lvid=vid
lstart=start
lncount=ncount
call START_LOG(ncgetv_begin)
#ifdef usenc3
ncstatus=nf_get_vara_double(lncid,lvid,lstart,lncount,vdat)
#else
ncstatus=nf90_get_var(lncid,lvid,vdat,start=lstart,count=lncount)
#endif
call END_LOG(ncgetv_end)
call ncmsg("get_vara_double4d",ncstatus)

return
end subroutine ccnf_get_vara_double4d
#endif

subroutine ccnf_get_att_text(ncid,vid,aname,atext,ierr)

use cc_mpi

implicit none

integer, intent(in) :: ncid, vid
integer, intent(out), optional :: ierr
integer ncstatus
integer(kind=4) lncid, lvid
character(len=*), intent(in) :: aname
character(len=*), intent(out) :: atext

atext=''
lncid=ncid
lvid=vid
#ifdef usenc3
ncstatus = nf_get_att_text(lncid,lvid,aname,atext)
#else
ncstatus = nf90_get_att(lncid,lvid,aname,atext)
#endif
if (present(ierr)) then
  ierr=ncstatus
else
  call ncmsg(aname,ncstatus)
end if

return
end subroutine ccnf_get_att_text

subroutine ccnf_get_att_real(ncid,vid,aname,vdat,ierr)

use cc_mpi

implicit none

integer, intent(in) :: ncid, vid
integer, intent(out), optional :: ierr
integer ncstatus
integer(kind=4) lncid, lvid
character(len=*), intent(in) :: aname
real, intent(out) :: vdat

lncid=ncid
lvid=vid
#ifdef usenc3
#ifdef i8r8
ncstatus = nf_get_att_double(lncid,lvid,aname,vdat)
#else
ncstatus = nf_get_att_real(lncid,lvid,aname,vdat)
#endif
#else
ncstatus = nf90_get_att(lncid,lvid,aname,vdat)
#endif
if (present(ierr)) then
  ierr=ncstatus
else
  call ncmsg(aname,ncstatus)
end if


return
end subroutine ccnf_get_att_real

subroutine ccnf_get_att_realg1r(ncid,aname,vdat,ierr)

use cc_mpi

implicit none

integer, intent(in) :: ncid
integer, intent(out), optional :: ierr
integer ncstatus
integer(kind=4) lncid
character(len=*), intent(in) :: aname
real, intent(out) :: vdat

lncid=ncid
#ifdef usenc3
#ifdef i8r8
ncstatus = nf_get_att_double(lncid,nf_global,aname,vdat)
#else
ncstatus = nf_get_att_real(lncid,nf_global,aname,vdat)
#endif
#else
ncstatus = nf90_get_att(lncid,nf90_global,aname,vdat)
#endif
if (present(ierr)) then
  ierr=ncstatus
else
  call ncmsg("get_attg",ncstatus)
end if

return
end subroutine ccnf_get_att_realg1r

subroutine ccnf_get_att_realg2r(ncid,aname,vdat)

use cc_mpi

implicit none

integer, intent(in) :: ncid
integer ncstatus
integer(kind=4) lncid
character(len=*), intent(in) :: aname
real, dimension(:), intent(out) :: vdat
#ifdef i8r8
real(kind=8), dimension(size(vdat)) :: rvals
#else
real(kind=4), dimension(size(vdat)) :: rvals
#endif

lncid=ncid
#ifdef usenc3
#ifdef i8r8
ncstatus = nf_get_att_double(lncid,nf_global,aname,rvals(:))
#else
ncstatus = nf_get_att_real(lncid,nf_global,aname,rvals(:))
#endif
vdat(:) = rvals(:)
#else
ncstatus = nf90_get_att(lncid,nf90_global,aname,vdat(:))
#endif
call ncmsg("get_attg",ncstatus)

return
end subroutine ccnf_get_att_realg2r

subroutine ccnf_get_att_intg1i(ncid,aname,vdat)

use cc_mpi

implicit none

integer, intent(in) :: ncid
integer, intent(out) :: vdat
integer ncstatus
integer(kind=4) :: lncid
#ifdef usenc3
integer(kind=4), dimension(1) :: ivals
#else
integer(kind=4) :: lvdat
#endif
character(len=*), intent(in) :: aname

lncid=ncid
#ifdef usenc3
ncstatus = nf_get_att_int(lncid,nf_global,aname,ivals)
vdat = ivals(1)
#else
ncstatus = nf90_get_att(lncid,nf90_global,aname,lvdat)
vdat=lvdat
#endif
call ncmsg("get_attg",ncstatus)

return
end subroutine ccnf_get_att_intg1i

subroutine ccnf_get_att_intg2i(ncid,aname,vdat)

use cc_mpi

implicit none

integer, intent(in) :: ncid
integer, dimension(:), intent(out) :: vdat
integer ncstatus
integer(kind=4) lncid
integer(kind=4), dimension(size(vdat)) :: lvdat
character(len=*), intent(in) :: aname

lncid=ncid
#ifdef usenc3
ncstatus = nf_get_att_int(lncid,nf_global,aname,lvdat)
#else
ncstatus = nf90_get_att(lncid,nf90_global,aname,lvdat)
#endif
vdat=lvdat
call ncmsg("get_attg",ncstatus)

return
end subroutine ccnf_get_att_intg2i

subroutine ccnf_read(fname,vname,vdat)

use cc_mpi

implicit none

include 'newmpar.h'

integer ncstatus
integer(kind=4) lncid, lvid
real, dimension(ifull), intent(out) :: vdat
real, dimension(ifull_g) :: vdat_g
character(len=*), intent(in) :: fname
character(len=*), intent(in) :: vname

if (myid==0) then
#ifdef usenc3
  ncstatus = nf_open(fname,nf_nowrite,lncid)
  call ncmsg(fname,ncstatus)
  ncstatus = nf_inq_varid(lncid,vname,lvid)
  call ncmsg(fname,ncstatus)
  call START_LOG(ncgetv_begin)
#ifdef i8r8
  ncstatus = nf_get_var_double(lncid,lvid,vdat_g)
#else
  ncstatus = nf_get_var_real(lncid,lvid,vdat_g)
#endif
  call END_LOG(ncgetv_end)
  call ncmsg(fname,ncstatus)
  ncstatus = nf_close(lncid)
  call ncmsg(fname,ncstatus)
#else
  ncstatus = nf90_open(fname,nf90_nowrite,lncid)
  call ncmsg(fname,ncstatus)
  ncstatus = nf90_inq_varid(lncid,vname,lvid)
  call ncmsg(fname,ncstatus)
  call START_LOG(ncgetv_begin)
  ncstatus = nf90_get_var(lncid,lvid,vdat_g)
  call END_LOG(ncgetv_end)
  call ncmsg(fname,ncstatus)
  ncstatus = nf90_close(lncid)
  call ncmsg(fname,ncstatus)
#endif
  call ccmpi_distribute(vdat,vdat_g)
else
  call ccmpi_distribute(vdat)
end if

return
end subroutine ccnf_read

subroutine ccnf_put_var_text2r(ncid,vid,vtxt)

use cc_mpi

implicit none

integer, intent(in) :: ncid, vid
integer ncstatus, i
integer(kind=4) lncid, lvid
integer(kind=4), dimension(1) :: lstart
character(len=*), dimension(:), intent(in) :: vtxt

lncid=ncid
lvid=vid
#ifdef usenc3
do i = 1,size(vtxt)
  lstart = i
  ncstatus = nf_put_var1_text(lncid,lvid,lstart,vtxt(i))
end do
#else
ncstatus = nf90_put_var(lncid,lvid,vtxt)
#endif
call ncmsg("put_var",ncstatus)

return
end subroutine ccnf_put_var_text2r

subroutine ccnf_put_var_int2i(ncid,vid,vdat)

use cc_mpi

implicit none

integer, intent(in) :: ncid,vid
integer ncstatus
integer, dimension(:), intent(in) :: vdat
integer(kind=4) lncid, lvid
integer(kind=4), dimension(size(vdat)) :: lvdat
integer(kind=4), dimension(1) :: lstart, lncount

lncid = ncid
lvid = vid
#ifdef usenc3
lvdat = vdat
lstart = 1
lncount = size(vdat)
ncstatus = nf_put_vara_int(lncid,lvid,lstart,lncount,lvdat)
#else
ncstatus = nf90_put_var(lncid,lvid,vdat)
#endif
call ncmsg("put_var",ncstatus)

return
end subroutine ccnf_put_var_int2i

subroutine ccnf_put_var_int3i(ncid,vid,vdat)

use cc_mpi

implicit none

integer, intent(in) :: ncid, vid
integer ncstatus
integer, dimension(:,:), intent(in) :: vdat
integer(kind=4) lncid, lvid
integer(kind=4), dimension(size(vdat,1),size(vdat,2)) :: lvdat
integer(kind=4), dimension(2) :: lstart, lncount

lncid = ncid
lvid = vid
#ifdef usenc3
lvdat = vdat
lstart(:) = 1
lncount(1) = size(vdat,1)
lncount(2) = size(vdat,2)
ncstatus = nf_put_vara_int(lncid,lvid,lstart,lncount,lvdat)
#else
ncstatus = nf90_put_var(lncid,lvid,vdat)
#endif
call ncmsg("put_var",ncstatus)

return
end subroutine ccnf_put_var_int3i

subroutine ccnf_put_var1_int(ncid,vid,start,vdat)

use cc_mpi

implicit none

integer, intent(in) :: ncid, vid, start
integer, intent(in) :: vdat
integer ncstatus
integer(kind=4) lncid, lvid
#ifdef usenc3
integer(kind=4), dimension(1) :: lstart, ldat
#else
integer(kind=4), dimension(1) :: lstart, lncount
integer, dimension(1) :: ldat
#endif

lncid=ncid
lvid=vid
#ifdef usenc3
lstart=start
ldat=vdat
ncstatus=nf_put_var1_int(lncid,lvid,lstart,ldat(1))
#else
lstart=start
lncount=1
ldat(1)=vdat
ncstatus=nf90_put_var(lncid,lvid,ldat,start=lstart,count=lncount)
#endif
call ncmsg("put_var1",ncstatus)

return
end subroutine ccnf_put_var1_int

subroutine ccnf_put_var1_real(ncid,vid,start,vdat)

use cc_mpi

implicit none

integer, intent(in) :: ncid, vid, start
integer ncstatus
integer(kind=4) lncid, lvid
#ifdef usenc3
integer(kind=4), dimension(1) :: lstart
#else
integer(kind=4), dimension(1) :: lstart, lcount
#endif
real, dimension(1) :: ldat
real, intent(in) :: vdat

lncid = ncid
lvid = vid
ldat = vdat
#ifdef usenc3
lstart = start
#ifdef i8r8
ncstatus = nf_put_var1_double(lncid,lvid,lstart,ldat(1))
#else
ncstatus = nf_put_var1_real(lncid,lvid,lstart,ldat(1))
#endif
#else
lstart = start
lcount = 1
ncstatus = nf90_put_var(lncid,lvid,ldat,start=lstart,count=lcount)
#endif
call ncmsg("put_var1",ncstatus)

return
end subroutine ccnf_put_var1_real

#ifndef i8r8
subroutine ccnf_put_var1_double(ncid,vid,start,vdat)

use cc_mpi

implicit none

integer, intent(in) :: ncid, vid, start
integer ncstatus
integer(kind=4) lncid, lvid
#ifdef usenc3
integer(kind=4), dimension(1) :: lstart
#else
integer(kind=4), dimension(1) :: lstart, lcount
#endif
real(kind=8), dimension(1) :: ldat
real(kind=8), intent(in) :: vdat

lncid = ncid
lvid = vid
ldat = vdat
#ifdef usenc3
lstart = start
ncstatus = nf_put_var1_double(lncid,lvid,lstart,ldat(1))
#else
lstart = start
lcount = 1
ncstatus = nf90_put_var(lncid,lvid,ldat,start=lstart,count=lcount)
#endif
call ncmsg("put_var1",ncstatus)

return
end subroutine ccnf_put_var1_double
#endif

subroutine ccnf_put_vara_real2r(ncid,vid,start,ncount,vdat)

use cc_mpi

implicit none

integer, intent(in) :: ncid, vid
integer ncstatus
integer, dimension(:), intent(in) :: start, ncount
integer(kind=4) lncid, lvid
integer(kind=4), dimension(size(start)) :: lstart
integer(kind=4), dimension(size(ncount)) :: lncount
real, dimension(:), intent(in) :: vdat

lncid=ncid
lvid=vid
lstart=start
lncount=ncount
#ifdef usenc3
#ifdef i8r8
ncstatus=nf_put_vara_double(lncid,lvid,lstart,lncount,vdat)
#else
ncstatus=nf_put_vara_real(lncid,lvid,lstart,lncount,vdat)
#endif
#else
ncstatus=nf90_put_var(lncid,lvid,vdat,start=lstart,count=lncount)
#endif
call ncmsg("put_vara_real2r",ncstatus)

return
end subroutine ccnf_put_vara_real2r

subroutine ccnf_put_vara_real3r(ncid,vid,start,ncount,vdat)

use cc_mpi

implicit none

integer, intent(in) :: ncid, vid
integer ncstatus
integer, dimension(:), intent(in) :: start, ncount
integer(kind=4) :: lncid, lvid
integer(kind=4), dimension(size(start)) :: lstart
integer(kind=4), dimension(size(ncount)) :: lncount
real, dimension(:,:), intent(in) :: vdat

lncid=ncid
lvid=vid
lstart=start
lncount=ncount
#ifdef usenc3
#ifdef i8r8
ncstatus=nf_put_vara_double(lncid,lvid,lstart,lncount,vdat)
#else
ncstatus=nf_put_vara_real(lncid,lvid,lstart,lncount,vdat)
#endif
#else
ncstatus=nf90_put_var(lncid,lvid,vdat,start=lstart,count=lncount)
#endif
call ncmsg("put_vara",ncstatus)

return
end subroutine ccnf_put_vara_real3r

#ifndef i8r8
subroutine ccnf_put_vara_double2r(ncid,vid,start,ncount,vdat)

use cc_mpi

implicit none

integer, intent(in) :: ncid, vid
integer ncstatus
integer, dimension(:), intent(in) :: start, ncount
integer(kind=4) lncid, lvid
integer(kind=4), dimension(size(start)) :: lstart
integer(kind=4), dimension(size(ncount)) :: lncount
real(kind=8), dimension(:), intent(in) :: vdat

lncid=ncid
lvid=vid
lstart=start
lncount=ncount
#ifdef usenc3
ncstatus=nf_put_vara_double(lncid,lvid,lstart,lncount,vdat)
#else
ncstatus=nf90_put_var(lncid,lvid,vdat,start=lstart,count=lncount)
#endif
call ncmsg("put_vara",ncstatus)

return
end subroutine ccnf_put_vara_double2r
#endif

subroutine ccnf_put_vara_int2i(ncid,vid,start,ncount,vdat)

use cc_mpi

implicit none

integer, intent(in) :: ncid, vid
integer ncstatus
integer, dimension(:), intent(in) :: start, ncount
integer, dimension(:), intent(in) :: vdat
integer(kind=4) lncid, lvid
integer(kind=4), dimension(size(start)) :: lstart
integer(kind=4), dimension(size(ncount)) :: lncount
integer(kind=4), dimension(size(vdat)) :: lvdat

lncid=ncid
lvid=vid
lstart=start
lncount=ncount
#ifdef usenc3
lvdat=vdat
ncstatus=nf_put_vara_int(lncid,lvid,lstart,lncount,lvdat)
#else
ncstatus=nf90_put_var(lncid,lvid,vdat,start=lstart,count=lncount)
#endif
call ncmsg("put_vara",ncstatus)

return
end subroutine ccnf_put_vara_int2i

subroutine ccnf_put_att_text(ncid,vid,aname,asize,atext)

use cc_mpi

implicit none

integer, intent(in) :: ncid,vid,asize
integer ncstatus
integer(kind=4) lncid, lvid, lasize
character(len=*), intent(in) :: aname
character(len=asize), intent(in) :: atext

lncid=ncid
lvid=vid
#ifdef usenc3
lasize=asize
ncstatus=nf_put_att_text(lncid,lvid,aname,lasize,atext)
#else
ncstatus=nf90_put_att(lncid,lvid,aname,atext)
#endif
call ncmsg("put_att",ncstatus)

return
end subroutine ccnf_put_att_text

subroutine ccnf_put_att_textg(ncid,aname,atext)

use cc_mpi

implicit none

integer, intent(in) :: ncid
integer ncstatus
integer(kind=4) lsize, lncid
character(len=*), intent(in) :: aname
character(len=*), intent(in) :: atext

lncid=ncid
#ifdef usenc3
lsize=len_trim(atext)
ncstatus=nf_put_att_text(lncid,nf_global,aname,lsize,atext)
#else
ncstatus=nf90_put_att(lncid,nf90_global,aname,atext)
#endif
call ncmsg("put_attg",ncstatus)

return
end subroutine ccnf_put_att_textg

subroutine ccnf_put_att_intg1(ncid,aname,vdat)

use cc_mpi

implicit none

integer, intent(in) :: ncid
integer ncstatus
integer, intent(in) :: vdat
integer(kind=4) lncid
integer(kind=4), dimension(1) :: ldat
character(len=*), intent(in) :: aname

lncid=ncid
ldat(1)=vdat
#ifdef usenc3
ncstatus=nf_put_att_int(lncid,nf_global,aname,nf_int,1_4,ldat)
#else
ncstatus=nf90_put_att(lncid,nf90_global,aname,ldat)
#endif
call ncmsg("put_attg",ncstatus)

return
end subroutine ccnf_put_att_intg1

subroutine ccnf_put_att_intg2(ncid,aname,vsize,vdat)

use cc_mpi

implicit none

integer, intent(in) :: ncid,vsize
integer ncstatus
integer, dimension(:), intent(in) :: vdat
integer(kind=4) lncid, lvsize
integer(kind=4), dimension(size(vdat)) :: lvdat
character(len=*), intent(in) :: aname

lncid=ncid
lvdat=vdat
#ifdef usenc3
lvsize=vsize
ncstatus=nf_put_att_int(lncid,nf_global,aname,nf_int,lvsize,lvdat)
#else
ncstatus=nf90_put_att(lncid,nf90_global,aname,lvdat)
#endif
call ncmsg("put_att",ncstatus)

return
end subroutine ccnf_put_att_intg2

subroutine ccnf_put_att_realg1(ncid,aname,vdat)

use cc_mpi

implicit none

integer, intent(in) :: ncid
integer ncstatus
integer(kind=4) lncid
real, intent(in) :: vdat
real, dimension(1) :: ldat
character(len=*), intent(in) :: aname

lncid=ncid
ldat(1)=vdat
#ifdef usenc3
#ifdef i8r8
ncstatus=nf_put_att_double(lncid,nf_global,aname,nf_float,1_4,ldat)
#else
ncstatus=nf_put_att_real(lncid,nf_global,aname,nf_float,1_4,ldat)
#endif
#else
ncstatus=nf90_put_att(lncid,nf90_global,aname,ldat)
#endif
call ncmsg("put_attg",ncstatus)

return
end subroutine ccnf_put_att_realg1

subroutine ccnf_put_att_realg2(ncid,aname,vsize,vdat)

use cc_mpi

implicit none

integer, intent(in) :: ncid, vsize
integer ncstatus
integer(kind=4) lncid, lvsize
real, dimension(:), intent(in) :: vdat
character(len=*), intent(in) :: aname

lncid=ncid
#ifdef usenc3
lvsize=vsize
#ifdef i8r8
ncstatus=nf_put_att_double(lncid,nf_global,aname,nf_float,lvsize,vdat)
#else
ncstatus=nf_put_att_real(lncid,nf_global,aname,nf_float,lvsize,vdat)
#endif
#else
ncstatus=nf90_put_att(lncid,nf90_global,aname,vdat)
#endif
call ncmsg("put_attg",ncstatus)

return
end subroutine ccnf_put_att_realg2

subroutine surfread(dat,varname,netcdfid,filename)

use cc_mpi

implicit none

include 'newmpar.h'

integer, intent(in), optional :: netcdfid
integer ifull_l
character(len=*), intent(in), optional :: filename
character(len=*), intent(in) :: varname
real, dimension(:), intent(out) :: dat

ifull_l=size(dat)
if (myid==0) then
  if (present(filename)) then
    call surfreadglob(dat,varname,filename=filename)  
  else if (present(netcdfid)) then
    call surfreadglob(dat,varname,netcdfid=netcdfid)  
  else
    write(6,*) 'Failed to specify input file'
    call ccmpi_abort(-1)
  end if
else
  if (ifull_l==ifull) then
    call ccmpi_distribute(dat)
  end if
end if

return
end subroutine surfread

! Read surface data and distribute over processors
! This version suports both netcdf and text file formats
subroutine surfreadglob(dat,vname,netcdfid,filename)

use cc_mpi

implicit none

include 'newmpar.h'   ! Grid parameters
include 'parm.h'      ! Model configuration
include 'parmgeom.h'  ! Coordinate data

integer, intent(in), optional :: netcdfid
integer, dimension(3) :: spos, npos
integer ifull_l, ncidx, iernc, varid, ierr
integer ilx, jlx, ndims
character(len=*), intent(in), optional :: filename
character(len=*), intent(in) :: vname
character(len=47) header
real, dimension(:), intent(out) :: dat
real, dimension(ifull_g) :: glob2d
real rlong0x, rlat0x, schmidtx, dsx

ifull_l=size(dat)

if (present(filename)) then
  call ccnf_open(filename,ncidx,iernc)
else if (present(netcdfid)) then
  ncidx=netcdfid
  iernc=0
end if

if (iernc==0) then ! Netcdf file

  call ccnf_inq_dimlen(ncidx,'longitude',ilx)
  call ccnf_inq_dimlen(ncidx,'latitude',jlx)
  call ccnf_get_attg(ncidx,'lon0',rlong0x)
  call ccnf_get_attg(ncidx,'lat0',rlat0x)
  call ccnf_get_attg(ncidx,'schmidt',schmidtx)
  if(ilx/=il_g.or.jlx/=jl_g.or.rlong0x/=rlong0.or.rlat0x/=rlat0.or.schmidtx/=schmidt) then
    write(6,*) 'wrong data file supplied'
    call ccmpi_abort(-1)
  end if

  spos(1:3)=1
  npos(1)=il_g
  npos(2)=6*il_g
  npos(3)=1
  call ccnf_inq_varid(ncidx,vname,varid)
  call ccnf_inq_varndims(ncidx,varid,ndims)
  call ccnf_get_vara(ncidx,varid,spos(1:ndims),npos(1:ndims),glob2d)
  if (present(filename)) then
    call ccnf_close(ncidx)
  end if
  
else ! ASCII file

  open(87,file=filename,status='old')
  read(87,*,iostat=ierr) ilx,jlx,rlong0x,rlat0x,schmidtx,dsx,header
  if ( ierr == 0 ) then
    write(6,*) ilx,jlx,rlong0x,rlat0x,schmidtx,dsx,header
    if(ilx/=il_g.or.jlx/=jl_g.or.rlong0x/=rlong0.or.rlat0x/=rlat0.or.schmidtx/=schmidt) then
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
if (ifull_l==ifull) then
  call ccmpi_distribute(dat, glob2d)
else if (ifull_l==ifull_g) then
  dat=glob2d
else
  write(6,*) "ERROR: Invalid array size in surfreadglob"
  call ccmpi_abort(-1)
end if
  
return
end subroutine surfreadglob

!--------------------------------------------------------------
! Trap netcdf error messages
subroutine ncmsg(txt,ierr)

use cc_mpi

implicit none

integer, intent(in) :: ierr
integer(kind=4) :: lierr
character(len=*), intent(in) :: txt

#ifdef usenc3
if (ierr/=nf_noerr) then
  lierr=ierr
  write(6,*) "ERROR: Netcdf error = ",ierr
  write(6,*) txt," ",nf_strerror(lierr)
#else
if (ierr/=nf90_noerr) then
  lierr=ierr
  write(6,*) "ERROR: Netcdf error = ",ierr
  write(6,*) txt," ",nf90_strerror(lierr)
#endif
  call ccmpi_abort(-1)
end if

return
end subroutine ncmsg

end module infile
