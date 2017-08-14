! Conformal Cubic Atmospheric Model
    
! Copyright 2015-2017 Commonwealth Scientific Industrial Research Organisation (CSIRO)
    
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
! supplied in parallel and then makes this data avaliable to 
! all processors for interpolation. The code can also identify
! restart files, in which case no additional message passing
! is required.

#ifdef usenc_mod
! use netcdf.mod interface
use netcdf
#else
! use netcdf.inc interface (default) or C interface (-Dncclib)
use netcdf_m
#endif

implicit none
            
private
public vertint, datefix, getzinp, ncmsg
public ptest, pfall, ncidold
public histopen, histclose, histrd3, histrd4, histrd5, surfread
public attrib, histwrt3, histwrt4, histwrt5, freqwrite
public ccnf_open, ccnf_create, ccnf_close, ccnf_sync, ccnf_enddef
public ccnf_redef, ccnf_nofill, ccnf_inq_varid, ccnf_inq_dimid
public ccnf_inq_dimlen, ccnf_inq_varndims, ccnf_def_dim, ccnf_def_dimu
public ccnf_def_var, ccnf_get_vara, ccnf_get_att, ccnf_get_attg
public ccnf_read, ccnf_put_vara, ccnf_put_att, ccnf_put_attg
public file_distribute
public pil_g, pjl_g, pka_g, pko_g, mynproc
public comm_ip

interface histrd3
  module procedure histrd3r4
#ifndef i8r8  
  module procedure histrd3r8
#endif
end interface histrd3
interface histrd4
  module procedure histrd4r4
#ifndef i8r8
  module procedure histrd4r8
#endif
end interface histrd4
interface histrd5
  module procedure histrd5r4
#ifndef i8r8
  module procedure histrd5r8
#endif
end interface histrd5
interface histwrt3
  module procedure histwrt3r4
#ifndef i8r8
  module procedure histwrt3r8
#endif
end interface histwrt3
interface histwrt4
  module procedure histwrt4r4
#ifndef i8r8
  module procedure histwrt4r8
#endif
end interface histwrt4
interface histwrt5
  module procedure histwrt5r4
#ifndef i8r8
  module procedure histwrt5r8
#endif
end interface histwrt5
interface ccnf_get_att
  module procedure ccnf_get_att_text, ccnf_get_att_real
end interface ccnf_get_att
interface ccnf_get_attg
  module procedure ccnf_get_att_intg1i, ccnf_get_att_intg2i
  module procedure ccnf_get_att_realg1r, ccnf_get_att_realg2r
end interface ccnf_get_attg
interface ccnf_get_vara
  module procedure ccnf_get_var_real, ccnf_get_var_int
  module procedure ccnf_get_vara_text_t
  module procedure ccnf_get_vara_real1r_s, ccnf_get_vara_real1r_t
  module procedure ccnf_get_vara_real2r_s, ccnf_get_vara_real2r_t, ccnf_get_vara_real2r
  module procedure ccnf_get_vara_real3r_t, ccnf_get_vara_real3r, ccnf_get_vara_real4r 
  module procedure ccnf_get_vara_int1i_s, ccnf_get_vara_int2i_t, ccnf_get_vara_int2i
  module procedure ccnf_get_vara_int3i_t
#ifndef i8r8
  module procedure ccnf_get_vara_double2r_t, ccnf_get_vara_double4d
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
  module procedure ccnf_put_var_text2r
  module procedure ccnf_put_var_int2i, ccnf_put_var_int3i
  module procedure ccnf_put_vara_real1r_s, ccnf_put_vara_real1r_t
  module procedure ccnf_put_vara_real2r_s, ccnf_put_vara_real2r_t
  module procedure ccnf_put_vara_real3r_t
  module procedure ccnf_put_vara_real2r, ccnf_put_vara_real3r
  module procedure ccnf_put_vara_int1i_s, ccnf_put_vara_int1i_t
  module procedure ccnf_put_vara_int2i
#ifndef i8r8
  module procedure ccnf_put_vara_double1r_s, ccnf_put_vara_double2r
#endif
end interface ccnf_put_vara
interface file_distribute
  module procedure file_distribute2
end interface file_distribute

integer(kind=4), dimension(:), allocatable, save :: pncid
integer, dimension(:), allocatable, save :: pprid
integer, save :: ncidold = -1
integer, save :: mynproc
integer, save :: pil_g, pjl_g, pka_g, pko_g
integer, save :: comm_ip
logical, dimension(:), allocatable, save :: pfown
logical, save :: ptest, pfall, resprocformat

integer(kind=2), parameter :: minv = -32500
integer(kind=2), parameter :: maxv =  32500
integer(kind=2), parameter :: missval = -32501

contains

!--------------------------------------------------------------
! Interface for reading 2D+time fields
subroutine histrd3r4(iarchi,ier,name,ik,var,ifull,nogather)
      
use cc_mpi, only : myid, ccmpi_reduce, histrd3_begin, histrd3_end, fnresid, &
                   start_log, end_log, ccmpi_distribute
use parm_m

implicit none
      
integer, intent(in) :: iarchi,ik,ifull
integer, intent(out) :: ier
real vmax, vmin, vmax_g, vmin_g
real, dimension(:), intent(inout) :: var ! may be dummy argument from myid/=0
logical, intent(in), optional :: nogather
logical ngflag
character(len=*), intent(in) :: name

call START_LOG(histrd3_begin)

ngflag = .false.
if ( present(nogather) ) then
  ngflag = nogather
end if

if ( (ifull/=6*ik*ik.and.ptest) .or. ngflag ) then
  ! read local arrays without gather and distribute (e.g., restart file)
  call hr3p(iarchi,ier,name,.true.,var)
  if ( ier==0 .and. nmaxpr==1 .and. myid<fnresid ) then
    vmax = maxval(var)
    vmin = minval(var) 
    call ccmpi_reduce(vmax,vmax_g,"max",0,comm_ip)
    call ccmpi_reduce(vmin,vmin_g,"min",0,comm_ip)
    if ( myid==0 ) then
      write(6,'(" done histrd3 ",a8,i4,i3,2e14.6)') trim(name),ier,iarchi,vmin,vmax
    end if
  else if ( ier==0 .and. mod(ktau,nmaxpr)==0 .and. myid==0 ) then
    write(6,'(" done histrd3 ",a48,i4,i3)') trim(name),ier,iarchi
  end if

else if ( ifull/=6*ik*ik ) then  
  ! read local arrays with gather and distribute (i.e., change in number of processors) 
  if ( myid==0 ) then
    call hr3a(iarchi,ier,name,ik,var,ifull)
  else
    call hr3p(iarchi,ier,name,.false.)
    call ccmpi_distribute(var)
  end if

else
  ! read global arrays and gather on myid==0
  if ( myid==0 ) then
    call hr3a(iarchi,ier,name,ik,var,ifull)
  else
    call hr3p(iarchi,ier,name,.false.)
  end if

end if

call END_LOG(histrd3_end)

return
end subroutine histrd3r4   

!--------------------------------------------------------------------
! Gather 2D+time fields on myid==0
subroutine hr3a(iarchi,ier,name,ik,var,ifull)
      
use cc_mpi, only : ccmpi_distribute
use parm_m
      
implicit none
      
integer, intent(in) :: iarchi, ik, ifull
integer, intent(out) :: ier
integer iq
real, dimension(:), intent(inout) :: var
real, dimension(:), allocatable :: globvar
real vmax, vmin
character(len=*), intent(in) :: name

allocate( globvar(6*ik*ik) )
globvar(:) = 0.

call hr3p(iarchi,ier,name,.false.,globvar)

if ( ier==0 .and. mod(ktau,nmaxpr)==0 ) then
  vmax = maxval(globvar)
  vmin = minval(globvar)
  iq = id + (jd-1)*ik
  if ( iq<=size(globvar) ) then
    write(6,'(" done histrd3 ",a8,i4,i3,3e14.6)') trim(name),ier,iarchi,vmin,vmax,globvar(iq)
  else
    write(6,'(" done histrd3 ",a20,i4,i3,2e14.6)') trim(name),ier,iarchi,vmin,vmax
  end if
end if

if ( ifull==6*ik*ik ) then
  ! read global arrays for myid==0
  var(1:6*ik*ik) = globvar(:) ! really ifull_g
else
  ! read local arrays with gather and distribute (no longer used)
  call ccmpi_distribute(var,globvar)
end if

deallocate( globvar )

return
end subroutine hr3a  

!--------------------------------------------------------------
! This subroutine reads 2D+time input files
      
! when qtest=.true. the input grid decomposition should
! match the current processor decomposition.  We can then
! skip the MPI gather and distribute steps.
subroutine hr3p(iarchi,ier,name,qtest,var)

implicit none

integer, intent(in) :: iarchi
integer, intent(out) :: ier
real, dimension(:), intent(inout), optional :: var
logical, intent(in) :: qtest
character(len=*), intent(in) :: name

if ( resprocformat .and. present(var) ) then
  call hr3p_procformat(iarchi,ier,name,qtest,var)
else if ( resprocformat ) then
  call hr3p_procformat(iarchi,ier,name,qtest)  
else if ( present(var) ) then
  call hr3p_para(iarchi,ier,name,qtest,var)
else
  call hr3p_para(iarchi,ier,name,qtest)  
end if

return
end subroutine hr3p

subroutine hr3p_procformat(iarchi,ier,name,qtest,var)

use cc_mpi
use newmpar_m
      
implicit none

integer, intent(in) :: iarchi
integer, intent(out) :: ier
integer(kind=4), dimension(4) :: start, ncount
integer ipf, ca
integer(kind=4) idv, ndims
real, dimension(:), intent(inout), optional :: var
real, dimension(pil*pjl*pnpan) :: rvar
real(kind=4) laddoff, lsf
logical, intent(in) :: qtest
character(len=*), intent(in) :: name

ier = 0
      
do ipf = 0,mynproc-1

  start  = (/ 1, 1, pprid(ipf), iarchi /)
  ncount = (/ pil, pjl*pnpan, 1, 1 /)
    
  rvar(:) = 0. ! default for missing field
  
  ! get variable idv
  ier=nf90_inq_varid(pncid(ipf),name,idv)
  if ( ier/=nf90_noerr ) then
    if ( myid==0 .and. ipf==0 ) then
      write(6,*) '***absent field for ncid,name,ier: ',pncid(0),name,ier
    end if
  else
    ! obtain scaling factors and offsets from attributes
    ier=nf90_get_att(pncid(ipf),idv,'add_offset',laddoff)
    if (ier/=nf90_noerr) laddoff=0.
    ier=nf90_get_att(pncid(ipf),idv,'scale_factor',lsf)
    if (ier/=nf90_noerr) lsf=1.
    ier=nf90_inquire_variable(pncid(ipf),idv,ndims=ndims)
    ier=nf90_get_var(pncid(ipf),idv,rvar,start=start(1:ndims),count=ncount(1:ndims))
    call ncmsg(name,ier)
    ! unpack compressed data
    rvar(:) = rvar(:)*real(lsf) + real(laddoff)
  end if ! ier
      
  if ( qtest ) then
    ! e.g., restart file or nogather=.true.
    ca = pil*pjl*pnpan*ipf
    var(1+ca:pil*pjl*pnpan+ca) = rvar(:)
  else
    ! e.g., mesonest file or nogather=.false.
    if ( myid==0 ) then
      call host_hr3p(ipf,rvar,var)
    else
      call proc_hr3p(rvar)
    end if
  end if ! qtest

end do ! ipf

return
end subroutine hr3p_procformat

subroutine hr3p_para(iarchi,ier,name,qtest,var)
      
use cc_mpi
use newmpar_m
      
implicit none

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
      
do ipf = 0,mynproc-1

  rvar(:) = 0. ! default for missing field
  
  ! get variable idv
  ier = nf90_inq_varid(pncid(ipf),name,idv)
  if ( ier/=nf90_noerr ) then
    if ( myid==0 .and. ipf==0 ) then
      write(6,*) '***absent field for ncid,name,ier: ',pncid(0),name,ier
    end if
  else
    ! obtain scaling factors and offsets from attributes
    ier = nf90_get_att(pncid(ipf),idv,'add_offset',laddoff)
    if ( ier/=nf90_noerr ) laddoff = 0.
    ier = nf90_get_att(pncid(ipf),idv,'scale_factor',lsf)
    if ( ier/=nf90_noerr ) lsf = 1.
    ier = nf90_inquire_variable(pncid(ipf),idv,ndims=ndims)
    ier = nf90_get_var(pncid(ipf),idv,rvar,start=start(1:ndims),count=ncount(1:ndims))
    call ncmsg(name,ier)
    ! unpack compressed data
    rvar(:) = rvar(:)*real(lsf)+real(laddoff)
  end if ! ier
      
  if ( qtest ) then
    ! e.g., restart file or nogather=.true.
    ca = pil*pjl*pnpan*ipf
    var(1+ca:pil*pjl*pnpan+ca) = rvar(:)
  else
    ! e.g., mesonest file or nogather=.false.
    if ( myid==0 ) then
      call host_hr3p(ipf,rvar,var)
    else
      call proc_hr3p(rvar)
    end if
  end if ! qtest

end do ! ipf

return
end subroutine hr3p_para

subroutine host_hr3p(ipf,rvar,var)

use cc_mpi
use newmpar_m

implicit none

integer, intent(in) :: ipf
integer jpf, ip, n, no, ca, cc, j
real, dimension(:), intent(inout) :: var
real, dimension(pil*pjl*pnpan), intent(in) :: rvar
real, dimension(:,:), allocatable :: globvar 

if ( fnresid==1 ) then
    
  ip = ipf*fnresid
  do n = 0,pnpan-1
    no = n - pnoff(ip) + 1
    ca = pioff(ip,no) + (pjoff(ip,no)-1)*pil_g + no*pil_g*pil_g
    cc = n*pil*pjl - pil
    do j = 1,pjl
      var(1+j*pil_g+ca:pil+j*pil_g+ca) = rvar(1+j*pil+cc:pil+j*pil+cc)
    end do
  end do
  
else

  allocate( globvar(pil*pjl*pnpan,fnresid) )
  call ccmpi_gatherx(globvar,rvar,0,comm_ip)
  do jpf = 1,fnresid
    ip = ipf*fnresid + jpf - 1
    do n = 0,pnpan-1
      no = n - pnoff(ip) + 1
      ca = pioff(ip,no) + (pjoff(ip,no)-1)*pil_g + no*pil_g*pil_g
      cc = n*pil*pjl - pil
      do j = 1,pjl
        var(1+j*pil_g+ca:pil+j*pil_g+ca) = globvar(1+j*pil+cc:pil+j*pil+cc,jpf)
      end do
    end do
  end do
  deallocate( globvar )
  
end if

return
end subroutine host_hr3p

subroutine proc_hr3p(rvar)

use cc_mpi

implicit none

real, dimension(pil*pjl*pnpan), intent(in) :: rvar
real, dimension(0,0) :: globvar 

call ccmpi_gatherx(globvar,rvar,0,comm_ip)

return
end subroutine proc_hr3p

!--------------------------------------------------------------
! Interface for reading 2D+time fields (double precision version)
subroutine histrd3r8(iarchi,ier,name,ik,var,ifull,nogather)
      
use cc_mpi, only : myid, ccmpi_reducer8, histrd3_begin, histrd3_end, fnresid, &
                   start_log, end_log, ccmpi_distributer8
use parm_m

implicit none
      
integer, intent(in) :: iarchi,ik,ifull
integer, intent(out) :: ier
real(kind=8) vmax, vmin, vmax_g, vmin_g
real(kind=8), dimension(:), intent(inout) :: var ! may be dummy argument from myid/=0
logical, intent(in), optional :: nogather
logical ngflag
character(len=*), intent(in) :: name

call START_LOG(histrd3_begin)

ngflag = .false.
if ( present(nogather) ) then
  ngflag = nogather
end if

if ( (ifull/=6*ik*ik.and.ptest) .or. ngflag ) then
  ! read local arrays without gather and distribute (e.g., restart file)
  call hr3pr8(iarchi,ier,name,.true.,var)
  if ( ier==0 .and. nmaxpr==1 .and. myid<fnresid ) then
    vmax = maxval(var)
    vmin = minval(var) 
    call ccmpi_reducer8(vmax,vmax_g,"max",0,comm_ip)
    call ccmpi_reducer8(vmin,vmin_g,"min",0,comm_ip)
    if ( myid==0 ) then
      write(6,'(" done histrd3r8 ",a8,i4,i3,2e14.6)') trim(name),ier,iarchi,vmin,vmax
    end if
  else if ( ier==0 .and. mod(ktau,nmaxpr)==0 .and. myid==0 ) then
    write(6,'(" done histrd3r8 ",a46,i4,i3)') trim(name),ier,iarchi
  end if

else if ( ifull/=6*ik*ik ) then  
  ! read local arrays with gather and distribute (i.e., change in number of processors) 
  if ( myid==0 ) then
    call hr3ar8(iarchi,ier,name,ik,var,ifull)
  else
    call hr3pr8(iarchi,ier,name,.false.)
    call ccmpi_distributer8(var)
  end if

else
  ! read global arrays and gather on myid==0
  if ( myid==0 ) then
    call hr3ar8(iarchi,ier,name,ik,var,ifull)
  else
    call hr3pr8(iarchi,ier,name,.false.)
  end if

end if

call END_LOG(histrd3_end)

return
end subroutine histrd3r8   

!--------------------------------------------------------------------
! Gather 2D+time fields on myid==0
subroutine hr3ar8(iarchi,ier,name,ik,var,ifull)
      
use cc_mpi, only : ccmpi_distributer8
use parm_m
      
implicit none
      
integer, intent(in) :: iarchi, ik, ifull
integer, intent(out) :: ier
integer iq
character(len=*), intent(in) :: name
real(kind=8), dimension(:), intent(inout) :: var
real(kind=8), dimension(:), allocatable :: globvar
real(kind=8) vmax, vmin

allocate( globvar(6*ik*ik) )
globvar(:) = 0.

call hr3pr8(iarchi,ier,name,.false.,globvar)

if ( ier==0 .and. mod(ktau,nmaxpr)==0 ) then
  vmax = maxval(globvar)
  vmin = minval(globvar)
  iq = id+(jd-1)*ik
  if ( iq<=size(globvar) ) then
    write(6,'(" done histrd3r8 ",a8,i4,i3,3e14.6)') trim(name),ier,iarchi,vmin,vmax,globvar(iq)
  else
    write(6,'(" done histrd3r8 ",a18,i4,i3,2e14.6)') trim(name),ier,iarchi,vmin,vmax
  end if
end if

if ( ifull==6*ik*ik ) then
  ! read global arrays for myid==0
  var(1:6*ik*ik)=globvar(:) ! really ifull_g
else
  ! read local arrays with gather and distribute (i.e., change in number of processors)
  call ccmpi_distributer8(var,globvar)
endif

deallocate( globvar )

return
end subroutine hr3ar8 

!--------------------------------------------------------------
! This subroutine reads 2D+time input files (double precision version)
      
! when qtest=.true. the input grid decomposition should
! match the current processor decomposition.  We can then
! skip the MPI gather and distribute steps.
subroutine hr3pr8(iarchi,ier,name,qtest,var)

implicit none

integer, intent(in) :: iarchi
integer, intent(out) :: ier
logical, intent(in) :: qtest
real(kind=8), dimension(:), intent(inout), optional :: var
character(len=*), intent(in) :: name

if ( resprocformat .and. present(var) ) then
  call hr3p_procformatr8(iarchi,ier,name,qtest,var)
else if ( resprocformat ) then
  call hr3p_procformatr8(iarchi,ier,name,qtest)  
else if ( present(var) ) then
  call hr3p_parar8(iarchi,ier,name,qtest,var)
else
  call hr3p_parar8(iarchi,ier,name,qtest)  
end if

return
end subroutine hr3pr8

subroutine hr3p_procformatr8(iarchi,ier,name,qtest,var)

use cc_mpi
use newmpar_m
      
implicit none

integer, intent(in) :: iarchi
integer, intent(out) :: ier
integer(kind=4), dimension(4) :: start, ncount
integer ipf, ca
integer(kind=4) idv, ndims
real(kind=8), dimension(:), intent(inout), optional :: var
real(kind=8), dimension(pil*pjl*pnpan) :: rvar
real(kind=4) laddoff, lsf
logical, intent(in) :: qtest
character(len=*), intent(in) :: name

ier = 0
      
do ipf = 0,mynproc-1

  start  = (/ 1, 1, pprid(ipf), iarchi /)
  ncount = (/ pil, pjl*pnpan, 1, 1 /)
    
  rvar(:) = 0._8 ! default for missing field
  
  ! get variable idv
  ier=nf90_inq_varid(pncid(ipf),name,idv)
  if ( ier/=nf90_noerr ) then
    if ( myid==0 .and. ipf==0 ) then
      write(6,*) '***absent field for ncid,name,ier: ',pncid(0),name,ier
    end if
  else
    ! obtain scaling factors and offsets from attributes
    ier=nf90_get_att(pncid(ipf),idv,'add_offset',laddoff)
    if (ier/=nf90_noerr) laddoff=0.
    ier=nf90_get_att(pncid(ipf),idv,'scale_factor',lsf)
    if (ier/=nf90_noerr) lsf=1.
    ier=nf90_inquire_variable(pncid(ipf),idv,ndims=ndims)
    call ncmsg(name,ier)
    ier=nf90_get_var(pncid(ipf),idv,rvar,start=start(1:ndims),count=ncount(1:ndims))
    call ncmsg(name,ier)
    ! unpack compressed data
    rvar(:) = rvar(:)*real(lsf,8) + real(laddoff,8)
  end if ! ier
      
  if ( qtest ) then
    ! e.g., restart file or nogather=.true.
    ca = pil*pjl*pnpan*ipf
    var(1+ca:pil*pjl*pnpan+ca) = rvar(:)
  else
    ! e.g., mesonest file or nogather=.false.
    if ( myid==0 ) then
      call host_hr3pr8(ipf,rvar,var)
    else
      call proc_hr3pr8(rvar)
    end if
  end if ! qtest

end do ! ipf

return
end subroutine hr3p_procformatr8

subroutine hr3p_parar8(iarchi,ier,name,qtest,var)
      
use cc_mpi
use newmpar_m
      
implicit none

integer, intent(in) :: iarchi
integer, intent(out) :: ier
integer(kind=4), dimension(3) :: start, ncount
integer ipf, ca
integer(kind=4) idv, ndims
real(kind=8), dimension(:), intent(inout), optional :: var
real(kind=8), dimension(pil*pjl*pnpan) :: rvar
real(kind=4) laddoff, lsf
logical, intent(in) :: qtest
character(len=*), intent(in) :: name

start  = (/ 1, 1, iarchi /)
ncount = (/ pil, pjl*pnpan, 1 /)
ier = 0
      
do ipf = 0,mynproc-1

  rvar(:)=0._8 ! default for missing field
  
  ! get variable idv
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
    ier=nf90_get_var(pncid(ipf),idv,rvar,start=start(1:ndims),count=ncount(1:ndims))
    call ncmsg(name,ier)
    ! unpack compressed data
    rvar(:)=rvar(:)*real(lsf,8)+real(laddoff,8)
  end if ! ier
      
  if (qtest) then
    ! e.g., restart file or nogather=.true.
    ca = pil*pjl*pnpan*ipf
    var(1+ca:pil*pjl*pnpan+ca)=rvar(:)
  else
    ! e.g., mesonest file or nogather=.false.
    if ( myid==0 ) then
      call host_hr3pr8(ipf,rvar,var)
    else
      call proc_hr3pr8(rvar)
    end if
  end if ! qtest

end do ! ipf

return
end subroutine hr3p_parar8

subroutine host_hr3pr8(ipf,rvar,var)

use cc_mpi
use newmpar_m

implicit none


integer, intent(in) :: ipf
integer jpf, ip, n, no, ca, cc, j
real(kind=8), dimension(:), intent(inout) :: var
real(kind=8), dimension(pil*pjl*pnpan), intent(in) :: rvar
real(kind=8), dimension(:,:), allocatable :: globvar 

if ( fnresid==1 ) then
    
  ip = ipf*fnresid
  do n = 0,pnpan-1
    no = n - pnoff(ip) + 1
    ca = pioff(ip,no) + (pjoff(ip,no)-1)*pil_g + no*pil_g*pil_g
    cc = n*pil*pjl - pil
    do j = 1,pjl
      var(1+j*pil_g+ca:pil+j*pil_g+ca) = rvar(1+j*pil+cc:pil+j*pil+cc)
    end do
  end do
  
else

  allocate( globvar(pil*pjl*pnpan,fnresid) )
  call ccmpi_gatherxr8(globvar,rvar,0,comm_ip)
  do jpf = 1,fnresid
    ip = ipf*fnresid + jpf - 1
    do n = 0,pnpan-1
      no = n - pnoff(ip) + 1
      ca = pioff(ip,no) + (pjoff(ip,no)-1)*pil_g + no*pil_g*pil_g
      cc = n*pil*pjl - pil
      do j = 1,pjl
        var(1+j*pil_g+ca:pil+j*pil_g+ca) = globvar(1+j*pil+cc:pil+j*pil+cc,jpf)
      end do
    end do
  end do
  deallocate( globvar )
  
end if

return
end subroutine host_hr3pr8

subroutine proc_hr3pr8(rvar)

use cc_mpi

implicit none

real(kind=8), dimension(pil*pjl*pnpan), intent(in) :: rvar
real(kind=8), dimension(0,0) :: globvar 

call ccmpi_gatherxr8(globvar,rvar,0,comm_ip)

return
end subroutine proc_hr3pr8


!--------------------------------------------------------------   
! Interface for reading 3D+time fields
subroutine histrd4r4(iarchi,ier,name,ik,kk,var,ifull,nogather)
      
use cc_mpi, only : myid, ccmpi_reduce, histrd4_begin, histrd4_end, fnresid, &
                   start_log, end_log, ccmpi_distribute
use parm_m
      
implicit none
      
integer, intent(in) :: iarchi, ik, kk, ifull
integer, intent(out) :: ier
real, dimension(:,:), intent(inout) :: var ! may be dummy argument from myid/=0
real vmax, vmin, vmax_g, vmin_g
logical, intent(in), optional :: nogather
logical ngflag
character(len=*), intent(in) :: name

call START_LOG(histrd4_begin)

ngflag = .false.
if ( present(nogather) ) then
  ngflag = nogather
end if

if ( (ifull/=6*ik*ik.and.ptest) .or. ngflag ) then
  ! read local arrays without gather and distribute
  call hr4p(iarchi,ier,name,kk,.true.,var)
  if ( ier==0 .and. nmaxpr==1 .and. myid<fnresid ) then
    vmax = maxval(var)
    vmin = minval(var) 
    call ccmpi_reduce(vmax,vmax_g,"max",0,comm_ip)
    call ccmpi_reduce(vmin,vmin_g,"min",0,comm_ip)
    if ( myid==0 ) then
      write(6,'(" done histrd4 ",a6,i3,i4,i3,2f12.4)') trim(name),kk,ier,iarchi,vmin,vmax
    end if
  else if ( ier==0 .and. mod(ktau,nmaxpr)==0 .and. myid==0 ) then  
    write(6,'(" done histrd4 ",a48,i4,i3)') trim(name),ier,iarchi
  end if

else if ( ifull/=6*ik*ik ) then  
  ! read local arrays with gather and distribute
  if ( myid==0 ) then
    call hr4sa(iarchi,ier,name,ik,kk,var,ifull)
  else 
    call hr4p(iarchi,ier,name,kk,.false.)
    call ccmpi_distribute(var)
  end if
  
else
  ! read global arrays and gather on myid==0
  if ( myid==0 ) then
    call hr4sa(iarchi,ier,name,ik,kk,var,ifull)
  else 
    call hr4p(iarchi,ier,name,kk,.false.)
  end if

end if

call END_LOG(histrd4_end)

return
end subroutine histrd4r4

!--------------------------------------------------------------------
! Gather 3D+time fields on myid==0

subroutine hr4sa(iarchi,ier,name,ik,kk,var,ifull)
      
use cc_mpi, only : ccmpi_distribute
use parm_m
      
implicit none
      
integer, intent(in) :: iarchi, ik, kk, ifull
integer, intent(out) :: ier
integer iq
character(len=*), intent(in) :: name
real, dimension(:,:) :: var
real, dimension(:,:), allocatable :: globvar
real vmax, vmin
      
allocate( globvar(6*ik*ik,kk) )
globvar(:,:) = 0.

call hr4p(iarchi,ier,name,kk,.false.,globvar)     

if( ier==0 .and. mod(ktau,nmaxpr)==0 ) then
  vmax = maxval(globvar)
  vmin = minval(globvar)
  iq = id+(jd-1)*ik
  if ( iq<=size(globvar,1) .and. nlv<=size(globvar,2) ) then
    write(6,'(" done histrd4 ",a6,i3,i4,i3,3f12.4)') trim(name),kk,ier,iarchi,vmin,vmax,globvar(id+(jd-1)*ik,nlv)
  else
    write(6,'(" done histrd4 ",a18,i3,i4,i3,2f12.4)') trim(name),kk,ier,iarchi,vmin,vmax
  end if
end if

! Have to return correct value of ier on all processes because it's 
! used for initialisation in calling routine
if ( ifull==6*ik*ik ) then
  ! read global arrays for myid==0
  var(1:6*ik*ik,1:kk) = globvar(:,:)
else
  ! read local arrays with gather and distribute (i.e., change in number of processors)
  call ccmpi_distribute(var,globvar)
endif

deallocate( globvar )

return
end subroutine hr4sa      

!--------------------------------------------------------------
! This subroutine reads 3D+time input files

! when qtest=.true. the input grid decomposition should
! match the current processor decomposition.  We can then
! skip the MPI gather and distribute steps.
subroutine hr4p(iarchi,ier,name,kk,qtest,var)

implicit none

integer, intent(in) :: iarchi, kk
integer, intent(out) :: ier
logical, intent(in) :: qtest
real, dimension(:,:), intent(inout), optional :: var
character(len=*), intent(in) :: name

if ( resprocformat .and. present(var) ) then
  call hr4p_procformat(iarchi,ier,name,kk,qtest,var)
else if ( resprocformat ) then
  call hr4p_procformat(iarchi,ier,name,kk,qtest)  
else if ( present(var) ) then
  call hr4p_para(iarchi,ier,name,kk,qtest,var)
else
  call hr4p_para(iarchi,ier,name,kk,qtest)  
end if

return
end subroutine hr4p

subroutine hr4p_procformat(iarchi,ier,name,kk,qtest,var)


use cc_mpi
use newmpar_m
      
implicit none

integer, intent(in) :: iarchi, kk
integer, intent(out) :: ier
integer(kind=4), dimension(5) :: start, ncount
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
  ier = nf90_inq_varid(pncid(ipf),name,idv)
  if ( ier==nf90_noerr ) then
    start(:)  = (/ 1, 1, 1, pprid(ipf), iarchi /)
    ncount(:) = (/ pil, pjl*pnpan, kk, 1, 1 /)   
    ! obtain scaling factors and offsets from attributes
    ier = nf90_get_att(pncid(ipf),idv,'add_offset',laddoff)
    if ( ier/=nf90_noerr ) laddoff = 0.
    ier = nf90_get_att(pncid(ipf),idv,'scale_factor',lsf)
    if ( ier/=nf90_noerr ) lsf = 1.
    ier = nf90_inquire_variable(pncid(ipf),idv,ndims=ndims)
    ier = nf90_get_var(pncid(ipf),idv,rvar,start=start(1:ndims),count=ncount(1:ndims))
    call ncmsg(name,ier)
    ! unpack data
    rvar(:,:) = rvar(:,:)*real(lsf) + real(laddoff)
  else
    start(1:4) = (/ 1, 1, pprid(ipf), iarchi /)
    ncount(1:4) = (/ pil, pjl*pnpan, 1, 1 /)
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
        if ( myid==0 .and. ipf==0 ) then
          write(6,*) '***absent field for ncid,name,ier: ',pncid(0),name,ier
        end if
        rvar(:,:) = 0. ! default value for missing field
        exit
      end if
      ! obtain scaling factors and offsets from attributes
      ier = nf90_get_att(pncid(ipf),idv,'add_offset',laddoff)
      if ( ier/=nf90_noerr ) laddoff = 0.
      ier = nf90_get_att(pncid(ipf),idv,'scale_factor',lsf)
      if ( ier/=nf90_noerr ) lsf = 1.
      ier = nf90_inquire_variable(pncid(ipf),idv,ndims=ndims)
      ier = nf90_get_var(pncid(ipf),idv,rvar(:,k),start=start(1:ndims),count=ncount(1:ndims))
      call ncmsg(name,ier)
      ! unpack data
      rvar(:,k) = rvar(:,k)*real(lsf) + real(laddoff)      
    end do
  end if ! ier

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
end subroutine hr4p_procformat

subroutine hr4p_para(iarchi,ier,name,kk,qtest,var)


use cc_mpi
use newmpar_m
      
implicit none

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
  ier = nf90_inq_varid(pncid(ipf),name,idv)
  if ( ier==nf90_noerr ) then
    start  = (/ 1, 1, 1, iarchi /)
    ncount = (/ pil, pjl*pnpan, kk, 1 /)   
    ! obtain scaling factors and offsets from attributes
    ier = nf90_get_att(pncid(ipf),idv,'add_offset',laddoff)
    if ( ier/=nf90_noerr ) laddoff = 0.
    ier = nf90_get_att(pncid(ipf),idv,'scale_factor',lsf)
    if ( ier/=nf90_noerr ) lsf = 1.
    ier = nf90_inquire_variable(pncid(ipf),idv,ndims=ndims)
    ier = nf90_get_var(pncid(ipf),idv,rvar,start=start(1:ndims),count=ncount(1:ndims))
    call ncmsg(name,ier)
    ! unpack data
    rvar(:,:) = rvar(:,:)*real(lsf) + real(laddoff)
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
      if ( ier/=nf90_noerr ) then
        if ( myid==0 .and. ipf==0 ) then
          write(6,*) '***absent field for ncid,name,ier: ',pncid(0),name,ier
        end if
        rvar(:,:) = 0. ! default value for missing field
        exit
      end if
      ! obtain scaling factors and offsets from attributes
      ier = nf90_get_att(pncid(ipf),idv,'add_offset',laddoff)
      if ( ier/=nf90_noerr ) laddoff = 0.
      ier = nf90_get_att(pncid(ipf),idv,'scale_factor',lsf)
      if ( ier/=nf90_noerr ) lsf = 1.
      ier = nf90_inquire_variable(pncid(ipf),idv,ndims=ndims)
      ier = nf90_get_var(pncid(ipf),idv,rvar(:,k),start=start(1:ndims),count=ncount(1:ndims))
      call ncmsg(name,ier)
      ! unpack data
      rvar(:,k) = rvar(:,k)*real(lsf) + real(laddoff)      
    end do
  end if ! ier

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
end subroutine hr4p_para


subroutine host_hr4p(ipf,kk,rvar,var)

use cc_mpi
use newmpar_m

implicit none

integer, intent(in) :: ipf, kk
integer jpf, ip, n, no, ca, cc, j, k
real, dimension(:,:), intent(inout) :: var
real, dimension(pil*pjl*pnpan,kk), intent(in) :: rvar
real, dimension(:,:,:), allocatable :: globvar 

if ( fnresid==1 ) then

  ip = ipf*fnresid             ! local file number
  do k = 1,kk
    do n = 0,pnpan-1
      no = n - pnoff(ip) + 1   ! global panel number of local file
      ca = pioff(ip,no) + pjoff(ip,no)*pil_g + no*pil_g*pil_g - pil_g
      cc = n*pil*pjl - pil
      do j = 1,pjl
        var(1+j*pil_g+ca:pil+j*pil_g+ca,k) = rvar(1+j*pil+cc:pil+j*pil+cc,k)
      end do
    end do
  end do
  
else

  allocate( globvar(pil*pjl*pnpan,kk,fnresid) )
  call ccmpi_gatherx(globvar,rvar,0,comm_ip)
  do jpf = 1,fnresid
    ip = ipf*fnresid + jpf - 1   ! local file number
    do k = 1,kk
      do n = 0,pnpan-1
        no = n - pnoff(ip) + 1   ! global panel number of local file
        ca = pioff(ip,no) + pjoff(ip,no)*pil_g + no*pil_g*pil_g - pil_g
        cc = n*pil*pjl - pil
        do j = 1,pjl
          var(1+j*pil_g+ca:pil+j*pil_g+ca,k) = globvar(1+j*pil+cc:pil+j*pil+cc,k,jpf)
        end do
      end do
    end do
  end do
  deallocate( globvar )

end if
  
return
end subroutine host_hr4p

subroutine proc_hr4p(kk,rvar)

use cc_mpi

implicit none

integer, intent(in) :: kk
real, dimension(pil*pjl*pnpan,kk), intent(in) :: rvar
real, dimension(0,0,0) :: globvar 

call ccmpi_gatherx(globvar,rvar,0,comm_ip)

return
end subroutine proc_hr4p

!--------------------------------------------------------------   
! Interface for reading 3D+time fields (double precision version)
subroutine histrd4r8(iarchi,ier,name,ik,kk,var,ifull,nogather)
      
use cc_mpi, only : myid, ccmpi_reducer8, histrd4_begin, histrd4_end, fnresid, &
                   start_log, end_log, ccmpi_distributer8
use parm_m
      
implicit none
      
integer, intent(in) :: iarchi, ik, kk, ifull
integer, intent(out) :: ier
real(kind=8), dimension(:,:), intent(inout) :: var ! may be dummy argument from myid/=0
real(kind=8) vmax, vmin, vmax_g, vmin_g
logical, intent(in), optional :: nogather
logical ngflag
character(len=*), intent(in) :: name

call START_LOG(histrd4_begin)

ngflag = .false.
if ( present(nogather) ) then
  ngflag = nogather
end if

if ( (ifull/=6*ik*ik.and.ptest) .or. ngflag ) then
  ! read local arrays without gather and distribute
  call hr4pr8(iarchi,ier,name,kk,.true.,var)
  if ( ier==0 .and. nmaxpr==1 .and. myid<fnresid ) then
    vmax = maxval(var)
    vmin = minval(var) 
    call ccmpi_reducer8(vmax,vmax_g,"max",0,comm_ip)
    call ccmpi_reducer8(vmin,vmin_g,"min",0,comm_ip)
    if ( myid==0 ) then
      write(6,'(" done histrd4r8 ",a6,i3,i4,i3,2f12.4)') trim(name),kk,ier,iarchi,vmin,vmax
    end if
  else if ( ier==0 .and. mod(ktau,nmaxpr)==0 .and. myid==0 ) then  
    write(6,'(" done histrd4r8 ",a46,i4,i3)') trim(name),ier,iarchi
  end if

else if ( ifull/=6*ik*ik ) then  
  ! read local arrays with gather and distribute
  if ( myid==0 ) then
    call hr4sar8(iarchi,ier,name,ik,kk,var,ifull)
  else 
    call hr4pr8(iarchi,ier,name,kk,.false.)
    call ccmpi_distributer8(var)
  end if
  
else
  ! read global arrays and gather on myid==0
  if ( myid==0 ) then
    call hr4sar8(iarchi,ier,name,ik,kk,var,ifull)
  else 
    call hr4pr8(iarchi,ier,name,kk,.false.)
  end if

end if

call END_LOG(histrd4_end)

return
end subroutine histrd4r8

!--------------------------------------------------------------------
! Gather 3D+time fields on myid==0

subroutine hr4sar8(iarchi,ier,name,ik,kk,var,ifull)
      
use cc_mpi, only : ccmpi_distributer8
use parm_m
      
implicit none
      
integer, intent(in) :: iarchi, ik, kk, ifull
integer, intent(out) :: ier
integer iq
character(len=*), intent(in) :: name
real(kind=8), dimension(:,:) :: var
real(kind=8), dimension(:,:), allocatable :: globvar
real(kind=8) vmax, vmin
      
allocate( globvar(6*ik*ik,kk) )
globvar(:,:) = 0._8

call hr4pr8(iarchi,ier,name,kk,.false.,globvar)     

if( ier==0 .and. mod(ktau,nmaxpr)==0 ) then
  vmax = maxval(globvar)
  vmin = minval(globvar)
  iq = id+(jd-1)*ik
  if ( iq<=size(globvar,1) .and. nlv<=size(globvar,2) ) then
    write(6,'(" done histrd4r8 ",a6,i3,i4,i3,3f12.4)') trim(name),kk,ier,iarchi,vmin,vmax,globvar(id+(jd-1)*ik,nlv)
  else
    write(6,'(" done histrd4r8 ",a16,i3,i4,i3,2f12.4)') trim(name),kk,ier,iarchi,vmin,vmax
  end if
end if

! Have to return correct value of ier on all processes because it's 
! used for initialisation in calling routine
if ( ifull==6*ik*ik ) then
  ! read global arrays for myid==0
  var(1:6*ik*ik,1:kk) = globvar(:,:)
else
  ! read local arrays with gather and distribute (i.e., change in number of processors)
  call ccmpi_distributer8(var,globvar)
endif

deallocate( globvar )

return
end subroutine hr4sar8      

!--------------------------------------------------------------
! This subroutine reads 3D+time input files (Double precision version)

! when qtest=.true. the input grid decomposition should
! match the current processor decomposition.  We can then
! skip the MPI gather and distribute steps.
subroutine hr4pr8(iarchi,ier,name,kk,qtest,var)

implicit none

integer, intent(in) :: iarchi, kk
integer, intent(out) :: ier
logical, intent(in) :: qtest
real(kind=8), dimension(:,:), intent(inout), optional :: var
character(len=*), intent(in) :: name

if ( resprocformat .and. present(var) ) then
  call hr4p_procformatr8(iarchi,ier,name,kk,qtest,var)
else if ( resprocformat ) then
  call hr4p_procformatr8(iarchi,ier,name,kk,qtest)  
else if ( present(var) ) then
  call hr4p_parar8(iarchi,ier,name,kk,qtest,var)
else
  call hr4p_parar8(iarchi,ier,name,kk,qtest)  
end if

return
end subroutine hr4pr8

subroutine hr4p_procformatr8(iarchi,ier,name,kk,qtest,var)


use cc_mpi
use newmpar_m
      
implicit none

integer, intent(in) :: iarchi, kk
integer, intent(out) :: ier
integer(kind=4), dimension(5) :: start, ncount
integer ipf, k, ca
integer(kind=4) idv, ndims
real(kind=8), dimension(:,:), intent(inout), optional :: var
real(kind=8), dimension(pil*pjl*pnpan,kk) :: rvar
real(kind=4) laddoff, lsf
logical, intent(in) :: qtest
character(len=*), intent(in) :: name
character(len=80) :: newname

ier = 0
      
do ipf = 0,mynproc-1
    
  ! get variable idv
  ier = nf90_inq_varid(pncid(ipf),name,idv)
  if ( ier==nf90_noerr ) then
    start(:)  = (/ 1, 1, 1, pprid(ipf), iarchi /)
    ncount(:) = (/ pil, pjl*pnpan, kk, 1, 1 /)   
    ! obtain scaling factors and offsets from attributes
    ier = nf90_get_att(pncid(ipf),idv,'add_offset',laddoff)
    if ( ier/=nf90_noerr ) laddoff = 0.
    ier = nf90_get_att(pncid(ipf),idv,'scale_factor',lsf)
    if ( ier/=nf90_noerr ) lsf = 1.
    ier = nf90_inquire_variable(pncid(ipf),idv,ndims=ndims)
    ier = nf90_get_var(pncid(ipf),idv,rvar,start=start(1:ndims),count=ncount(1:ndims))
    call ncmsg(name,ier)
    ! unpack data
    rvar(:,:) = rvar(:,:)*real(lsf,8) + real(laddoff,8)
  else
    start(1:4) = (/ 1, 1, pprid(ipf), iarchi /)
    ncount(1:4) = (/ pil, pjl*pnpan, 1, 1 /)
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
        if ( myid==0 .and. ipf==0 ) then
          write(6,*) '***absent field for ncid,name,ier: ',pncid(0),name,ier
        end if
        rvar(:,:) = 0._8 ! default value for missing field
        exit
      end if
      ! obtain scaling factors and offsets from attributes
      ier = nf90_get_att(pncid(ipf),idv,'add_offset',laddoff)
      if ( ier/=nf90_noerr ) laddoff = 0.
      ier = nf90_get_att(pncid(ipf),idv,'scale_factor',lsf)
      if ( ier/=nf90_noerr ) lsf = 1.
      ier = nf90_inquire_variable(pncid(ipf),idv,ndims=ndims)
      ier = nf90_get_var(pncid(ipf),idv,rvar(:,k),start=start(1:ndims),count=ncount(1:ndims))
      call ncmsg(name,ier)
      ! unpack data
      rvar(:,k) = rvar(:,k)*real(lsf,8) + real(laddoff,8)      
    end do
  end if ! ier

  if ( qtest ) then
    ! e.g., restart file or nogather=.true.
    ca = pil*pjl*pnpan*ipf
    var(1+ca:pil*pjl*pnpan+ca,1:kk) = rvar(:,:)
  else
    ! e.g., mesonest file
    if ( myid==0 ) then
      call host_hr4pr8(ipf,kk,rvar,var)
    else
      call proc_hr4pr8(kk,rvar)
    end if
  end if ! qtest

end do ! ipf

return
end subroutine hr4p_procformatr8

subroutine hr4p_parar8(iarchi,ier,name,kk,qtest,var)


use cc_mpi
use newmpar_m
      
implicit none

integer, intent(in) :: iarchi, kk
integer, intent(out) :: ier
integer(kind=4), dimension(4) :: start, ncount
integer ipf, k, ca
integer(kind=4) idv, ndims
real(kind=8), dimension(:,:), intent(inout), optional :: var
real(kind=8), dimension(pil*pjl*pnpan,kk) :: rvar
real(kind=4) laddoff, lsf
logical, intent(in) :: qtest
character(len=*), intent(in) :: name
character(len=80) :: newname

ier = 0
      
do ipf = 0,mynproc-1

  ! get variable idv
  ier = nf90_inq_varid(pncid(ipf),name,idv)
  if ( ier==nf90_noerr ) then
    start  = (/ 1, 1, 1, iarchi /)
    ncount = (/ pil, pjl*pnpan, kk, 1 /)   
    ! obtain scaling factors and offsets from attributes
    ier = nf90_get_att(pncid(ipf),idv,'add_offset',laddoff)
    if ( ier/=nf90_noerr ) laddoff = 0.
    ier = nf90_get_att(pncid(ipf),idv,'scale_factor',lsf)
    if ( ier/=nf90_noerr ) lsf = 1.
    ier = nf90_inquire_variable(pncid(ipf),idv,ndims=ndims)
    ier = nf90_get_var(pncid(ipf),idv,rvar,start=start(1:ndims),count=ncount(1:ndims))
    call ncmsg(name,ier)
    ! unpack data
    rvar(:,:) = rvar(:,:)*real(lsf,8) + real(laddoff,8)
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
      if ( ier/=nf90_noerr ) then
        if ( myid==0 .and. ipf==0 ) then
          write(6,*) '***absent field for ncid,name,ier: ',pncid(0),name,ier
        end if
        rvar(:,:) = 0._8 ! default value for missing field
        exit
      end if
      ! obtain scaling factors and offsets from attributes
      ier = nf90_get_att(pncid(ipf),idv,'add_offset',laddoff)
      if ( ier/=nf90_noerr ) laddoff = 0.
      ier = nf90_get_att(pncid(ipf),idv,'scale_factor',lsf)
      if ( ier/=nf90_noerr ) lsf = 1.
      ier = nf90_inquire_variable(pncid(ipf),idv,ndims=ndims)
      ier = nf90_get_var(pncid(ipf),idv,rvar(:,k),start=start(1:ndims),count=ncount(1:ndims))
      call ncmsg(name,ier)
      ! unpack data
      rvar(:,k) = rvar(:,k)*real(lsf,8) + real(laddoff,8)      
    end do
  end if ! ier

  if ( qtest ) then
    ! e.g., restart file or nogather=.true.
    ca = pil*pjl*pnpan*ipf
    var(1+ca:pil*pjl*pnpan+ca,1:kk) = rvar(:,:)
  else
    ! e.g., mesonest file
    if ( myid==0 ) then
      call host_hr4pr8(ipf,kk,rvar,var)
    else
      call proc_hr4pr8(kk,rvar)
    end if
  end if ! qtest

end do ! ipf

return
end subroutine hr4p_parar8


subroutine host_hr4pr8(ipf,kk,rvar,var)

use cc_mpi
use newmpar_m

implicit none

integer, intent(in) :: ipf, kk
integer jpf, ip, n, no, ca, cc, j, k
real(kind=8), dimension(:,:), intent(inout) :: var
real(kind=8), dimension(pil*pjl*pnpan,kk), intent(in) :: rvar
real(kind=8), dimension(:,:,:), allocatable :: globvar

if ( fnresid==1 ) then
    
  ip = ipf*fnresid             ! local file number
  do k = 1,kk
    do n = 0,pnpan-1
      no = n - pnoff(ip) + 1   ! global panel number of local file
      ca = pioff(ip,no) + pjoff(ip,no)*pil_g + no*pil_g*pil_g - pil_g
      cc = n*pil*pjl - pil
      do j = 1,pjl
        var(1+j*pil_g+ca:pil+j*pil_g+ca,k) = rvar(1+j*pil+cc:pil+j*pil+cc,k)
      end do
    end do
  end do
  
else

  allocate( globvar(pil*pjl*pnpan,kk,fnresid) )
  call ccmpi_gatherxr8(globvar,rvar,0,comm_ip)
  do jpf = 1,fnresid
    ip = ipf*fnresid + jpf - 1   ! local file number
    do k = 1,kk
      do n = 0,pnpan-1
        no = n - pnoff(ip) + 1   ! global panel number of local file
        ca = pioff(ip,no) + pjoff(ip,no)*pil_g + no*pil_g*pil_g - pil_g
        cc = n*pil*pjl - pil
        do j = 1,pjl
          var(1+j*pil_g+ca:pil+j*pil_g+ca,k) = globvar(1+j*pil+cc:pil+j*pil+cc,k,jpf)
        end do
      end do
    end do
  end do
  deallocate( globvar )
  
end if

return
end subroutine host_hr4pr8

subroutine proc_hr4pr8(kk,rvar)

use cc_mpi

implicit none

integer, intent(in) :: kk
real(kind=8), dimension(pil*pjl*pnpan,kk), intent(in) :: rvar
real(kind=8), dimension(0,0,0) :: globvar 

call ccmpi_gatherxr8(globvar,rvar,0,comm_ip)

return
end subroutine proc_hr4pr8

subroutine histrd5r4(iarchi,ier,name,ik,kk,ll,var,ifull,nogather)
      
use cc_mpi, only : myid, ccmpi_reduce, histrd5_begin, histrd5_end, &
                   start_log, end_log, ccmpi_distribute
use parm_m
      
implicit none
      
integer, intent(in) :: iarchi, ik, kk, ll, ifull
integer, intent(out) :: ier
real, dimension(:,:,:), intent(inout) :: var ! may be dummy argument from myid/=0
logical, intent(in), optional :: nogather
logical ngflag
character(len=*), intent(in) :: name

call START_LOG(histrd5_begin)

ngflag = .false.
if ( present(nogather) ) then
  ngflag = nogather
end if

if ( (ifull/=6*ik*ik.and.ptest) .or. ngflag ) then
  ! read local arrays without gather and distribute
  call hr5p(iarchi,ier,name,kk,ll,.true.,var)
  if ( ier==0 .and. mod(ktau,nmaxpr)==0 .and. myid==0 ) then  
    write(6,'(" done histrd5 ",a48,i4,i3)') trim(name),ier,iarchi
  end if

else if ( ifull/=6*ik*ik ) then  
  ! read local arrays with gather and distribute
  if ( myid==0 ) then
    call hr5sa(iarchi,ier,name,ik,kk,ll,var,ifull)
  else 
    call hr5p(iarchi,ier,name,kk,ll,.false.)
    call ccmpi_distribute(var)
  end if
  
else
  ! read global arrays and gather on myid==0
  if ( myid==0 ) then
    call hr5sa(iarchi,ier,name,ik,kk,ll,var,ifull)
  else 
    call hr5p(iarchi,ier,name,kk,ll,.false.)
  end if

end if

call END_LOG(histrd5_end)

return
end subroutine histrd5r4

!--------------------------------------------------------------------
! Gather 3D+time fields on myid==0

subroutine hr5sa(iarchi,ier,name,ik,kk,ll,var,ifull)
      
use cc_mpi, only : ccmpi_distribute
use parm_m
      
implicit none
      
integer, intent(in) :: iarchi, ik, kk, ll, ifull
integer, intent(out) :: ier
character(len=*), intent(in) :: name
real, dimension(:,:,:) :: var
real, dimension(:,:,:), allocatable :: globvar
real vmax, vmin
      
allocate( globvar(6*ik*ik,kk,ll) )
globvar(:,:,:) = 0.

call hr5p(iarchi,ier,name,kk,ll,.false.,globvar)     

if( ier==0 .and. mod(ktau,nmaxpr)==0 ) then
  vmax = maxval(globvar)
  vmin = minval(globvar)
  write(6,'(" done histrd5 ",a18,2i3,i4,i3,2f12.4)') trim(name),kk,ll,ier,iarchi,vmin,vmax
end if

! Have to return correct value of ier on all processes because it's 
! used for initialisation in calling routine
if ( ifull==6*ik*ik ) then
  ! read global arrays for myid==0
  var(1:6*ik*ik,1:kk,1:ll) = globvar(:,:,:)
else
  ! read local arrays with gather and distribute (i.e., change in number of processors)
  call ccmpi_distribute(var,globvar)
endif

deallocate( globvar )

return
end subroutine hr5sa      

!--------------------------------------------------------------
! This subroutine reads 3D+time input files

! when qtest=.true. the input grid decomposition should
! match the current processor decomposition.  We can then
! skip the MPI gather and distribute steps.
subroutine hr5p(iarchi,ier,name,kk,ll,qtest,var)

implicit none

integer, intent(in) :: iarchi, kk, ll
integer, intent(out) :: ier
logical, intent(in) :: qtest
real, dimension(:,:,:), intent(inout), optional :: var
character(len=*), intent(in) :: name

if ( resprocformat .and. present(var) ) then
  call hr5p_procformat(iarchi,ier,name,kk,ll,qtest,var)
else if ( resprocformat ) then
  call hr5p_procformat(iarchi,ier,name,kk,ll,qtest)  
else if ( present(var) ) then
  call hr5p_para(iarchi,ier,name,kk,ll,qtest,var)
else
  call hr5p_para(iarchi,ier,name,kk,ll,qtest)  
end if

return
end subroutine hr5p

subroutine hr5p_procformat(iarchi,ier,name,kk,ll,qtest,var)


use cc_mpi
use newmpar_m
      
implicit none

integer, intent(in) :: iarchi, kk, ll
integer, intent(out) :: ier
integer(kind=4), dimension(6) :: start, ncount
integer ipf, ca
integer(kind=4) idv, ndims
real, dimension(:,:,:), intent(inout), optional :: var
real, dimension(pil*pjl*pnpan,kk,ll) :: rvar
real(kind=4) laddoff, lsf
logical, intent(in) :: qtest
character(len=*), intent(in) :: name

ier = 0
      
do ipf = 0,mynproc-1
    
  ! get variable idv
  ier = nf90_inq_varid(pncid(ipf),name,idv)
  start(:)  = (/ 1, 1, 1, 1, pprid(ipf), iarchi /)
  ncount(:) = (/ pil, pjl*pnpan, kk, ll, 1, 1 /)   
  ! obtain scaling factors and offsets from attributes
  ier = nf90_get_att(pncid(ipf),idv,'add_offset',laddoff)
  if ( ier/=nf90_noerr ) laddoff = 0.
  ier = nf90_get_att(pncid(ipf),idv,'scale_factor',lsf)
  if ( ier/=nf90_noerr ) lsf = 1.
  ier = nf90_inquire_variable(pncid(ipf),idv,ndims=ndims)
  ier = nf90_get_var(pncid(ipf),idv,rvar,start=start(1:ndims),count=ncount(1:ndims))
  call ncmsg(name,ier)
   ! unpack data
  rvar(:,:,:) = rvar(:,:,:)*real(lsf) + real(laddoff)

  if ( qtest ) then
    ! e.g., restart file or nogather=.true.
    ca = pil*pjl*pnpan*ipf
    var(1+ca:pil*pjl*pnpan+ca,1:kk,1:ll) = rvar(:,:,:)
  else
    ! e.g., mesonest file
    if ( myid==0 ) then
      call host_hr5p(ipf,kk,ll,rvar,var)
    else
      call proc_hr5p(kk,ll,rvar)
    end if
  end if ! qtest

end do ! ipf

return
end subroutine hr5p_procformat

subroutine hr5p_para(iarchi,ier,name,kk,ll,qtest,var)

use cc_mpi
use newmpar_m
      
implicit none

integer, intent(in) :: iarchi, kk, ll
integer, intent(out) :: ier
integer(kind=4), dimension(5) :: start, ncount
integer ipf, ca
integer(kind=4) idv, ndims
real, dimension(:,:,:), intent(inout), optional :: var
real, dimension(pil*pjl*pnpan,kk,ll) :: rvar
real(kind=4) laddoff, lsf
logical, intent(in) :: qtest
character(len=*), intent(in) :: name

ier = 0
      
do ipf = 0,mynproc-1

  ! get variable idv
  ier = nf90_inq_varid(pncid(ipf),name,idv)
  start  = (/ 1, 1, 1, 1, iarchi /)
  ncount = (/ pil, pjl*pnpan, kk, ll, 1 /)   
  ! obtain scaling factors and offsets from attributes
  ier = nf90_get_att(pncid(ipf),idv,'add_offset',laddoff)
  if ( ier/=nf90_noerr ) laddoff = 0.
  ier = nf90_get_att(pncid(ipf),idv,'scale_factor',lsf)
  if ( ier/=nf90_noerr ) lsf = 1.
  ier = nf90_inquire_variable(pncid(ipf),idv,ndims=ndims)
  ier = nf90_get_var(pncid(ipf),idv,rvar,start=start(1:ndims),count=ncount(1:ndims))
  call ncmsg(name,ier)
  ! unpack data
  rvar(:,:,:) = rvar(:,:,:)*real(lsf) + real(laddoff)

  if ( qtest ) then
    ! e.g., restart file or nogather=.true.
    ca = pil*pjl*pnpan*ipf
    var(1+ca:pil*pjl*pnpan+ca,1:kk,1:ll) = rvar(:,:,:)
  else
    ! e.g., mesonest file
    if ( myid==0 ) then
      call host_hr5p(ipf,kk,ll,rvar,var)
    else
      call proc_hr5p(kk,ll,rvar)
    end if
  end if ! qtest

end do ! ipf

return
end subroutine hr5p_para

subroutine host_hr5p(ipf,kk,ll,rvar,var)

use cc_mpi
use newmpar_m

implicit none

integer, intent(in) :: ipf, kk, ll
integer jpf, ip, n, no, ca, cc, j, k, l
real, dimension(:,:,:), intent(inout) :: var
real, dimension(pil*pjl*pnpan,kk,ll), intent(in) :: rvar
real, dimension(:,:,:,:), allocatable :: globvar 

if ( fnresid==1 ) then
    
  ip = ipf*fnresid               ! local file number
  do l = 1,ll
    do k = 1,kk
      do n = 0,pnpan-1
        no = n - pnoff(ip) + 1   ! global panel number of local file
        ca = pioff(ip,no) + pjoff(ip,no)*pil_g + no*pil_g*pil_g - pil_g
        cc = n*pil*pjl - pil
        do j = 1,pjl
          var(1+j*pil_g+ca:pil+j*pil_g+ca,k,l) = rvar(1+j*pil+cc:pil+j*pil+cc,k,l)
        end do
      end do
    end do
  end do
  
else

  allocate( globvar(pil*pjl*pnpan,kk,ll,fnresid) )
  call ccmpi_gatherx(globvar,rvar,0,comm_ip)
  do jpf = 1,fnresid
    ip = ipf*fnresid + jpf - 1   ! local file number
    do l = 1,ll
      do k = 1,kk
        do n = 0,pnpan-1
          no = n - pnoff(ip) + 1   ! global panel number of local file
          ca = pioff(ip,no) + pjoff(ip,no)*pil_g + no*pil_g*pil_g - pil_g
          cc = n*pil*pjl - pil
          do j = 1,pjl
            var(1+j*pil_g+ca:pil+j*pil_g+ca,k,l) = globvar(1+j*pil+cc:pil+j*pil+cc,k,l,jpf)
          end do
        end do
      end do
    end do
  end do
  deallocate( globvar )
  
end if

return
end subroutine host_hr5p

subroutine proc_hr5p(kk,ll,rvar)

use cc_mpi

implicit none

integer, intent(in) :: kk, ll
real, dimension(pil*pjl*pnpan,kk,ll), intent(in) :: rvar
real, dimension(0,0,0) :: globvar 

call ccmpi_gatherx(globvar,rvar,0,comm_ip)

return
end subroutine proc_hr5p

!--------------------------------------------------------------   
! Interface for reading 3D+time fields (double precision version)
subroutine histrd5r8(iarchi,ier,name,ik,kk,ll,var,ifull,nogather)
      
use cc_mpi, only : myid, ccmpi_reducer8, histrd5_begin, histrd5_end, &
                   start_log, end_log, ccmpi_distributer8
use parm_m
      
implicit none
      
integer, intent(in) :: iarchi, ik, kk, ll, ifull
integer, intent(out) :: ier
real(kind=8), dimension(:,:,:), intent(inout) :: var ! may be dummy argument from myid/=0
logical, intent(in), optional :: nogather
logical ngflag
character(len=*), intent(in) :: name

call START_LOG(histrd5_begin)

ngflag = .false.
if ( present(nogather) ) then
  ngflag = nogather
end if

if ( (ifull/=6*ik*ik.and.ptest) .or. ngflag ) then
  ! read local arrays without gather and distribute
  call hr5pr8(iarchi,ier,name,kk,ll,.true.,var)
  if ( ier==0 .and. mod(ktau,nmaxpr)==0 .and. myid==0 ) then  
    write(6,'(" done histrd5r8 ",a46,i4,i3)') trim(name),ier,iarchi
  end if

else if ( ifull/=6*ik*ik ) then  
  ! read local arrays with gather and distribute
  if ( myid==0 ) then
    call hr5sar8(iarchi,ier,name,ik,kk,ll,var,ifull)
  else 
    call hr5pr8(iarchi,ier,name,kk,ll,.false.)
    call ccmpi_distributer8(var)
  end if
  
else
  ! read global arrays and gather on myid==0
  if ( myid==0 ) then
    call hr5sar8(iarchi,ier,name,ik,kk,ll,var,ifull)
  else 
    call hr5pr8(iarchi,ier,name,kk,ll,.false.)
  end if

end if

call END_LOG(histrd5_end)

return
end subroutine histrd5r8

!--------------------------------------------------------------------
! Gather 3D+time fields on myid==0

subroutine hr5sar8(iarchi,ier,name,ik,kk,ll,var,ifull)
      
use cc_mpi, only : ccmpi_distributer8
use parm_m
      
implicit none
      
integer, intent(in) :: iarchi, ik, kk, ll, ifull
integer, intent(out) :: ier
character(len=*), intent(in) :: name
real(kind=8), dimension(:,:,:) :: var
real(kind=8), dimension(:,:,:), allocatable :: globvar
real(kind=8) vmax, vmin
      
allocate( globvar(6*ik*ik,kk,ll) )
globvar(:,:,:) = 0._8

call hr5pr8(iarchi,ier,name,kk,ll,.false.,globvar)     

if( ier==0 .and. mod(ktau,nmaxpr)==0 ) then
  vmax = maxval(globvar)
  vmin = minval(globvar)
  write(6,'(" done histrd5r8 ",a16,2i3,i4,i3,2f12.4)') trim(name),kk,ll,ier,iarchi,vmin,vmax
end if

! Have to return correct value of ier on all processes because it's 
! used for initialisation in calling routine
if ( ifull==6*ik*ik ) then
  ! read global arrays for myid==0
  var(1:6*ik*ik,1:kk,1:ll) = globvar(:,:,:)
else
  ! read local arrays with gather and distribute (i.e., change in number of processors
  ! or single input file)
  call ccmpi_distributer8(var,globvar)
endif

deallocate( globvar )

return
end subroutine hr5sar8      

!--------------------------------------------------------------
! This subroutine reads 3D+time input files (Double precision version)

! when qtest=.true. the input grid decomposition should
! match the current processor decomposition.  We can then
! skip the MPI gather and distribute steps.
subroutine hr5pr8(iarchi,ier,name,kk,ll,qtest,var)

implicit none

integer, intent(in) :: iarchi, kk, ll
integer, intent(out) :: ier
logical, intent(in) :: qtest
real(kind=8), dimension(:,:,:), intent(inout), optional :: var
character(len=*), intent(in) :: name

if ( resprocformat .and. present(var) ) then
  call hr5p_procformatr8(iarchi,ier,name,kk,ll,qtest,var)
else if ( resprocformat ) then
  call hr5p_procformatr8(iarchi,ier,name,kk,ll,qtest)  
else if ( present(var) ) then
  call hr5p_parar8(iarchi,ier,name,kk,ll,qtest,var)
else
  call hr5p_parar8(iarchi,ier,name,kk,ll,qtest)  
end if

return
end subroutine hr5pr8

subroutine hr5p_procformatr8(iarchi,ier,name,kk,ll,qtest,var)

use cc_mpi
use newmpar_m
      
implicit none

integer, intent(in) :: iarchi, kk, ll
integer, intent(out) :: ier
integer(kind=4), dimension(6) :: start, ncount
integer ipf, ca
integer(kind=4) idv, ndims
real(kind=8), dimension(:,:,:), intent(inout), optional :: var
real(kind=8), dimension(pil*pjl*pnpan,kk,ll) :: rvar
real(kind=4) laddoff, lsf
logical, intent(in) :: qtest
character(len=*), intent(in) :: name

ier = 0
      
do ipf = 0,mynproc-1
    
  ! get variable idv
  ier = nf90_inq_varid(pncid(ipf),name,idv)
  start(:)  = (/ 1, 1, 1, 1, pprid(ipf), iarchi /)
  ncount(:) = (/ pil, pjl*pnpan, kk, ll, 1, 1 /)   
  ! obtain scaling factors and offsets from attributes
  ier = nf90_get_att(pncid(ipf),idv,'add_offset',laddoff)
  if ( ier/=nf90_noerr ) laddoff = 0.
  ier = nf90_get_att(pncid(ipf),idv,'scale_factor',lsf)
  if ( ier/=nf90_noerr ) lsf = 1.
  ier = nf90_inquire_variable(pncid(ipf),idv,ndims=ndims)
  ier = nf90_get_var(pncid(ipf),idv,rvar,start=start(1:ndims),count=ncount(1:ndims))
  call ncmsg(name,ier)
  ! unpack data
  rvar(:,:,:) = rvar(:,:,:)*real(lsf,8) + real(laddoff,8)

  if ( qtest ) then
    ! e.g., restart file or nogather=.true.
    ca = pil*pjl*pnpan*ipf
    var(1+ca:pil*pjl*pnpan+ca,1:kk,1:ll) = rvar(:,:,:)
  else
    ! e.g., mesonest file
    if ( myid==0 ) then
      call host_hr5pr8(ipf,kk,ll,rvar,var)
    else
      call proc_hr5pr8(kk,ll,rvar)
    end if
  end if ! qtest

end do ! ipf

return
end subroutine hr5p_procformatr8

subroutine hr5p_parar8(iarchi,ier,name,kk,ll,qtest,var)


use cc_mpi
use newmpar_m
      
implicit none

integer, intent(in) :: iarchi, kk, ll
integer, intent(out) :: ier
integer(kind=4), dimension(5) :: start, ncount
integer ipf, ca
integer(kind=4) idv, ndims
real(kind=8), dimension(:,:,:), intent(inout), optional :: var
real(kind=8), dimension(pil*pjl*pnpan,kk,ll) :: rvar
real(kind=4) laddoff, lsf
logical, intent(in) :: qtest
character(len=*), intent(in) :: name

ier = 0
      
do ipf = 0,mynproc-1

  ! get variable idv
  ier = nf90_inq_varid(pncid(ipf),name,idv)
  start  = (/ 1, 1, 1, 1, iarchi /)
  ncount = (/ pil, pjl*pnpan, kk, ll, 1 /)   
  ! obtain scaling factors and offsets from attributes
  ier = nf90_get_att(pncid(ipf),idv,'add_offset',laddoff)
  if ( ier/=nf90_noerr ) laddoff = 0.
  ier = nf90_get_att(pncid(ipf),idv,'scale_factor',lsf)
  if ( ier/=nf90_noerr ) lsf = 1.
  ier = nf90_inquire_variable(pncid(ipf),idv,ndims=ndims)
  ier = nf90_get_var(pncid(ipf),idv,rvar,start=start(1:ndims),count=ncount(1:ndims))
  call ncmsg(name,ier)
  ! unpack data
  rvar(:,:,:) = rvar(:,:,:)*real(lsf,8) + real(laddoff,8)

  if ( qtest ) then
    ! e.g., restart file or nogather=.true.
    ca = pil*pjl*pnpan*ipf
    var(1+ca:pil*pjl*pnpan+ca,1:kk,1:ll) = rvar(:,:,:)
  else
    ! e.g., mesonest file
    if ( myid==0 ) then
      call host_hr5pr8(ipf,kk,ll,rvar,var)
    else
      call proc_hr5pr8(kk,ll,rvar)
    end if
  end if ! qtest

end do ! ipf

return
end subroutine hr5p_parar8

subroutine host_hr5pr8(ipf,kk,ll,rvar,var)

use cc_mpi
use newmpar_m

implicit none

integer, intent(in) :: ipf, kk, ll
integer jpf, ip, n, no, ca, cc, j, k, l
real(kind=8), dimension(:,:,:), intent(inout) :: var
real(kind=8), dimension(pil*pjl*pnpan,kk,ll), intent(in) :: rvar
real(kind=8), dimension(:,:,:,:), allocatable :: globvar 

if ( fnresid==1 ) then
    
  ip = ipf*fnresid               ! local file number
  do l = 1,ll
    do k = 1,kk
      do n = 0,pnpan-1
        no = n - pnoff(ip) + 1   ! global panel number of local file
        ca = pioff(ip,no) + pjoff(ip,no)*pil_g + no*pil_g*pil_g - pil_g
        cc = n*pil*pjl - pil
        do j = 1,pjl
          var(1+j*pil_g+ca:pil+j*pil_g+ca,k,l) = rvar(1+j*pil+cc:pil+j*pil+cc,k,l)
        end do
      end do  
    end do
  end do
  
else

  allocate( globvar(pil*pjl*pnpan,kk,ll,fnresid) )
  call ccmpi_gatherxr8(globvar,rvar,0,comm_ip)
  do jpf = 1,fnresid
    ip = ipf*fnresid + jpf - 1   ! local file number
    do l = 1,ll
      do k = 1,kk
        do n = 0,pnpan-1
          no = n - pnoff(ip) + 1   ! global panel number of local file
          ca = pioff(ip,no) + pjoff(ip,no)*pil_g + no*pil_g*pil_g - pil_g
          cc = n*pil*pjl - pil
          do j = 1,pjl
            var(1+j*pil_g+ca:pil+j*pil_g+ca,k,l) = globvar(1+j*pil+cc:pil+j*pil+cc,k,l,jpf)
          end do
        end do  
      end do
    end do
  end do
  deallocate( globvar )
  
end if

return
end subroutine host_hr5pr8

subroutine proc_hr5pr8(kk,ll,rvar)

use cc_mpi

implicit none

integer, intent(in) :: kk,ll
real(kind=8), dimension(pil*pjl*pnpan,kk,ll), intent(in) :: rvar
real(kind=8), dimension(0,0,0,0) :: globvar 

call ccmpi_gatherxr8(globvar,rvar,0,comm_ip)

return
end subroutine proc_hr5pr8

!--------------------------------------------------------------
! This subroutine opens parallel input files
subroutine histopen(ncid,ifile,ier)
      
use cc_mpi
use newmpar_m
use parm_m
      
implicit none
      
integer, parameter :: nihead = 54
      
integer, dimension(0:5) :: duma, dumb
integer, dimension(12) :: idum
integer, intent(out) :: ncid, ier
integer is, ipf, dmode
integer ipin, ipin_f, ipin_new, nxpr, nypr
integer ltst, der, myrank
integer resprocmode
integer, dimension(:,:), allocatable, save :: dum_off
integer, dimension(:), allocatable, save :: resprocmap_inv
integer, dimension(:), allocatable, save :: procfileowner
integer(kind=4), dimension(1) :: start, ncount
integer(kind=4), dimension(nihead) :: lahead
integer(kind=4) lncid, lidum, ldid, lvid, llen
character(len=*), intent(in) :: ifile
character(len=1031) pfile
character(len=8) fdecomp

if ( myid==0 ) then
  ! attempt to open single file with myid==0
  ier = nf90_open(ifile,nf90_nowrite,lncid)
  ncid = lncid
  fnproc = 1              ! number of files to be read over all processors
  dmode = 0               ! Single file (dmode=0), Face decomposition (dmode=1),
                          ! Depreciated (dmode=2) or Uniform decomposition (dmode=3)
  pil = 0                 ! Number of X grid points within a file panel
  pjl = 0                 ! Number of Y grid points within a file panel
  pnpan = 0               ! Number of panels in file
  ptest = .false.         ! Files match current processor (e.g., Restart file), allowing MPI gather/scatter to be avoided
  pfall = .false.         ! Every processor has been assigned at least one file, no need to Bcast metadata data
  resprocformat = .false. ! procformat file format with multiple processes per file
  resprocmode = 0         ! base number of processes in procformat files
      
  ! attempt to open parallel files
  if ( ier/=nf90_noerr ) then
    write(pfile,"(a,'.',i6.6)") trim(ifile), 0
    ier = nf90_open(pfile,nf90_nowrite,lncid)
    ncid = lncid
    if ( ier/=nf90_noerr ) then
      write(6,*) "WARN: Cannot open ",trim(pfile)
      write(6,*) "WARN: Cannot open ",trim(ifile)
    else  
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
        case('uniform')  ! old uniform
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
        
    if ( allocated(pioff) ) then
      deallocate( pioff, pjoff, pnoff )
    end if
    allocate( pioff(0:fnproc-1,0:5), pjoff(0:fnproc-1,0:5) )
    allocate( pnoff(0:fnproc-1) )
        
    select case(dmode)
      case(0) ! no decomposition
        pnpan=6
        pnoff=1
        pioff=0
        pjoff=0
        pil=pil_g
        pjl=pil_g
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

    ptest = .false.
    if ( nproc==fnproc ) then
      if ( pil_g==il_g .and. pjl_g==jl_g ) then
        if ( uniform_decomp ) then
          if ( dmode==3 ) then
            ptest = .true.
          end if
        else
          if ( dmode==1 ) then
            ptest = .true.
          end if
        end if
      end if
    end if

    write(6,*) "Found pil_g,pjl_g,fnproc        ",pil_g,pjl_g,fnproc
    write(6,*) "Found dmode,ptest               ",dmode,ptest
    write(6,*) "Found resprocformat,resprocmode ",resprocformat,resprocmode
    
  end if

  idum(1) = fnproc
  if ( resprocformat ) then
    idum(2) = 1
  else
    idum(2) = 0
  end if
  idum(3) = resprocmode
  idum(4) = pil
  idum(5) = pjl
  idum(6) = pnpan
  if (ptest) then
    idum(7) = 1
  else
    idum(7) = 0      
  end if
  idum(8) = ier
  idum(9) = pka_g
  idum(10) = pko_g
  idum(11) = pil_g
  idum(12) = pjl_g
  
  if ( resprocformat ) then
    allocate( resprocmap_inv(0:fnproc-1) )
    der = nf90_inq_varid(lncid,'gprocessor',lvid)
    start(1) = 1
    ncount(1) = fnproc
    der = nf90_get_var(lncid,lvid,resprocmap_inv,start,ncount)
  end if
  
  write(6,*) "Broadcasting file metadata"
end if

! Broadcast file metadata
call ccmpi_bcast(idum(1:12),0,comm_world)
fnproc        = idum(1)      ! number of files to be read
resprocformat = (idum(2)==1) ! test for procformat file format
resprocmode   = idum(3)      ! base number of 'files' per input in procformat
pil           = idum(4)      ! width of panel in each file
pjl           = idum(5)      ! length of panel in each file
pnpan         = idum(6)      ! number of panels in each file
ptest         = (idum(7)==1) ! test for match between files and processes
ier           = idum(8)      ! file error flag
pka_g         = idum(9)      ! number of atmosphere levels
pko_g         = idum(10)     ! number of ocean levels
pil_g         = idum(11)     ! global grid size
pjl_g         = idum(12)     ! global grid size

if ( ier/=nf90_noerr ) return

if ( myid==0 ) then
  write(6,*) "Opening data files"
end if

! calculate number of files to be read on this processor
fnresid = min( fnproc, nproc ) 
do while ( mod(fnproc,fnresid)/=0 )
  fnresid = fnresid - 1     ! limit on processor ranks that will read files    
end do
fncount = fnproc/fnresid
if ( myid<fnresid) then
  mynproc = fncount  ! calculate the number of files to be read per process
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
  allocate(pfown(0:mynproc-1))
  pfown(:) = .false.
end if
! Rank 0 can start with the second file, because the first file has already been opened
if ( myid==0 ) then 
  is=1
  pncid(0)=ncid
else
  is=0
end if


! distribute comms
if ( myid==0 ) then
  write(6,*) "Splitting comms for distributing file data with fnresid ",fnresid
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

    ! procformat  
      
    allocate( pprid(0:mynproc-1) )
      
    ! copy process map
    allocate( procfileowner(0:fnproc-1) )
    procfileowner(:) = -1
    if ( .not.allocated(resprocmap_inv) ) then
      allocate(resprocmap_inv(0:fnproc-1))
    end if
    call ccmpi_bcast(resprocmap_inv,0,comm_ip)  

    ! update required process and load files
    do ipf = 0,mynproc-1
      ipin = ipf*fnresid + myid                   ! parallel file number
      ipin_new = resprocmap_inv(ipin)             ! remap files
      pprid(ipf) = mod(ipin_new, resprocmode) + 1 ! procformat process
    end do
    
    if ( myid==0 ) then
      procfileowner(0) = 0
      pfown(0) = .true.
    end if
    
    do ipf = is,mynproc-1
      ipin = ipf*fnresid + myid               ! parallel file number
      ipin_new = resprocmap_inv(ipin)         ! remap files
      ipin_f = ipin_new/resprocmode           ! procformat file
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
    deallocate(procfileowner)
    deallocate(resprocmap_inv)    
   
  else
    
    ! usual parallel file
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
   
  end if ! resprocformat ..else..
end if   ! mynproc>0



pfall=fnresid==nproc  ! are all processes associated with a file?
                      ! this means we do not need to Bcst file metadata
if ( mynproc>0 ) then
  ncid=pncid(0)       ! set ncid to the first file handle as onthefly
                      ! assumes changes in ncid reflect a new file
                      ! and hence updates the metadata
end if


allocate( dum_off(0:fnproc-1,1:13) )
if ( myid==0 ) then
  write(6,*) "Broadcast file coordinate data"
  dum_off(0:fnproc-1,1:6)  = pioff(0:fnproc-1,0:5)
  dum_off(0:fnproc-1,7:12) = pjoff(0:fnproc-1,0:5)
  dum_off(0:fnproc-1,13)   = pnoff(0:fnproc-1)
else
  allocate( pioff(0:fnproc-1,0:5), pjoff(0:fnproc-1,0:5) )
  allocate( pnoff(0:fnproc-1))
end if
call ccmpi_bcast(dum_off,0,comm_world)
pioff(0:fnproc-1,0:5) = dum_off(0:fnproc-1,1:6)
pjoff(0:fnproc-1,0:5) = dum_off(0:fnproc-1,7:12)
pnoff(0:fnproc-1)     = dum_off(0:fnproc-1,13)
deallocate( dum_off )


if ( myid==0 ) then
  write(6,*) "Create file RMA windows with kblock = ",kblock
end if
call ccmpi_filewincreate(kblock)

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

if ( allocated(pprid) ) then
  deallocate(pprid)
end if

if ( allocated(pncid) ) then
  plen = size(pncid)
  do ipf = 0,plen-1
    if ( pfown(ipf) ) then
      ierr = nf90_close(pncid(ipf))
    end if
  end do
  deallocate(pfown)
  deallocate(pncid)
end if

if ( myid==0 ) then
  write(6,*) 'Removing file communication group'
end if
call ccmpi_commfree(comm_ip)

if (allocated(pioff)) then
  deallocate(pioff,pjoff,pnoff)
end if

if ( myid==0 ) then
  write(6,*) 'Removing file RMA windows'
end if
call ccmpi_filewinfree

ncidold = -1 ! flag onthefly to load metadata

return
end subroutine histclose

!--------------------------------------------------------------
! transforms 3d array from dimension kk in vertical to kl

! This version of vertint can interpolate from host models with
! a greater number of vertical levels than the nested model.
subroutine vertint(told,t,n,kk,sigin)

use sigs_m
use newmpar_m
use parm_m

implicit none

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
    t(1:ifull,1:kl) = told(:,:)
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
  do k = 1,kl
    if ( sig(k)>=sigin(1) ) then
      ka(k) = 2
      wta(k) = 0.
      klapse = k   ! i.e. T lapse correction for k<=klapse
    else if ( sig(k)<=sigin(kk) ) then   ! at top
      ka(k) = kk
      wta(k) = 1.
    else
      do while ( sig(k)<=sigin(kin) .and. kin<kk )
        kin = kin + 1
      end do
      ka(k) = kin
      wta(k) = (sigin(kin-1)-sig(k))/(sigin(kin-1)-sigin(kin))
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

do k = 1,kl
  ! N.B. "a" denotes "above", "b" denotes "below"
  t(1:ifull,k) = wta(k)*told(:,ka(k)) + (1.-wta(k))*told(:,ka(k)-1)
enddo    ! k loop

if ( n==1 .and. klapse/=0 ) then  ! for T lapse correction
  do k = 1,klapse
    ! assume 6.5 deg/km, with dsig=.1 corresponding to 1 km
    t(1:ifull,k) = t(1:ifull,k) + (sig(k)-sigin(1))*6.5/.1
  enddo    ! k loop
else if ( n==2 ) then  ! for qg do a -ve fix
  t(1:ifull,:) = max(t(1:ifull,:), qgmin)
else if ( n==5 ) then  ! for qfg, qlg do a -ve fix
  t(1:ifull,:) = max(t(1:ifull,:), 0.)
endif
      
return
end subroutine vertint

!--------------------------------------------------------------------
! This subroutine advances input date by the amount of time defined by mtimer_r
subroutine datefix(kdate_r,ktime_r,mtimer_r,allleap)

use cc_mpi
use newmpar_m
use parm_m

implicit none

integer, intent(inout) :: kdate_r,ktime_r,mtimer_r
integer, intent(in), optional :: allleap
integer, dimension(12) :: mdays = (/31,28,31,30,31,30,31,31,30,31,30,31/)
integer iyr,imo,iday,ihr,imins
integer mtimerh,mtimerm,mtimer
integer mdays_save, leap_l
integer, parameter :: minsday = 1440

if ( present(allleap) ) then
  leap_l = allleap
else
  leap_l = leap
end if

if ( kdate_r>=00600000 .and. kdate_r<=00991231 ) then   ! old 1960-1999
  kdate_r=kdate_r+19000000
  if ( myid==0 ) write(6,*) 'For Y2K kdate_r altered to: ',kdate_r
endif

iyr=kdate_r/10000
imo=(kdate_r-10000*iyr)/100
iday=kdate_r-10000*iyr-100*imo
ihr=ktime_r/100
imins=ktime_r-100*ihr
if ( myid==0 ) then
  write(6,*) 'entering datefix'
  write(6,*) 'iyr,imo,iday:       ',iyr,imo,iday
  write(6,*) 'ihr,imins,mtimer_r: ',ihr,imins,mtimer_r
end if

mdays(2)=28
if ( leap_l==1 ) then
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
    if ( leap_l==1 ) then
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
subroutine getzinp(jyear,jmonth,jday,jhour,jmin,mins,allleap)

use cc_mpi
use dates_m
use parm_m

implicit none

integer, intent(out) :: jyear,jmonth,jday,jhour,jmin ! start date of run
integer, intent(out) :: mins                         ! elapsed time from start of year
integer mstart
integer, dimension(12) :: ndoy
! days from beginning of year (1st Jan is 0)
integer, dimension(12), parameter :: odoy=(/ 0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334 /) 
logical, intent(in), optional :: allleap ! force use of leap days even if leap=0
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

mstart = 1440*(ndoy(jmonth)+jday-1) + 60*jhour + jmin ! mins from start of year
! mtimer contains number of minutes since the start of the run.
mins = mtimer + mstart

return
end subroutine getzinp

!--------------------------------------------------------------------
! DEFINE ATTRIBUTES
subroutine attrib(cdfid,dim,ndim,name,lname,units,xmin,xmax,daily,itype)

use cc_mpi
use newmpar_m
use parm_m

implicit none

integer, intent(in) :: cdfid, itype, ndim
integer, intent(in) :: daily
integer, dimension(ndim), intent(in) :: dim
integer ier
integer(kind=4) vtype, idv, lcdfid, lsize, lcompression
integer(kind=4), dimension(ndim) :: ldim
integer(kind=4), dimension(ndim) :: chunks
real, intent(in) :: xmin, xmax
real(kind=4) lscalef, laddoff
character(len=*), intent(in) :: name
character(len=*), intent(in) :: lname
character(len=*), intent(in) :: units

if ( itype==1 ) then
  vtype = nf90_short
else if ( itype==2 ) then
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
#ifdef usenc3
ier = nf90_def_var(lcdfid, name, vtype, ldim, idv)
#else
if ( procformat .and. ndim>3 ) then
  ! MJT notes - PR identified (/il, jl, kl,vnode_nproc, min(10, tlen)/) as optimal.
  ! However, here we simplify the code and PR reports that the performance is
  ! similar
  select case(ndim)
    case(6)
      chunks = (/ il, jl, 1, 1, vnode_nproc, 1 /)
    case(5)
      chunks = (/ il, jl, 1, vnode_nproc, 1 /)
    case(4)
      chunks = (/ il, jl, vnode_nproc, 1 /)
    case default
      write(6,*) "ERROR: Invalid ndim in attrib ",ndim
      call ccmpi_abort(-1)
  end select
  ! MJT notes - compression can degrade model runtime when using a subset of processes
  ! for file IO.  However, file sizes can increase significantly without compression.
  lcompression = compression   
  ier = nf90_def_var(lcdfid, name, vtype, ldim, idv)
  ier = nf90_def_var_deflate(lcdfid, idv, 1_4, 1_4, lcompression)
  ier = nf90_def_var_chunking(lcdfid, idv, NF90_CHUNKED, chunks )
else
  lcompression = compression
  ier = nf90_def_var(lcdfid, name, vtype, ldim, idv)
  ier = nf90_def_var_deflate(lcdfid, idv, 1_4, 1_4, lcompression)
end if
#endif
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
real, dimension(ifull,1) :: wvar
character(len=*), intent(in) :: sname
logical, intent(in) :: local, lwrite

if (.not.lwrite) then
  wvar(:,1)=real(nf90_fill_float)
else
  wvar(:,1)=var(:)
end if

if ( local .and. procformat ) then
  call fw3lp(wvar,sname,idnc,iarch,1)  
else if ( procformat ) then
  call fw3p(wvar,1)  
else if ( local ) then
  call fw3l(wvar,sname,idnc,iarch,1)
else if ( myid==0 ) then
  call fw3a(wvar,sname,idnc,iarch,1)
else
  call ccmpi_gather(wvar(1:ifull,1:1))
end if

return
end subroutine histwrt3r4

!--------------------------------------------------------------------
! 3D NETCDF WRITE ARRAY ROUTINES (double precision version)
subroutine histwrt3r8(var,sname,idnc,iarch,local,lwrite)

use cc_mpi              ! CC MPI routines
use newmpar_m           ! Grid parameters
use parm_m              ! Model configuration

implicit none

integer, intent(in) :: idnc, iarch
real(kind=8), dimension(ifull), intent(in) :: var
real(kind=8), dimension(ifull,1) :: wvar
character(len=*), intent(in) :: sname
logical, intent(in) :: local, lwrite

if (.not.lwrite) then
  wvar(:,1)=real(nf90_fill_float,8)
else
  wvar(:,1)=var(:)
end if

if ( local .and. procformat ) then
  call fw3lpr8(wvar,sname,idnc,iarch,1)  
else if ( procformat ) then
  call fw3pr8(wvar,1)  
else if ( local ) then
  call fw3lr8(wvar,sname,idnc,iarch,1)
else if ( myid==0 ) then
  call fw3ar8(wvar,sname,idnc,iarch,1)
else
  call ccmpi_gatherr8(wvar(1:ifull,1:1))
end if

return
end subroutine histwrt3r8


subroutine freqwrite(fncid,cname,fiarch,istep,local,datain)

use cc_mpi               ! CC MPI routines
use newmpar_m            ! Grid parameters
use parm_m               ! Model configuration

implicit none
      
integer, intent(in) :: fncid, fiarch, istep
real, dimension(ifull,istep), intent(in) :: datain
logical, intent(in) :: local
character(len=*), intent(in) :: cname
      
if ( local .and. procformat ) then
  call fw3lp(datain,cname,fncid,fiarch,istep)   
else if ( procformat ) then
  call fw3p(datain,istep)  
else if ( local ) then
  call fw3l(datain,cname,fncid,fiarch,istep)
else if ( myid==0 ) then
  call fw3a(datain,cname,fncid,fiarch,istep)
else
  call ccmpi_gather(datain(1:ifull,1:istep))
endif
     
return
end subroutine freqwrite

! procformat and local(write)
subroutine fw3lp(var,sname,idnc,iarch,istep)

use cc_mpi               ! CC MPI routines
use newmpar_m            ! Grid parameters
use parm_m               ! Model configuration
      
implicit none
      
integer, intent(in) :: idnc, iarch, istep
integer ier, i, v
integer(kind=4) :: lidnc, mid, vtype, ndims
integer(kind=4), dimension(4) :: start, ncount
integer(kind=2), dimension(ifull,vnode_nproc,istep) :: ipack_g
real, dimension(ifull,istep), intent(in) :: var
real, dimension(ifull,vnode_nproc,istep) :: var_g
real, dimension(ifull,istep,vnode_nproc) :: var_t
real(kind=4) laddoff, lscale_f
character(len=*), intent(in) :: sname

start = (/ 1, 1, 1, iarch /)
ncount = (/ il, jl, vnode_nproc, istep /)

!if ( useiobuffer ) then
!  ! MJT notes - move this to its own subroutine ...  
!  !call add_iobuffer(idnc,mid,ndims,ifull,istep,vnode_nproc,start,ncount,var)
!  write(6,*) "ERROR: iobuffer not yet implemented"
!  call ccmpi_abort(-1)
!  return
!end if

call ccmpi_gatherx(var_t,var,0,comm_vnode)
do v = 1,vnode_nproc
  do i = 1,istep    
    var_g(1:ifull,v,i) = var_t(1:ifull,i,v)
  end do
end do

lidnc = idnc
ier = nf90_inq_varid(lidnc,sname,mid)
call ncmsg(sname,ier)
ier = nf90_inquire_variable(lidnc,mid,xtype=vtype,ndims=ndims)
if ( vtype==nf90_short ) then
  if ( all(var_g>9.8E36) ) then
    ipack_g(:,:,:) = missval
  else
    ier = nf90_get_att(lidnc,mid,'add_offset',laddoff)
    ier = nf90_get_att(lidnc,mid,'scale_factor',lscale_f)
    do i = 1,istep
      do v = 1,vnode_nproc        
        ipack_g(:,v,i) = nint(max(min((var_g(:,v,i)-real(laddoff))/real(lscale_f),real(maxv)),real(minv)),2)
      end do
    end do
  end if
  ier = nf90_put_var(lidnc,mid,ipack_g,start=start(1:ndims),count=ncount(1:ndims))
else
  ier = nf90_put_var(lidnc,mid,var_g,start=start(1:ndims),count=ncount(1:ndims))
end if
call ncmsg(sname,ier)

if ( mod(ktau,nmaxpr)==0 .and. myid==0 ) then
  if ( any(abs(var-real(nf90_fill_float))<1.e-20) ) then
    write(6,'(" histwrt3 ",a20,i8,a7)') sname,iarch,"missing"
  else
    write(6,'(" histwrt3 ",a20,i8)') sname,iarch
  end if
end if

return
end subroutine fw3lp

! procformat without local(write)
subroutine fw3p(var,istep)

use cc_mpi               ! CC MPI routines
use newmpar_m            ! Grid parameters
use parm_m               ! Model configuration
      
implicit none
      
integer, intent(in) :: istep
real, dimension(ifull,istep), intent(in) :: var
real, dimension(0,0,0) :: var_t

!if ( useiobuffer ) then
!  ! MJT notes - move this to its own subroutine ...  
!  !call add_iobuffer(idnc,mid,ndims,ifull,istep,vnode_nproc,start,ncount,var)
!  write(6,*) "ERROR: iobuffer not yet implemented"
!  call ccmpi_abort(-1)
!end if

call ccmpi_gatherx(var_t,var,0,comm_vnode)

return
end subroutine fw3p

! pure local(write)
subroutine fw3l(var,sname,idnc,iarch,istep)

use cc_mpi               ! CC MPI routines
use newmpar_m            ! Grid parameters
use parm_m               ! Model configuration
      
implicit none
      
integer, intent(in) :: idnc, iarch, istep
integer ier, i
integer(kind=4) :: lidnc, mid, vtype, ndims
integer(kind=4), dimension(3) :: start, ncount
integer(kind=2), dimension(ifull,istep) :: ipack
real, dimension(ifull,istep), intent(in) :: var
real(kind=4) laddoff, lscale_f
character(len=*), intent(in) :: sname

start = (/ 1, 1, iarch /)
ncount = (/ il, jl, istep /)

lidnc = idnc
ier = nf90_inq_varid(lidnc,sname,mid)
call ncmsg(sname,ier)
ier = nf90_inquire_variable(lidnc,mid,xtype=vtype,ndims=ndims)
if ( vtype==nf90_short ) then
  if ( all(var>9.8E36) ) then
    ipack = missval
  else
    ier = nf90_get_att(lidnc,mid,'add_offset',laddoff)
    ier = nf90_get_att(lidnc,mid,'scale_factor',lscale_f)
    do i = 1,istep
      ipack(:,i) = nint(max(min((var(:,i)-real(laddoff))/real(lscale_f),real(maxv)),real(minv)),2)
    end do
  end if
  ier = nf90_put_var(lidnc,mid,ipack,start=start(1:ndims),count=ncount(1:ndims))
else
  ier = nf90_put_var(lidnc,mid,var,start=start(1:ndims),count=ncount(1:ndims))
end if
call ncmsg(sname,ier)
if ( mod(ktau,nmaxpr)==0 .and. myid==0 ) then
  if ( any(abs(var-real(nf90_fill_float))<1.e-20) ) then
    write(6,'(" histwrt3 ",a20,i8,a7)') sname,iarch,"missing"
  else
    write(6,'(" histwrt3 ",a20,i8)') sname,iarch
  end if
end if

return
end subroutine fw3l
      
! global(write) with single file
subroutine fw3a(var,sname,idnc,iarch,istep)

use cc_mpi               ! CC MPI routines
use newmpar_m            ! Grid parameters
use parm_m               ! Model configuration

implicit none

integer, intent(in) :: idnc, iarch, istep
integer ier, imn, imx, jmn, jmx, iq, i
integer(kind=4) lidnc, mid, vtype, ndims
integer(kind=4), dimension(3) :: start, ncount
integer(kind=2), dimension(:,:), allocatable :: ipack
real, dimension(ifull,istep), intent(in) :: var
real, dimension(:,:), allocatable :: globvar
real varn, varx
real(kind=4) laddoff, lscale_f
character(len=*), intent(in) :: sname
      
allocate( globvar(ifull_g,istep) )
call ccmpi_gather(var(1:ifull,1:istep), globvar(1:ifull_g,1:istep))

start = (/ 1, 1, iarch /)
ncount = (/ il_g, jl_g, istep /)

!     find variable index
lidnc = idnc
ier = nf90_inq_varid(lidnc,sname,mid)
call ncmsg(sname,ier)
ier = nf90_inquire_variable(lidnc,mid,xtype=vtype,ndims=ndims)
if ( vtype==nf90_short ) then
  allocate( ipack(ifull_g,istep) )  
  if ( all(globvar>9.8e36) ) then
    ipack = missval
  else
    ier = nf90_get_att(lidnc,mid,'add_offset',laddoff)
    ier = nf90_get_att(lidnc,mid,'scale_factor',lscale_f)
    do i = 1,istep
      ipack(:,i) = nint(max(min((globvar(:,i)-real(laddoff))/real(lscale_f),real(maxv)),real(minv)),2)
    end do
  endif
  ier = nf90_put_var(lidnc,mid,ipack,start=start(1:ndims),count=ncount(1:ndims))
  deallocate( ipack )
else
  ier = nf90_put_var(lidnc,mid,globvar,start=start(1:ndims),count=ncount(1:ndims))
endif
call ncmsg(sname,ier)

if ( mod(ktau,nmaxpr)==0 ) then
  if ( any(abs(globvar-real(nf90_fill_float))<1.e-20) ) then
    write(6,'(" histwrt3 ",a20,i4,a7)') sname,iarch,"missing"
  else
    varn = minval(globvar(:,1))
    varx = maxval(globvar(:,1))
    ! This should work ???? but sum trick is more portable???
    ! iq = minloc(globvar,dim=1)
    iq = sum(minloc(globvar(:,1)))
    ! Convert this 1D index to 2D
    imn = 1 + modulo(iq-1,il_g)
    jmn = 1 + (iq-1)/il_g
    iq = sum(maxloc(globvar(:,1)))
    ! Convert this 1D index to 2D
    imx = 1 + modulo(iq-1,il_g)
    jmx = 1 + (iq-1)/il_g
    write(6,'(" histwrt3 ",a20,i4,f12.4,2i4,f12.4,2i4,f12.4)') &
                    sname,iarch,varn,imn,jmn,varx,imx,jmx,    &
                    globvar(id+(jd-1)*il_g,1)
  end if
end if

deallocate( globvar )

return
end subroutine fw3a

! procformat and local(write) (double precision version)
subroutine fw3lpr8(var,sname,idnc,iarch,istep)

use cc_mpi               ! CC MPI routines
use newmpar_m            ! Grid parameters
use parm_m               ! Model configuration
      
implicit none
      
integer, intent(in) :: idnc, iarch, istep
integer ier, i, v
integer(kind=4) :: lidnc, mid, vtype, ndims
integer(kind=4), dimension(4) :: start, ncount
integer(kind=2), dimension(ifull,vnode_nproc,istep) :: ipack_g
real(kind=8), dimension(ifull,istep), intent(in) :: var
real(kind=8), dimension(ifull,vnode_nproc,istep) :: var_g
real(kind=8), dimension(ifull,istep,vnode_nproc) :: var_t
real(kind=4) laddoff, lscale_f
character(len=*), intent(in) :: sname

start = (/ 1, 1, 1, iarch /)
ncount = (/ il, jl, vnode_nproc, istep /)

!if ( useiobuffer ) then
!  ! MJT notes - move this to its own subroutine ...  
!  !call add_iobuffer(idnc,mid,ndims,ifull,istep,vnode_nproc,start,ncount,var)
!  write(6,*) "ERROR: iobuffer not yet implemented"
!  call ccmpi_abort(-1)
!  return
!end if

call ccmpi_gatherxr8(var_t,var,0,comm_vnode)
do v = 1,vnode_nproc
  do i = 1,istep    
    var_g(1:ifull,v,i) = var_t(1:ifull,i,v)
  end do
end do

lidnc = idnc
ier = nf90_inq_varid(lidnc,sname,mid)
call ncmsg(sname,ier)
ier = nf90_inquire_variable(lidnc,mid,xtype=vtype,ndims=ndims)
if ( vtype==nf90_short ) then
  if ( all(var_g>9.8E36_8) ) then
    ipack_g(:,:,:) = missval
  else
    ier = nf90_get_att(lidnc,mid,'add_offset',laddoff)
    ier = nf90_get_att(lidnc,mid,'scale_factor',lscale_f)
    do i = 1,istep
      do v = 1,vnode_nproc        
        ipack_g(:,v,i) = nint(max(min((var_g(:,v,i)-real(laddoff,8))/real(lscale_f,8),real(maxv,8)),real(minv,8)),2)
      end do
    end do
  end if
  ier = nf90_put_var(lidnc,mid,ipack_g,start=start(1:ndims),count=ncount(1:ndims))
else
  ier = nf90_put_var(lidnc,mid,var_g,start=start(1:ndims),count=ncount(1:ndims))
end if
call ncmsg(sname,ier)

if ( mod(ktau,nmaxpr)==0 .and. myid==0 ) then
  if ( any(abs(var-real(nf90_fill_float,8))<1.e-20_8) ) then
    write(6,'(" histwrt3r8 ",a20,i8,a7)') sname,iarch,"missing"
  else
    write(6,'(" histwrt3r8 ",a20,i8)') sname,iarch
  end if
end if

return
end subroutine fw3lpr8

! procformat without local(write) (double precision version)
subroutine fw3pr8(var,istep)

use cc_mpi               ! CC MPI routines
use newmpar_m            ! Grid parameters
use parm_m               ! Model configuration
      
implicit none
      
integer, intent(in) :: istep
real(kind=8), dimension(ifull,istep), intent(in) :: var
real(kind=8), dimension(0,0,0) :: var_t

!if ( useiobuffer ) then
!  ! MJT notes - move this to its own subroutine ...  
!  !call add_iobuffer(idnc,mid,ndims,ifull,istep,vnode_nproc,start,ncount,var)
!  write(6,*) "ERROR: iobuffer not yet implemented"
!  call ccmpi_abort(-1)
!end if

call ccmpi_gatherxr8(var_t,var,0,comm_vnode)

return
end subroutine fw3pr8

! pure local(write) (Double precision version)
subroutine fw3lr8(var,sname,idnc,iarch,istep)

use cc_mpi               ! CC MPI routines
use newmpar_m            ! Grid parameters
use parm_m               ! Model configuration
      
implicit none
      
integer, intent(in) :: idnc, iarch, istep
integer ier, i
integer(kind=4) :: lidnc, mid, vtype, ndims
integer(kind=4), dimension(3) :: start, ncount
integer(kind=2), dimension(ifull,istep) :: ipack
real(kind=8), dimension(ifull,istep), intent(in) :: var
real(kind=4) laddoff, lscale_f
character(len=*), intent(in) :: sname

start = (/ 1, 1, iarch /)
ncount = (/ il, jl, istep /)

lidnc = idnc
ier = nf90_inq_varid(lidnc,sname,mid)
call ncmsg(sname,ier)
ier = nf90_inquire_variable(lidnc,mid,xtype=vtype,ndims=ndims)
if ( vtype==nf90_short ) then
  if ( all(var>9.8E36_8) ) then
    ipack = missval
  else
    ier = nf90_get_att(lidnc,mid,'add_offset',laddoff)
    ier = nf90_get_att(lidnc,mid,'scale_factor',lscale_f)
    do i = 1,istep
      ipack(:,i) = nint(max(min((var(:,i)-real(laddoff,8))/real(lscale_f,8),real(maxv,8)),real(minv,8)),2)
    end do
  end if
  ier = nf90_put_var(lidnc,mid,ipack,start=start(1:ndims),count=ncount(1:ndims))
else
  ier = nf90_put_var(lidnc,mid,var,start=start(1:ndims),count=ncount(1:ndims))
end if
call ncmsg(sname,ier)
if ( mod(ktau,nmaxpr)==0 .and. myid==0 ) then
  if ( any(abs(var-real(nf90_fill_float,8))<1.e-20_8) ) then
    write(6,'(" histwrt3r8 ",a20,i8,a7)') sname,iarch,"missing"
  else
    write(6,'(" histwrt3r8 ",a20,i8)') sname,iarch
  end if
end if

return
end subroutine fw3lr8
      
! global(write) with single file (double precision version)
subroutine fw3ar8(var,sname,idnc,iarch,istep)

use cc_mpi               ! CC MPI routines
use newmpar_m            ! Grid parameters
use parm_m               ! Model configuration

implicit none

integer, intent(in) :: idnc, iarch, istep
integer ier, imn, imx, jmn, jmx, iq, i
integer(kind=4) lidnc, mid, vtype, ndims
integer(kind=4), dimension(3) :: start, ncount
integer(kind=2), dimension(:,:), allocatable :: ipack
real(kind=8), dimension(ifull,istep), intent(in) :: var
real(kind=8), dimension(:,:), allocatable :: globvar
real(kind=8) varn, varx
real(kind=4) laddoff, lscale_f
character(len=*), intent(in) :: sname
      
allocate( globvar(ifull_g,istep) )
call ccmpi_gatherr8(var(1:ifull,1:istep), globvar(1:ifull_g,1:istep))

start = (/ 1, 1, iarch /)
ncount = (/ il_g, jl_g, istep /)

!     find variable index
lidnc = idnc
ier = nf90_inq_varid(lidnc,sname,mid)
call ncmsg(sname,ier)
ier = nf90_inquire_variable(lidnc,mid,xtype=vtype,ndims=ndims)
if ( vtype==nf90_short ) then
  allocate( ipack(ifull_g,istep) )
  if ( all(globvar>9.8e36_8) ) then
    ipack = missval
  else
    ier = nf90_get_att(lidnc,mid,'add_offset',laddoff)
    ier = nf90_get_att(lidnc,mid,'scale_factor',lscale_f)
    do i = 1,istep
      ipack(:,i) = nint(max(min((globvar(:,i)-real(laddoff,8))/real(lscale_f,8),real(maxv,8)),real(minv,8)),2)
    end do
  endif
  ier = nf90_put_var(lidnc,mid,ipack,start=start(1:ndims),count=ncount(1:ndims))
  deallocate( ipack )
else
  ier = nf90_put_var(lidnc,mid,globvar,start=start(1:ndims),count=ncount(1:ndims))
endif
call ncmsg(sname,ier)

if ( mod(ktau,nmaxpr)==0 ) then
  if ( any(abs(globvar-real(nf90_fill_float,8))<1.e-20_8) ) then
    write(6,'(" histwrt3r8 ",a20,i4,a7)') sname,iarch,"missing"
  else
    varn = minval(globvar(:,1))
    varx = maxval(globvar(:,1))
    ! This should work ???? but sum trick is more portable???
    ! iq = minloc(globvar,dim=1)
    iq = sum(minloc(globvar(:,1)))
    ! Convert this 1D index to 2D
    imn = 1 + modulo(iq-1,il_g)
    jmn = 1 + (iq-1)/il_g
    iq = sum(maxloc(globvar(:,1)))
    ! Convert this 1D index to 2D
    imx = 1 + modulo(iq-1,il_g)
    jmx = 1 + (iq-1)/il_g
    write(6,'(" histwrt3r8 ",a20,i4,f12.4,2i4,f12.4,2i4,f12.4)') &
                    sname,iarch,varn,imn,jmn,varx,imx,jmx,    &
                    globvar(id+(jd-1)*il_g,1)
  end if
end if

deallocate( globvar )

return
end subroutine fw3ar8

!--------------------------------------------------------------------
! 4D NETCDF WRITE ARRAY ROUTINES
subroutine histwrt4r4(var,sname,idnc,iarch,local,lwrite)

use cc_mpi              ! CC MPI routines
use newmpar_m           ! Grid parameters
use parm_m              ! Model configuration

implicit none

integer, intent(in) :: idnc, iarch
integer ll
real, dimension(:,:), intent(in) :: var
real, dimension(ifull,size(var,2)) :: wvar
character(len=*), intent(in) :: sname
logical, intent(in) :: local, lwrite

ll = size(var,2)

if ( .not.lwrite ) then
  wvar=real(nf90_fill_float)
else
  wvar(:,:)=var(1:ifull,1:ll)
endif

if ( local .and. procformat ) then
  call hw4lp(wvar,sname,idnc,iarch)  
else if ( procformat ) then
  call hw4p(wvar)  
else if ( local ) then
  call hw4l(wvar,sname,idnc,iarch)
else if ( myid==0 ) then
  call hw4a(wvar,sname,idnc,iarch)
else
  call ccmpi_gather(wvar(1:ifull,1:ll))
endif

return
end subroutine histwrt4r4

subroutine histwrt4r8(var,sname,idnc,iarch,local,lwrite)

use cc_mpi              ! CC MPI routines
use newmpar_m           ! Grid parameters
use parm_m              ! Model configuration

implicit none

integer, intent(in) :: idnc, iarch
integer ll
real(kind=8), dimension(:,:), intent(in) :: var
real(kind=8), dimension(ifull,size(var,2)) :: wvar
character(len=*), intent(in) :: sname
logical, intent(in) :: local, lwrite

ll = size(var,2)

if ( .not.lwrite ) then
  wvar=real(nf90_fill_float)
else
  wvar(:,:)=var(1:ifull,1:ll)
endif

if ( local .and. procformat ) then
  call hw4lpr8(wvar,sname,idnc,iarch)  
else if ( procformat ) then
  call hw4pr8(wvar)  
else if ( local ) then
  call hw4lr8(wvar,sname,idnc,iarch)
else if ( myid==0 ) then
  call hw4ar8(wvar,sname,idnc,iarch)
else
  call ccmpi_gatherr8(wvar(1:ifull,1:ll))
endif

return
end subroutine histwrt4r8

! procformat and local(write)
subroutine hw4lp(var,sname,idnc,iarch)

use cc_mpi               ! CC MPI routines
use newmpar_m            ! Grid parameters
use parm_m               ! Model configuration

implicit none

integer, intent(in) :: idnc, iarch
integer iq, k, ier, v, ll
integer(kind=4) mid, vtype, lidnc, ndims
integer(kind=4), dimension(5) :: start, ncount
real, dimension(:,:), intent(in) :: var
real, dimension(ifull,size(var,2),vnode_nproc) :: var_g
real(kind=4) laddoff, lscale_f
character(len=*), intent(in) :: sname
integer(kind=2), dimension(ifull,size(var,2),vnode_nproc) :: ipack_g

ll = size(var,2)
start = (/ 1, 1, 1, 1, iarch /)
ncount = (/ il, jl, ll, vnode_nproc, 1 /)

!if ( useiobuffer ) then
!  ! MJT notes - move this to its own subroutine ...  
!  !call add_iobuffer(idnc,mid,ndims,ifull,istep,vnode_nproc,start,ncount,var)
!  write(6,*) "ERROR: iobuffer not yet implemented"
!  call ccmpi_abort(-1)
!  return
!end if

call ccmpi_gatherx(var_g,var,0,comm_vnode)

lidnc = idnc
ier = nf90_inq_varid(lidnc,sname,mid)
call ncmsg(sname,ier)
ier = nf90_inquire_variable(lidnc,mid,xtype=vtype,ndims=ndims)

if ( vtype==nf90_short ) then
  if ( all(var>9.8e36) ) then
    ipack_g(:,:,:) = missval
  else
    ier = nf90_get_att(lidnc,mid,'add_offset',laddoff)
    ier = nf90_get_att(lidnc,mid,'scale_factor',lscale_f)
    do v = 1,vnode_nproc
      do k = 1,ll
        do iq = 1,ifull
          ipack_g(iq,k,v) = nint(max(min((var_g(iq,k,v)-real(laddoff))/real(lscale_f),real(maxv)),real(minv)),2)
        end do
      end do
    end do
  end if
  ier = nf90_put_var(lidnc,mid,ipack_g,start=start(1:ndims),count=ncount(1:ndims))
else
  ier = nf90_put_var(lidnc,mid,var_g,start=start(1:ndims),count=ncount(1:ndims))
end if
call ncmsg(sname,ier)

if ( mod(ktau,nmaxpr)==0 .and. myid==0 ) then
  if ( any(abs(var-real(nf90_fill_float))<1.e-20) ) then
    write(6,'(" histwrt4 ",a20,i4,a7)') sname,iarch,"missing"
  else
    write(6,'(" histwrt4 ",a20,i4)') sname,iarch
  end if
end if

return
end subroutine hw4lp      

! procformat without local(write)
subroutine hw4p(var)

use cc_mpi               ! CC MPI routines
use newmpar_m            ! Grid parameters
use parm_m               ! Model configuration

implicit none

real, dimension(:,:), intent(in) :: var
real, dimension(0,0,0) :: var_g

!if ( useiobuffer ) then
!  ! MJT notes - move this to its own subroutine ...  
!  !call add_iobuffer(idnc,mid,ndims,ifull,istep,vnode_nproc,start,ncount,var)
!  write(6,*) "ERROR: iobuffer not yet implemented"
!  call ccmpi_abort(-1)
!end if

call ccmpi_gatherx(var_g,var,0,comm_vnode)

return
end subroutine hw4p      

! pure local(write)
subroutine hw4l(var,sname,idnc,iarch)

use cc_mpi               ! CC MPI routines
use newmpar_m            ! Grid parameters
use parm_m               ! Model configuration

implicit none

integer, intent(in) :: idnc, iarch
integer iq, k, ier, ll
integer(kind=4) mid, vtype, lidnc, ndims
integer(kind=4), dimension(4) :: start, ncount
real, dimension(:,:), intent(in) :: var
real(kind=4) laddoff, lscale_f
character(len=*), intent(in) :: sname
integer(kind=2), dimension(ifull,size(var,2)) :: ipack

ll = size(var,2)
start = (/ 1, 1, 1, iarch /)
ncount = (/ il, jl, ll, 1 /)

lidnc = idnc
ier = nf90_inq_varid(lidnc,sname,mid)
call ncmsg(sname,ier)
ier = nf90_inquire_variable(lidnc,mid,xtype=vtype,ndims=ndims)
if ( vtype==nf90_short ) then
  if ( all(var>9.8e36) ) then
    ipack = missval
  else
    ier = nf90_get_att(lidnc,mid,'add_offset',laddoff)
    ier = nf90_get_att(lidnc,mid,'scale_factor',lscale_f)
    do k = 1,ll
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

if ( mod(ktau,nmaxpr)==0 .and. myid==0 ) then
  if ( any(abs(var-real(nf90_fill_float))<1.e-20) ) then
    write(6,'(" histwrt4 ",a20,i4,a7)') sname,iarch,"missing"
  else
    write(6,'(" histwrt4 ",a20,i4)') sname,iarch
  end if
end if

return
end subroutine hw4l      

! global write with single file
subroutine hw4a(var,sname,idnc,iarch)

use cc_mpi              ! CC MPI routines
use newmpar_m           ! Grid parameters
use parm_m              ! Model configuration

implicit none

integer, intent(in) :: idnc, iarch
integer ier, imx, jmx, kmx, iq, k, ll
integer, dimension(2) :: max_result
integer(kind=4) mid, vtype, lidnc, ndims
integer(kind=4), dimension(4) :: start, ncount
real varn, varx
real, dimension(:,:), intent(in) :: var
real, dimension(:,:), allocatable :: globvar
real(kind=4) laddoff, lscale_f
character(len=*), intent(in) :: sname
integer(kind=2), dimension(:,:), allocatable :: ipack
      
ll = size(var,2)
allocate( globvar(ifull_g,ll) )
call ccmpi_gather(var(1:ifull,1:ll), globvar(1:ifull_g,1:ll))
start = (/ 1, 1, 1, iarch /)
ncount = (/ il_g, jl_g, ll, 1 /)

!     find variable index
lidnc = idnc
ier = nf90_inq_varid(lidnc,sname,mid)
call ncmsg(sname,ier)
ier = nf90_inquire_variable(lidnc, mid, xtype=vtype, ndims=ndims)
if ( vtype==nf90_short ) then
  allocate( ipack(ifull_g,ll) )  
  if ( all(globvar>9.8e36) )then
    ipack = missval
  else
    ier = nf90_get_att(lidnc,mid,'add_offset',laddoff)
    ier = nf90_get_att(lidnc,mid,'scale_factor',lscale_f)
    do k = 1,ll
      do iq = 1,ifull_g
        ipack(iq,k) = nint(max(min((globvar(iq,k)-real(laddoff))/real(lscale_f),real(maxv)),real(minv)),2)
      end do
    end do
  end if
  ier = nf90_put_var(lidnc,mid,ipack,start=start(1:ndims),count=ncount(1:ndims))
  deallocate( ipack )
else
  ier = nf90_put_var(lidnc,mid,globvar,start=start(1:ndims),count=ncount(1:ndims))
endif
call ncmsg(sname,ier)

if ( mod(ktau,nmaxpr)==0 ) then
  if ( any(abs(globvar-real(nf90_fill_float))<1.e-20) ) then
    write(6,'(" histwrt4 ",a7,i4,a7)') sname,iarch,"missing"
  else
    varn = minval(globvar)
    varx = maxval(globvar)
    max_result = maxloc(globvar)
    kmx = max_result(2)
    iq = max_result(1)
    ! Convert this 1D index to 2D
    imx = 1 + modulo(iq-1,il_g)
    jmx = 1 + (iq-1)/il_g
    write(6,'(" histwrt4 ",a20,i4,2f12.4,3i4,f12.4)') sname,iarch,varn,varx,imx,jmx,kmx,globvar(id+(jd-1)*il_g,nlv)
  end if
end if

deallocate( globvar )

return
end subroutine hw4a      

! procformat and local(write)
subroutine hw4lpr8(var,sname,idnc,iarch)

use cc_mpi               ! CC MPI routines
use newmpar_m            ! Grid parameters
use parm_m               ! Model configuration

implicit none

integer, intent(in) :: idnc, iarch
integer iq, k, ier, v, ll
integer(kind=4) mid, vtype, lidnc, ndims
integer(kind=4), dimension(5) :: start, ncount
real(kind=8), dimension(:,:), intent(in) :: var
real(kind=8), dimension(ifull,size(var,2),vnode_nproc) :: var_g
real(kind=4) laddoff, lscale_f
character(len=*), intent(in) :: sname
integer(kind=2), dimension(ifull,size(var,2),vnode_nproc) :: ipack_g

ll = size(var,2)
start = (/ 1, 1, 1, 1, iarch /)
ncount = (/ il, jl, ll, vnode_nproc, 1 /)

!if ( useiobuffer ) then
!  ! MJT notes - move this to its own subroutine ...  
!  !call add_iobuffer(idnc,mid,ndims,ifull,istep,vnode_nproc,start,ncount,var)
!  write(6,*) "ERROR: iobuffer not yet implemented"
!  call ccmpi_abort(-1)
!  return
!end if

call ccmpi_gatherxr8(var_g,var,0,comm_vnode)

lidnc = idnc
ier = nf90_inq_varid(lidnc,sname,mid)
call ncmsg(sname,ier)
ier = nf90_inquire_variable(lidnc,mid,xtype=vtype,ndims=ndims)

if ( vtype==nf90_short ) then
  if ( all(var>9.8e36_8) ) then
    ipack_g(:,:,:) = missval
  else
    ier = nf90_get_att(lidnc,mid,'add_offset',laddoff)
    ier = nf90_get_att(lidnc,mid,'scale_factor',lscale_f)
    do v = 1,vnode_nproc
      do k = 1,ll
        do iq = 1,ifull
          ipack_g(iq,k,v) = nint(max(min((var_g(iq,k,v)-real(laddoff,8))/real(lscale_f,8),real(maxv,8)),real(minv,8)),2)
        end do
      end do
    end do
  end if
  ier = nf90_put_var(lidnc,mid,ipack_g,start=start(1:ndims),count=ncount(1:ndims))
else
  ier = nf90_put_var(lidnc,mid,var_g,start=start(1:ndims),count=ncount(1:ndims))
end if
call ncmsg(sname,ier)

if ( mod(ktau,nmaxpr)==0 .and. myid==0 ) then
  if ( any(abs(var-real(nf90_fill_float,8))<1.e-20_8) ) then
    write(6,'(" histwrt4r8 ",a20,i4,a7)') sname,iarch,"missing"
  else
    write(6,'(" histwrt4r8 ",a20,i4)') sname,iarch
  end if
end if

return
end subroutine hw4lpr8    

! procformat without local(write)
subroutine hw4pr8(var)

use cc_mpi               ! CC MPI routines
use newmpar_m            ! Grid parameters
use parm_m               ! Model configuration

implicit none

real(kind=8), dimension(:,:), intent(in) :: var
real(kind=8), dimension(0,0,0) :: var_g

!if ( useiobuffer ) then
!  ! MJT notes - move this to its own subroutine ...  
!  !call add_iobuffer(idnc,mid,ndims,ifull,istep,vnode_nproc,start,ncount,var)
!  write(6,*) "ERROR: iobuffer not yet implemented"
!  call ccmpi_abort(-1)
!end if

call ccmpi_gatherxr8(var_g,var,0,comm_vnode)

return
end subroutine hw4pr8      

! pure local(write)
subroutine hw4lr8(var,sname,idnc,iarch)

use cc_mpi               ! CC MPI routines
use newmpar_m            ! Grid parameters
use parm_m               ! Model configuration

implicit none

integer, intent(in) :: idnc, iarch
integer iq, k, ier, ll
integer(kind=4) mid, vtype, lidnc, ndims
integer(kind=4), dimension(4) :: start, ncount
real(kind=8), dimension(:,:), intent(in) :: var
real(kind=4) laddoff, lscale_f
character(len=*), intent(in) :: sname
integer(kind=2), dimension(ifull,size(var,2)) :: ipack

ll = size(var,2)
start = (/ 1, 1, 1, iarch /)
ncount = (/ il, jl, ll, 1 /)

lidnc = idnc
ier = nf90_inq_varid(lidnc,sname,mid)
call ncmsg(sname,ier)
ier = nf90_inquire_variable(lidnc,mid,xtype=vtype,ndims=ndims)
if ( vtype==nf90_short ) then
  if ( all(var>9.8e36_8) ) then
    ipack = missval
  else
    ier = nf90_get_att(lidnc,mid,'add_offset',laddoff)
    ier = nf90_get_att(lidnc,mid,'scale_factor',lscale_f)
    do k = 1,ll
      do iq = 1,ifull
        ipack(iq,k) = nint(max(min((var(iq,k)-real(laddoff,8))/real(lscale_f,8),real(maxv,8)),real(minv,8)),2)
      end do
    end do
  end if
  ier = nf90_put_var(lidnc,mid,ipack,start=start(1:ndims),count=ncount(1:ndims))
else
  ier = nf90_put_var(lidnc,mid,var,start=start(1:ndims),count=ncount(1:ndims))
endif
call ncmsg(sname,ier)

if ( mod(ktau,nmaxpr)==0 .and. myid==0 ) then
  if ( any(abs(var-real(nf90_fill_float,8))<1.e-20_8) ) then
    write(6,'(" histwrt4r8 ",a20,i4,a7)') sname,iarch,"missing"
  else
    write(6,'(" histwrt4r8 ",a20,i4)') sname,iarch
  end if
end if

return
end subroutine hw4lr8      

! global write with single file
subroutine hw4ar8(var,sname,idnc,iarch)

use cc_mpi              ! CC MPI routines
use newmpar_m           ! Grid parameters
use parm_m              ! Model configuration

implicit none

integer, intent(in) :: idnc, iarch
integer ier, imx, jmx, kmx, iq, k, ll
integer, dimension(2) :: max_result
integer(kind=4) mid, vtype, lidnc, ndims
integer(kind=4), dimension(4) :: start, ncount
real(kind=8) varn, varx
real(kind=8), dimension(:,:), intent(in) :: var
real(kind=8), dimension(:,:), allocatable :: globvar
real(kind=4) laddoff, lscale_f
character(len=*), intent(in) :: sname
integer(kind=2), dimension(:,:), allocatable :: ipack
      
ll = size(var,2)
allocate( globvar(ifull_g,ll) )
call ccmpi_gatherr8(var(1:ifull,1:ll), globvar(1:ifull_g,1:ll))
start = (/ 1, 1, 1, iarch /)
ncount = (/ il_g, jl_g, ll, 1 /)

!     find variable index
lidnc = idnc
ier = nf90_inq_varid(lidnc,sname,mid)
call ncmsg(sname,ier)
ier = nf90_inquire_variable(lidnc, mid, xtype=vtype, ndims=ndims)
if ( vtype==nf90_short ) then
  allocate( ipack(ifull_g,ll) )  
  if ( all(globvar>9.8e36_8) )then
    ipack = missval
  else
    ier = nf90_get_att(lidnc,mid,'add_offset',laddoff)
    ier = nf90_get_att(lidnc,mid,'scale_factor',lscale_f)
    do k = 1,ll
      do iq = 1,ifull_g
        ipack(iq,k) = nint(max(min((globvar(iq,k)-real(laddoff,8))/real(lscale_f,8),real(maxv,8)),real(minv,8)),2)
      end do
    end do
  end if
  ier = nf90_put_var(lidnc,mid,ipack,start=start(1:ndims),count=ncount(1:ndims))
  deallocate( ipack )
else
  ier = nf90_put_var(lidnc,mid,globvar,start=start(1:ndims),count=ncount(1:ndims))
endif
call ncmsg(sname,ier)

if ( mod(ktau,nmaxpr)==0 ) then
  if ( any(abs(globvar-real(nf90_fill_float,8))<1.e-20_8) ) then
    write(6,'(" histwrt4r8 ",a20,i4,a7)') sname,iarch,"missing"
  else
    varn = minval(globvar)
    varx = maxval(globvar)
    max_result = maxloc(globvar)
    kmx = max_result(2)
    iq = max_result(1)
    ! Convert this 1D index to 2D
    imx = 1 + modulo(iq-1,il_g)
    jmx = 1 + (iq-1)/il_g
    write(6,'(" histwrt4r8 ",a20,i4,2f12.4,3i4,f12.4)') sname,iarch,varn,varx,imx,jmx,kmx,globvar(id+(jd-1)*il_g,nlv)
  end if
end if

deallocate( globvar )

return
end subroutine hw4ar8

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
character(len=*), intent(in) :: sname
logical, intent(in) :: local, lwrite

kk = size(var,2)
ll = size(var,3)

if ( .not.lwrite ) then
  wvar=real(nf90_fill_float)
else
  wvar(:,:,:)=var(1:ifull,1:kk,1:ll)
endif

if ( local .and. procformat ) then
  call hw5lp(wvar,sname,idnc,iarch)  
else if ( procformat ) then
  call hw5p(wvar)  
else if ( local ) then
  call hw5l(wvar,sname,idnc,iarch)
else if ( myid==0 ) then
  call hw5a(wvar,sname,idnc,iarch)
else
  call ccmpi_gather(wvar(1:ifull,1:kk,1:ll))
endif

return
end subroutine histwrt5r4

subroutine histwrt5r8(var,sname,idnc,iarch,local,lwrite)

use cc_mpi              ! CC MPI routines
use newmpar_m           ! Grid parameters
use parm_m              ! Model configuration

implicit none

integer, intent(in) :: idnc, iarch
integer kk, ll
real(kind=8), dimension(:,:,:), intent(in) :: var
real(kind=8), dimension(ifull,size(var,2),size(var,3)) :: wvar
character(len=*), intent(in) :: sname
logical, intent(in) :: local, lwrite

kk = size(var,2)
ll = size(var,3)

if ( .not.lwrite ) then
  wvar=real(nf90_fill_float)
else
  wvar(:,:,:)=var(1:ifull,1:kk,1:ll)
endif

if ( local .and. procformat ) then
  call hw5lpr8(wvar,sname,idnc,iarch)  
else if ( procformat ) then
  call hw5pr8(wvar)  
else if ( local ) then
  call hw5lr8(wvar,sname,idnc,iarch)
else if ( myid==0 ) then
  call hw5ar8(wvar,sname,idnc,iarch)
else
  call ccmpi_gatherr8(wvar(1:ifull,1:kk,1:ll))
endif

return
end subroutine histwrt5r8

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
real, dimension(:,:,:), intent(in) :: var
real, dimension(ifull,size(var,2),size(var,3),vnode_nproc) :: var_g
real(kind=4) laddoff, lscale_f
character(len=*), intent(in) :: sname
integer(kind=2), dimension(ifull,size(var,2),size(var,3),vnode_nproc) :: ipack_g

kk = size(var,2)
ll = size(var,3)
start = (/ 1, 1, 1, 1, 1, iarch /)
ncount = (/ il, jl, kk, ll, vnode_nproc, 1 /)

!if ( useiobuffer ) then
!  ! MJT notes - move this to its own subroutine ...  
!  !call add_iobuffer(idnc,mid,ndims,ifull,istep,vnode_nproc,start,ncount,var)
!  write(6,*) "ERROR: iobuffer not yet implemented"
!  call ccmpi_abort(-1)
!  return
!end if

call ccmpi_gatherx(var_g,var,0,comm_vnode)

lidnc = idnc
ier = nf90_inq_varid(lidnc,sname,mid)
call ncmsg(sname,ier)
ier = nf90_inquire_variable(lidnc,mid,xtype=vtype,ndims=ndims)

if ( vtype==nf90_short ) then
  if ( all(var>9.8e36) ) then
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
else
  ier = nf90_put_var(lidnc,mid,var_g,start=start(1:ndims),count=ncount(1:ndims))
end if
call ncmsg(sname,ier)

if ( mod(ktau,nmaxpr)==0 .and. myid==0 ) then
  if ( any(abs(var-real(nf90_fill_float))<1.e-20) ) then
    write(6,'(" histwrt5 ",a20,i4,a7)') sname,iarch,"missing"
  else
    write(6,'(" histwrt5 ",a20,i4)') sname,iarch
  end if
end if

return
end subroutine hw5lp      

! procformat without local(write)
subroutine hw5p(var)

use cc_mpi               ! CC MPI routines
use newmpar_m            ! Grid parameters
use parm_m               ! Model configuration

implicit none

real, dimension(:,:,:), intent(in) :: var
real, dimension(0,0,0,0) :: var_g

!if ( useiobuffer ) then
!  ! MJT notes - move this to its own subroutine ...  
!  !call add_iobuffer(idnc,mid,ndims,ifull,istep,vnode_nproc,start,ncount,var)
!  write(6,*) "ERROR: iobuffer not yet implemented"
!  call ccmpi_abort(-1)
!end if

call ccmpi_gatherx(var_g,var,0,comm_vnode)

return
end subroutine hw5p      

! pure local(write)
subroutine hw5l(var,sname,idnc,iarch)

use cc_mpi               ! CC MPI routines
use newmpar_m            ! Grid parameters
use parm_m               ! Model configuration

implicit none

integer, intent(in) :: idnc, iarch
integer iq, k, ier, kk, l, ll
integer(kind=4) mid, vtype, lidnc, ndims
integer(kind=4), dimension(5) :: start, ncount
real, dimension(:,:,:), intent(in) :: var
real(kind=4) laddoff, lscale_f
character(len=*), intent(in) :: sname
integer(kind=2), dimension(ifull,size(var,2),size(var,3)) :: ipack

kk = size(var,2)
ll = size(var,3)
start = (/ 1, 1, 1, 1, iarch /)
ncount = (/ il, jl, kk, ll, 1 /)

lidnc = idnc
ier = nf90_inq_varid(lidnc,sname,mid)
call ncmsg(sname,ier)
ier = nf90_inquire_variable(lidnc,mid,xtype=vtype,ndims=ndims)
if ( vtype==nf90_short ) then
  if ( all(var>9.8e36) ) then
    ipack = missval
  else
    ier = nf90_get_att(lidnc,mid,'add_offset',laddoff)
    ier = nf90_get_att(lidnc,mid,'scale_factor',lscale_f)
    do l = 1,ll
      do k = 1,kk
        do iq = 1,ifull
          ipack(iq,k,l) = nint(max(min((var(iq,k,l)-real(laddoff))/real(lscale_f),real(maxv)),real(minv)),2)
        end do  
      end do
    end do
  end if
  ier = nf90_put_var(lidnc,mid,ipack,start=start(1:ndims),count=ncount(1:ndims))
else
  ier = nf90_put_var(lidnc,mid,var,start=start(1:ndims),count=ncount(1:ndims))
endif
call ncmsg(sname,ier)

if ( mod(ktau,nmaxpr)==0 .and. myid==0 ) then
  if ( any(abs(var-real(nf90_fill_float))<1.e-20) ) then
    write(6,'(" histwrt5 ",a20,i4,a7)') sname,iarch,"missing"
  else
    write(6,'(" histwrt5 ",a20,i4)') sname,iarch
  end if
end if

return
end subroutine hw5l      

! global write with single file
subroutine hw5a(var,sname,idnc,iarch)

use cc_mpi              ! CC MPI routines
use newmpar_m           ! Grid parameters
use parm_m              ! Model configuration

implicit none

integer, intent(in) :: idnc, iarch
integer ier, iq, k, kk, l, ll
integer(kind=4) mid, vtype, lidnc, ndims
integer(kind=4), dimension(5) :: start, ncount
real varn, varx
real, dimension(:,:,:), intent(in) :: var
real, dimension(:,:,:), allocatable :: globvar
real(kind=4) laddoff, lscale_f
character(len=*), intent(in) :: sname
integer(kind=2), dimension(:,:,:), allocatable :: ipack
      
kk = size(var,2)
ll = size(var,3)
allocate( globvar(ifull_g,kk,ll) )
call ccmpi_gather(var(1:ifull,1:kk,1:ll), globvar(1:ifull_g,1:kk,1:ll))
start = (/ 1, 1, 1, 1, iarch /)
ncount = (/ il_g, jl_g, kk, ll, 1 /)

!     find variable index
lidnc = idnc
ier = nf90_inq_varid(lidnc,sname,mid)
call ncmsg(sname,ier)
ier = nf90_inquire_variable(lidnc, mid, xtype=vtype, ndims=ndims)
if ( vtype==nf90_short ) then
  allocate( ipack(ifull_g,kk,ll) )
  if ( all(globvar>9.8e36) )then
    ipack = missval
  else
    ier = nf90_get_att(lidnc,mid,'add_offset',laddoff)
    ier = nf90_get_att(lidnc,mid,'scale_factor',lscale_f)
    do l = 1,ll
      do k = 1,kk
        do iq = 1,ifull_g
          ipack(iq,k,l) = nint(max(min((globvar(iq,k,l)-real(laddoff))/real(lscale_f),real(maxv)),real(minv)),2)
        end do  
      end do
    end do
  end if
  ier = nf90_put_var(lidnc,mid,ipack,start=start(1:ndims),count=ncount(1:ndims))
  deallocate( ipack )
else
  ier = nf90_put_var(lidnc,mid,globvar,start=start(1:ndims),count=ncount(1:ndims))
endif
call ncmsg(sname,ier)

if ( mod(ktau,nmaxpr)==0 ) then
  if ( any(abs(globvar-real(nf90_fill_float))<1.e-20) ) then
    write(6,'(" histwrt5 ",a7,i4,a7)') sname,iarch,"missing"
  else
    varn = minval(globvar)
    varx = maxval(globvar)
    write(6,'(" histwrt5 ",a20,i4,2f12.4)') sname,iarch,varn,varx
  end if
end if

deallocate( globvar )

return
end subroutine hw5a      

! procformat and local(write)
subroutine hw5lpr8(var,sname,idnc,iarch)

use cc_mpi               ! CC MPI routines
use newmpar_m            ! Grid parameters
use parm_m               ! Model configuration

implicit none

integer, intent(in) :: idnc, iarch
integer iq, k, ier, v, kk, l, ll
integer(kind=4) mid, vtype, lidnc, ndims
integer(kind=4), dimension(6) :: start, ncount
real(kind=8), dimension(:,:,:), intent(in) :: var
real(kind=8), dimension(ifull,size(var,2),size(var,3),vnode_nproc) :: var_g
real(kind=4) laddoff, lscale_f
character(len=*), intent(in) :: sname
integer(kind=2), dimension(ifull,size(var,2),size(var,3),vnode_nproc) :: ipack_g

kk = size(var,2)
ll = size(var,3)
start = (/ 1, 1, 1, 1, 1, iarch /)
ncount = (/ il, jl, kk, ll, vnode_nproc, 1 /)

!if ( useiobuffer ) then
!  ! MJT notes - move this to its own subroutine ...  
!  !call add_iobuffer(idnc,mid,ndims,ifull,istep,vnode_nproc,start,ncount,var)
!  write(6,*) "ERROR: iobuffer not yet implemented"
!  call ccmpi_abort(-1)
!  return
!end if

call ccmpi_gatherxr8(var_g,var,0,comm_vnode)

lidnc = idnc
ier = nf90_inq_varid(lidnc,sname,mid)
call ncmsg(sname,ier)
ier = nf90_inquire_variable(lidnc,mid,xtype=vtype,ndims=ndims)

if ( vtype==nf90_short ) then
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
else
  ier = nf90_put_var(lidnc,mid,var_g,start=start(1:ndims),count=ncount(1:ndims))
end if
call ncmsg(sname,ier)

if ( mod(ktau,nmaxpr)==0 .and. myid==0 ) then
  if ( any(abs(var-real(nf90_fill_float,8))<1.e-20_8) ) then
    write(6,'(" histwrt5r8 ",a20,i4,a7)') sname,iarch,"missing"
  else
    write(6,'(" histwrt5r8 ",a20,i4)') sname,iarch
  end if
end if

return
end subroutine hw5lpr8    

! procformat without local(write)
subroutine hw5pr8(var)

use cc_mpi               ! CC MPI routines
use newmpar_m            ! Grid parameters
use parm_m               ! Model configuration

implicit none

real(kind=8), dimension(:,:,:), intent(in) :: var
real(kind=8), dimension(0,0,0,0) :: var_g

!if ( useiobuffer ) then
!  ! MJT notes - move this to its own subroutine ...  
!  !call add_iobuffer(idnc,mid,ndims,ifull,istep,vnode_nproc,start,ncount,var)
!  write(6,*) "ERROR: iobuffer not yet implemented"
!  call ccmpi_abort(-1)
!end if

call ccmpi_gatherxr8(var_g,var,0,comm_vnode)

return
end subroutine hw5pr8      

! pure local(write)
subroutine hw5lr8(var,sname,idnc,iarch)

use cc_mpi               ! CC MPI routines
use newmpar_m            ! Grid parameters
use parm_m               ! Model configuration

implicit none

integer, intent(in) :: idnc, iarch
integer iq, k, ier, kk, l, ll
integer(kind=4) mid, vtype, lidnc, ndims
integer(kind=4), dimension(5) :: start, ncount
real(kind=8), dimension(:,:,:), intent(in) :: var
real(kind=4) laddoff, lscale_f
character(len=*), intent(in) :: sname
integer(kind=2), dimension(ifull,size(var,2),size(var,3)) :: ipack

kk = size(var,2)
ll = size(var,3)
start = (/ 1, 1, 1, 1, iarch /)
ncount = (/ il, jl, kk, ll, 1 /)

lidnc = idnc
ier = nf90_inq_varid(lidnc,sname,mid)
call ncmsg(sname,ier)
ier = nf90_inquire_variable(lidnc,mid,xtype=vtype,ndims=ndims)
if ( vtype==nf90_short ) then
  if ( all(var>9.8e36_8) ) then
    ipack = missval
  else
    ier = nf90_get_att(lidnc,mid,'add_offset',laddoff)
    ier = nf90_get_att(lidnc,mid,'scale_factor',lscale_f)
    do l = 1,ll
      do k = 1,kk
        do iq = 1,ifull
          ipack(iq,k,l) = nint(max(min((var(iq,k,l)-real(laddoff,8))/real(lscale_f,8),real(maxv,8)),real(minv,8)),2)
        end do  
      end do
    end do
  end if
  ier = nf90_put_var(lidnc,mid,ipack,start=start(1:ndims),count=ncount(1:ndims))
else
  ier = nf90_put_var(lidnc,mid,var,start=start(1:ndims),count=ncount(1:ndims))
endif
call ncmsg(sname,ier)

if ( mod(ktau,nmaxpr)==0 .and. myid==0 ) then
  if ( any(abs(var-real(nf90_fill_float,8))<1.e-20_8) ) then
    write(6,'(" histwrt5r8 ",a20,i4,a7)') sname,iarch,"missing"
  else
    write(6,'(" histwrt5r8 ",a20,i4)') sname,iarch
  end if
end if

return
end subroutine hw5lr8      

! global write with single file
subroutine hw5ar8(var,sname,idnc,iarch)

use cc_mpi              ! CC MPI routines
use newmpar_m           ! Grid parameters
use parm_m              ! Model configuration

implicit none

integer, intent(in) :: idnc, iarch
integer ier, iq, k, kk, l, ll
integer(kind=4) mid, vtype, lidnc, ndims
integer(kind=4), dimension(5) :: start, ncount
real(kind=8) varn, varx
real(kind=8), dimension(:,:,:), intent(in) :: var
real(kind=8), dimension(:,:,:), allocatable :: globvar
real(kind=4) laddoff, lscale_f
character(len=*), intent(in) :: sname
integer(kind=2), dimension(:,:,:), allocatable :: ipack
      
kk = size(var,2)
ll = size(var,3)
allocate( globvar(ifull_g,kk,ll) )
call ccmpi_gatherr8(var(1:ifull,1:kk,1:ll), globvar(1:ifull_g,1:kk,1:ll))
start = (/ 1, 1, 1, 1, iarch /)
ncount = (/ il_g, jl_g, kk, ll, 1 /)

!     find variable index
lidnc = idnc
ier = nf90_inq_varid(lidnc,sname,mid)
call ncmsg(sname,ier)
ier = nf90_inquire_variable(lidnc, mid, xtype=vtype, ndims=ndims)
if ( vtype==nf90_short ) then
  allocate( ipack(ifull_g,kk,ll) )
  if ( all(globvar>9.8e36_8) )then
    ipack = missval
  else
    ier = nf90_get_att(lidnc,mid,'add_offset',laddoff)
    ier = nf90_get_att(lidnc,mid,'scale_factor',lscale_f)
    do l = 1,ll
      do k = 1,kk
        do iq = 1,ifull_g
          ipack(iq,k,l) = nint(max(min((globvar(iq,k,l)-real(laddoff,8))/real(lscale_f,8),real(maxv,8)),real(minv,8)),2)
        end do  
      end do
    end do
  end if
  ier = nf90_put_var(lidnc,mid,ipack,start=start(1:ndims),count=ncount(1:ndims))
  deallocate( ipack )
else
  ier = nf90_put_var(lidnc,mid,globvar,start=start(1:ndims),count=ncount(1:ndims))
endif
call ncmsg(sname,ier)

if ( mod(ktau,nmaxpr)==0 ) then
  if ( any(abs(globvar-real(nf90_fill_float,8))<1.e-20_8) ) then
    write(6,'(" histwrt5r8 ",a20,i4,a7)') sname,iarch,"missing"
  else
    varn = minval(globvar)
    varx = maxval(globvar)
    write(6,'(" histwrt5r8 ",a20,i4,2f12.4)') sname,iarch,varn,varx
  end if
end if

deallocate( globvar )

return
end subroutine hw5ar8

!--------------------------------------------------------------------

subroutine ccnf_open(fname,ncid,status)

use cc_mpi

implicit none

integer, intent(out) :: ncid
integer, intent(out), optional :: status
integer ncstatus
integer(kind=4) lncid
character(len=*), intent(in) :: fname

ncstatus = nf90_open(fname,nf90_nowrite,lncid)
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
use parm_m

implicit none

integer, intent(out) :: ncid
integer ncstatus
integer(kind=4) :: lncid, lcmode
character(len=*), intent(in) :: fname

#ifdef usenc3
if ( synchist ) then
  lcmode = ior(nf90_64bit_offset, nf90_share)
else
  lcmode = nf90_64bit_offset
end if
ncstatus = nf90_create(fname,lcmode,lncid)
#else
if ( synchist ) then
  lcmode = ior(nf90_netcdf4, nf90_share)
else
  lcmode = nf90_netcdf4
end if
ncstatus = nf90_create(fname,lcmode,lncid)
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
ncstatus = nf90_close(lncid)
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
ncstatus = nf90_set_fill(lncid,nf90_nofill,lomode)
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
ncstatus=nf90_inq_dimid(lncid,dname,ldid)
ltst=(ncstatus/=nf90_noerr)
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
ncstatus = nf90_inq_varid(lncid,vname,lvid)
ltst=(ncstatus/=nf90_noerr)
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
ncstatus = nf90_inquire_variable(lncid,lvid,ndims=lndims)
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
ncstatus=nf90_def_dim(lncid,dname,lnsize,ldid)
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
ncstatus=nf90_def_dim(lncid,dname,nf90_unlimited,ldid)
did=ldid
call ncmsg("def_dimu",ncstatus)

return
end subroutine ccnf_def_dimu

subroutine ccnf_def_var_v(ncid,vname,vtype,vndim,dims,vid)

use cc_mpi

implicit none

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
    ltype=nf90_int
  case('float')
    ltype=nf90_float
  case('double')
    ltype=nf90_double
  case('char')
    ltype=nf90_char
  case default
    write(6,*) "ERROR: Unknown option for ccnf_def_var ",vtype
    call ccmpi_abort(-1)
end select

lncid=ncid
ldims=dims
#ifdef usenc3
ncstatus = nf90_def_var(lncid,vname,ltype,ldims,lvid)
#else
ncstatus = nf90_def_var(lncid,vname,ltype,ldims,lvid,deflate_level=1_4)
#endif
vid=lvid
call ncmsg("def_var",ncstatus)

return
end subroutine ccnf_def_var_v

subroutine ccnf_def_var_s(ncid,vname,vtype,vid)

use cc_mpi

implicit none

integer, intent(in) :: ncid
integer, intent(out) :: vid
integer ncstatus
integer(kind=4) lncid, lvid, ltype
character(len=*), intent(in) :: vname
character(len=*), intent(in) :: vtype

select case(vtype)
  case('int')
    ltype=nf90_int
  case('float')
    ltype=nf90_float
  case('double')
    ltype=nf90_double
  case('char')
    ltype=nf90_char
  case default
    write(6,*) "ERROR: Unknown option for ccnf_def_var ",vtype
    call ccmpi_abort(-1)
end select

lncid=ncid
ncstatus = nf90_def_var(lncid,vname,ltype,lvid)
vid=lvid
call ncmsg("def_var0",ncstatus)

return
end subroutine ccnf_def_var_s

subroutine ccnf_enddef(ncid)

use cc_mpi

implicit none

integer, intent(in) :: ncid
integer ncstatus
integer(kind=4) lncid

lncid=ncid
ncstatus=nf90_enddef(lncid)
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
ncstatus=nf90_redef(lncid)
call ncmsg("redef",ncstatus)

return
end subroutine ccnf_redef

subroutine ccnf_sync(ncid)

use cc_mpi

implicit none

integer, intent(in) :: ncid
integer ncstatus
integer(kind=4) lncid

lncid = ncid
ncstatus = nf90_sync(lncid)
call ncmsg("sync",ncstatus)

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
ncstatus = nf90_inq_varid(lncid,vname,lvid)
call ncmsg("get_var_varid",ncstatus)
ncstatus = nf90_get_var(lncid,lvid,vdat)
call ncmsg("get_var",ncstatus)

return
end subroutine ccnf_get_var_real

subroutine ccnf_get_vara_text_t(ncid,name,start,ncount,vdat)

use cc_mpi

implicit none

integer, intent(in) :: ncid
integer, dimension(:), intent(in) :: start, ncount
integer ncstatus
integer(kind=4) lncid, lvid
integer(kind=4), dimension(size(start)) :: lstart
integer(kind=4), dimension(size(ncount)) :: lncount
character(len=*), intent(out) :: vdat
character(len=*), intent(in) :: name

lncid=ncid
ncstatus=nf90_inq_varid(lncid,name,lvid)
call ncmsg(name,ncstatus)
lstart(:)=start(:)
lncount(:)=ncount(:)
ncstatus=nf90_get_var(lncid,lvid,vdat,start=lstart,count=lncount)
call ncmsg("get_vara_text_t",ncstatus)

return
end subroutine ccnf_get_vara_text_t 

subroutine ccnf_get_var_int(ncid,vname,vdat)

use cc_mpi

implicit none

integer, intent(in) :: ncid
integer ncstatus
integer, dimension(:), intent(out) :: vdat
integer(kind=4) :: lncid, lvid
character(len=*), intent(in) :: vname

lncid=ncid
ncstatus = nf90_inq_varid(lncid,vname,lvid)
call ncmsg("get_var_varid",ncstatus)
ncstatus = nf90_get_var(lncid,lvid,vdat)
call ncmsg("get_var",ncstatus)

return
end subroutine ccnf_get_var_int

subroutine ccnf_get_vara_real1r_s(ncid,vid,start,vdat)

use cc_mpi

implicit none

integer, intent(in) :: ncid, vid, start
integer ncstatus
integer(kind=4) lncid, lvid
integer(kind=4), dimension(1) :: lstart
integer(kind=4), dimension(1) :: lncount
real, intent(out) :: vdat
real, dimension(1) :: lvdat

lncid=ncid
lvid=vid
lstart=start
lncount=1
ncstatus=nf90_get_var(lncid,lvid,lvdat,start=lstart,count=lncount)
call ncmsg("get_vara_real1r",ncstatus)
vdat=lvdat(1)

return
end subroutine ccnf_get_vara_real1r_s

subroutine ccnf_get_vara_real1r_t(ncid,name,start,vdat)

use cc_mpi

implicit none

integer, intent(in) :: ncid, start
integer ncstatus
integer(kind=4) lncid, lvid
integer(kind=4), dimension(1) :: lstart
integer(kind=4), dimension(1) :: lncount
real, intent(out) :: vdat
real, dimension(1) :: lvdat
character(len=*), intent(in) :: name

lncid=ncid
ncstatus=nf90_inq_varid(lncid,name,lvid)
call ncmsg(name,ncstatus)
lstart=start
lncount=1
ncstatus=nf90_get_var(lncid,lvid,lvdat,start=lstart,count=lncount)
call ncmsg("get_vara_real1r",ncstatus)
vdat=lvdat(1)

return
end subroutine ccnf_get_vara_real1r_t

subroutine ccnf_get_vara_real2r_s(ncid,vid,start,ncount,vdat)

use cc_mpi

implicit none

integer, intent(in) :: ncid, vid
integer, intent(in) :: start, ncount
integer ncstatus
integer(kind=4) lncid, lvid
integer(kind=4), dimension(1) :: lstart
integer(kind=4), dimension(1) :: lncount
real, dimension(:), intent(out) :: vdat

lncid=ncid
lvid=vid
lstart=start
lncount=ncount
ncstatus=nf90_get_var(lncid,lvid,vdat,start=lstart,count=lncount)
call ncmsg("get_vara_real2r",ncstatus)

return
end subroutine ccnf_get_vara_real2r_s 

subroutine ccnf_get_vara_real2r_t(ncid,name,start,ncount,vdat)

use cc_mpi

implicit none

integer, intent(in) :: ncid
integer, dimension(:), intent(in) :: start, ncount
integer ncstatus
integer(kind=4) lncid, lvid
integer(kind=4), dimension(size(start)) :: lstart
integer(kind=4), dimension(size(ncount)) :: lncount
real, dimension(:), intent(out) :: vdat
character(len=*), intent(in) :: name

lncid=ncid
ncstatus=nf90_inq_varid(lncid,name,lvid)
call ncmsg(name,ncstatus)
lstart(:)=start(:)
lncount(:)=ncount(:)
ncstatus=nf90_get_var(lncid,lvid,vdat,start=lstart,count=lncount)
call ncmsg("get_vara_real2r_t",ncstatus)

return
end subroutine ccnf_get_vara_real2r_t 

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
ncstatus=nf90_get_var(lncid,lvid,vdat,start=lstart,count=lncount)
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
ncstatus=nf90_get_var(lncid,lvid,vdat,start=lstart,count=lncount)
call ncmsg("get_vara_real3r",ncstatus)

return
end subroutine ccnf_get_vara_real3r

subroutine ccnf_get_vara_real3r_t(ncid,name,start,ncount,vdat)

use cc_mpi

implicit none

integer, intent(in) :: ncid
integer, dimension(:), intent(in) :: start, ncount
integer ncstatus
integer(kind=4) lncid, lvid
integer(kind=4), dimension(size(start)) :: lstart
integer(kind=4), dimension(size(ncount)) :: lncount
real, dimension(:,:), intent(out) :: vdat
character(len=*), intent(in) :: name

lncid=ncid
ncstatus=nf90_inq_varid(lncid,name,lvid)
call ncmsg(name,ncstatus)
lstart(:)=start(:)
lncount(:)=ncount(:)
ncstatus=nf90_get_var(lncid,lvid,vdat,start=lstart,count=lncount)
call ncmsg("get_vara_real3r_t",ncstatus)

return
end subroutine ccnf_get_vara_real3r_t 

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
ncstatus=nf90_get_var(lncid,lvid,vdat,start=lstart,count=lncount)
call ncmsg("get_vara_real4r",ncstatus)

return
end subroutine ccnf_get_vara_real4r

subroutine ccnf_get_vara_int1i_s(ncid,vid,start,vdat)

use cc_mpi

implicit none

integer, intent(in) :: ncid, vid, start
integer ncstatus
integer, intent(out) :: vdat
integer(kind=4) :: lncid, lvid
integer(kind=4), dimension(1) :: lstart
integer(kind=4), dimension(1) :: lncount
integer(kind=4), dimension(1) :: lvdat

lncid=ncid
lvid=vid
lstart=start
lncount=1
ncstatus=nf90_get_var(lncid,lvid,lvdat,start=lstart,count=lncount)
call ncmsg("get_vara_int1i",ncstatus)
vdat=lvdat(1)

return
end subroutine ccnf_get_vara_int1i_s

subroutine ccnf_get_vara_int2i_t(ncid,name,start,ncount,vdat)

use cc_mpi

implicit none

integer, intent(in) :: ncid
integer, dimension(:) :: start, ncount
integer ncstatus
integer, dimension(:), intent(out) :: vdat
integer(kind=4) :: lncid, lvid
integer(kind=4), dimension(size(start)) :: lstart
integer(kind=4), dimension(size(ncount)) :: lncount
character(len=*), intent(in) :: name

lncid=ncid
ncstatus=nf90_inq_varid(lncid,name,lvid)
lstart(:)=start(:)
lncount(:)=ncount(:)
ncstatus=nf90_get_var(lncid,lvid,vdat,start=lstart,count=lncount)
call ncmsg("get_vara_int2i_t",ncstatus)

return
end subroutine ccnf_get_vara_int2i_t

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

lncid=ncid
lvid=vid
lstart=start
lncount=ncount
ncstatus=nf90_get_var(lncid,lvid,vdat,start=lstart,count=lncount)
call ncmsg("get_vara_int2i",ncstatus)

return
end subroutine ccnf_get_vara_int2i

subroutine ccnf_get_vara_int3i_t(ncid,name,start,ncount,vdat)

use cc_mpi

implicit none

integer, intent(in) :: ncid
integer, dimension(:) :: start, ncount
integer ncstatus
integer, dimension(:,:), intent(out) :: vdat
integer(kind=4) :: lncid, lvid
integer(kind=4), dimension(size(start)) :: lstart
integer(kind=4), dimension(size(ncount)) :: lncount
character(len=*), intent(in) :: name

lncid=ncid
ncstatus=nf90_inq_varid(lncid,name,lvid)
lstart(:)=start(:)
lncount(:)=ncount(:)
ncstatus=nf90_get_var(lncid,lvid,vdat,start=lstart,count=lncount)
call ncmsg("get_vara_int3i_t",ncstatus)

return
end subroutine ccnf_get_vara_int3i_t

#ifndef i8r8
subroutine ccnf_get_vara_double2r_t(ncid,name,start,ncount,vdat)

use cc_mpi

implicit none

integer, intent(in) :: ncid
integer, dimension(:), intent(in) :: start, ncount
integer ncstatus
integer(kind=4) lncid, lvid
integer(kind=4), dimension(size(start)) :: lstart
integer(kind=4), dimension(size(ncount)) :: lncount
real(kind=8), dimension(:), intent(out) :: vdat
character(len=*), intent(in) :: name

lncid=ncid
ncstatus=nf90_inq_varid(lncid,name,lvid)
call ncmsg(name,ncstatus)
lstart(:)=start(:)
lncount(:)=ncount(:)
ncstatus=nf90_get_var(lncid,lvid,vdat,start=lstart,count=lncount)
call ncmsg("get_vara_double2r_t",ncstatus)

return
end subroutine ccnf_get_vara_double2r_t 

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
ncstatus=nf90_get_var(lncid,lvid,vdat,start=lstart,count=lncount)
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
ncstatus = nf90_get_att(lncid,lvid,aname,atext)
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
ncstatus = nf90_get_att(lncid,lvid,aname,vdat)
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

lncid = ncid
ncstatus = nf90_get_att(lncid,nf90_global,aname,vdat)
if (present(ierr)) then
  ierr = ncstatus
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

lncid=ncid
ncstatus = nf90_get_att(lncid,nf90_global,aname,vdat(:))
call ncmsg("get_attg",ncstatus)

return
end subroutine ccnf_get_att_realg2r

subroutine ccnf_get_att_intg1i(ncid,aname,vdat,tst)

use cc_mpi

implicit none

integer, intent(in) :: ncid
integer, intent(out) :: vdat
integer ncstatus
integer(kind=4) :: lncid
integer(kind=4) :: lvdat
logical, intent(out), optional :: tst
character(len=*), intent(in) :: aname

lncid=ncid
ncstatus = nf90_get_att(lncid,nf90_global,aname,lvdat)
vdat=lvdat
if (present(tst)) then
  tst=ncstatus/=nf90_noerr
else
  call ncmsg("get_attg",ncstatus)
end if

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
ncstatus = nf90_get_att(lncid,nf90_global,aname,lvdat)
vdat=lvdat
call ncmsg("get_attg",ncstatus)

return
end subroutine ccnf_get_att_intg2i

subroutine ccnf_read(fname,vname,vdat)

use cc_mpi
use newmpar_m

implicit none

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

subroutine ccnf_put_var_text2r(ncid,vid,vtxt)

use cc_mpi

implicit none

integer, intent(in) :: ncid, vid
integer ncstatus
integer(kind=4) lncid, lvid
character(len=*), dimension(:), intent(in) :: vtxt

lncid=ncid
lvid=vid
ncstatus = nf90_put_var(lncid,lvid,vtxt)
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

lncid = ncid
lvid = vid
ncstatus = nf90_put_var(lncid,lvid,vdat)
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

lncid = ncid
lvid = vid
ncstatus = nf90_put_var(lncid,lvid,vdat)
call ncmsg("put_var",ncstatus)

return
end subroutine ccnf_put_var_int3i

subroutine ccnf_put_vara_real1r_t(ncid,name,start,vdat)

use cc_mpi

implicit none

integer, intent(in) :: ncid
integer ncstatus
integer, intent(in) :: start
integer(kind=4) lncid, lvid
integer(kind=4), dimension(1) :: lstart
integer(kind=4), dimension(1) :: lncount
real, intent(in) :: vdat
real, dimension(1) :: lvdat
character(len=*), intent(in) :: name

lncid=ncid
ncstatus=nf90_inq_varid(lncid,name,lvid)
call ncmsg(name,ncstatus)
lstart=start
lncount=1
lvdat=vdat
ncstatus=nf90_put_var(lncid,lvid,lvdat,start=lstart,count=lncount)
call ncmsg("put_vara_real1r",ncstatus)

return
end subroutine ccnf_put_vara_real1r_t

subroutine ccnf_put_vara_real1r_s(ncid,vid,start,vdat)

use cc_mpi

implicit none

integer, intent(in) :: ncid, vid
integer ncstatus
integer, intent(in) :: start
integer(kind=4) lncid, lvid
integer(kind=4), dimension(1) :: lstart
integer(kind=4), dimension(1) :: lncount
real, intent(in) :: vdat
real, dimension(1) :: lvdat

lncid=ncid
lvid=vid
lstart=start
lncount=1
lvdat=vdat
ncstatus=nf90_put_var(lncid,lvid,lvdat,start=lstart,count=lncount)
call ncmsg("put_vara_real1r",ncstatus)

return
end subroutine ccnf_put_vara_real1r_s

subroutine ccnf_put_vara_real2r_t(ncid,name,start,ncount,vdat)

use cc_mpi

implicit none

integer, intent(in) :: ncid
integer ncstatus
integer, intent(in) :: start, ncount
integer(kind=4) lncid, lvid
integer(kind=4), dimension(1) :: lstart
integer(kind=4), dimension(1) :: lncount
real, dimension(:), intent(in) :: vdat
character(len=*), intent(in) :: name

lncid=ncid
ncstatus=nf90_inq_varid(lncid,name,lvid)
call ncmsg("put_vara_real2r",ncstatus)
lstart=start
lncount=ncount
ncstatus=nf90_put_var(lncid,lvid,vdat,start=lstart,count=lncount)
call ncmsg("put_vara_real2r",ncstatus)

return
end subroutine ccnf_put_vara_real2r_t

subroutine ccnf_put_vara_real2r_s(ncid,vid,start,ncount,vdat)

use cc_mpi

implicit none

integer, intent(in) :: ncid, vid
integer ncstatus
integer, intent(in) :: start, ncount
integer(kind=4) lncid, lvid
integer(kind=4), dimension(1) :: lstart
integer(kind=4), dimension(1) :: lncount
real, dimension(:), intent(in) :: vdat

lncid=ncid
lvid=vid
lstart=start
lncount=ncount
ncstatus=nf90_put_var(lncid,lvid,vdat,start=lstart,count=lncount)
call ncmsg("put_vara_real2r",ncstatus)

return
end subroutine ccnf_put_vara_real2r_s

subroutine ccnf_put_vara_real3r_t(ncid,name,start,ncount,vdat)

use cc_mpi

implicit none

integer, intent(in) :: ncid
integer ncstatus
integer, dimension(2), intent(in) :: start, ncount
integer(kind=4) lncid, lvid
integer(kind=4), dimension(2) :: lstart
integer(kind=4), dimension(2) :: lncount
real, dimension(:,:), intent(in) :: vdat
character(len=*), intent(in) :: name

lncid=ncid
ncstatus=nf90_inq_varid(lncid,name,lvid)
call ncmsg("put_vara_real2r",ncstatus)
lstart(:)=start(:)
lncount(:)=ncount(:)
ncstatus=nf90_put_var(lncid,lvid,vdat,start=lstart,count=lncount)
call ncmsg("put_vara_real2r",ncstatus)

return
end subroutine ccnf_put_vara_real3r_t

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
ncstatus=nf90_put_var(lncid,lvid,vdat,start=lstart,count=lncount)
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
ncstatus=nf90_put_var(lncid,lvid,vdat,start=lstart,count=lncount)
call ncmsg("put_vara",ncstatus)

return
end subroutine ccnf_put_vara_real3r

#ifndef i8r8
subroutine ccnf_put_vara_double1r_s(ncid,vid,start,vdat)

use cc_mpi

implicit none

integer, intent(in) :: ncid, vid, start
integer ncstatus
integer(kind=4) lncid, lvid
integer(kind=4), dimension(1) :: lstart
integer(kind=4), dimension(1) :: lncount
real(kind=8), intent(in) :: vdat
real(kind=8), dimension(1) :: lvdat

lncid=ncid
lvid=vid
lstart=start
lncount=1
lvdat=vdat
ncstatus=nf90_put_var(lncid,lvid,lvdat,start=lstart,count=lncount)
call ncmsg("put_vara",ncstatus)

return
end subroutine ccnf_put_vara_double1r_s
#endif

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
ncstatus=nf90_put_var(lncid,lvid,vdat,start=lstart,count=lncount)
call ncmsg("put_vara",ncstatus)

return
end subroutine ccnf_put_vara_double2r
#endif

subroutine ccnf_put_vara_int1i_t(ncid,name,start,vdat)

use cc_mpi

implicit none

integer, intent(in) :: ncid
integer ncstatus
integer, intent(in) :: start
integer(kind=4) lncid, lvid
integer(kind=4), dimension(1) :: lstart
integer(kind=4), dimension(1) :: lncount
integer, intent(in) :: vdat
integer, dimension(1) :: lvdat
character(len=*), intent(in) :: name

lncid=ncid
ncstatus=nf90_inq_varid(lncid,name,lvid)
call ncmsg("put_vara_int1i",ncstatus)
lstart=start
lncount=1
lvdat=vdat
ncstatus=nf90_put_var(lncid,lvid,lvdat,start=lstart,count=lncount)
call ncmsg("put_vara_int1i",ncstatus)

return
end subroutine ccnf_put_vara_int1i_t

subroutine ccnf_put_vara_int1i_s(ncid,vid,start,vdat)

use cc_mpi

implicit none

integer, intent(in) :: ncid, vid
integer ncstatus
integer, intent(in) :: start
integer(kind=4) lncid, lvid
integer(kind=4), dimension(1) :: lstart
integer(kind=4), dimension(1) :: lncount
integer, intent(in) :: vdat
integer, dimension(1) :: lvdat

lncid=ncid
lvid=vid
lstart=start
lncount=1
lvdat=vdat
ncstatus=nf90_put_var(lncid,lvid,lvdat,start=lstart,count=lncount)
call ncmsg("put_vara_int1i",ncstatus)

return
end subroutine ccnf_put_vara_int1i_s

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

lncid=ncid
lvid=vid
lstart=start
lncount=ncount
ncstatus=nf90_put_var(lncid,lvid,vdat,start=lstart,count=lncount)
call ncmsg("put_vara",ncstatus)

return
end subroutine ccnf_put_vara_int2i

subroutine ccnf_put_att_text(ncid,vid,aname,atext)

use cc_mpi

implicit none

integer, intent(in) :: ncid, vid
integer ncstatus
integer(kind=4) lncid, lvid
character(len=*), intent(in) :: aname
character(len=*), intent(in) :: atext

lncid=ncid
lvid=vid
ncstatus=nf90_put_att(lncid,lvid,aname,atext)
call ncmsg("put_att",ncstatus)

return
end subroutine ccnf_put_att_text

subroutine ccnf_put_att_rs(ncid,vid,aname,rval)

use cc_mpi

implicit none

integer, intent(in) :: ncid, vid
integer ncstatus
integer(kind=4) lncid, lvid
real, intent(in) :: rval
character(len=*), intent(in) :: aname

lncid=ncid
lvid=vid
ncstatus=nf90_put_att(lncid,lvid,aname,rval)
call ncmsg("put_att",ncstatus)

return
end subroutine ccnf_put_att_rs

subroutine ccnf_put_att_textg(ncid,aname,atext)

use cc_mpi

implicit none

integer, intent(in) :: ncid
integer ncstatus
integer(kind=4) lncid
character(len=*), intent(in) :: aname
character(len=*), intent(in) :: atext

lncid=ncid
ncstatus=nf90_put_att(lncid,nf90_global,aname,atext)
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
ncstatus=nf90_put_att(lncid,nf90_global,aname,ldat)
call ncmsg("put_attg",ncstatus)

return
end subroutine ccnf_put_att_intg1

subroutine ccnf_put_att_intg2(ncid,aname,vdat)

use cc_mpi

implicit none

integer, intent(in) :: ncid
integer ncstatus
integer, dimension(:), intent(in) :: vdat
integer(kind=4) lncid
integer(kind=4), dimension(size(vdat)) :: lvdat
character(len=*), intent(in) :: aname

lncid=ncid
lvdat=vdat
ncstatus=nf90_put_att(lncid,nf90_global,aname,lvdat)
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
ncstatus=nf90_put_att(lncid,nf90_global,aname,ldat)
call ncmsg("put_attg",ncstatus)

return
end subroutine ccnf_put_att_realg1

subroutine ccnf_put_att_realg2(ncid,aname,vdat)

use cc_mpi

implicit none

integer, intent(in) :: ncid
integer ncstatus
integer(kind=4) lncid
real, dimension(:), intent(in) :: vdat
character(len=*), intent(in) :: aname

lncid=ncid
ncstatus=nf90_put_att(lncid,nf90_global,aname,vdat)
call ncmsg("put_attg",ncstatus)

return
end subroutine ccnf_put_att_realg2

subroutine surfread(dat,varname,netcdfid,filename)

use cc_mpi
use newmpar_m

implicit none

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

implicit none

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

ifull_l=size(dat)

allocate( glob2d(ifull_g) )

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
  if(ilx/=il_g.or.jlx/=jl_g.or.abs(rlong0x-rlong0)>=1.e-20.or.abs(rlat0x-rlat0)>=1.e-20.or. &
      abs(schmidtx-schmidt)>=1.e-20) then
    write(6,*) 'wrong data file supplied'
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

implicit none

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

subroutine file_distribute2(rvar,gvar)

use cc_mpi

implicit none

! Convert standard 1D arrays to face form and distribute to processors
real, dimension(:), intent(out) :: rvar
real, dimension(:), intent(in), optional :: gvar

call START_LOG(distribute_begin)

! Copy internal region
if ( myid==0 ) then
  if ( .not.present(gvar) ) then
    write(6,*) "Error: file_distribute argument required on proc 0"
    call ccmpi_abort(-1)
  end if
  call host_filedistribute2(rvar,gvar)
else
  call proc_filedistribute2(rvar)
end if

call END_LOG(distribute_end)

return
end subroutine file_distribute2

subroutine host_filedistribute2(rvar,gvar)

use cc_mpi

implicit none

integer ipf, jpf, ip, n, fsize, no, ca, cc, j
real, dimension(:), intent(out) :: rvar
real, dimension(:), intent(in) :: gvar
real, dimension(:,:), allocatable :: bufvar

fsize = pil*pjl*pnpan

allocate( bufvar(fsize,fnresid) )

! map array in order of processor rank
do ipf = 0,mynproc-1
  do jpf = 1,fnresid
    ip = ipf*fnresid + jpf - 1
    do n = 0,pnpan-1
      no = n - pnoff(ip) + 1
      ca = pioff(ip,no) + (pjoff(ip,no)-1)*pil_g + no*pil_g*pil_g
      cc = n*pil*pjl - pil
      do j = 1,pjl
        bufvar(1+j*pil+cc:pil+j*pil+cc,jpf) = gvar(1+j*pil_g+ca:pil+j*pil_g+ca)
      end do
    end do
  end do    
  ca = ipf*fsize
  call ccmpi_scatterx(bufvar,rvar(1+ca:fsize+ca),0,comm_ip)
end do

deallocate( bufvar )

return
end subroutine host_filedistribute2

subroutine proc_filedistribute2(rvar)
   
use cc_mpi

implicit none

integer ipf, fsize, ca
real, dimension(:), intent(out) :: rvar
real, dimension(0,0) :: bufvar

fsize = pil*pjl*pnpan

do ipf = 0,mynproc-1
  ca = ipf*fsize
  call ccmpi_scatterx(bufvar,rvar(1+ca:fsize+ca),0,comm_ip)
end do

return
end subroutine proc_filedistribute2


end module infile
