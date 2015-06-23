#ifdef ncclib
! C interface
module netcdf_m

use, intrinsic :: ISO_C_BINDING, only: C_SHORT, C_INT, C_FLOAT, C_DOUBLE, C_SIZE_T, C_LOC, C_NULL_CHAR, C_PTR, &
                                       C_F_POINTER

implicit none

private
public nf_unlimited
public nf_noerr, nf_enameinuse
public nf_nowrite, nf_64bit_offset, nf_clobber, nf_nofill, nf_write
public nf_global
public nf_short, nf_int2, nf_int, nf_float, nf_real, nf_double, nf_char
public nf_fill_float
public nf_open, nf_close, nf_create, nf_enddef, nf_set_fill, nf_redef, nf_sync, nf_strerror
public nf__open, nf__create, nf__enddef, nf_abort
public nf_inq_varndims, nf_inq_vardimid, nf_inq_dimlen, nf_inq_varid, nf_inq_dimid, nf_inq_vartype
public nf_inq, nf_inq_varname, nf_inq_ndims, nf_inq_nvars, nf_inq_libvers, nf_inq_dim
public nf_get_att_text, nf_get_att_real, nf_get_att_int
public nf_get_vara_real, nf_get_vara_int, nf_get_vara_int2, nf_get_vara_double, nf_get_var_real
public nf_get_var1_real, nf_get_var1_int
public nf_def_dim, nf_def_var
public nf_rename_dim
public nf_put_att_text, nf_put_att_int2, nf_put_att_real, nf_put_att_int
public nf_put_vara_real, nf_put_vara_int, nf_put_vara_int2, nf_put_vara_double, nf_put_var_real
public nf_put_var1_int, nf_put_var1_real, nf_put_var1_double, nf_put_var1_text
public nf_copy_att

interface

integer (C_INT) function nc_open(path,omode,ncidp) bind(C, name='nc_open')
  use, intrinsic :: ISO_C_BINDING
  implicit none
  integer (C_INT), value :: omode
  type (C_PTR), value :: ncidp
  character, dimension(*) :: path
end function nc_open

integer (C_INT) function nc_close(ncid) bind(C, name='nc_close')
  use, intrinsic :: ISO_C_BINDING
  implicit none
  integer (C_INT), value :: ncid
end function nc_close

integer (C_INT) function nc_create(path,omode,ncidp) bind(C, name='nc_create')
  use, intrinsic :: ISO_C_BINDING
  implicit none
  integer (C_INT), value :: omode
  type (C_PTR), value :: ncidp
  character, dimension(*) :: path
end function nc_create

integer (C_INT) function nc_enddef(ncid) bind(C, name='nc_enddef')
  use, intrinsic :: ISO_C_BINDING
  implicit none
  integer (C_INT), value :: ncid
end function nc_enddef

integer (C_INT) function nc_set_fill(ncid,fillmode,o_modep) bind(C, name='nc_set_fill')
  use, intrinsic :: ISO_C_BINDING
  implicit none
  integer (C_INT), value :: ncid, fillmode
  type (C_PTR), value :: o_modep
end function nc_set_fill

integer (C_INT) function nc_redef(ncid) bind(C, name='nc_redef')
  use, intrinsic :: ISO_C_BINDING
  implicit none
  integer (C_INT), value :: ncid
end function nc_redef

integer (C_INT) function nc_sync(ncid) bind(C, name='nc_sync')
  use, intrinsic :: ISO_C_BINDING
  implicit none
  integer (C_INT), value :: ncid
end function nc_sync

type (C_PTR) function nc_strerror(ncerr) bind(C, name='nc_strerror')
  use, intrinsic :: ISO_C_BINDING
  implicit none
  integer (C_INT), value :: ncerr
end function nc_strerror

integer (C_INT) function nc__open(path,mode,bufrsizehintp,ncidp) bind(C, name='nc__open')
  use, intrinsic :: ISO_C_BINDING
  implicit none
  integer (C_INT), value :: mode
  type (C_PTR), value :: bufrsizehintp
  type (C_PTR), value :: ncidp
  character, dimension(*) :: path
end function nc__open

integer (C_INT) function nc__create(path,cmode,initialsz,bufrsizehintp,ncidp) bind(C, name='nc__create')
  use, intrinsic :: ISO_C_BINDING
  implicit none
  integer (C_INT), value :: cmode
  integer (C_SIZE_T), value :: initialsz
  type (C_PTR), value :: bufrsizehintp
  type (C_PTR), value :: ncidp
  character, dimension(*) :: path
end function nc__create

integer (C_INT) function nc__enddef(ncid,h_minfree,v_align,v_minfree,r_align) bind(C, name='nc__enddef')
  use, intrinsic :: ISO_C_BINDING
  implicit none
  integer (C_INT), value :: ncid
  integer (C_SIZE_T), value :: h_minfree, v_align, v_minfree, r_align
end function nc__enddef

integer (C_INT) function nc_abort(ncid) bind(C, name='nc_abort')
  use, intrinsic :: ISO_C_BINDING
  implicit none
  integer (C_INT), value :: ncid
end function nc_abort

integer (C_INT) function nc_inq_varndims(ncid,varid,ndimsp) bind(C, name='nc_inq_varndims')
  use, intrinsic :: ISO_C_BINDING
  implicit none
  integer (C_INT), value :: ncid, varid
  type (C_PTR), value :: ndimsp
end function nc_inq_varndims

integer (C_INT) function nc_inq_vardimid(ncid,varid,dimidsp) bind(C, name='nc_inq_vardimid')
  use, intrinsic :: ISO_C_BINDING
  implicit none
  integer (C_INT), value :: ncid, varid
  type (C_PTR), value :: dimidsp
end function nc_inq_vardimid

integer (C_INT) function nc_inq_dimlen(ncid,dimid,lengthp) bind(C, name='nc_inq_dimlen')
  use, intrinsic :: ISO_C_BINDING
  implicit none
  integer (C_INT), value :: ncid, dimid
  type (C_PTR), value :: lengthp
end function nc_inq_dimlen

integer (C_INT) function nc_inq_varid(ncid,name,varidp) bind(C, name='nc_inq_varid')
  use, intrinsic :: ISO_C_BINDING
  implicit none
  integer (C_INT), value :: ncid
  character, dimension(*) :: name
  type (C_PTR), value :: varidp
end function nc_inq_varid

integer (C_INT) function nc_inq_dimid(ncid,name,dimidp) bind(C, name='nc_inq_dimid')
  use, intrinsic :: ISO_C_BINDING
  implicit none
  integer (C_INT), value :: ncid
  character, dimension(*) :: name
  type (C_PTR), value :: dimidp
end function nc_inq_dimid

integer (C_INT) function nc_inq_vartype(ncid,varid,xtypep) bind(C, name='nc_inq_vartype')
  use, intrinsic :: ISO_C_BINDING
  implicit none
  integer (C_INT), value :: ncid, varid
  type (C_PTR), value :: xtypep
end function nc_inq_vartype

integer (C_INT) function nc_inq(ncid,ndimsp,nvarsp,ngattsp,unlimdimidp) bind(C, name='nc_inq')
  use, intrinsic :: ISO_C_BINDING
  implicit none
  integer (C_INT), value :: ncid
  type (C_PTR), value :: ndimsp, nvarsp, ngattsp, unlimdimidp
end function nc_inq

integer (C_INT) function nc_inq_varname(ncid,varid,tp) bind(C, name='nc_inq_varname')
  use, intrinsic :: ISO_C_BINDING
  implicit none
  integer (C_INT), value :: ncid, varid
  character, dimension(*) :: tp
end function nc_inq_varname

integer (C_INT) function nc_inq_ndims(ncid,ndimsp) bind(C, name='nc_inq_ndims')
  use, intrinsic :: ISO_C_BINDING
  implicit none
  integer (C_INT), value :: ncid
  type (C_PTR), value :: ndimsp
end function nc_inq_ndims

integer (C_INT) function nc_inq_nvars(ncid,nvarsp) bind(C, name='nc_inq_nvars')
  use, intrinsic :: ISO_C_BINDING
  implicit none
  integer (C_INT), value :: ncid
  type (C_PTR), value :: nvarsp
end function nc_inq_nvars

type (C_PTR) function nc_inq_libvers() bind(C, name='nc_inq_libvers')
  use, intrinsic :: ISO_C_BINDING
  implicit none
end function nc_inq_libvers

integer (C_INT) function nc_inq_dim(ncid,dimid,name,lengthp) bind(C, name='nc_inq_dim')
  use, intrinsic :: ISO_C_BINDING
  implicit none
  integer (C_INT), value :: ncid, dimid
  character, dimension(*) :: name
  type (C_PTR), value :: lengthp
end function nc_inq_dim

integer (C_INT) function nc_get_att_text(ncid,varid,name,tp) bind(C, name='nc_get_att_text')
  use, intrinsic :: ISO_C_BINDING
  implicit none
  integer (C_INT), value :: ncid, varid
  character, dimension(*) :: name
  type (C_PTR), value :: tp
end function nc_get_att_text

integer (C_INT) function nc_get_att_float(ncid,varid,name,rp) bind(C, name='nc_get_att_float')
  use, intrinsic :: ISO_C_BINDING
  implicit none
  integer (C_INT), value :: ncid, varid
  type (C_PTR), value :: rp
  character, dimension(*) :: name
end function nc_get_att_float

integer (C_INT) function nc_get_att_int(ncid,varid,name,ip) bind(C, name='nc_get_att_int')
  use, intrinsic :: ISO_C_BINDING
  implicit none
  integer (C_INT), value :: ncid, varid
  type (C_PTR), value :: ip
  character, dimension(*) :: name
end function nc_get_att_int

integer (C_INT) function nc_get_vara_float(ncid,varid,start,count,rp) bind(C, name='nc_get_vara_float')
  use, intrinsic :: ISO_C_BINDING
  implicit none
  integer (C_INT), value :: ncid, varid
  integer (C_SIZE_T), dimension(*) :: start, count
  type (C_PTR), value :: rp
end function nc_get_vara_float

integer (C_INT) function nc_get_vara_int(ncid,varid,start,count,ip) bind(C, name='nc_get_vara_int')
  use, intrinsic :: ISO_C_BINDING
  implicit none
  integer (C_INT), value :: ncid, varid
  integer (C_SIZE_T), dimension(*) :: start, count
  type (C_PTR), value :: ip
end function nc_get_vara_int

integer (C_INT) function nc_get_vara_short(ncid,varid,start,count,sp) bind(C, name='nc_get_vara_short')
  use, intrinsic :: ISO_C_BINDING
  implicit none
  integer (C_INT), value :: ncid, varid
  integer (C_SIZE_T), dimension(*) :: start, count
  type (C_PTR), value :: sp
end function nc_get_vara_short

integer (C_INT) function nc_get_vara_double(ncid,varid,start,count,dp) bind(C, name='nc_get_vara_double')
  use, intrinsic :: ISO_C_BINDING
  implicit none
  integer (C_INT), value :: ncid, varid
  integer (C_SIZE_T), dimension(*) :: start, count
  type (C_PTR), value :: dp
end function nc_get_vara_double

integer (C_INT) function nc_get_var_float(ncid,varid,rp) bind(C, name='nc_get_var_float')
  use, intrinsic :: ISO_C_BINDING
  implicit none
  integer (C_INT), value :: ncid, varid
  type (C_PTR), value :: rp
end function nc_get_var_float

integer (C_INT) function nc_get_var1_float(ncid,varid,start,rp) bind(C, name='nc_get_var1_float')
  use, intrinsic :: ISO_C_BINDING
  implicit none
  integer (C_INT), value :: ncid, varid
  integer (C_SIZE_T), dimension(*) :: start
  type (C_PTR), value :: rp
end function nc_get_var1_float

integer (C_INT) function nc_get_var1_int(ncid,varid,start,ip) bind(C, name='nc_get_var1_int')
  use, intrinsic :: ISO_C_BINDING
  implicit none
  integer (C_INT), value :: ncid, varid
  integer (C_SIZE_T), dimension(*) :: start
  type (C_PTR), value :: ip
end function nc_get_var1_int

integer (C_INT) function nc_def_dim(ncid,name,size,dimidp) bind(C, name='nc_def_dim')
  use, intrinsic :: ISO_C_BINDING
  implicit none
  integer (C_INT), value :: ncid
  integer (C_SIZE_T), value :: size
  type (C_PTR), value :: dimidp
  character, dimension(*) :: name
end function nc_def_dim

integer (C_INT) function nc_def_var(ncid,name,xtype,ndims,dimids,varidp) bind(C, name='nc_def_var')
  use, intrinsic :: ISO_C_BINDING
  implicit none
  integer (C_INT), value :: ncid, xtype, ndims
  integer (C_INT), dimension(*) :: dimids
  type (C_PTR), value :: varidp
  character, dimension(*) :: name
end function nc_def_var

integer (C_INT) function nc_rename_dim(ncid,dimid,name) bind(C, name='nc_rename_dim')
  use, intrinsic :: ISO_C_BINDING
  implicit none
  integer (C_INT), value :: ncid, dimid
  character, dimension(*) :: name
end function nc_rename_dim
    
integer (C_INT) function nc_put_att_text(ncid,varid,name,size,tp) bind(C, name='nc_put_att_text')
  use, intrinsic :: ISO_C_BINDING
  implicit none
  integer (C_INT), value :: ncid, varid
  integer (C_SIZE_T), value :: size
  character, dimension(*) :: name
  character, dimension(*) :: tp
end function nc_put_att_text

integer (C_INT) function nc_put_att_short(ncid,varid,name,xtype,size,sp) bind(C, name='nc_put_att_short')
  use, intrinsic :: ISO_C_BINDING
  implicit none
  integer (C_INT), value :: ncid, varid, xtype
  integer (C_SIZE_T), value :: size
  type (C_PTR), value :: sp
  character, dimension(*) :: name
end function nc_put_att_short

integer (C_INT) function nc_put_att_float(ncid,varid,name,xtype,size,fp) bind(C, name='nc_put_att_float')
  use, intrinsic :: ISO_C_BINDING
  implicit none
  integer (C_INT), value :: ncid, varid, xtype
  integer (C_SIZE_T), value :: size
  type (C_PTR), value :: fp
  character, dimension(*) :: name
end function nc_put_att_float

integer (C_INT) function nc_put_att_int(ncid,varid,name,xtype,size,ip) bind(C, name='nc_put_att_int')
  use, intrinsic :: ISO_C_BINDING
  implicit none
  integer (C_INT), value :: ncid, varid, xtype
  integer (C_SIZE_T), value :: size
  type (C_PTR), value :: ip
  character, dimension(*) :: name
end function nc_put_att_int

integer (C_INT) function nc_put_vara_float(ncid,varid,start,count,rp) bind(C, name='nc_put_vara_float')
  use, intrinsic :: ISO_C_BINDING
  implicit none
  integer (C_INT), value :: ncid, varid
  integer (C_SIZE_T), dimension(*) :: start, count
  type (C_PTR), value :: rp
end function nc_put_vara_float

integer (C_INT) function nc_put_vara_int(ncid,varid,start,count,ip) bind(C, name='nc_put_vara_int')
  use, intrinsic :: ISO_C_BINDING
  implicit none
  integer (C_INT), value :: ncid, varid
  integer (C_SIZE_T), dimension(*) :: start, count
  type (C_PTR), value :: ip
end function nc_put_vara_int

integer (C_INT) function nc_put_vara_short(ncid,varid,start,count,sp) bind(C, name='nc_put_vara_short')
  use, intrinsic :: ISO_C_BINDING
  implicit none
  integer (C_INT), value :: ncid, varid
  integer (C_SIZE_T), dimension(*) :: start, count
  type (C_PTR), value :: sp
end function nc_put_vara_short

integer (C_INT) function nc_put_vara_double(ncid,varid,start,count,dp) bind(C, name='nc_put_vara_double')
  use, intrinsic :: ISO_C_BINDING
  implicit none
  integer (C_INT), value :: ncid, varid
  integer (C_SIZE_T), dimension(*) :: start, count
  type (C_PTR), value :: dp
end function nc_put_vara_double

integer (C_INT) function nc_put_var_float(ncid,varid,rp) bind(C, name='nc_put_var_float')
  use, intrinsic :: ISO_C_BINDING
  implicit none
  integer (C_INT), value :: ncid, varid
  type (C_PTR), value :: rp
end function nc_put_var_float

integer (C_INT) function nc_put_var1_int(ncid,varid,start,ip) bind(C, name='nc_put_var1_int')
  use, intrinsic :: ISO_C_BINDING
  implicit none
  integer (C_INT), value :: ncid, varid
  integer (C_SIZE_T), dimension(*) :: start
  type (C_PTR), value :: ip
end function nc_put_var1_int

integer (C_INT) function nc_put_var1_float(ncid,varid,start,rp) bind(C, name='nc_put_var1_float')
  use, intrinsic :: ISO_C_BINDING
  implicit none
  integer (C_INT), value :: ncid, varid
  integer (C_SIZE_T), dimension(*) :: start
  type (C_PTR), value :: rp
end function nc_put_var1_float

integer (C_INT) function nc_put_var1_double(ncid,varid,start,dp) bind(C, name='nc_put_var1_double')
  use, intrinsic :: ISO_C_BINDING
  implicit none
  integer (C_INT), value :: ncid, varid
  integer (C_SIZE_T), dimension(*) :: start
  type (C_PTR), value :: dp
end function nc_put_var1_double

integer (C_INT) function nc_put_var1_text(ncid,varid,start,tp) bind(C, name='nc_put_var1_text')
  use, intrinsic :: ISO_C_BINDING
  implicit none
  integer (C_INT), value :: ncid, varid
  integer (C_SIZE_T), dimension(*) :: start
  character, dimension(*) :: tp
end function nc_put_var1_text

integer (C_INT) function nc_copy_att(ncidin,varidin,name,ncidout,varidout) bind(C, name='nc_copy_att')
  use, intrinsic :: ISO_C_BINDING
  implicit none
  integer (C_INT), value :: ncidin, varidin, ncidout, varidout
  character, dimension(*) :: name
end function nc_copy_att
    
end interface

interface nf_get_att_real
  module procedure nf_get_att_real_s, nf_get_att_real_v
end interface nf_get_att_real

interface nf_get_att_int
  module procedure nf_get_att_int_s, nf_get_att_int_v
end interface nf_get_att_int

interface nf_get_vara_real
  module procedure nf_get_vara_real_d1, nf_get_vara_real_d2, nf_get_vara_real_d3, nf_get_vara_real_d4
end interface nf_get_vara_real

interface nf_get_vara_int
  module procedure nf_get_vara_int_d1, nf_get_vara_int_d2, nf_get_vara_int_d3, nf_get_vara_int_d4
end interface nf_get_vara_int

interface nf_get_vara_int2
  module procedure nf_get_vara_int2_d1, nf_get_vara_int2_d2, nf_get_vara_int2_d3, nf_get_vara_int2_d4
end interface nf_get_vara_int2

interface nf_get_vara_double
  module procedure nf_get_vara_double_d1, nf_get_vara_double_d2, nf_get_vara_double_d3, nf_get_vara_double_d4
end interface nf_get_vara_double

interface nf_get_var_real
  module procedure nf_get_var_real_d1, nf_get_var_real_d2, nf_get_var_real_d3, nf_get_var_real_d4
end interface nf_get_var_real

interface nf_get_var1_real
  module procedure nf_get_var1_real_s, nf_get_var1_real_v
end interface nf_get_var1_real

interface nf_get_var1_int
  module procedure nf_get_var1_int_s, nf_get_var1_int_v
end interface nf_get_var1_int
    
interface nf_put_vara_real
  module procedure nf_put_vara_real_d1, nf_put_vara_real_d2, nf_put_vara_real_d3, nf_put_vara_real_d4
end interface nf_put_vara_real

interface nf_put_vara_int
  module procedure nf_put_vara_int_d1, nf_put_vara_int_d2, nf_put_vara_int_d3, nf_put_vara_int_d4
end interface nf_put_vara_int

interface nf_put_vara_int2
  module procedure nf_put_vara_int2_d1, nf_put_vara_int2_d2, nf_put_vara_int2_d3, nf_put_vara_int2_d4
end interface nf_put_vara_int2
    
interface nf_put_vara_double
  module procedure nf_put_vara_double_d1, nf_put_vara_double_d2, nf_put_vara_double_d3, nf_put_vara_double_d4
end interface nf_put_vara_double

interface nf_put_var_real
  module procedure nf_put_var_real_d1, nf_put_var_real_d2, nf_put_var_real_d3, nf_put_var_real_d4
end interface nf_put_var_real

integer, parameter :: nf_unlimited = 0

integer, parameter :: nf_noerr = 0
integer, parameter :: nf_enameinuse = -42

integer, parameter :: nf_nowrite = 0
integer, parameter :: nf_write = 1
integer, parameter :: nf_clobber = 0
integer, parameter :: nf_nofill = 256
integer, parameter :: nf_64bit_offset = 512

integer, parameter :: nf_global = -1

integer, parameter :: nf_char = 2
integer, parameter :: nf_short = 3
integer, parameter :: nf_int2 = nf_short
integer, parameter :: nf_int = 4
integer, parameter :: nf_float = 5
integer, parameter :: nf_real = nf_float
integer, parameter :: nf_double = 6

real, parameter :: nf_fill_float = 9.9692099683868690e+36

integer, parameter :: charsize = 255

contains
    
integer function nf_open(name,mode,ncid) result(ierr)
  implicit none
  integer, intent(in) :: mode
  integer, intent(out) :: ncid
  character(len=*), intent(in) :: name
  integer (C_INT) :: c_mode
  integer (C_INT), target :: c_ncid
  character, dimension(charsize) :: c_name
  integer i, ix
  if (len(name)>charsize) then
    write(6,*) "ERROR: charsize in netcdf_m is too small"
    write(6,*) "Req ",len(name)," Actual ",charsize
    stop
  end if
  c_mode = mode
  ix = len_trim(name)
  do i=1,len_trim(name)
    c_name(i) = name(i:i)
  end do
  c_name(ix+1:ix+1) = C_NULL_CHAR
  ierr = nc_open(c_name,c_mode,C_LOC(c_ncid))
  ncid = c_ncid
end function nf_open

integer function nf_close(ncid) result(ierr)
  implicit none
  integer, intent(in) :: ncid
  integer (C_INT) :: c_ncid
  c_ncid = ncid
  ierr = nc_close(c_ncid)
end function nf_close

integer function nf_create(name,mode,ncid) result(ierr)
  implicit none
  integer, intent(in) :: mode
  integer, intent(out) :: ncid
  character(len=*), intent(in) :: name
  integer (C_INT) :: c_mode
  integer (C_INT), target :: c_ncid
  character, dimension(charsize) :: c_name
  integer i, ix
  if (len(name)>charsize) then
    write(6,*) "ERROR: charsize in netcdf_m is too small"
    write(6,*) "Req ",len(name)," Actual ",charsize
    stop
  end if
  c_mode = mode
  ix = len_trim(name)
  do i=1,len_trim(name)
    c_name(i) = name(i:i)
  end do
  c_name(ix+1:ix+1) = C_NULL_CHAR
  ierr = nc_create(c_name,c_mode,C_LOC(c_ncid))
  ncid = c_ncid
end function nf_create

integer function nf_enddef(ncid) result(ierr)
  implicit none
  integer, intent(in) :: ncid
  integer (C_INT) :: c_ncid
  c_ncid = ncid
  ierr = nc_enddef(c_ncid)
end function nf_enddef

integer function nf_set_fill(ncid,fillmode,omode) result(ierr)
  implicit none
  integer, intent(in) :: ncid, fillmode
  integer, intent(out) :: omode
  integer (C_INT) :: c_ncid, c_fillmode
  integer (C_INT), target :: c_omode
  c_ncid = ncid
  c_fillmode = fillmode
  ierr = nc_set_fill(c_ncid,c_fillmode,C_LOC(c_omode))
  omode = c_omode
end function nf_set_fill

integer function nf_redef(ncid) result(ierr)
  implicit none
  integer, intent(in) :: ncid
  integer (C_INT) :: c_ncid
  c_ncid = ncid
  ierr = nc_redef(c_ncid)
end function nf_redef

integer function nf_sync(ncid) result(ierr)
  implicit none
  integer, intent(in) :: ncid
  integer (C_INT) :: c_ncid
  c_ncid = ncid
  ierr = nc_sync(c_ncid)
end function nf_sync

character(len=charsize) function nf_strerror(ncerr) result(msg)
  implicit none
  integer, intent(in) :: ncerr
  integer (C_INT) :: c_ncerr
  character, dimension(:), pointer :: c_tp
  character(len=charsize) :: temp
  type (C_PTR) :: c_msg
  integer ix
  integer, dimension(1) :: string_shape
  c_ncerr = ncerr
  c_msg = nc_strerror(c_ncerr)
  string_shape = charsize
  call c_f_pointer(c_msg,c_tp,string_shape)
  temp = ''
  temp = transfer(c_tp,temp)
  msg = ''
  ix = index(temp,C_NULL_CHAR) - 1
  msg = temp(1:ix)
end function nf_strerror

integer function nf__open(name,mode,bufrsizehint,ncid) result(ierr)
  implicit none
  integer, intent(in) :: mode
  integer, intent(inout) :: bufrsizehint
  integer, intent(out) :: ncid
  character(len=*), intent(in) :: name
  integer (C_INT) :: c_mode
  integer (C_INT), target :: c_ncid, c_bufrsizehint
  character, dimension(charsize) :: c_name
  integer i, ix
  if (len(name)>charsize) then
    write(6,*) "ERROR: charsize in netcdf_m is too small"
    write(6,*) "Req ",len(name)," Actual ",charsize
    stop
  end if
  c_mode = mode
  ix = len_trim(name)
  do i=1,len_trim(name)
    c_name(i) = name(i:i)
  end do
  c_name(ix+1:ix+1) = C_NULL_CHAR
  c_bufrsizehint = bufrsizehint
  ierr = nc__open(c_name,c_mode,C_LOC(c_bufrsizehint),C_LOC(c_ncid))
  bufrsizehint = c_bufrsizehint
  ncid = c_ncid
end function nf__open

integer function nf__create(name,mode,initialsz,bufrsizehint,ncid) result(ierr)
  implicit none
  integer, intent(in) :: mode, initialsz
  integer, intent(inout) :: bufrsizehint
  integer, intent(out) :: ncid
  character(len=*), intent(in) :: name
  integer (C_INT) :: c_mode
  integer (C_SIZE_T) :: c_initialsz
  integer (C_SIZE_T), target :: c_bufrsizehint
  integer (C_INT), target :: c_ncid
  character, dimension(charsize) :: c_name
  integer i, ix
  if (len(name)>charsize) then
    write(6,*) "ERROR: charsize in netcdf_m is too small"
    write(6,*) "Req ",len(name)," Actual ",charsize
    stop
  end if
  c_mode = mode
  ix = len_trim(name)
  do i=1,len_trim(name)
    c_name(i) = name(i:i)
  end do
  c_name(ix+1:ix+1) = C_NULL_CHAR
  c_initialsz = initialsz
  c_bufrsizehint = bufrsizehint
  ierr = nc__create(c_name,c_mode,c_initialsz,C_LOC(c_bufrsizehint),C_LOC(c_ncid))
  bufrsizehint = c_bufrsizehint
  ncid = c_ncid
end function nf__create

integer function nf__enddef(ncid,h_minfree,v_align,v_minfree,r_align) result(ierr)
  implicit none
  integer, intent(in) :: ncid, h_minfree, v_align, v_minfree, r_align
  integer (C_INT) :: c_ncid
  integer (C_SIZE_T) :: c_h_minfree, c_v_align, c_v_minfree, c_r_align
  c_ncid = ncid
  c_h_minfree = h_minfree
  c_v_align = v_align
  c_v_minfree = v_minfree
  c_r_align = r_align
  ierr = nc__enddef(c_ncid,c_h_minfree,c_v_align,c_v_minfree,c_r_align)
end function nf__enddef

integer function nf_abort(ncid) result(ierr)
  implicit none
  integer, intent(in) :: ncid
  integer (C_INT) :: c_ncid
  c_ncid = ncid
  ierr = nc_abort(c_ncid)
end function nf_abort

integer function nf_inq_varndims(ncid,varid,ndims) result(ierr)
  implicit none
  integer, intent(in) :: ncid, varid
  integer, intent(out) :: ndims
  integer (C_INT) :: c_ncid, c_varid
  integer (C_INT), target :: c_ndims
  c_ncid = ncid
  c_varid = varid
  ierr = nc_inq_varndims(c_ncid,c_varid,C_LOC(c_ndims))
  ndims = c_ndims
end function nf_inq_varndims

integer function nf_inq_vardimid(ncid,varid,dimids) result(ierr)
  implicit none
  integer, intent(in) :: ncid, varid
  integer, dimension(:), intent(out) :: dimids
  integer (C_INT) :: c_ncid, c_varid
  integer (C_INT), dimension(size(dimids)), target :: c_dimids
  integer i, ix
  c_ncid = ncid
  c_varid = varid
  ierr = nc_inq_vardimid(c_ncid,c_varid,C_LOC(c_dimids))
  ix = size(dimids)
  do i=1,ix
    dimids(i) = c_dimids(ix-i+1)
  end do
end function nf_inq_vardimid

integer function nf_inq_dimlen(ncid,dimid,length) result(ierr)
  implicit none
  integer, intent(in) :: ncid, dimid
  integer, intent(out) :: length
  integer (C_INT) :: c_ncid, c_dimid
  integer (C_SIZE_T), target :: c_length
  c_ncid = ncid
  c_dimid = dimid
  ierr = nc_inq_dimlen(c_ncid,c_dimid,C_LOC(c_length))
  length = c_length
end function nf_inq_dimlen

integer function nf_inq_varid(ncid,name,varid) result(ierr)
  implicit none
  integer, intent(in) :: ncid
  integer, intent(out) :: varid
  character(len=*), intent(in) :: name
  integer (C_INT) :: c_ncid
  integer (C_INT), target :: c_varid
  character, dimension(charsize) :: c_name
  integer i, ix
  if (len(name)>charsize) then
    write(6,*) "ERROR: charsize in netcdf_m is too small"
    write(6,*) "Req ",len(name)," Actual ",charsize
    stop
  end if
  c_ncid = ncid
  ix = len_trim(name)
  do i=1,len_trim(name)
    c_name(i) = name(i:i)
  end do
  c_name(ix+1:ix+1) = C_NULL_CHAR
  ierr = nc_inq_varid(c_ncid,c_name,C_LOC(c_varid))
  varid = c_varid
end function nf_inq_varid

integer function nf_inq_dimid(ncid,name,dimid) result(ierr)
  implicit none
  integer, intent(in) :: ncid
  integer, intent(out) :: dimid
  character(len=*), intent(in) :: name
  integer (C_INT) :: c_ncid
  integer (C_INT), target :: c_dimid
  character, dimension(charsize) :: c_name
  integer i, ix
  if (len(name)>charsize) then
    write(6,*) "ERROR: charsize in netcdf_m is too small"
    write(6,*) "Req ",len(name)," Actual ",charsize
    stop
  end if
  c_ncid = ncid
  ix = len_trim(name)
  do i=1,len_trim(name)
    c_name(i) = name(i:i)
  end do
  c_name(ix+1:ix+1) = C_NULL_CHAR
  ierr = nc_inq_dimid(c_ncid,c_name,C_LOC(c_dimid))
  dimid = c_dimid
end function nf_inq_dimid

integer function nf_inq_vartype(ncid,varid,xtype) result(ierr)
  implicit none
  integer, intent(in) :: ncid, varid
  integer, intent(out) :: xtype
  integer (C_INT) :: c_ncid, c_varid
  integer (C_INT), target :: c_xtype
  c_ncid = ncid
  c_varid = varid
  ierr = nc_inq_vartype(c_ncid,c_varid,C_LOC(c_xtype))
  xtype = c_xtype
end function nf_inq_vartype

integer function nf_inq(ncid,ndims,nvars,ngatts,unlimdimid) result(ierr)
  implicit none
  integer, intent(in) :: ncid
  integer, intent(out) :: ndims, nvars, ngatts, unlimdimid
  integer (C_INT) :: c_ncid
  integer (C_INT), target :: c_ndims, c_nvars, c_ngatts, c_unlimdimid
  c_ncid = ncid
  ierr = nc_inq(c_ncid,C_LOC(c_ndims),C_LOC(c_nvars),C_LOC(c_ngatts),C_LOC(c_unlimdimid))
  ndims = c_ndims
  nvars = c_nvars
  ngatts = c_ngatts
  unlimdimid = c_unlimdimid
end function nf_inq

integer function nf_inq_varname(ncid,varid,tp) result(ierr)
  implicit none
  integer, intent(in) :: ncid, varid
  character(len=*), intent(out) :: tp
  character(len=charsize) :: temp
  integer (C_INT) :: c_ncid, c_varid
  character, dimension(charsize) :: c_tp
  integer i, ix
  if (len(tp)>charsize) then
    write(6,*) "ERROR: charsize in netcdf_m is too small"
    write(6,*) "Req ",len(tp)," Actual ",charsize
    stop
  end if
  c_ncid = ncid
  c_varid = varid
  c_tp(:) = ''  
  ierr = nc_inq_varname(c_ncid,c_varid,c_tp)
  temp = ''
  temp = transfer(c_tp,temp)
  !tp = ''
  !ix = index(temp,C_NULL_CHAR) - 1
  !tp = temp(1:ix)
  tp = temp
end function nf_inq_varname

integer function nf_inq_ndims(ncid,ndims) result(ierr)
  implicit none
  integer, intent(in) :: ncid
  integer, intent(out) :: ndims
  integer (C_INT) :: c_ncid
  integer (C_INT), target :: c_ndims
  c_ncid = ncid
  ierr = nc_inq_ndims(c_ncid,C_LOC(c_ndims))
  ndims = c_ndims
end function nf_inq_ndims

integer function nf_inq_nvars(ncid,nvars) result(ierr)
  implicit none
  integer, intent(in) :: ncid
  integer, intent(out) :: nvars
  integer (C_INT) :: c_ncid
  integer (C_INT), target :: c_nvars
  c_ncid = ncid
  ierr = nc_inq_nvars(c_ncid,C_LOC(c_nvars))
  nvars = c_nvars
end function nf_inq_nvars

character(len=charsize) function nf_inq_libvers() result(msg)
  implicit none
  character, dimension(:), pointer :: c_tp
  character(len=charsize) :: temp
  type (C_PTR) :: c_msg
  integer ix
  integer, dimension(1) :: string_shape
  c_msg = nc_inq_libvers()
  string_shape = charsize
  call c_f_pointer(c_msg,c_tp,string_shape)
  temp = ''
  temp = transfer(c_tp,temp)
  msg = ''
  ix = index(temp,C_NULL_CHAR) - 1
  msg = temp(1:ix)
end function nf_inq_libvers

integer function nf_inq_dim(ncid,dimid,name,length) result(ierr)
  implicit none
  integer, intent(in) :: ncid, dimid
  integer, intent(out) :: length
  character(len=*), intent(out) :: name
  character(len=charsize) :: temp
  integer (C_INT) :: c_ncid, c_dimid  
  integer (C_SIZE_T), target :: c_length
  character, dimension(charsize) :: c_name
  integer ix
  c_ncid = ncid
  c_dimid = dimid
  ierr = nc_inq_dim(c_ncid,c_dimid,c_name,C_LOC(c_length))
  temp = ''
  temp = transfer(c_name,temp)
  !name = ''
  !ix = index(temp,C_NULL_CHAR) - 1
  !name = temp(1:ix)
  name = temp
  length = c_length
end function nf_inq_dim

integer function nf_get_att_text(ncid,varid,name,tp) result(ierr)
  implicit none
  integer, intent(in) :: ncid, varid
  character(len=*), intent(in) :: name
  character(len=*), intent(out) :: tp
  character(len=charsize) :: temp
  integer (C_INT) :: c_ncid, c_varid
  character, dimension(charsize) :: c_name
  character, dimension(charsize), target :: c_tp
  integer i, ix
  if (len(name)>charsize) then
    write(6,*) "ERROR: charsize in netcdf_m is too small"
    write(6,*) "Req ",len(name)," Actual ",charsize
    stop
  end if
  if (len(tp)>charsize) then
    write(6,*) "ERROR: charsize in netcdf_m is too small"
    write(6,*) "Req ",len(tp)," Actual ",charsize
    stop
  end if  
  c_ncid = ncid
  c_varid = varid
  ix = len_trim(name)
  do i=1,len_trim(name)
    c_name(i) = name(i:i)
  end do
  c_name(ix+1:ix+1) = C_NULL_CHAR
  c_tp = ''
  ierr = nc_get_att_text(c_ncid,c_varid,c_name,C_LOC(c_tp))
  temp = ''
  temp = transfer(c_tp,temp)
  !tp = ''
  !ix = index(trim(temp,C_NULL_CHAR) - 1
  !tp = temp(1:ix)
  tp = temp
end function nf_get_att_text

integer function nf_get_att_real_s(ncid,varid,name,rp) result(ierr)
  implicit none
  integer, intent(in) :: ncid, varid
  real, intent(out) :: rp
  character(len=*), intent(in) :: name
  real (C_FLOAT), target :: c_rp
  integer (C_INT) :: c_ncid, c_varid
  character, dimension(charsize) :: c_name
  integer i, ix
  if (len(name)>charsize) then
    write(6,*) "ERROR: charsize in netcdf_m is too small"
    write(6,*) "Req ",len(name)," Actual ",charsize
    stop
  end if
  c_ncid = ncid
  c_varid = varid
  ix = len_trim(name)
  do i=1,len_trim(name)
    c_name(i) = name(i:i)
  end do
  c_name(ix+1:ix+1) = C_NULL_CHAR
  ierr = nc_get_att_float(c_ncid,c_varid,c_name,C_LOC(c_rp))
  rp = c_rp  
end function nf_get_att_real_s

integer function nf_get_att_real_v(ncid,varid,name,rp) result(ierr)
  implicit none
  integer, intent(in) :: ncid, varid
  real, dimension(:), intent(out) :: rp
  character(len=*), intent(in) :: name
  real (C_FLOAT), dimension(size(rp)), target :: c_rp
  integer (C_INT) :: c_ncid, c_varid
  character, dimension(charsize) :: c_name
  integer i, ix
  if (len(name)>charsize) then
    write(6,*) "ERROR: charsize in netcdf_m is too small"
    write(6,*) "Req ",len(name)," Actual ",charsize
    stop
  end if
  c_ncid = ncid
  c_varid = varid
  ix = len_trim(name)
  do i=1,len_trim(name)
    c_name(i) = name(i:i)
  end do
  c_name(ix+1:ix+1) = C_NULL_CHAR
  ierr = nc_get_att_float(c_ncid,c_varid,c_name,C_LOC(c_rp))
  rp = c_rp  
end function nf_get_att_real_v

integer function nf_get_att_int_s(ncid,varid,name,ip) result(ierr)
  implicit none
  integer, intent(in) :: ncid, varid
  integer, intent(out) :: ip
  character(len=*), intent(in) :: name
  integer (C_INT), target :: c_ip
  integer (C_INT) :: c_ncid, c_varid
  character, dimension(charsize) :: c_name
  integer i, ix
  if (len(name)>charsize) then
    write(6,*) "ERROR: charsize in netcdf_m is too small"
    write(6,*) "Req ",len(name)," Actual ",charsize
    stop
  end if
  c_ncid = ncid
  c_varid = varid
  ix = len_trim(name)
  do i=1,len_trim(name)
    c_name(i) = name(i:i)
  end do
  c_name(ix+1:ix+1) = C_NULL_CHAR
  ierr = nc_get_att_int(c_ncid,c_varid,c_name,C_LOC(c_ip))
  ip = c_ip  
end function nf_get_att_int_s

integer function nf_get_att_int_v(ncid,varid,name,ip) result(ierr)
  implicit none
  integer, intent(in) :: ncid, varid
  integer, dimension(:), intent(out) :: ip
  character(len=*), intent(in) :: name
  integer (C_INT), dimension(size(ip)), target :: c_ip
  integer (C_INT) :: c_ncid, c_varid
  character, dimension(charsize) :: c_name
  integer i, ix
  if (len(name)>charsize) then
    write(6,*) "ERROR: charsize in netcdf_m is too small"
    write(6,*) "Req ",len(name)," Actual ",charsize
    stop
  end if
  c_ncid = ncid
  c_varid = varid
  ix = len_trim(name)
  do i=1,len_trim(name)
    c_name(i) = name(i:i)
  end do
  c_name(ix+1:ix+1) = C_NULL_CHAR
  ierr = nc_get_att_int(c_ncid,c_varid,c_name,C_LOC(c_ip))
  ip = c_ip  
end function nf_get_att_int_v

integer function nf_get_vara_real_d1(ncid,varid,start,ncount,fp) result(ierr)
  implicit none
  integer, intent(in) :: ncid, varid
  integer, dimension(:), intent(in) :: start, ncount
  real, dimension(:), intent(out) :: fp
  integer (C_INT) :: c_ncid, c_varid
  integer (C_SIZE_T), dimension(size(start)) :: c_start, c_ncount
  real (C_FLOAT), dimension(size(fp)), target :: c_fp
  integer i, ix
  c_ncid = ncid
  c_varid = varid
  ix = size(start)
  do i = 1,ix
    c_start(ix-i+1) = start(i) - 1
    c_ncount(ix-i+1) = ncount(i)
  end do
  ierr = nc_get_vara_float(c_ncid,c_varid,c_start,c_ncount,C_LOC(c_fp))
  fp = c_fp
end function nf_get_vara_real_d1

integer function nf_get_vara_real_d2(ncid,varid,start,ncount,fp) result(ierr)
  implicit none
  integer, intent(in) :: ncid, varid
  integer, dimension(:), intent(in) :: start, ncount
  real, dimension(:,:), intent(out) :: fp
  integer (C_INT) :: c_ncid, c_varid
  integer (C_SIZE_T), dimension(size(start)) :: c_start, c_ncount
  real (C_FLOAT), dimension(size(fp,1),size(fp,2)), target :: c_fp
  integer i, ix
  c_ncid = ncid
  c_varid = varid
  ix = size(start)
  do i = 1,ix
    c_start(ix-i+1) = start(i) - 1
    c_ncount(ix-i+1) = ncount(i)
  end do
  ierr = nc_get_vara_float(c_ncid,c_varid,c_start,c_ncount,C_LOC(c_fp))
  fp = c_fp
end function nf_get_vara_real_d2

integer function nf_get_vara_real_d3(ncid,varid,start,ncount,fp) result(ierr)
  implicit none
  integer, intent(in) :: ncid, varid
  integer, dimension(:), intent(in) :: start, ncount
  real, dimension(:,:,:), intent(out) :: fp
  integer (C_INT) :: c_ncid, c_varid
  integer (C_SIZE_T), dimension(size(start)) :: c_start, c_ncount
  real (C_FLOAT), dimension(size(fp,1),size(fp,2),size(fp,3)), target :: c_fp
  integer i, ix
  c_ncid = ncid
  c_varid = varid
  ix = size(start)
  do i = 1,ix
    c_start(ix-i+1) = start(i) - 1
    c_ncount(ix-i+1) = ncount(i)
  end do
  ierr = nc_get_vara_float(c_ncid,c_varid,c_start,c_ncount,C_LOC(c_fp))
  fp = c_fp
end function nf_get_vara_real_d3

integer function nf_get_vara_real_d4(ncid,varid,start,ncount,fp) result(ierr)
  implicit none
  integer, intent(in) :: ncid, varid
  integer, dimension(:), intent(in) :: start, ncount
  real, dimension(:,:,:,:), intent(out) :: fp
  integer (C_INT) :: c_ncid, c_varid
  integer (C_SIZE_T), dimension(size(start)) :: c_start, c_ncount
  real (C_FLOAT), dimension(size(fp,1),size(fp,2),size(fp,3),size(fp,4)), target :: c_fp
  integer i, ix
  c_ncid = ncid
  c_varid = varid
  ix = size(start)
  do i = 1,ix
    c_start(ix-i+1) = start(i) - 1
    c_ncount(ix-i+1) = ncount(i)
  end do
  ierr = nc_get_vara_float(c_ncid,c_varid,c_start,c_ncount,C_LOC(c_fp))
  fp = c_fp
end function nf_get_vara_real_d4

integer function nf_get_vara_int_d1(ncid,varid,start,ncount,ip) result(ierr)
  implicit none
  integer, intent(in) :: ncid, varid
  integer, dimension(:), intent(in) :: start, ncount
  integer, dimension(:), intent(out) :: ip
  integer (C_INT) :: c_ncid, c_varid
  integer (C_SIZE_T), dimension(size(start)) :: c_start, c_ncount
  integer (C_INT), dimension(size(ip)), target :: c_ip
  integer i, ix
  c_ncid = ncid
  c_varid = varid
  ix = size(start)
  do i = 1,ix
    c_start(ix-i+1) = start(i) - 1
    c_ncount(ix-i+1) = ncount(i)
  end do
  ierr = nc_get_vara_int(c_ncid,c_varid,c_start,c_ncount,C_LOC(c_ip))
  ip = c_ip
end function nf_get_vara_int_d1

integer function nf_get_vara_int_d2(ncid,varid,start,ncount,ip) result(ierr)
  implicit none
  integer, intent(in) :: ncid, varid
  integer, dimension(:), intent(in) :: start, ncount
  integer, dimension(:,:), intent(out) :: ip
  integer (C_INT) :: c_ncid, c_varid
  integer (C_SIZE_T), dimension(size(start)) :: c_start, c_ncount
  integer (C_INT), dimension(size(ip,1),size(ip,2)), target :: c_ip
  integer i, ix
  c_ncid = ncid
  c_varid = varid
  ix = size(start)
  do i = 1,ix
    c_start(ix-i+1) = start(i) - 1
    c_ncount(ix-i+1) = ncount(i)
  end do
  ierr = nc_get_vara_int(c_ncid,c_varid,c_start,c_ncount,C_LOC(c_ip))
  ip = c_ip
end function nf_get_vara_int_d2

integer function nf_get_vara_int_d3(ncid,varid,start,ncount,ip) result(ierr)
  implicit none
  integer, intent(in) :: ncid, varid
  integer, dimension(:), intent(in) :: start, ncount
  integer, dimension(:,:,:), intent(out) :: ip
  integer (C_INT) :: c_ncid, c_varid
  integer (C_SIZE_T), dimension(size(start)) :: c_start, c_ncount
  integer (C_INT), dimension(size(ip,1),size(ip,2),size(ip,3)), target :: c_ip
  integer i, ix
  c_ncid = ncid
  c_varid = varid
  ix = size(start)
  do i = 1,ix
    c_start(ix-i+1) = start(i) - 1
    c_ncount(ix-i+1) = ncount(i)
  end do
  ierr = nc_get_vara_int(c_ncid,c_varid,c_start,c_ncount,C_LOC(c_ip))
  ip = c_ip
end function nf_get_vara_int_d3

integer function nf_get_vara_int_d4(ncid,varid,start,ncount,ip) result(ierr)
  implicit none
  integer, intent(in) :: ncid, varid
  integer, dimension(:), intent(in) :: start, ncount
  integer, dimension(:,:,:,:), intent(out) :: ip
  integer (C_INT) :: c_ncid, c_varid
  integer (C_SIZE_T), dimension(size(start)) :: c_start, c_ncount
  integer (C_INT), dimension(size(ip,1),size(ip,2),size(ip,3),size(ip,4)), target :: c_ip
  integer i, ix
  c_ncid = ncid
  c_varid = varid
  ix = size(start)
  do i = 1,ix
    c_start(ix-i+1) = start(i) - 1
    c_ncount(ix-i+1) = ncount(i)
  end do
  ierr = nc_get_vara_int(c_ncid,c_varid,c_start,c_ncount,C_LOC(c_ip))
  ip = c_ip
end function nf_get_vara_int_d4

integer function nf_get_vara_int2_d1(ncid,varid,start,ncount,sp) result(ierr)
  implicit none
  integer, intent(in) :: ncid, varid
  integer, dimension(:), intent(in) :: start, ncount
  integer(kind=2), dimension(:), intent(out) :: sp
  integer (C_INT) :: c_ncid, c_varid
  integer (C_SIZE_T), dimension(size(start)) :: c_start, c_ncount
  integer (C_SHORT), dimension(size(sp)), target :: c_sp
  integer i, ix
  c_ncid = ncid
  c_varid = varid
  ix = size(start)
  do i = 1,ix
    c_start(ix-i+1) = start(i) - 1
    c_ncount(ix-i+1) = ncount(i)
  end do
  ierr = nc_get_vara_short(c_ncid,c_varid,c_start,c_ncount,C_LOC(c_sp))
  sp = c_sp
end function nf_get_vara_int2_d1

integer function nf_get_vara_int2_d2(ncid,varid,start,ncount,sp) result(ierr)
  implicit none
  integer, intent(in) :: ncid, varid
  integer, dimension(:), intent(in) :: start, ncount
  integer(kind=2), dimension(:,:), intent(out) :: sp
  integer (C_INT) :: c_ncid, c_varid
  integer (C_SIZE_T), dimension(size(start)) :: c_start, c_ncount
  integer (C_SHORT), dimension(size(sp,1),size(sp,2)), target :: c_sp
  integer i, ix
  c_ncid = ncid
  c_varid = varid
  ix = size(start)
  do i = 1,ix
    c_start(ix-i+1) = start(i) - 1
    c_ncount(ix-i+1) = ncount(i)
  end do
  ierr = nc_get_vara_short(c_ncid,c_varid,c_start,c_ncount,C_LOC(c_sp))
  sp = c_sp
end function nf_get_vara_int2_d2

integer function nf_get_vara_int2_d3(ncid,varid,start,ncount,sp) result(ierr)
  implicit none
  integer, intent(in) :: ncid, varid
  integer, dimension(:), intent(in) :: start, ncount
  integer(kind=2), dimension(:,:,:), intent(out) :: sp
  integer (C_INT) :: c_ncid, c_varid
  integer (C_SIZE_T), dimension(size(start)) :: c_start, c_ncount
  integer (C_SHORT), dimension(size(sp,1),size(sp,2),size(sp,3)), target :: c_sp
  integer i, ix
  c_ncid = ncid
  c_varid = varid
  ix = size(start)
  do i = 1,ix
    c_start(ix-i+1) = start(i) - 1
    c_ncount(ix-i+1) = ncount(i)
  end do
  ierr = nc_get_vara_short(c_ncid,c_varid,c_start,c_ncount,C_LOC(c_sp))
  sp = c_sp
end function nf_get_vara_int2_d3

integer function nf_get_vara_int2_d4(ncid,varid,start,ncount,sp) result(ierr)
  implicit none
  integer, intent(in) :: ncid, varid
  integer, dimension(:), intent(in) :: start, ncount
  integer(kind=2), dimension(:,:,:,:), intent(out) :: sp
  integer (C_INT) :: c_ncid, c_varid
  integer (C_SIZE_T), dimension(size(start)) :: c_start, c_ncount
  integer (C_SHORT), dimension(size(sp,1),size(sp,2),size(sp,3),size(sp,4)), target :: c_sp
  integer i, ix
  c_ncid = ncid
  c_varid = varid
  ix = size(start)
  do i = 1,ix
    c_start(ix-i+1) = start(i) - 1
    c_ncount(ix-i+1) = ncount(i)
  end do
  ierr = nc_get_vara_short(c_ncid,c_varid,c_start,c_ncount,C_LOC(c_sp))
  sp = c_sp
end function nf_get_vara_int2_d4

integer function nf_get_vara_double_d1(ncid,varid,start,ncount,dp) result(ierr)
  implicit none
  integer, intent(in) :: ncid, varid
  integer, dimension(:), intent(in) :: start, ncount
  real(kind=8), dimension(:), intent(out) :: dp
  integer (C_INT) :: c_ncid, c_varid
  integer (C_SIZE_T), dimension(size(start)) :: c_start, c_ncount
  real (C_DOUBLE), dimension(size(dp)), target :: c_dp
  integer i, ix
  c_ncid = ncid
  c_varid = varid
  ix = size(start)
  do i = 1,ix
    c_start(ix-i+1) = start(i) - 1
    c_ncount(ix-i+1) = ncount(i)
  end do
  ierr = nc_get_vara_double(c_ncid,c_varid,c_start,c_ncount,C_LOC(c_dp))
  dp = c_dp
end function nf_get_vara_double_d1

integer function nf_get_vara_double_d2(ncid,varid,start,ncount,dp) result(ierr)
  implicit none
  integer, intent(in) :: ncid, varid
  integer, dimension(:), intent(in) :: start, ncount
  real(kind=8), dimension(:,:), intent(out) :: dp
  integer (C_INT) :: c_ncid, c_varid
  integer (C_SIZE_T), dimension(size(start)) :: c_start, c_ncount
  real (C_DOUBLE), dimension(size(dp,1),size(dp,2)), target :: c_dp
  integer i, ix
  c_ncid = ncid
  c_varid = varid
  ix = size(start)
  do i = 1,ix
    c_start(ix-i+1) = start(i) - 1
    c_ncount(ix-i+1) = ncount(i)
  end do
  ierr = nc_get_vara_double(c_ncid,c_varid,c_start,c_ncount,C_LOC(c_dp))
  dp = c_dp
end function nf_get_vara_double_d2

integer function nf_get_vara_double_d3(ncid,varid,start,ncount,dp) result(ierr)
  implicit none
  integer, intent(in) :: ncid, varid
  integer, dimension(:), intent(in) :: start, ncount
  real(kind=8), dimension(:,:,:), intent(out) :: dp
  integer (C_INT) :: c_ncid, c_varid
  integer (C_SIZE_T), dimension(size(start)) :: c_start, c_ncount
  real (C_DOUBLE), dimension(size(dp,1),size(dp,2),size(dp,3)), target :: c_dp
  integer i, ix
  c_ncid = ncid
  c_varid = varid
  ix = size(start)
  do i = 1,ix
    c_start(ix-i+1) = start(i) - 1
    c_ncount(ix-i+1) = ncount(i)
  end do
  ierr = nc_get_vara_double(c_ncid,c_varid,c_start,c_ncount,C_LOC(c_dp))
  dp = c_dp
end function nf_get_vara_double_d3

integer function nf_get_vara_double_d4(ncid,varid,start,ncount,dp) result(ierr)
  implicit none
  integer, intent(in) :: ncid, varid
  integer, dimension(:), intent(in) :: start, ncount
  real(kind=8), dimension(:,:,:,:), intent(out) :: dp
  integer (C_INT) :: c_ncid, c_varid
  integer (C_SIZE_T), dimension(size(start)) :: c_start, c_ncount
  real (C_DOUBLE), dimension(size(dp,1),size(dp,2),size(dp,3),size(dp,4)), target :: c_dp
  integer i, ix
  c_ncid = ncid
  c_varid = varid
  ix = size(start)
  do i = 1,ix
    c_start(ix-i+1) = start(i) - 1
    c_ncount(ix-i+1) = ncount(i)
  end do
  ierr = nc_get_vara_double(c_ncid,c_varid,c_start,c_ncount,C_LOC(c_dp))
  dp = c_dp
end function nf_get_vara_double_d4

integer function nf_get_var_real_d1(ncid,varid,fp) result(ierr)
  implicit none
  integer, intent(in) :: ncid, varid
  real, dimension(:), intent(out) :: fp
  integer (C_INT) :: c_ncid, c_varid
  real (C_FLOAT), dimension(size(fp)), target :: c_fp
  c_ncid = ncid
  c_varid = varid
  ierr = nc_get_var_float(c_ncid,c_varid,C_LOC(c_fp))
  fp = c_fp
end function nf_get_var_real_d1

integer function nf_get_var_real_d2(ncid,varid,fp) result(ierr)
  implicit none
  integer, intent(in) :: ncid, varid
  real, dimension(:,:), intent(out) :: fp
  integer (C_INT) :: c_ncid, c_varid
  real (C_FLOAT), dimension(size(fp,1),size(fp,2)), target :: c_fp
  c_ncid = ncid
  c_varid = varid
  ierr = nc_get_var_float(c_ncid,c_varid,C_LOC(c_fp))
  fp = c_fp
end function nf_get_var_real_d2

integer function nf_get_var_real_d3(ncid,varid,fp) result(ierr)
  implicit none
  integer, intent(in) :: ncid, varid
  real, dimension(:,:,:), intent(out) :: fp
  integer (C_INT) :: c_ncid, c_varid
  real (C_FLOAT), dimension(size(fp,1),size(fp,2),size(fp,3)), target :: c_fp
  c_ncid = ncid
  c_varid = varid
  ierr = nc_get_var_float(c_ncid,c_varid,C_LOC(c_fp))
  fp = c_fp
end function nf_get_var_real_d3

integer function nf_get_var_real_d4(ncid,varid,fp) result(ierr)
  implicit none
  integer, intent(in) :: ncid, varid
  real, dimension(:,:,:,:), intent(out) :: fp
  integer (C_INT) :: c_ncid, c_varid
  real (C_FLOAT), dimension(size(fp,1),size(fp,2),size(fp,3),size(fp,4)), target :: c_fp
  c_ncid = ncid
  c_varid = varid
  ierr = nc_get_var_float(c_ncid,c_varid,C_LOC(c_fp))
  fp = c_fp
end function nf_get_var_real_d4

integer function nf_get_var1_real_s(ncid,varid,start,fp) result(ierr)
  implicit none
  integer, intent(in) :: ncid, varid
  integer, intent(in) :: start
  real, intent(out) :: fp
  integer (C_INT) :: c_ncid, c_varid
  integer (C_SIZE_T), dimension(1) :: c_start
  real (C_FLOAT), target :: c_fp
  c_ncid = ncid
  c_varid = varid
  c_start = start - 1
  ierr = nc_get_var1_float(c_ncid,c_varid,c_start,C_LOC(c_fp))
  fp = c_fp
end function nf_get_var1_real_s

integer function nf_get_var1_real_v(ncid,varid,start,fp) result(ierr)
  implicit none
  integer, intent(in) :: ncid, varid
  integer, dimension(:), intent(in) :: start
  real, intent(out) :: fp
  integer (C_INT) :: c_ncid, c_varid
  integer (C_SIZE_T), dimension(size(start)) :: c_start
  real (C_FLOAT), target :: c_fp
  integer i, ix
  c_ncid = ncid
  c_varid = varid
  ix = size(start)
  do i = 1,ix
    c_start(ix-i+1) = start(i) - 1
  end do
  ierr = nc_get_var1_float(c_ncid,c_varid,c_start,C_LOC(c_fp))
  fp = c_fp
end function nf_get_var1_real_v

integer function nf_get_var1_int_s(ncid,varid,start,ip) result(ierr)
  implicit none
  integer, intent(in) :: ncid, varid
  integer, intent(in) :: start
  integer, intent(out) :: ip
  integer (C_INT) :: c_ncid, c_varid
  integer (C_SIZE_T), dimension(1) :: c_start
  integer (C_INT), target :: c_ip
  c_ncid = ncid
  c_varid = varid
  c_start = start - 1
  ierr = nc_get_var1_int(c_ncid,c_varid,c_start,C_LOC(c_ip))
  ip = c_ip
end function nf_get_var1_int_s

integer function nf_get_var1_int_v(ncid,varid,start,ip) result(ierr)
  implicit none
  integer, intent(in) :: ncid, varid
  integer, dimension(:), intent(in) :: start
  integer, intent(out) :: ip
  integer (C_INT) :: c_ncid, c_varid
  integer (C_SIZE_T), dimension(size(start)) :: c_start
  integer (C_INT), target :: c_ip
  integer i, ix
  c_ncid = ncid
  c_varid = varid
  ix = size(start)
  do i = 1,ix
    c_start(ix-i+1) = start(i) - 1
  end do
  ierr = nc_get_var1_int(c_ncid,c_varid,c_start,C_LOC(c_ip))
  ip = c_ip
end function nf_get_var1_int_v

integer function nf_def_dim(ncid,name,size,dimid) result(ierr)
  implicit none
  integer, intent(in) :: ncid, size
  integer, intent(out) :: dimid
  character(len=*), intent(in) :: name
  integer (C_INT) :: c_ncid
  integer (C_SIZE_T) :: c_size
  integer (C_INT), target :: c_dimid
  character, dimension(charsize) :: c_name
  integer i, ix
  if (len(name)>charsize) then
    write(6,*) "ERROR: charsize in netcdf_m is too small"
    write(6,*) "Req ",len(name)," Actual ",charsize
    stop
  end if
  c_ncid = ncid
  ix = len_trim(name)
  do i=1,len_trim(name)
    c_name(i) = name(i:i)
  end do
  c_name(ix+1:ix+1) = C_NULL_CHAR
  c_size = size
  ierr = nc_def_dim(c_ncid,c_name,c_size,C_LOC(c_dimid))
  dimid = c_dimid
end function nf_def_dim

integer function nf_def_var(ncid,name,xtype,ndims,dimids,varid) result(ierr)
  implicit none
  integer, intent(in) :: ncid, xtype, ndims
  integer, intent(out) :: varid
  integer, dimension(ndims), intent(in) :: dimids
  character(len=*), intent(in) :: name
  integer (C_INT) :: c_ncid, c_xtype, c_ndims
  integer (C_INT), target :: c_varid
  integer (C_INT), dimension(ndims) :: c_dimids
  character, dimension(charsize) :: c_name
  integer i, ix
  if (len(name)>charsize) then
    write(6,*) "ERROR: charsize in netcdf_m is too small"
    write(6,*) "Req ",len(name)," Actual ",charsize
    stop
  end if
  c_ncid = ncid
  ix = len_trim(name)
  do i=1,len_trim(name)
    c_name(i) = name(i:i)
  end do
  c_name(ix+1:ix+1) = C_NULL_CHAR
  c_xtype = xtype
  c_ndims = ndims
  do i=1,ndims
    c_dimids(ndims-i+1) = dimids(i)
  end do
  ierr = nc_def_var(c_ncid,c_name,c_xtype,c_ndims,c_dimids,C_LOC(c_varid))
  varid = c_varid
end function nf_def_var

integer function nf_rename_dim(ncid,dimid,name) result(ierr)
  implicit none
  integer, intent(in) :: ncid, dimid
  character(len=*), intent(in) :: name
  integer (C_INT) :: c_ncid, c_dimid
  character, dimension(charsize) :: c_name
  integer i, ix
  if (len(name)>charsize) then
    write(6,*) "ERROR: charsize in netcdf_m is too small"
    write(6,*) "Req ",len(name)," Actual ",charsize
    stop
  end if
  c_ncid = ncid
  c_dimid = dimid
  ix = len_trim(name)
  do i=1,len_trim(name)
    c_name(i) = name(i:i)
  end do
  c_name(ix+1:ix+1) = C_NULL_CHAR
  ierr = nc_rename_dim(c_ncid,c_dimid,c_name)
end function nf_rename_dim

integer function nf_put_att_text(ncid,varid,name,clen,tp) result(ierr)
  implicit none
  integer, intent(in) :: ncid, varid, clen
  character(len=*), intent(in) :: name
  character(len=*), intent(in) :: tp
  integer (C_INT) :: c_ncid, c_varid
  integer (C_SIZE_T) :: c_clen
  character, dimension(charsize) :: c_name
  character, dimension(charsize) :: c_tp
  integer i, ix
  if (len(name)>charsize) then
    write(6,*) "ERROR: charsize in netcdf_m is too small"
    write(6,*) "Req ",len(name)," Actual ",charsize
    stop
  end if
  if (len(tp)>charsize) then
    write(6,*) "ERROR: charsize in netcdf_m is too small"
    write(6,*) "Req ",len(tp)," Actual ",charsize
    stop
  end if  
  c_ncid = ncid
  c_varid = varid
  ix = len_trim(name)
  do i=1,len_trim(name)
    c_name(i) = name(i:i)
  end do
  c_name(ix+1:ix+1) = C_NULL_CHAR
  c_clen = clen
  ix = len_trim(tp)
  do i=1,ix
    c_tp(i) = tp(i:i)
  end do
  c_tp(ix+1:ix+1) = C_NULL_CHAR
  ierr = nc_put_att_text(c_ncid,c_varid,c_name,c_clen,c_tp)
end function nf_put_att_text

integer function nf_put_att_int2(ncid,varid,name,xtype,slen,sp) result(ierr)
  implicit none
  integer, intent(in) :: ncid, varid, xtype, slen
  integer(kind=2), dimension(slen), intent(in) :: sp
  character(len=*), intent(in) :: name
  integer (C_INT) :: c_ncid, c_varid, c_xtype
  integer (C_SIZE_T) :: c_slen
  integer (C_SHORT), dimension(slen), target :: c_sp  
  character, dimension(charsize) :: c_name
  integer i, ix
  if (len(name)>charsize) then
    write(6,*) "ERROR: charsize in netcdf_m is too small"
    write(6,*) "Req ",len(name)," Actual ",charsize
    stop
  end if
  c_ncid = ncid
  c_varid = varid
  ix = len_trim(name)
  do i=1,len_trim(name)
    c_name(i) = name(i:i)
  end do
  c_name(ix+1:ix+1) = C_NULL_CHAR
  c_xtype = xtype
  c_slen = slen
  c_sp = sp
  ierr = nc_put_att_short(c_ncid,c_varid,c_name,c_xtype,c_slen,C_LOC(c_sp))
end function nf_put_att_int2

integer function nf_put_att_real(ncid,varid,name,xtype,rlen,fp) result(ierr)
  implicit none
  integer, intent(in) :: ncid, varid, xtype, rlen
  real, dimension(rlen), intent(in) :: fp
  character(len=*), intent(in) :: name
  integer (C_INT) :: c_ncid, c_varid, c_xtype
  integer (C_SIZE_T) :: c_rlen
  real (C_FLOAT), dimension(rlen), target :: c_fp  
  character, dimension(charsize) :: c_name
  integer i, ix
  if (len(name)>charsize) then
    write(6,*) "ERROR: charsize in netcdf_m is too small"
    write(6,*) "Req ",len(name)," Actual ",charsize
    stop
  end if
  c_ncid = ncid
  c_varid = varid
  ix = len_trim(name)
  do i=1,len_trim(name)
    c_name(i) = name(i:i)
  end do
  c_name(ix+1:ix+1) = C_NULL_CHAR
  c_xtype = xtype
  c_rlen = rlen
  c_fp = fp
  ierr = nc_put_att_float(c_ncid,c_varid,c_name,c_xtype,c_rlen,C_LOC(c_fp))
end function nf_put_att_real

integer function nf_put_att_int(ncid,varid,name,xtype,slen,ip) result(ierr)
  implicit none
  integer, intent(in) :: ncid, varid, xtype, slen
  integer, dimension(slen), intent(in) :: ip
  character(len=*), intent(in) :: name
  integer (C_INT) :: c_ncid, c_varid, c_xtype
  integer (C_SIZE_T) :: c_slen
  integer (C_INT), dimension(slen), target :: c_ip  
  character, dimension(charsize) :: c_name
  integer i, ix
  if (len(name)>charsize) then
    write(6,*) "ERROR: charsize in netcdf_m is too small"
    write(6,*) "Req ",len(name)," Actual ",charsize
    stop
  end if
  c_ncid = ncid
  c_varid = varid
  ix = len_trim(name)
  do i=1,len_trim(name)
    c_name(i) = name(i:i)
  end do
  c_name(ix+1:ix+1) = C_NULL_CHAR
  c_xtype = xtype
  c_slen = slen
  c_ip = ip
  ierr = nc_put_att_int(c_ncid,c_varid,c_name,c_xtype,c_slen,C_LOC(c_ip))
end function nf_put_att_int

integer function nf_put_vara_real_d1(ncid,varid,start,ncount,fp) result(ierr)
  implicit none
  integer, intent(in) :: ncid, varid
  integer, dimension(:), intent(in) :: start, ncount
  real, dimension(:), intent(in) :: fp
  integer (C_INT) :: c_ncid, c_varid
  integer (C_SIZE_T), dimension(size(start)) :: c_start, c_ncount
  real (C_FLOAT), dimension(size(fp)), target :: c_fp
  integer i, ix
  c_ncid = ncid
  c_varid = varid
  ix = size(start)
  do i = 1,ix
    c_start(ix-i+1) = start(i) - 1
    c_ncount(ix-i+1) = ncount(i)
  end do
  c_fp = fp
  ierr = nc_put_vara_float(c_ncid,c_varid,c_start,c_ncount,C_LOC(c_fp))
end function nf_put_vara_real_d1

integer function nf_put_vara_real_d2(ncid,varid,start,ncount,fp) result(ierr)
  implicit none
  integer, intent(in) :: ncid, varid
  integer, dimension(:), intent(in) :: start, ncount
  real, dimension(:,:), intent(in) :: fp
  integer (C_INT) :: c_ncid, c_varid
  integer (C_SIZE_T), dimension(size(start)) :: c_start, c_ncount
  real (C_FLOAT), dimension(size(fp,1),size(fp,2)), target :: c_fp
  integer i, ix
  c_ncid = ncid
  c_varid = varid
  ix = size(start)
  do i = 1,ix
    c_start(ix-i+1) = start(i) - 1
    c_ncount(ix-i+1) = ncount(i)
  end do
  c_fp = fp
  ierr = nc_put_vara_float(c_ncid,c_varid,c_start,c_ncount,C_LOC(c_fp))
end function nf_put_vara_real_d2

integer function nf_put_vara_real_d3(ncid,varid,start,ncount,fp) result(ierr)
  implicit none
  integer, intent(in) :: ncid, varid
  integer, dimension(:), intent(in) :: start, ncount
  real, dimension(:,:,:), intent(in) :: fp
  integer (C_INT) :: c_ncid, c_varid
  integer (C_SIZE_T), dimension(size(start)) :: c_start, c_ncount
  real (C_FLOAT), dimension(size(fp,1),size(fp,2),size(fp,3)), target :: c_fp
  integer i, ix
  c_ncid = ncid
  c_varid = varid
  ix = size(start)
  do i = 1,ix
    c_start(ix-i+1) = start(i) - 1
    c_ncount(ix-i+1) = ncount(i)
  end do
  c_fp = fp
  ierr = nc_put_vara_float(c_ncid,c_varid,c_start,c_ncount,C_LOC(c_fp))
end function nf_put_vara_real_d3

integer function nf_put_vara_real_d4(ncid,varid,start,ncount,fp) result(ierr)
  implicit none
  integer, intent(in) :: ncid, varid
  integer, dimension(:), intent(in) :: start, ncount
  real, dimension(:,:,:,:), intent(in) :: fp
  integer (C_INT) :: c_ncid, c_varid
  integer (C_SIZE_T), dimension(size(start)) :: c_start, c_ncount
  real (C_FLOAT), dimension(size(fp,1),size(fp,2),size(fp,3),size(fp,4)), target :: c_fp
  integer i, ix
  c_ncid = ncid
  c_varid = varid
  ix = size(start)
  do i = 1,ix
    c_start(ix-i+1) = start(i) - 1
    c_ncount(ix-i+1) = ncount(i)
  end do
  c_fp = fp
  ierr = nc_put_vara_float(c_ncid,c_varid,c_start,c_ncount,C_LOC(c_fp))
end function nf_put_vara_real_d4

integer function nf_put_vara_int_d1(ncid,varid,start,ncount,ip) result(ierr)
  implicit none
  integer, intent(in) :: ncid, varid
  integer, dimension(:), intent(in) :: start, ncount
  integer, dimension(:), intent(in) :: ip
  integer (C_INT) :: c_ncid, c_varid
  integer (C_SIZE_T), dimension(size(start)) :: c_start, c_ncount
  integer (C_INT), dimension(size(ip)), target :: c_ip
  integer i, ix
  c_ncid = ncid
  c_varid = varid
  ix = size(start)
  do i = 1,ix
    c_start(ix-i+1) = start(i) - 1
    c_ncount(ix-i+1) = ncount(i)
  end do
  c_ip = ip
  ierr = nc_put_vara_int(c_ncid,c_varid,c_start,c_ncount,C_LOC(c_ip))
end function nf_put_vara_int_d1

integer function nf_put_vara_int_d2(ncid,varid,start,ncount,ip) result(ierr)
  implicit none
  integer, intent(in) :: ncid, varid
  integer, dimension(:), intent(in) :: start, ncount
  integer, dimension(:,:), intent(in) :: ip
  integer (C_INT) :: c_ncid, c_varid
  integer (C_SIZE_T), dimension(size(start)) :: c_start, c_ncount
  integer (C_INT), dimension(size(ip,1),size(ip,2)), target :: c_ip
  integer i, ix
  c_ncid = ncid
  c_varid = varid
  ix = size(start)
  do i = 1,ix
    c_start(ix-i+1) = start(i) - 1
    c_ncount(ix-i+1) = ncount(i)
  end do
  c_ip = ip
  ierr = nc_put_vara_int(c_ncid,c_varid,c_start,c_ncount,C_LOC(c_ip))
end function nf_put_vara_int_d2

integer function nf_put_vara_int_d3(ncid,varid,start,ncount,ip) result(ierr)
  implicit none
  integer, intent(in) :: ncid, varid
  integer, dimension(:), intent(in) :: start, ncount
  integer, dimension(:,:,:), intent(in) :: ip
  integer (C_INT) :: c_ncid, c_varid
  integer (C_SIZE_T), dimension(size(start)) :: c_start, c_ncount
  integer (C_INT), dimension(size(ip,1),size(ip,2),size(ip,3)), target :: c_ip
  integer i, ix
  c_ncid = ncid
  c_varid = varid
  ix = size(start)
  do i = 1,ix
    c_start(ix-i+1) = start(i) - 1
    c_ncount(ix-i+1) = ncount(i)
  end do
  c_ip = ip
  ierr = nc_put_vara_int(c_ncid,c_varid,c_start,c_ncount,C_LOC(c_ip))
end function nf_put_vara_int_d3

integer function nf_put_vara_int_d4(ncid,varid,start,ncount,ip) result(ierr)
  implicit none
  integer, intent(in) :: ncid, varid
  integer, dimension(:), intent(in) :: start, ncount
  integer, dimension(:,:,:,:), intent(in) :: ip
  integer (C_INT) :: c_ncid, c_varid
  integer (C_SIZE_T), dimension(size(start)) :: c_start, c_ncount
  integer (C_INT), dimension(size(ip,1),size(ip,2),size(ip,3),size(ip,4)), target :: c_ip
  integer i, ix
  c_ncid = ncid
  c_varid = varid
  ix = size(start)
  do i = 1,ix
    c_start(ix-i+1) = start(i) - 1
    c_ncount(ix-i+1) = ncount(i)
  end do
  c_ip = ip
  ierr = nc_put_vara_int(c_ncid,c_varid,c_start,c_ncount,C_LOC(c_ip))
end function nf_put_vara_int_d4

integer function nf_put_vara_int2_d1(ncid,varid,start,ncount,sp) result(ierr)
  implicit none
  integer, intent(in) :: ncid, varid
  integer, dimension(:), intent(in) :: start, ncount
  integer(kind=2), dimension(:), intent(in) :: sp
  integer (C_INT) :: c_ncid, c_varid
  integer (C_SIZE_T), dimension(size(start)) :: c_start, c_ncount
  integer (C_SHORT), dimension(size(sp)), target :: c_sp
  integer i, ix
  c_ncid = ncid
  c_varid = varid
  ix = size(start)
  do i = 1,ix
    c_start(ix-i+1) = start(i) - 1
    c_ncount(ix-i+1) = ncount(i)
  end do
  c_sp = sp
  ierr = nc_put_vara_short(c_ncid,c_varid,c_start,c_ncount,C_LOC(c_sp))
end function nf_put_vara_int2_d1

integer function nf_put_vara_int2_d2(ncid,varid,start,ncount,sp) result(ierr)
  implicit none
  integer, intent(in) :: ncid, varid
  integer, dimension(:), intent(in) :: start, ncount
  integer(kind=2), dimension(:,:), intent(in) :: sp
  integer (C_INT) :: c_ncid, c_varid
  integer (C_SIZE_T), dimension(size(start)) :: c_start, c_ncount
  integer (C_SHORT), dimension(size(sp,1),size(sp,2)), target :: c_sp
  integer i, ix
  c_ncid = ncid
  c_varid = varid
  ix = size(start)
  do i = 1,ix
    c_start(ix-i+1) = start(i) - 1
    c_ncount(ix-i+1) = ncount(i)
  end do
  c_sp = sp
  ierr = nc_put_vara_short(c_ncid,c_varid,c_start,c_ncount,C_LOC(c_sp))
end function nf_put_vara_int2_d2

integer function nf_put_vara_int2_d3(ncid,varid,start,ncount,sp) result(ierr)
  implicit none
  integer, intent(in) :: ncid, varid
  integer, dimension(:), intent(in) :: start, ncount
  integer(kind=2), dimension(:,:,:), intent(in) :: sp
  integer (C_INT) :: c_ncid, c_varid
  integer (C_SIZE_T), dimension(size(start)) :: c_start, c_ncount
  integer (C_SHORT), dimension(size(sp,1),size(sp,2),size(sp,3)), target :: c_sp
  integer i, ix
  c_ncid = ncid
  c_varid = varid
  ix = size(start)
  do i = 1,ix
    c_start(ix-i+1) = start(i) - 1
    c_ncount(ix-i+1) = ncount(i)
  end do
  c_sp = sp
  ierr = nc_put_vara_short(c_ncid,c_varid,c_start,c_ncount,C_LOC(c_sp))
end function nf_put_vara_int2_d3

integer function nf_put_vara_int2_d4(ncid,varid,start,ncount,sp) result(ierr)
  implicit none
  integer, intent(in) :: ncid, varid
  integer, dimension(:), intent(in) :: start, ncount
  integer(kind=2), dimension(:,:,:,:), intent(in) :: sp
  integer (C_INT) :: c_ncid, c_varid
  integer (C_SIZE_T), dimension(size(start)) :: c_start, c_ncount
  integer (C_SHORT), dimension(size(sp,1),size(sp,2),size(sp,3),size(sp,4)), target :: c_sp
  integer i, ix
  c_ncid = ncid
  c_varid = varid
  ix = size(start)
  do i = 1,ix
    c_start(ix-i+1) = start(i) - 1
    c_ncount(ix-i+1) = ncount(i)
  end do
  c_sp = sp
  ierr = nc_put_vara_short(c_ncid,c_varid,c_start,c_ncount,C_LOC(c_sp))
end function nf_put_vara_int2_d4

integer function nf_put_vara_double_d1(ncid,varid,start,ncount,dp) result(ierr)
  implicit none
  integer, intent(in) :: ncid, varid
  integer, dimension(:), intent(in) :: start, ncount
  real(kind=8), dimension(:), intent(in) :: dp
  integer (C_INT) :: c_ncid, c_varid
  integer (C_SIZE_T), dimension(size(start)) :: c_start, c_ncount
  real (C_DOUBLE), dimension(size(dp)), target :: c_dp
  integer i, ix
  c_ncid = ncid
  c_varid = varid
  ix = size(start)
  do i = 1,ix
    c_start(ix-i+1) = start(i) - 1
    c_ncount(ix-i+1) = ncount(i)
  end do
  c_dp = dp
  ierr = nc_put_vara_double(c_ncid,c_varid,c_start,c_ncount,C_LOC(c_dp))
end function nf_put_vara_double_d1

integer function nf_put_vara_double_d2(ncid,varid,start,ncount,dp) result(ierr)
  implicit none
  integer, intent(in) :: ncid, varid
  integer, dimension(:), intent(in) :: start, ncount
  real(kind=8), dimension(:,:), intent(in) :: dp
  integer (C_INT) :: c_ncid, c_varid
  integer (C_SIZE_T), dimension(size(start)) :: c_start, c_ncount
  real (C_DOUBLE), dimension(size(dp,1),size(dp,2)), target :: c_dp
  integer i, ix
  c_ncid = ncid
  c_varid = varid
  ix = size(start)
  do i = 1,ix
    c_start(ix-i+1) = start(i) - 1
    c_ncount(ix-i+1) = ncount(i)
  end do
  c_dp = dp
  ierr = nc_put_vara_double(c_ncid,c_varid,c_start,c_ncount,C_LOC(c_dp))
end function nf_put_vara_double_d2

integer function nf_put_vara_double_d3(ncid,varid,start,ncount,dp) result(ierr)
  implicit none
  integer, intent(in) :: ncid, varid
  integer, dimension(:), intent(in) :: start, ncount
  real(kind=8), dimension(:,:,:), intent(in) :: dp
  integer (C_INT) :: c_ncid, c_varid
  integer (C_SIZE_T), dimension(size(start)) :: c_start, c_ncount
  real (C_DOUBLE), dimension(size(dp,1),size(dp,2),size(dp,3)), target :: c_dp
  integer i, ix
  c_ncid = ncid
  c_varid = varid
  ix = size(start)
  do i = 1,ix
    c_start(ix-i+1) = start(i) - 1
    c_ncount(ix-i+1) = ncount(i)
  end do
  c_dp = dp
  ierr = nc_put_vara_double(c_ncid,c_varid,c_start,c_ncount,C_LOC(c_dp))
end function nf_put_vara_double_d3

integer function nf_put_vara_double_d4(ncid,varid,start,ncount,dp) result(ierr)
  implicit none
  integer, intent(in) :: ncid, varid
  integer, dimension(:), intent(in) :: start, ncount
  real(kind=8), dimension(:,:,:,:), intent(in) :: dp
  integer (C_INT) :: c_ncid, c_varid
  integer (C_SIZE_T), dimension(size(start)) :: c_start, c_ncount
  real (C_DOUBLE), dimension(size(dp,1),size(dp,2),size(dp,3),size(dp,4)), target :: c_dp
  integer i, ix
  c_ncid = ncid
  c_varid = varid
  ix = size(start)
  do i = 1,ix
    c_start(ix-i+1) = start(i) - 1
    c_ncount(ix-i+1) = ncount(i)
  end do
  c_dp = dp
  ierr = nc_put_vara_double(c_ncid,c_varid,c_start,c_ncount,C_LOC(c_dp))
end function nf_put_vara_double_d4

integer function nf_put_var_real_d1(ncid,varid,fp) result(ierr)
  implicit none
  integer, intent(in) :: ncid, varid
  real, dimension(:), intent(in) :: fp
  integer (C_INT) :: c_ncid, c_varid
  real (C_FLOAT), dimension(size(fp)), target :: c_fp
  c_ncid = ncid
  c_varid = varid
  c_fp = fp
  ierr = nc_put_var_float(c_ncid,c_varid,C_LOC(c_fp))
end function nf_put_var_real_d1

integer function nf_put_var_real_d2(ncid,varid,fp) result(ierr)
  implicit none
  integer, intent(in) :: ncid, varid
  real, dimension(:,:), intent(in) :: fp
  integer (C_INT) :: c_ncid, c_varid
  real (C_FLOAT), dimension(size(fp,1),size(fp,2)), target :: c_fp
  c_ncid = ncid
  c_varid = varid
  c_fp = fp
  ierr = nc_put_var_float(c_ncid,c_varid,C_LOC(c_fp))
end function nf_put_var_real_d2

integer function nf_put_var_real_d3(ncid,varid,fp) result(ierr)
  implicit none
  integer, intent(in) :: ncid, varid
  real, dimension(:,:,:), intent(in) :: fp
  integer (C_INT) :: c_ncid, c_varid
  real (C_FLOAT), dimension(size(fp,1),size(fp,2),size(fp,3)), target :: c_fp
  c_ncid = ncid
  c_varid = varid
  c_fp = fp
  ierr = nc_put_var_float(c_ncid,c_varid,C_LOC(c_fp))
end function nf_put_var_real_d3

integer function nf_put_var_real_d4(ncid,varid,fp) result(ierr)
  implicit none
  integer, intent(in) :: ncid, varid
  real, dimension(:,:,:,:), intent(in) :: fp
  integer (C_INT) :: c_ncid, c_varid
  real (C_FLOAT), dimension(size(fp,1),size(fp,2),size(fp,3),size(fp,4)), target :: c_fp
  c_ncid = ncid
  c_varid = varid
  c_fp = fp
  ierr = nc_put_var_float(c_ncid,c_varid,C_LOC(c_fp))
end function nf_put_var_real_d4

integer function nf_put_var1_int(ncid,varid,start,ip) result(ierr)
  implicit none
  integer, intent(in) :: ncid, varid
  integer, dimension(:), intent(in) :: start
  integer, intent(in) :: ip
  integer (C_INT) :: c_ncid, c_varid
  integer (C_SIZE_T), dimension(size(start)) :: c_start
  integer (C_INT), target :: c_ip
  integer i, ix
  c_ncid = ncid
  c_varid = varid
  ix = size(start)
  do i = 1,ix
    c_start(ix-i+1) = start(i) - 1
  end do
  c_ip = ip
  ierr = nc_put_var1_int(c_ncid,c_varid,c_start,C_LOC(c_ip))
end function nf_put_var1_int

integer function nf_put_var1_real(ncid,varid,start,rp) result(ierr)
  implicit none
  integer, intent(in) :: ncid, varid
  integer, dimension(:), intent(in) :: start
  real, intent(in) :: rp
  integer (C_INT) :: c_ncid, c_varid
  integer (C_SIZE_T), dimension(size(start)) :: c_start
  real (C_FLOAT), target :: c_rp
  integer i, ix
  c_ncid = ncid
  c_varid = varid
  ix = size(start)
  do i = 1,ix
    c_start(ix-i+1) = start(i) - 1
  end do
  c_rp = rp
  ierr = nc_put_var1_float(c_ncid,c_varid,c_start,C_LOC(c_rp))
end function nf_put_var1_real

integer function nf_put_var1_double(ncid,varid,start,dp) result(ierr)
  implicit none
  integer, intent(in) :: ncid, varid
  integer, dimension(:), intent(in) :: start
  real(kind=8), intent(in) :: dp
  integer (C_INT) :: c_ncid, c_varid
  integer (C_SIZE_T), dimension(size(start)) :: c_start
  real (C_DOUBLE), target :: c_dp
  integer i, ix
  c_ncid = ncid
  c_varid = varid
  ix = size(start)
  do i = 1,ix
    c_start(ix-i+1) = start(i) - 1
  end do
  c_dp = dp
  ierr = nc_put_var1_double(c_ncid,c_varid,c_start,C_LOC(c_dp))
end function nf_put_var1_double

integer function nf_put_var1_text(ncid,varid,start,tp) result(ierr)
  implicit none
  integer, intent(in) :: ncid, varid
  integer, dimension(:), intent(in) :: start  
  character(len=*), intent(in) :: tp
  integer (C_INT) :: c_ncid, c_varid
  integer (C_SIZE_T), dimension(size(start)) :: c_start
  character, dimension(charsize) :: c_tp
  integer i, ix
  if (len(tp)>charsize) then
    write(6,*) "ERROR: charsize in netcdf_m is too small"
    write(6,*) "Req ",len(tp)," Actual ",charsize
    stop
  end if  
  c_ncid = ncid
  c_varid = varid
  ix = size(start)
  do i = 1,ix
    c_start(ix-i+1) = start(i) - 1
  end do
  ix = len_trim(tp)
  do i=1,ix
    c_tp(i) = tp(i:i)
  end do
  c_tp(ix+1:ix+1) = C_NULL_CHAR
  ierr = nc_put_var1_text(c_ncid,c_varid,c_start,c_tp)
end function nf_put_var1_text

integer function nf_copy_att(ncidin,varidin,name,ncidout,varidout) result(ierr)
  implicit none
  integer, intent(in) :: ncidin, varidin, ncidout, varidout
  character(len=*), intent(in) :: name
  integer (C_INT) :: c_ncidin, c_ncidout, c_varidin, c_varidout
  character, dimension(charsize) :: c_name
  integer i, ix
  if (len(name)>charsize) then
    write(6,*) "ERROR: charsize in netcdf_m is too small"
    write(6,*) "Req ",len(name)," Actual ",charsize
    stop
  end if
  c_ncidin = ncidin
  c_varidin = varidin
  c_ncidout = ncidout
  c_varidout = varidout
  ix = len_trim(name)
  do i=1,len_trim(name)
    c_name(i) = name(i:i)
  end do
  c_name(ix+1:ix+1) = C_NULL_CHAR
  ierr = nc_copy_att(c_ncidin,c_varidin,c_name,c_ncidout,c_varidout)
end function nf_copy_att

end module netcdf_m
#else
! Fortran 77 interface
module netcdf_m
public
include 'netcdf.inc'
end module netcdf_m
#endif
