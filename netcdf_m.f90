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

#ifdef ncclib
! C interface
module netcdf_m

use, intrinsic :: ISO_C_BINDING, only: C_SHORT, C_INT, C_FLOAT, C_DOUBLE, C_SIZE_T, C_LOC, C_NULL_CHAR, C_PTR, &
                                       C_SIGNED_CHAR, C_F_POINTER, C_INTPTR_T

implicit none

private
public nf_unlimited
public nf_noerr, nf_ebadid, nf_eexist, nf_einval, nf_enotindefine, nf_eindefine, nf_einvalcoords
public nf_emaxdims, nf_enameinuse, nf_enotatt, nf_emaxatts, nf_ebadtype, nf_ebaddim, nf_eunlimpos
public nf_emaxvars, nf_enotvar, nf_eglobal, nf_enotnc, nf_ests, nf_emaxname, nf_eunlimit
public nf_enorecvars, nf_echar, nf_eedge, nf_estride, nf_ebadname, nf_erange, nf_enomem
public nf_evarsize, nf_edimsize
public nf_nowrite, nf_write, nf_clobber, nf_noclobber, nf_fill, nf_nofill, nf_lock, nf_share
public nf_64bit_offset, nf_sizehint_default, nf_align_chunk, nf_format_classic, nf_format_64bit
public nf_global
public nf_byte, nf_int1, nf_char, nf_short, nf_int2, nf_int, nf_float, nf_real, nf_double
public nf_fill_byte, nf_fill_int1, nf_fill_char, nf_fill_short, nf_fill_int2, nf_fill_int
public nf_fill_float, nf_fill_real, nf_fill_double
public nf_max_dims, nf_max_attrs, nf_max_vars, nf_max_name, nf_max_var_dims
public nf_fatal, nf_verbose
public nf_open, nf_close, nf_create, nf_enddef, nf_set_fill, nf_redef, nf_sync, nf_strerror
public nf__open, nf__create, nf__enddef, nf_abort
public nf_inq_varndims, nf_inq_vardimid, nf_inq_dimlen, nf_inq_varid, nf_inq_dimid, nf_inq_vartype
public nf_inq, nf_inq_varname, nf_inq_ndims, nf_inq_nvars, nf_inq_libvers, nf_inq_dim
public nf_inq_attname, nf_inq_attid, nf_inq_att, nf_inq_var, nf_inq_attlen
public nf_get_att_text, nf_get_att_real, nf_get_att_int, nf_get_att_int1, nf_get_att_int2
public nf_get_att_double
public nf_get_vara_real, nf_get_vara_int, nf_get_vara_int2, nf_get_vara_double
public nf_get_var_real, nf_get_var_double
public nf_get_var1_real, nf_get_var1_int, nf_get_var1_int1, nf_get_var1_int2, nf_get_var1_double
public nf_get_vars_text, nf_get_vars_int1, nf_get_vars_int2, nf_get_vars_int, nf_get_vars_real
public nf_get_vars_double
public nf_get_varm_int1, nf_get_varm_int2, nf_get_varm_int, nf_get_varm_real, nf_get_varm_double
public nf_def_dim, nf_def_var
public nf_rename_dim, nf_rename_att, nf_rename_var
public nf_put_att_text, nf_put_att_int2, nf_put_att_real, nf_put_att_int, nf_put_att_int1
public nf_put_att_double
public nf_put_vara_real, nf_put_vara_int, nf_put_vara_int2, nf_put_vara_double
public nf_put_var_real
public nf_put_var1_int, nf_put_var1_real, nf_put_var1_double, nf_put_var1_text, nf_put_var1_int1
public nf_put_var1_int2
public nf_put_vars_text, nf_put_vars_int1, nf_put_vars_int2, nf_put_vars_int, nf_put_vars_real
public nf_put_vars_double
public nf_put_varm_int1, nf_put_varm_int2, nf_put_varm_int, nf_put_varm_real, nf_put_varm_double
public nf_copy_att, nf_del_att

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

integer (C_INT) function nc_inq_attname(ncid,varid,attnum,tp) bind(C, name='nc_inq_attname')
  use, intrinsic :: ISO_C_BINDING
  implicit none
  integer (C_INT), value :: ncid, varid, attnum
  character, dimension(*) :: tp
end function nc_inq_attname

integer (C_INT) function nc_inq_attid(ncid,varid,name,attnump) bind(C, name='nc_inq_attid')
  use, intrinsic :: ISO_C_BINDING
  implicit none
  integer (C_INT), value :: ncid, varid
  character, dimension(*) :: name
  type (C_PTR), value :: attnump
end function nc_inq_attid

integer (C_INT) function nc_inq_att(ncid,varid,name,xtypep,lenp) bind(C, name='nc_inq_att')
  use, intrinsic :: ISO_C_BINDING
  implicit none
  integer (C_INT), value :: ncid, varid
  character, dimension(*) :: name
  type (C_PTR), value :: xtypep
  type (C_PTR), value :: lenp
end function nc_inq_att

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

integer (C_INT) function nc_inq_var(ncid,varid,name,xtypep,ndimsp,dimidsp,nattsp) bind(C, name='nc_inq_var')
  use, intrinsic :: ISO_C_BINDING
  implicit none
  integer (C_INT), value :: ncid, varid
  character, dimension(*) :: name
  type (C_PTR), value :: xtypep, ndimsp, dimidsp, nattsp
end function nc_inq_var

integer (C_INT) function nc_inq_attlen(ncid,varid,name,lenp) bind(C, name='nc_inq_attlen')
  use, intrinsic :: ISO_C_BINDING
  implicit none
  integer (C_INT), value :: ncid, varid
  character, dimension(*) :: name
  type (C_PTR), value :: lenp
end function nc_inq_attlen
    
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

integer (C_INT) function nc_get_att_schar(ncid,varid,name,bp) bind(C, name='nc_get_att_schar')
  use, intrinsic :: ISO_C_BINDING
  implicit none
  integer (C_INT), value :: ncid, varid
  type (C_PTR), value :: bp
  character, dimension(*) :: name
end function nc_get_att_schar

integer (C_INT) function nc_get_att_short(ncid,varid,name,sp) bind(C, name='nc_get_att_short')
  use, intrinsic :: ISO_C_BINDING
  implicit none
  integer (C_INT), value :: ncid, varid
  type (C_PTR), value :: sp
  character, dimension(*) :: name
end function nc_get_att_short

integer (C_INT) function nc_get_att_double(ncid,varid,name,dp) bind(C, name='nc_get_att_double')
  use, intrinsic :: ISO_C_BINDING
  implicit none
  integer (C_INT), value :: ncid, varid
  type (C_PTR), value :: dp
  character, dimension(*) :: name
end function nc_get_att_double
    
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

integer (C_INT) function nc_get_var_double(ncid,varid,dp) bind(C, name='nc_get_var_double')
  use, intrinsic :: ISO_C_BINDING
  implicit none
  integer (C_INT), value :: ncid, varid
  type (C_PTR), value :: dp
end function nc_get_var_double    
    
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

integer (C_INT) function nc_get_var1_schar(ncid,varid,start,bp) bind(C, name='nc_get_var1_schar')
  use, intrinsic :: ISO_C_BINDING
  implicit none
  integer (C_INT), value :: ncid, varid
  integer (C_SIZE_T), dimension(*) :: start
  type (C_PTR), value :: bp
end function nc_get_var1_schar

integer (C_INT) function nc_get_var1_short(ncid,varid,start,sp) bind(C, name='nc_get_var1_short')
  use, intrinsic :: ISO_C_BINDING
  implicit none
  integer (C_INT), value :: ncid, varid
  integer (C_SIZE_T), dimension(*) :: start
  type (C_PTR), value :: sp
end function nc_get_var1_short

integer (C_INT) function nc_get_var1_double(ncid,varid,start,dp) bind(C, name='nc_get_var1_double')
  use, intrinsic :: ISO_C_BINDING
  implicit none
  integer (C_INT), value :: ncid, varid
  integer (C_SIZE_T), dimension(*) :: start
  type (C_PTR), value :: dp
end function nc_get_var1_double
    
integer (C_INT) function nc_get_vars_text(ncid,varid,start,ncount,stride,tp) bind(C, name='nc_get_vars_text')
  use, intrinsic :: ISO_C_BINDING
  implicit none
  integer (C_INT), value :: ncid, varid
  integer (C_SIZE_T), dimension(*) :: start, ncount
  integer (C_INTPTR_T), dimension(*) :: stride
  type (C_PTR), value :: tp
end function nc_get_vars_text

integer (C_INT) function nc_get_vars_schar(ncid,varid,start,ncount,stride,bp) bind(C, name='nc_get_vars_schar')
  use, intrinsic :: ISO_C_BINDING
  implicit none
  integer (C_INT), value :: ncid, varid
  integer (C_SIZE_T), dimension(*) :: start, ncount
  integer (C_INTPTR_T), dimension(*) :: stride
  type (C_PTR), value :: bp
end function nc_get_vars_schar

integer (C_INT) function nc_get_vars_short(ncid,varid,start,ncount,stride,sp) bind(C, name='nc_get_vars_short')
  use, intrinsic :: ISO_C_BINDING
  implicit none
  integer (C_INT), value :: ncid, varid
  integer (C_SIZE_T), dimension(*) :: start, ncount
  integer (C_INTPTR_T), dimension(*) :: stride
  type (C_PTR), value :: sp
end function nc_get_vars_short 

integer (C_INT) function nc_get_vars_int(ncid,varid,start,ncount,stride,ip) bind(C, name='nc_get_vars_int')
  use, intrinsic :: ISO_C_BINDING
  implicit none
  integer (C_INT), value :: ncid, varid
  integer (C_SIZE_T), dimension(*) :: start, ncount
  integer (C_INTPTR_T), dimension(*) :: stride
  type (C_PTR), value :: ip
end function nc_get_vars_int

integer (C_INT) function nc_get_vars_float(ncid,varid,start,ncount,stride,fp) bind(C, name='nc_get_vars_float')
  use, intrinsic :: ISO_C_BINDING
  implicit none
  integer (C_INT), value :: ncid, varid
  integer (C_SIZE_T), dimension(*) :: start, ncount
  integer (C_INTPTR_T), dimension(*) :: stride
  type (C_PTR), value :: fp
end function nc_get_vars_float    

integer (C_INT) function nc_get_vars_double(ncid,varid,start,ncount,stride,dp) bind(C, name='nc_get_vars_double')
  use, intrinsic :: ISO_C_BINDING
  implicit none
  integer (C_INT), value :: ncid, varid
  integer (C_SIZE_T), dimension(*) :: start, ncount
  integer (C_INTPTR_T), dimension(*) :: stride
  type (C_PTR), value :: dp
end function nc_get_vars_double
    
integer (C_INT) function nc_get_varm_schar(ncid,varid,start,ncount,stride,imap,bp) bind(C, name='nc_get_varm_schar')
  use, intrinsic :: ISO_C_BINDING
  implicit none
  integer (C_INT), value :: ncid, varid
  integer (C_SIZE_T), dimension(*) :: start, ncount
  integer (C_INTPTR_T), dimension(*) :: stride, imap
  type (C_PTR), value :: bp
end function nc_get_varm_schar

integer (C_INT) function nc_get_varm_short(ncid,varid,start,ncount,stride,imap,sp) bind(C, name='nc_get_varm_short')
  use, intrinsic :: ISO_C_BINDING
  implicit none
  integer (C_INT), value :: ncid, varid
  integer (C_SIZE_T), dimension(*) :: start, ncount
  integer (C_INTPTR_T), dimension(*) :: stride, imap
  type (C_PTR), value :: sp
end function nc_get_varm_short   

integer (C_INT) function nc_get_varm_int(ncid,varid,start,ncount,stride,imap,ip) bind(C, name='nc_get_varm_int')
  use, intrinsic :: ISO_C_BINDING
  implicit none
  integer (C_INT), value :: ncid, varid
  integer (C_SIZE_T), dimension(*) :: start, ncount
  integer (C_INTPTR_T), dimension(*) :: stride, imap
  type (C_PTR), value :: ip
end function nc_get_varm_int

integer (C_INT) function nc_get_varm_float(ncid,varid,start,ncount,stride,imap,fp) bind(C, name='nc_get_varm_float')
  use, intrinsic :: ISO_C_BINDING
  implicit none
  integer (C_INT), value :: ncid, varid
  integer (C_SIZE_T), dimension(*) :: start, ncount
  integer (C_INTPTR_T), dimension(*) :: stride, imap
  type (C_PTR), value :: fp
end function nc_get_varm_float   

integer (C_INT) function nc_get_varm_double(ncid,varid,start,ncount,stride,imap,dp) bind(C, name='nc_get_varm_double')
  use, intrinsic :: ISO_C_BINDING
  implicit none
  integer (C_INT), value :: ncid, varid
  integer (C_SIZE_T), dimension(*) :: start, ncount
  integer (C_INTPTR_T), dimension(*) :: stride, imap
  type (C_PTR), value :: dp
end function nc_get_varm_double
    
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

integer (C_INT) function nc_rename_att(ncid,varid,name,newname) bind(C, name='nc_rename_att')
  use, intrinsic :: ISO_C_BINDING
  implicit none
  integer (C_INT), value :: ncid, varid
  character, dimension(*) :: name
  character, dimension(*) :: newname
end function nc_rename_att  

integer (C_INT) function nc_rename_var(ncid,varid,name) bind(C, name='nc_rename_var')
  use, intrinsic :: ISO_C_BINDING
  implicit none
  integer (C_INT), value :: ncid, varid
  character, dimension(*) :: name
end function nc_rename_var 
    
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

integer (C_INT) function nc_put_att_schar(ncid,varid,name,xtype,size,bp) bind(C, name='nc_put_att_schar')
  use, intrinsic :: ISO_C_BINDING
  implicit none
  integer (C_INT), value :: ncid, varid, xtype
  integer (C_SIZE_T), value :: size
  character, dimension(*) :: name
  type (C_PTR), value :: bp
end function nc_put_att_schar

integer (C_INT) function nc_put_att_double(ncid,varid,name,xtype,size,dp) bind(C, name='nc_put_att_double')
  use, intrinsic :: ISO_C_BINDING
  implicit none
  integer (C_INT), value :: ncid, varid, xtype
  integer (C_SIZE_T), value :: size
  type (C_PTR), value :: dp
  character, dimension(*) :: name
end function nc_put_att_double
    
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

integer (C_INT) function nc_put_var1_schar(ncid,varid,start,bp) bind(C, name='nc_put_var1_schar')
  use, intrinsic :: ISO_C_BINDING
  implicit none
  integer (C_INT), value :: ncid, varid
  integer (C_SIZE_T), dimension(*) :: start
  type (C_PTR), value :: bp
end function nc_put_var1_schar

integer (C_INT) function nc_put_var1_short(ncid,varid,start,sp) bind(C, name='nc_put_var1_short')
  use, intrinsic :: ISO_C_BINDING
  implicit none
  integer (C_INT), value :: ncid, varid
  integer (C_SIZE_T), dimension(*) :: start
  type (C_PTR), value :: sp
end function nc_put_var1_short
    
integer (C_INT) function nc_put_vars_text(ncid,varid,start,ncount,stride,tp) bind(C, name='nc_put_vars_text')
  use, intrinsic :: ISO_C_BINDING
  implicit none
  integer (C_INT), value :: ncid, varid
  integer (C_SIZE_T), dimension(*) :: start, ncount
  integer (C_INTPTR_T), dimension(*) :: stride
  character, dimension(*) :: tp
end function nc_put_vars_text

integer (C_INT) function nc_put_vars_schar(ncid,varid,start,count,stride,bp) bind(C, name='nc_put_vars_schar')
  use, intrinsic :: ISO_C_BINDING
  implicit none
  integer (C_INT), value :: ncid, varid
  integer (C_SIZE_T), dimension(*) :: start, count
  integer (C_INTPTR_T), dimension(*) :: stride
  type (C_PTR), value :: bp
end function nc_put_vars_schar

integer (C_INT) function nc_put_vars_short(ncid,varid,start,count,stride,sp) bind(C, name='nc_put_vars_short')
  use, intrinsic :: ISO_C_BINDING
  implicit none
  integer (C_INT), value :: ncid, varid
  integer (C_SIZE_T), dimension(*) :: start, count
  integer (C_INTPTR_T), dimension(*) :: stride
  type (C_PTR), value :: sp
end function nc_put_vars_short

integer (C_INT) function nc_put_vars_int(ncid,varid,start,count,stride,ip) bind(C, name='nc_put_vars_int')
  use, intrinsic :: ISO_C_BINDING
  implicit none
  integer (C_INT), value :: ncid, varid
  integer (C_SIZE_T), dimension(*) :: start, count
  integer (C_INTPTR_T), dimension(*) :: stride
  type (C_PTR), value :: ip
end function nc_put_vars_int
   
integer (C_INT) function nc_put_vars_float(ncid,varid,start,count,stride,fp) bind(C, name='nc_put_vars_float')
  use, intrinsic :: ISO_C_BINDING
  implicit none
  integer (C_INT), value :: ncid, varid
  integer (C_SIZE_T), dimension(*) :: start, count
  integer (C_INTPTR_T), dimension(*) :: stride
  type (C_PTR), value :: fp
end function nc_put_vars_float

integer (C_INT) function nc_put_vars_double(ncid,varid,start,count,stride,dp) bind(C, name='nc_put_vars_double')
  use, intrinsic :: ISO_C_BINDING
  implicit none
  integer (C_INT), value :: ncid, varid
  integer (C_SIZE_T), dimension(*) :: start, count
  integer (C_INTPTR_T), dimension(*) :: stride
  type (C_PTR), value :: dp
end function nc_put_vars_double
    
integer (C_INT) function nc_put_varm_schar(ncid,varid,start,count,stride,imap,bp) bind(C, name='nc_put_varm_schar')
  use, intrinsic :: ISO_C_BINDING
  implicit none
  integer (C_INT), value :: ncid, varid
  integer (C_SIZE_T), dimension(*) :: start, count
  integer (C_INTPTR_T), dimension(*) :: stride, imap
  type (C_PTR), value :: bp
end function nc_put_varm_schar

integer (C_INT) function nc_put_varm_short(ncid,varid,start,count,stride,imap,sp) bind(C, name='nc_put_varm_short')
  use, intrinsic :: ISO_C_BINDING
  implicit none
  integer (C_INT), value :: ncid, varid
  integer (C_SIZE_T), dimension(*) :: start, count
  integer (C_INTPTR_T), dimension(*) :: stride, imap
  type (C_PTR), value :: sp
end function nc_put_varm_short

integer (C_INT) function nc_put_varm_int(ncid,varid,start,count,stride,imap,ip) bind(C, name='nc_put_varm_int')
  use, intrinsic :: ISO_C_BINDING
  implicit none
  integer (C_INT), value :: ncid, varid
  integer (C_SIZE_T), dimension(*) :: start, count
  integer (C_INTPTR_T), dimension(*) :: stride, imap
  type (C_PTR), value :: ip
end function nc_put_varm_int

integer (C_INT) function nc_put_varm_float(ncid,varid,start,count,stride,imap,fp) bind(C, name='nc_put_varm_float')
  use, intrinsic :: ISO_C_BINDING
  implicit none
  integer (C_INT), value :: ncid, varid
  integer (C_SIZE_T), dimension(*) :: start, count
  integer (C_INTPTR_T), dimension(*) :: stride, imap
  type (C_PTR), value :: fp
end function nc_put_varm_float

integer (C_INT) function nc_put_varm_double(ncid,varid,start,count,stride,imap,dp) bind(C, name='nc_put_varm_double')
  use, intrinsic :: ISO_C_BINDING
  implicit none
  integer (C_INT), value :: ncid, varid
  integer (C_SIZE_T), dimension(*) :: start, count
  integer (C_INTPTR_T), dimension(*) :: stride, imap
  type (C_PTR), value :: dp
end function nc_put_varm_double
    
integer (C_INT) function nc_copy_att(ncidin,varidin,name,ncidout,varidout) bind(C, name='nc_copy_att')
  use, intrinsic :: ISO_C_BINDING
  implicit none
  integer (C_INT), value :: ncidin, varidin, ncidout, varidout
  character, dimension(*) :: name
end function nc_copy_att

integer (C_INT) function nc_del_att(ncid,varid,name) bind(C, name='nc_del_att')
  use, intrinsic :: ISO_C_BINDING
  implicit none
  integer (C_INT), value :: ncid, varid
  character, dimension(*) :: name
end function nc_del_att
    
end interface

interface nf_get_att_real
  module procedure nf_get_att_real_s, nf_get_att_real_v
end interface nf_get_att_real

interface nf_get_att_int
  module procedure nf_get_att_int_s, nf_get_att_int_v
end interface nf_get_att_int

interface nf_get_att_int1
  module procedure nf_get_att_int1_s, nf_get_att_int1_v
end interface nf_get_att_int1

interface nf_get_att_int2
  module procedure nf_get_att_int2_s, nf_get_att_int2_v
end interface nf_get_att_int2

interface nf_get_att_double
  module procedure nf_get_att_double_s, nf_get_att_double_v
end interface nf_get_att_double
    
interface nf_get_vara_real
  module procedure nf_get_vara_real_d1, nf_get_vara_real_d2, nf_get_vara_real_d3, nf_get_vara_real_d4, &
                   nf_get_vara_real_d5
end interface nf_get_vara_real

interface nf_get_vara_int
  module procedure nf_get_vara_int_d1, nf_get_vara_int_d2, nf_get_vara_int_d3, nf_get_vara_int_d4, &
                   nf_get_vara_int_d5
end interface nf_get_vara_int

interface nf_get_vara_int2
  module procedure nf_get_vara_int2_d1, nf_get_vara_int2_d2, nf_get_vara_int2_d3, nf_get_vara_int2_d4, &
                   nf_get_vara_int2_d5
end interface nf_get_vara_int2

interface nf_get_vara_double
  module procedure nf_get_vara_double_d1, nf_get_vara_double_d2, nf_get_vara_double_d3, nf_get_vara_double_d4, &
                   nf_get_vara_double_d5
end interface nf_get_vara_double

interface nf_get_var_real
  module procedure nf_get_var_real_d1, nf_get_var_real_d2, nf_get_var_real_d3, nf_get_var_real_d4, &
                   nf_get_var_real_d5
end interface nf_get_var_real

interface nf_get_var_double
  module procedure nf_get_var_double_d1, nf_get_var_double_d2, nf_get_var_double_d3, nf_get_var_double_d4, &
                   nf_get_var_double_d5
end interface nf_get_var_double    
    
interface nf_get_var1_real
  module procedure nf_get_var1_real_s, nf_get_var1_real_v
end interface nf_get_var1_real

interface nf_get_var1_int
  module procedure nf_get_var1_int_s, nf_get_var1_int_v
end interface nf_get_var1_int

interface nf_get_var1_int1
  module procedure nf_get_var1_int1_s, nf_get_var1_int1_v
end interface nf_get_var1_int1
    
interface nf_get_var1_int2
  module procedure nf_get_var1_int2_s, nf_get_var1_int2_v
end interface nf_get_var1_int2

interface nf_get_var1_double
  module procedure nf_get_var1_double_s, nf_get_var1_double_v
end interface nf_get_var1_double

interface nf_get_vars_int1
  module procedure nf_get_vars_int1_d1, nf_get_vars_int1_d2, nf_get_vars_int1_d3, nf_get_vars_int1_d4, &
                   nf_get_vars_int1_d5, nf_get_vars_int1_d6, nf_get_vars_int1_d7
end interface nf_get_vars_int1

interface nf_get_vars_int2
  module procedure nf_get_vars_int2_d1, nf_get_vars_int2_d2, nf_get_vars_int2_d3, nf_get_vars_int2_d4, &
                   nf_get_vars_int2_d5, nf_get_vars_int2_d6, nf_get_vars_int2_d7
end interface nf_get_vars_int2

interface nf_get_vars_int
  module procedure nf_get_vars_int_d1, nf_get_vars_int_d2, nf_get_vars_int_d3, nf_get_vars_int_d4, &
                   nf_get_vars_int_d5, nf_get_vars_int_d6, nf_get_vars_int_d7
end interface nf_get_vars_int    

interface nf_get_vars_real
  module procedure nf_get_vars_real_d1, nf_get_vars_real_d2, nf_get_vars_real_d3, nf_get_vars_real_d4, &
                   nf_get_vars_real_d5, nf_get_vars_real_d6, nf_get_vars_real_d7
end interface nf_get_vars_real

interface nf_get_vars_double
  module procedure nf_get_vars_double_d1, nf_get_vars_double_d2, nf_get_vars_double_d3, nf_get_vars_double_d4, &
                   nf_get_vars_double_d5, nf_get_vars_double_d6, nf_get_vars_double_d7
end interface nf_get_vars_double

interface nf_get_varm_int1
  module procedure nf_get_varm_int1_d1, nf_get_varm_int1_d2, nf_get_varm_int1_d3, nf_get_varm_int1_d4, &
                   nf_get_varm_int1_d5, nf_get_varm_int1_d6, nf_get_varm_int1_d7
end interface nf_get_varm_int1

interface nf_get_varm_int2
  module procedure nf_get_varm_int2_d1, nf_get_varm_int2_d2, nf_get_varm_int2_d3, nf_get_varm_int2_d4, &
                   nf_get_varm_int2_d5, nf_get_varm_int2_d6, nf_get_varm_int2_d7
end interface nf_get_varm_int2

interface nf_get_varm_int
  module procedure nf_get_varm_int_d1, nf_get_varm_int_d2, nf_get_varm_int_d3, nf_get_varm_int_d4, &
                   nf_get_varm_int_d5, nf_get_varm_int_d6, nf_get_varm_int_d7
end interface nf_get_varm_int
    
interface nf_get_varm_real
  module procedure nf_get_varm_real_d1, nf_get_varm_real_d2, nf_get_varm_real_d3, nf_get_varm_real_d4, &
                   nf_get_varm_real_d5, nf_get_varm_real_d6, nf_get_varm_real_d7
end interface nf_get_varm_real

interface nf_get_varm_double
  module procedure nf_get_varm_double_d1, nf_get_varm_double_d2, nf_get_varm_double_d3, nf_get_varm_double_d4, &
                   nf_get_varm_double_d5, nf_get_varm_double_d6, nf_get_varm_double_d7
end interface nf_get_varm_double
    
interface nf_def_var
  module procedure nf_def_var_s, nf_def_var_v
end interface nf_def_var

interface nf_put_att_int2
  module procedure nf_put_att_int2_s, nf_put_att_int2_v
end interface nf_put_att_int2

interface nf_put_att_real
  module procedure nf_put_att_real_s, nf_put_att_real_v
end interface nf_put_att_real

interface nf_put_att_int
  module procedure nf_put_att_int_s, nf_put_att_int_v
end interface nf_put_att_int
    
interface nf_put_att_int1
  module procedure nf_put_att_int1_s, nf_put_att_int1_v
end interface nf_put_att_int1

interface nf_put_att_double
  module procedure nf_put_att_double_s, nf_put_att_double_v
end interface nf_put_att_double
    
interface nf_put_vara_real
  module procedure nf_put_vara_real_d1, nf_put_vara_real_d2, nf_put_vara_real_d3, nf_put_vara_real_d4, &
                   nf_put_vara_real_d5
end interface nf_put_vara_real

interface nf_put_vara_int
  module procedure nf_put_vara_int_d1, nf_put_vara_int_d2, nf_put_vara_int_d3, nf_put_vara_int_d4, &
                   nf_put_vara_int_d5
end interface nf_put_vara_int

interface nf_put_vara_int2
  module procedure nf_put_vara_int2_d1, nf_put_vara_int2_d2, nf_put_vara_int2_d3, nf_put_vara_int2_d4, &
                   nf_put_vara_int2_d5
end interface nf_put_vara_int2
    
interface nf_put_vara_double
  module procedure nf_put_vara_double_d1, nf_put_vara_double_d2, nf_put_vara_double_d3, nf_put_vara_double_d4, &
                   nf_put_vara_double_d5
end interface nf_put_vara_double

interface nf_put_var_real
  module procedure nf_put_var_real_d1, nf_put_var_real_d2, nf_put_var_real_d3, nf_put_var_real_d4, &
                   nf_put_var_real_d5
end interface nf_put_var_real

interface nf_put_vars_int1
  module procedure nf_put_vars_int1_d1, nf_put_vars_int1_d2, nf_put_vars_int1_d3, nf_put_vars_int1_d4, &
                   nf_put_vars_int1_d5, nf_put_vars_int1_d6, nf_put_vars_int1_d7
end interface nf_put_vars_int1

interface nf_put_vars_int2
  module procedure nf_put_vars_int2_d1, nf_put_vars_int2_d2, nf_put_vars_int2_d3, nf_put_vars_int2_d4, &
                   nf_put_vars_int2_d5, nf_put_vars_int2_d6, nf_put_vars_int2_d7
end interface nf_put_vars_int2

interface nf_put_vars_int
  module procedure nf_put_vars_int_d1, nf_put_vars_int_d2, nf_put_vars_int_d3, nf_put_vars_int_d4, &
                   nf_put_vars_int_d5, nf_put_vars_int_d6, nf_put_vars_int_d7
end interface nf_put_vars_int
    
interface nf_put_vars_real
  module procedure nf_put_vars_real_d1, nf_put_vars_real_d2, nf_put_vars_real_d3, nf_put_vars_real_d4, &
                   nf_put_vars_real_d5, nf_put_vars_real_d6, nf_put_vars_real_d7
end interface nf_put_vars_real

interface nf_put_vars_double
  module procedure nf_put_vars_double_d1, nf_put_vars_double_d2, nf_put_vars_double_d3, nf_put_vars_double_d4, &
                   nf_put_vars_double_d5, nf_put_vars_double_d6, nf_put_vars_double_d7
end interface nf_put_vars_double
    
interface nf_put_varm_int1
  module procedure nf_put_varm_int1_d1, nf_put_varm_int1_d2, nf_put_varm_int1_d3, nf_put_varm_int1_d4, &
                   nf_put_varm_int1_d5, nf_put_varm_int1_d6, nf_put_varm_int1_d7
end interface nf_put_varm_int1    

interface nf_put_varm_int2
  module procedure nf_put_varm_int2_d1, nf_put_varm_int2_d2, nf_put_varm_int2_d3, nf_put_varm_int2_d4, &
                   nf_put_varm_int2_d5, nf_put_varm_int2_d6, nf_put_varm_int2_d7
end interface nf_put_varm_int2

interface nf_put_varm_int
  module procedure nf_put_varm_int_d1, nf_put_varm_int_d2, nf_put_varm_int_d3, nf_put_varm_int_d4, &
                   nf_put_varm_int_d5, nf_put_varm_int_d6, nf_put_varm_int_d7
end interface nf_put_varm_int

interface nf_put_varm_real
  module procedure nf_put_varm_real_d1, nf_put_varm_real_d2, nf_put_varm_real_d3, nf_put_varm_real_d4, &
                   nf_put_varm_real_d5, nf_put_varm_real_d6, nf_put_varm_real_d7
end interface nf_put_varm_real

interface nf_put_varm_double
  module procedure nf_put_varm_double_d1, nf_put_varm_double_d2, nf_put_varm_double_d3, nf_put_varm_double_d4, &
                   nf_put_varm_double_d5, nf_put_varm_double_d6, nf_put_varm_double_d7
end interface nf_put_varm_double
    
integer, parameter :: nf_unlimited = 0

integer, parameter :: nf_noerr = 0
integer, parameter :: nf_ebadid = -33
integer, parameter :: nf_eexist = -35
integer, parameter :: nf_einval = -36
integer, parameter :: nf_eparm = -37
integer, parameter :: nf_enotindefine = -38
integer, parameter :: nf_eindefine = -39
integer, parameter :: nf_einvalcoords = -40
integer, parameter :: nf_emaxdims = -41
integer, parameter :: nf_enameinuse = -42
integer, parameter :: nf_enotatt = -43
integer, parameter :: nf_emaxatts = -44
integer, parameter :: nf_ebadtype = -45
integer, parameter :: nf_ebaddim = -46
integer, parameter :: nf_eunlimpos = -47
integer, parameter :: nf_emaxvars = -48
integer, parameter :: nf_enotvar = -49
integer, parameter :: nf_eglobal = -50
integer, parameter :: nf_enotnc = -51
integer, parameter :: nf_ests = - 52
integer, parameter :: nf_emaxname = -53
integer, parameter :: nf_eunlimit = -54
integer, parameter :: nf_enorecvars = -55
integer, parameter :: nf_echar = -56
integer, parameter :: nf_eedge = -57
integer, parameter :: nf_estride = - 58
integer, parameter :: nf_ebadname = -59
integer, parameter :: nf_erange = -60
integer, parameter :: nf_enomem = -61
integer, parameter :: nf_evarsize = -62
integer, parameter :: nf_edimsize = -63

integer, parameter :: nf_nowrite = 0
integer, parameter :: nf_write = 1
integer, parameter :: nf_clobber = 0
integer, parameter :: nf_noclobber = 4
integer, parameter :: nf_fill = 0
integer, parameter :: nf_nofill = 256
integer, parameter :: nf_lock = 1024
integer, parameter :: nf_share = 2048
integer, parameter :: nf_64bit_offset = 512
integer, parameter :: nf_sizehint_default = 0
integer, parameter :: nf_align_chunk = -1
integer, parameter :: nf_format_classic = 1
integer, parameter :: nf_format_64bit = 2

integer, parameter :: nf_global = 0

integer, parameter :: nf_byte = 1
integer, parameter :: nf_int1 = nf_byte
integer, parameter :: nf_char = 2
integer, parameter :: nf_short = 3
integer, parameter :: nf_int2 = nf_short
integer, parameter :: nf_int = 4
integer, parameter :: nf_float = 5
integer, parameter :: nf_real = nf_float
integer, parameter :: nf_double = 6

integer, parameter :: nf_fill_byte = -127
integer, parameter :: nf_fill_int1 = nf_fill_byte
integer, parameter :: nf_fill_char = 0
integer, parameter :: nf_fill_short = -32767
integer, parameter :: nf_fill_int2 = nf_fill_short
integer, parameter :: nf_fill_int = -2147483647
real, parameter :: nf_fill_float = 9.9692099683868690e+36
real, parameter :: nf_fill_real = nf_fill_float
real(kind=8), parameter :: nf_fill_double = 9.9692099683868690e+36

integer, parameter :: nf_max_dims = 512
integer, parameter :: nf_max_attrs = 4096
integer, parameter :: nf_max_vars = 4096
integer, parameter :: nf_max_name = 128
integer, parameter :: nf_max_var_dims = nf_max_dims

integer, parameter :: nf_fatal = 1
integer, parameter :: nf_verbose = 2

integer, parameter :: charsize = 80

contains
    
integer function nf_open(name,mode,ncid) result(ierr)
  implicit none
  integer, intent(in) :: mode
  integer, intent(out) :: ncid
  character(len=*), intent(in) :: name
  integer (C_INT) :: c_mode
  integer (C_INT), target :: c_ncid
  character, dimension(len(name)) :: c_name
  c_mode = mode
  call cf_strcopy(name,c_name)
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
  character, dimension(len(name)) :: c_name
  c_mode = mode
  call cf_strcopy(name,c_name)
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
  character, dimension(len(name)) :: c_name
  c_mode = mode
  call cf_strcopy(name,c_name)
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
  character, dimension(len(name)) :: c_name
  c_mode = mode
  call cf_strcopy(name,c_name)
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
  c_varid = varid - 1
  ierr = nc_inq_varndims(c_ncid,c_varid,C_LOC(c_ndims))
  ndims = c_ndims
end function nf_inq_varndims

integer function nf_inq_vardimid(ncid,varid,dimids) result(ierr)
  implicit none
  integer, intent(in) :: ncid, varid
  integer, dimension(:), intent(out) :: dimids
  integer (C_INT) :: c_ncid, c_varid
  integer (C_INT), target :: c_ndims
  integer (C_INT), dimension(size(dimids)), target :: c_dimids
  integer i, ndims
  c_ncid = ncid
  c_varid = varid - 1
  c_ndims = size(dimids)
  ierr = nc_inq_varndims(c_ncid,c_varid,C_LOC(c_ndims))
  ndims = c_ndims
  ierr = nc_inq_vardimid(c_ncid,c_varid,C_LOC(c_dimids))
  do i = 1,ndims
    dimids(ndims-i+1) = c_dimids(i) + 1
  end do
end function nf_inq_vardimid

integer function nf_inq_dimlen(ncid,dimid,length) result(ierr)
  implicit none
  integer, intent(in) :: ncid, dimid
  integer, intent(out) :: length
  integer (C_INT) :: c_ncid, c_dimid
  integer (C_SIZE_T), target :: c_length
  c_ncid = ncid
  c_dimid = dimid - 1
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
  character, dimension(len(name)) :: c_name
  c_ncid = ncid
  call cf_strcopy(name,c_name)
  ierr = nc_inq_varid(c_ncid,c_name,C_LOC(c_varid))
  varid = c_varid + 1
end function nf_inq_varid

integer function nf_inq_dimid(ncid,name,dimid) result(ierr)
  implicit none
  integer, intent(in) :: ncid
  integer, intent(out) :: dimid
  character(len=*), intent(in) :: name
  integer (C_INT) :: c_ncid
  integer (C_INT), target :: c_dimid
  character, dimension(len(name)) :: c_name
  c_ncid = ncid
  call cf_strcopy(name,c_name)
  ierr = nc_inq_dimid(c_ncid,c_name,C_LOC(c_dimid))
  dimid = c_dimid + 1
end function nf_inq_dimid

integer function nf_inq_vartype(ncid,varid,xtype) result(ierr)
  implicit none
  integer, intent(in) :: ncid, varid
  integer, intent(out) :: xtype
  integer (C_INT) :: c_ncid, c_varid
  integer (C_INT), target :: c_xtype
  c_ncid = ncid
  c_varid = varid - 1
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
  integer (C_INT) :: c_ncid, c_varid
  character, dimension(len(tp)) :: c_tp
  c_ncid = ncid
  c_varid = varid - 1
  c_tp(:) = ''  
  ierr = nc_inq_varname(c_ncid,c_varid,c_tp)
  call fc_strcopy(c_tp,tp)
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
  ix = index(temp,C_NULL_CHAR) - 1
  msg = ''  
  msg = temp(1:ix)
end function nf_inq_libvers

integer function nf_inq_dim(ncid,dimid,name,length) result(ierr)
  implicit none
  integer, intent(in) :: ncid, dimid
  integer, intent(out) :: length
  character(len=*), intent(out) :: name
    integer (C_INT) :: c_ncid, c_dimid  
  integer (C_SIZE_T), target :: c_length
  character, dimension(len(name)) :: c_name
  c_ncid = ncid
  c_dimid = dimid - 1
  ierr = nc_inq_dim(c_ncid,c_dimid,c_name,C_LOC(c_length))
  call fc_strcopy(c_name,name)
  length = c_length
end function nf_inq_dim

integer function nf_inq_attname(ncid,varid,attnum,tp) result(ierr)
  implicit none
  integer, intent(in) :: ncid, varid, attnum
  character(len=*), intent(out) :: tp
  integer (C_INT) :: c_ncid, c_varid, c_attnum
  character, dimension(len(tp)) :: c_tp
  integer i, ix
  c_ncid = ncid
  c_varid = varid - 1
  c_attnum = attnum
  c_tp(:) = ''  
  ierr = nc_inq_attname(c_ncid,c_varid,c_attnum,c_tp)
  call fc_strcopy(c_tp,tp)
end function nf_inq_attname

integer function nf_inq_attid(ncid,varid,name,attnum) result(ierr)
  implicit none
  integer, intent(in) :: ncid, varid
  integer, intent(out) :: attnum
  character(len=*), intent(in) :: name
  integer (C_INT) :: c_ncid, c_varid
  integer (C_INT), target :: c_attnum
  character, dimension(len(name)) :: c_name
  c_ncid = ncid
  c_varid = varid - 1
  call cf_strcopy(name,c_name)
  ierr = nc_inq_attid(c_ncid,c_varid,c_name,C_LOC(c_attnum))
  attnum = c_attnum
end function nf_inq_attid

integer function nf_inq_att(ncid,varid,name,xtypep,lenp) result(ierr)
  implicit none
  integer, intent(in) :: ncid, varid
  integer, intent(out) :: xtypep, lenp
  character(len=*), intent(in) :: name
  integer (C_INT) :: c_ncid, c_varid
  integer (C_INT), target :: c_xtypep
  integer (C_SIZE_T), target :: c_lenp
  character, dimension(len(name)) :: c_name
  c_ncid = ncid
  c_varid = varid - 1
  call cf_strcopy(name,c_name)
  ierr = nc_inq_att(c_ncid,c_varid,c_name,C_LOC(c_xtypep),C_LOC(c_lenp))
  xtypep = c_xtypep
  lenp = c_lenp
end function nf_inq_att

integer function nf_inq_var(ncid,varid,name,xtypep,ndimsp,dimidsp,nattsp) result(ierr)
  implicit none
  integer, intent(in) :: ncid, varid
  integer, intent(out) :: xtypep, ndimsp, nattsp
  integer, dimension(:), intent(out) :: dimidsp
  character(len=*), intent(out) :: name
  integer (C_INT) :: c_ncid, c_varid
  integer (C_INT), target :: c_xtypep, c_ndimsp, c_nattsp
  integer (C_INT), dimension(size(dimidsp)), target :: c_dimidsp
  character, dimension(len(name)) :: c_name
  integer i
  c_ncid = ncid
  c_varid = varid - 1
  c_name = ''
  ierr = nc_inq_var(c_ncid,c_varid,c_name,C_LOC(c_xtypep),C_LOC(c_ndimsp),C_LOC(c_dimidsp),C_LOC(c_nattsp))
  call fc_strcopy(c_name,name)
  xtypep = c_xtypep
  ndimsp = c_ndimsp
  do i = 1,ndimsp
    dimidsp(ndimsp-i+1) = c_dimidsp(i) + 1
  end do
  nattsp = c_nattsp
end function nf_inq_var

integer function nf_inq_attlen(ncid,varid,name,lenp) result(ierr)
  implicit none
  integer, intent(in) :: ncid, varid
  integer, intent(out) :: lenp
  character(len=*), intent(in) :: name
  integer (C_INT) :: c_ncid, c_varid
  integer (C_SIZE_T), target :: c_lenp
  character, dimension(len(name)) :: c_name
  c_ncid = ncid
  c_varid = varid - 1
  call cf_strcopy(name,c_name)
  ierr = nc_inq_attlen(c_ncid,c_varid,c_name,C_LOC(c_lenp))
  lenp = c_lenp
end function nf_inq_attlen

integer function nf_get_att_text(ncid,varid,name,tp) result(ierr)
  implicit none
  integer, intent(in) :: ncid, varid
  character(len=*), intent(in) :: name
  character(len=*), intent(out) :: tp
  integer (C_INT) :: c_ncid, c_varid
  character, dimension(len(name)) :: c_name
  character, dimension(len(tp)), target :: c_tp
  c_ncid = ncid
  c_varid = varid - 1
  call cf_strcopy(name,c_name)
  c_tp(:) = ''
  ierr = nc_get_att_text(c_ncid,c_varid,c_name,C_LOC(c_tp))
  call fc_strcopy(c_tp,tp)
end function nf_get_att_text

integer function nf_get_att_real_s(ncid,varid,name,rp) result(ierr)
  implicit none
  integer, intent(in) :: ncid, varid
  real(kind=4), intent(out) :: rp
  character(len=*), intent(in) :: name
  real (C_FLOAT), target :: c_rp
  integer (C_INT) :: c_ncid, c_varid
  character, dimension(len(name)) :: c_name
  c_ncid = ncid
  c_varid = varid - 1
  call cf_strcopy(name,c_name)
  ierr = nc_get_att_float(c_ncid,c_varid,c_name,C_LOC(c_rp))
  rp = c_rp  
end function nf_get_att_real_s

integer function nf_get_att_real_v(ncid,varid,name,rp) result(ierr)
  implicit none
  integer, intent(in) :: ncid, varid
  real(kind=4), dimension(:), intent(out) :: rp
  character(len=*), intent(in) :: name
  real (C_FLOAT), dimension(size(rp)), target :: c_rp
  integer (C_INT) :: c_ncid, c_varid
  character, dimension(len(name)) :: c_name
  c_ncid = ncid
  c_varid = varid - 1
  call cf_strcopy(name,c_name)
  ierr = nc_get_att_float(c_ncid,c_varid,c_name,C_LOC(c_rp))
  rp = c_rp  
end function nf_get_att_real_v

integer function nf_get_att_int_s(ncid,varid,name,ip) result(ierr)
  implicit none
  integer, intent(in) :: ncid, varid
  integer(kind=4), intent(out) :: ip
  character(len=*), intent(in) :: name
  integer (C_INT), target :: c_ip
  integer (C_INT) :: c_ncid, c_varid
  character, dimension(len(name)) :: c_name
  c_ncid = ncid
  c_varid = varid - 1
  call cf_strcopy(name,c_name)
  ierr = nc_get_att_int(c_ncid,c_varid,c_name,C_LOC(c_ip))
  ip = c_ip  
end function nf_get_att_int_s

integer function nf_get_att_int_v(ncid,varid,name,ip) result(ierr)
  implicit none
  integer, intent(in) :: ncid, varid
  integer(kind=4), dimension(:), intent(out) :: ip
  character(len=*), intent(in) :: name
  integer (C_INT), dimension(size(ip)), target :: c_ip
  integer (C_INT) :: c_ncid, c_varid
  character, dimension(len(name)) :: c_name
  c_ncid = ncid
  c_varid = varid - 1
  call cf_strcopy(name,c_name)
  ierr = nc_get_att_int(c_ncid,c_varid,c_name,C_LOC(c_ip))
  ip = c_ip  
end function nf_get_att_int_v

integer function nf_get_att_int1_s(ncid,varid,name,bp) result(ierr)
  implicit none
  integer, intent(in) :: ncid, varid
  integer(kind=1), intent(out) :: bp
  character(len=*), intent(in) :: name
  integer (C_SIGNED_CHAR), target :: c_bp
  integer (C_INT) :: c_ncid, c_varid
  character, dimension(len(name)) :: c_name
  c_ncid = ncid
  c_varid = varid - 1
  call cf_strcopy(name,c_name)
  ierr = nc_get_att_schar(c_ncid,c_varid,c_name,C_LOC(c_bp))
  bp = c_bp  
end function nf_get_att_int1_s

integer function nf_get_att_int1_v(ncid,varid,name,bp) result(ierr)
  implicit none
  integer, intent(in) :: ncid, varid
  integer(kind=1), dimension(:), intent(out) :: bp
  character(len=*), intent(in) :: name
  integer (C_SIGNED_CHAR), dimension(size(bp)), target :: c_bp
  integer (C_INT) :: c_ncid, c_varid
  character, dimension(len(name)) :: c_name
  c_ncid = ncid
  c_varid = varid - 1
  call cf_strcopy(name,c_name)
  ierr = nc_get_att_schar(c_ncid,c_varid,c_name,C_LOC(c_bp))
  bp = c_bp  
end function nf_get_att_int1_v

integer function nf_get_att_int2_s(ncid,varid,name,sp) result(ierr)
  implicit none
  integer, intent(in) :: ncid, varid
  integer(kind=2), intent(out) :: sp
  character(len=*), intent(in) :: name
  integer (C_SHORT), target :: c_sp
  integer (C_INT) :: c_ncid, c_varid
  character, dimension(len(name)) :: c_name
  c_ncid = ncid
  c_varid = varid - 1
  call cf_strcopy(name,c_name)
  ierr = nc_get_att_short(c_ncid,c_varid,c_name,C_LOC(c_sp))
  sp = c_sp  
end function nf_get_att_int2_s

integer function nf_get_att_int2_v(ncid,varid,name,sp) result(ierr)
  implicit none
  integer, intent(in) :: ncid, varid
  integer(kind=2), dimension(:), intent(out) :: sp
  character(len=*), intent(in) :: name
  integer (C_SHORT), dimension(size(sp)), target :: c_sp
  integer (C_INT) :: c_ncid, c_varid
  character, dimension(len(name)) :: c_name
  c_ncid = ncid
  c_varid = varid - 1
  call cf_strcopy(name,c_name)
  ierr = nc_get_att_short(c_ncid,c_varid,c_name,C_LOC(c_sp))
  sp = c_sp  
end function nf_get_att_int2_v

integer function nf_get_att_double_s(ncid,varid,name,dp) result(ierr)
  implicit none
  integer, intent(in) :: ncid, varid
  real(kind=8), intent(out) :: dp
  character(len=*), intent(in) :: name
  real (C_DOUBLE), target :: c_dp
  integer (C_INT) :: c_ncid, c_varid
  character, dimension(len(name)) :: c_name
  c_ncid = ncid
  c_varid = varid - 1
  call cf_strcopy(name,c_name)
  ierr = nc_get_att_double(c_ncid,c_varid,c_name,C_LOC(c_dp))
  dp = c_dp  
end function nf_get_att_double_s

integer function nf_get_att_double_v(ncid,varid,name,dp) result(ierr)
  implicit none
  integer, intent(in) :: ncid, varid
  real(kind=8), dimension(:), intent(out) :: dp
  character(len=*), intent(in) :: name
  real (C_DOUBLE), dimension(size(dp)), target :: c_dp
  integer (C_INT) :: c_ncid, c_varid
  character, dimension(len(name)) :: c_name
  c_ncid = ncid
  c_varid = varid - 1
  call cf_strcopy(name,c_name)
  ierr = nc_get_att_double(c_ncid,c_varid,c_name,C_LOC(c_dp))
  dp = c_dp  
end function nf_get_att_double_v

integer function nf_get_vara_real_d1(ncid,varid,start,ncount,fp) result(ierr)
  implicit none
  integer, intent(in) :: ncid, varid
  integer, dimension(:), intent(in) :: start, ncount
  real(kind=4), dimension(:), intent(out) :: fp
  integer (C_INT) :: c_ncid, c_varid
  integer (C_INT), target :: c_ndims
  integer (C_SIZE_T), dimension(size(start)) :: c_start, c_ncount
  real (C_FLOAT), dimension(size(fp)), target :: c_fp
  integer i, ndims
  c_ncid = ncid
  c_varid = varid - 1
  c_ndims = size(start)
  ierr = nc_inq_varndims(c_ncid,c_varid,C_LOC(c_ndims))
  ndims = c_ndims  
  do i = 1,ndims
    c_start(ndims-i+1) = start(i) - 1
    c_ncount(ndims-i+1) = ncount(i)
  end do
  ierr = nc_get_vara_float(c_ncid,c_varid,c_start,c_ncount,C_LOC(c_fp))
  fp = c_fp
end function nf_get_vara_real_d1

integer function nf_get_vara_real_d2(ncid,varid,start,ncount,fp) result(ierr)
  implicit none
  integer, intent(in) :: ncid, varid
  integer, dimension(:), intent(in) :: start, ncount
  real(kind=4), dimension(:,:), intent(out) :: fp
  integer (C_INT) :: c_ncid, c_varid
  integer (C_INT), target :: c_ndims
  integer (C_SIZE_T), dimension(size(start)) :: c_start, c_ncount
  real (C_FLOAT), dimension(size(fp,1),size(fp,2)), target :: c_fp
  integer i, ndims
  c_ncid = ncid
  c_varid = varid - 1
  c_ndims = size(start)
  ierr = nc_inq_varndims(c_ncid,c_varid,C_LOC(c_ndims))
  ndims = c_ndims  
  do i = 1,ndims
    c_start(ndims-i+1) = start(i) - 1
    c_ncount(ndims-i+1) = ncount(i)
  end do
  ierr = nc_get_vara_float(c_ncid,c_varid,c_start,c_ncount,C_LOC(c_fp))
  fp = c_fp
end function nf_get_vara_real_d2

integer function nf_get_vara_real_d3(ncid,varid,start,ncount,fp) result(ierr)
  implicit none
  integer, intent(in) :: ncid, varid
  integer, dimension(:), intent(in) :: start, ncount
  real(kind=4), dimension(:,:,:), intent(out) :: fp
  integer (C_INT) :: c_ncid, c_varid
  integer (C_INT), target :: c_ndims
  integer (C_SIZE_T), dimension(size(start)) :: c_start, c_ncount
  real (C_FLOAT), dimension(size(fp,1),size(fp,2),size(fp,3)), target :: c_fp
  integer i, ndims
  c_ncid = ncid
  c_varid = varid - 1
  c_ndims = size(start)
  ierr = nc_inq_varndims(c_ncid,c_varid,C_LOC(c_ndims))
  ndims = c_ndims
  do i = 1,ndims
    c_start(ndims-i+1) = start(i) - 1
    c_ncount(ndims-i+1) = ncount(i)
  end do
  ierr = nc_get_vara_float(c_ncid,c_varid,c_start,c_ncount,C_LOC(c_fp))
  fp = c_fp
end function nf_get_vara_real_d3

integer function nf_get_vara_real_d4(ncid,varid,start,ncount,fp) result(ierr)
  implicit none
  integer, intent(in) :: ncid, varid
  integer, dimension(:), intent(in) :: start, ncount
  real(kind=4), dimension(:,:,:,:), intent(out) :: fp
  integer (C_INT) :: c_ncid, c_varid
  integer (C_INT), target :: c_ndims
  integer (C_SIZE_T), dimension(size(start)) :: c_start, c_ncount
  real (C_FLOAT), dimension(size(fp,1),size(fp,2),size(fp,3),size(fp,4)), target :: c_fp
  integer i, ndims
  c_ncid = ncid
  c_varid = varid - 1
  c_ndims = size(start)
  ierr = nc_inq_varndims(c_ncid,c_varid,C_LOC(c_ndims))
  ndims = c_ndims
  do i = 1,ndims
    c_start(ndims-i+1) = start(i) - 1
    c_ncount(ndims-i+1) = ncount(i)
  end do
  ierr = nc_get_vara_float(c_ncid,c_varid,c_start,c_ncount,C_LOC(c_fp))
  fp = c_fp
end function nf_get_vara_real_d4

integer function nf_get_vara_real_d5(ncid,varid,start,ncount,fp) result(ierr)
  implicit none
  integer, intent(in) :: ncid, varid
  integer, dimension(:), intent(in) :: start, ncount
  real(kind=4), dimension(:,:,:,:,:), intent(out) :: fp
  integer (C_INT) :: c_ncid, c_varid
  integer (C_INT), target :: c_ndims
  integer (C_SIZE_T), dimension(size(start)) :: c_start, c_ncount
  real (C_FLOAT), dimension(size(fp,1),size(fp,2),size(fp,3),size(fp,4),size(fp,5)), &
      target :: c_fp
  integer i, ndims
  c_ncid = ncid
  c_varid = varid - 1
  c_ndims = size(start)
  ierr = nc_inq_varndims(c_ncid,c_varid,C_LOC(c_ndims))
  ndims = c_ndims
  do i = 1,ndims
    c_start(ndims-i+1) = start(i) - 1
    c_ncount(ndims-i+1) = ncount(i)
  end do
  ierr = nc_get_vara_float(c_ncid,c_varid,c_start,c_ncount,C_LOC(c_fp))
  fp = c_fp
end function nf_get_vara_real_d5

integer function nf_get_vara_int_d1(ncid,varid,start,ncount,ip) result(ierr)
  implicit none
  integer, intent(in) :: ncid, varid
  integer, dimension(:), intent(in) :: start, ncount
  integer(kind=4), dimension(:), intent(out) :: ip
  integer (C_INT) :: c_ncid, c_varid
  integer (C_INT), target :: c_ndims
  integer (C_SIZE_T), dimension(size(start)) :: c_start, c_ncount
  integer (C_INT), dimension(size(ip)), target :: c_ip
  integer i, ndims
  c_ncid = ncid
  c_varid = varid - 1
  c_ndims = size(start)
  ierr = nc_inq_varndims(c_ncid,c_varid,C_LOC(c_ndims))
  ndims = c_ndims
  do i = 1,ndims
    c_start(ndims-i+1) = start(i) - 1
    c_ncount(ndims-i+1) = ncount(i)
  end do
  ierr = nc_get_vara_int(c_ncid,c_varid,c_start,c_ncount,C_LOC(c_ip))
  ip = c_ip
end function nf_get_vara_int_d1

integer function nf_get_vara_int_d2(ncid,varid,start,ncount,ip) result(ierr)
  implicit none
  integer, intent(in) :: ncid, varid
  integer, dimension(:), intent(in) :: start, ncount
  integer(kind=4), dimension(:,:), intent(out) :: ip
  integer (C_INT) :: c_ncid, c_varid
  integer (C_INT), target :: c_ndims
  integer (C_SIZE_T), dimension(size(start)) :: c_start, c_ncount
  integer (C_INT), dimension(size(ip,1),size(ip,2)), target :: c_ip
  integer i, ndims
  c_ncid = ncid
  c_varid = varid - 1
  c_ndims = size(start)
  ierr = nc_inq_varndims(c_ncid,c_varid,C_LOC(c_ndims))
  ndims = c_ndims
  do i = 1,ndims
    c_start(ndims-i+1) = start(i) - 1
    c_ncount(ndims-i+1) = ncount(i)
  end do
  ierr = nc_get_vara_int(c_ncid,c_varid,c_start,c_ncount,C_LOC(c_ip))
  ip = c_ip
end function nf_get_vara_int_d2

integer function nf_get_vara_int_d3(ncid,varid,start,ncount,ip) result(ierr)
  implicit none
  integer, intent(in) :: ncid, varid
  integer, dimension(:), intent(in) :: start, ncount
  integer(kind=4), dimension(:,:,:), intent(out) :: ip
  integer (C_INT) :: c_ncid, c_varid
  integer (C_INT), target :: c_ndims
  integer (C_SIZE_T), dimension(size(start)) :: c_start, c_ncount
  integer (C_INT), dimension(size(ip,1),size(ip,2),size(ip,3)), target :: c_ip
  integer i, ndims
  c_ncid = ncid
  c_varid = varid - 1
  c_ndims = size(start)
  ierr = nc_inq_varndims(c_ncid,c_varid,C_LOC(c_ndims))
  ndims = c_ndims
  do i = 1,ndims
    c_start(ndims-i+1) = start(i) - 1
    c_ncount(ndims-i+1) = ncount(i)
  end do
  ierr = nc_get_vara_int(c_ncid,c_varid,c_start,c_ncount,C_LOC(c_ip))
  ip = c_ip
end function nf_get_vara_int_d3

integer function nf_get_vara_int_d4(ncid,varid,start,ncount,ip) result(ierr)
  implicit none
  integer, intent(in) :: ncid, varid
  integer, dimension(:), intent(in) :: start, ncount
  integer(kind=4), dimension(:,:,:,:), intent(out) :: ip
  integer (C_INT) :: c_ncid, c_varid
  integer (C_INT), target :: c_ndims
  integer (C_SIZE_T), dimension(size(start)) :: c_start, c_ncount
  integer (C_INT), dimension(size(ip,1),size(ip,2),size(ip,3),size(ip,4)), target :: c_ip
  integer i, ndims
  c_ncid = ncid
  c_varid = varid - 1
  c_ndims = size(start)
  ierr = nc_inq_varndims(c_ncid,c_varid,C_LOC(c_ndims))
  ndims = c_ndims
  do i = 1,ndims
    c_start(ndims-i+1) = start(i) - 1
    c_ncount(ndims-i+1) = ncount(i)
  end do
  ierr = nc_get_vara_int(c_ncid,c_varid,c_start,c_ncount,C_LOC(c_ip))
  ip = c_ip
end function nf_get_vara_int_d4

integer function nf_get_vara_int_d5(ncid,varid,start,ncount,ip) result(ierr)
  implicit none
  integer, intent(in) :: ncid, varid
  integer, dimension(:), intent(in) :: start, ncount
  integer(kind=4), dimension(:,:,:,:,:), intent(out) :: ip
  integer (C_INT) :: c_ncid, c_varid
  integer (C_INT), target :: c_ndims
  integer (C_SIZE_T), dimension(size(start)) :: c_start, c_ncount
  integer (C_INT), dimension(size(ip,1),size(ip,2),size(ip,3),size(ip,4),size(ip,5)), &
      target :: c_ip
  integer i, ndims
  c_ncid = ncid
  c_varid = varid - 1
  c_ndims = size(start)
  ierr = nc_inq_varndims(c_ncid,c_varid,C_LOC(c_ndims))
  ndims = c_ndims
  do i = 1,ndims
    c_start(ndims-i+1) = start(i) - 1
    c_ncount(ndims-i+1) = ncount(i)
  end do
  ierr = nc_get_vara_int(c_ncid,c_varid,c_start,c_ncount,C_LOC(c_ip))
  ip = c_ip
end function nf_get_vara_int_d5

integer function nf_get_vara_int2_d1(ncid,varid,start,ncount,sp) result(ierr)
  implicit none
  integer, intent(in) :: ncid, varid
  integer, dimension(:), intent(in) :: start, ncount
  integer(kind=2), dimension(:), intent(out) :: sp
  integer (C_INT) :: c_ncid, c_varid
  integer (C_INT), target :: c_ndims
  integer (C_SIZE_T), dimension(size(start)) :: c_start, c_ncount
  integer (C_SHORT), dimension(size(sp)), target :: c_sp
  integer i, ndims
  c_ncid = ncid
  c_varid = varid - 1
  c_ndims = size(start)
  ierr = nc_inq_varndims(c_ncid,c_varid,C_LOC(c_ndims))
  ndims = c_ndims
  do i = 1,ndims
    c_start(ndims-i+1) = start(i) - 1
    c_ncount(ndims-i+1) = ncount(i)
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
  integer (C_INT), target :: c_ndims
  integer (C_SIZE_T), dimension(size(start)) :: c_start, c_ncount
  integer (C_SHORT), dimension(size(sp,1),size(sp,2)), target :: c_sp
  integer i, ndims
  c_ncid = ncid
  c_varid = varid - 1
  c_ndims = size(start)
  ierr = nc_inq_varndims(c_ncid,c_varid,C_LOC(c_ndims))
  ndims = c_ndims
  do i = 1,ndims
    c_start(ndims-i+1) = start(i) - 1
    c_ncount(ndims-i+1) = ncount(i)
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
  integer (C_INT), target :: c_ndims
  integer (C_SIZE_T), dimension(size(start)) :: c_start, c_ncount
  integer (C_SHORT), dimension(size(sp,1),size(sp,2),size(sp,3)), target :: c_sp
  integer i, ndims
  c_ncid = ncid
  c_varid = varid - 1
  c_ndims = size(start)
  ierr = nc_inq_varndims(c_ncid,c_varid,C_LOC(c_ndims))
  ndims = c_ndims
  do i = 1,ndims
    c_start(ndims-i+1) = start(i) - 1
    c_ncount(ndims-i+1) = ncount(i)
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
  integer (C_INT), target :: c_ndims
  integer (C_SIZE_T), dimension(size(start)) :: c_start, c_ncount
  integer (C_SHORT), dimension(size(sp,1),size(sp,2),size(sp,3),size(sp,4)), target :: c_sp
  integer i, ndims
  c_ncid = ncid
  c_varid = varid - 1
  c_ndims = size(start)
  ierr = nc_inq_varndims(c_ncid,c_varid,C_LOC(c_ndims))
  ndims = c_ndims
  do i = 1,ndims
    c_start(ndims-i+1) = start(i) - 1
    c_ncount(ndims-i+1) = ncount(i)
  end do
  ierr = nc_get_vara_short(c_ncid,c_varid,c_start,c_ncount,C_LOC(c_sp))
  sp = c_sp
end function nf_get_vara_int2_d4

integer function nf_get_vara_int2_d5(ncid,varid,start,ncount,sp) result(ierr)
  implicit none
  integer, intent(in) :: ncid, varid
  integer, dimension(:), intent(in) :: start, ncount
  integer(kind=2), dimension(:,:,:,:,:), intent(out) :: sp
  integer (C_INT) :: c_ncid, c_varid
  integer (C_INT), target :: c_ndims
  integer (C_SIZE_T), dimension(size(start)) :: c_start, c_ncount
  integer (C_SHORT), dimension(size(sp,1),size(sp,2),size(sp,3),size(sp,4),size(sp,5)), &
      target :: c_sp
  integer i, ndims
  c_ncid = ncid
  c_varid = varid - 1
  c_ndims = size(start)
  ierr = nc_inq_varndims(c_ncid,c_varid,C_LOC(c_ndims))
  ndims = c_ndims
  do i = 1,ndims
    c_start(ndims-i+1) = start(i) - 1
    c_ncount(ndims-i+1) = ncount(i)
  end do
  ierr = nc_get_vara_short(c_ncid,c_varid,c_start,c_ncount,C_LOC(c_sp))
  sp = c_sp
end function nf_get_vara_int2_d5

integer function nf_get_vara_double_d1(ncid,varid,start,ncount,dp) result(ierr)
  implicit none
  integer, intent(in) :: ncid, varid
  integer, dimension(:), intent(in) :: start, ncount
  real(kind=8), dimension(:), intent(out) :: dp
  integer (C_INT) :: c_ncid, c_varid
  integer (C_INT), target :: c_ndims
  integer (C_SIZE_T), dimension(size(start)) :: c_start, c_ncount
  real (C_DOUBLE), dimension(size(dp)), target :: c_dp
  integer i, ndims
  c_ncid = ncid
  c_varid = varid - 1
  c_ndims = size(start)
  ierr = nc_inq_varndims(c_ncid,c_varid,C_LOC(c_ndims))
  ndims = c_ndims
  do i = 1,ndims
    c_start(ndims-i+1) = start(i) - 1
    c_ncount(ndims-i+1) = ncount(i)
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
  integer (C_INT), target :: c_ndims
  integer (C_SIZE_T), dimension(size(start)) :: c_start, c_ncount
  real (C_DOUBLE), dimension(size(dp,1),size(dp,2)), target :: c_dp
  integer i, ndims
  c_ncid = ncid
  c_varid = varid - 1
  c_ndims = size(start)
  ierr = nc_inq_varndims(c_ncid,c_varid,C_LOC(c_ndims))
  ndims = c_ndims
  do i = 1,ndims
    c_start(ndims-i+1) = start(i) - 1
    c_ncount(ndims-i+1) = ncount(i)
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
  integer (C_INT), target :: c_ndims
  integer (C_SIZE_T), dimension(size(start)) :: c_start, c_ncount
  real (C_DOUBLE), dimension(size(dp,1),size(dp,2),size(dp,3)), target :: c_dp
  integer i, ndims
  c_ncid = ncid
  c_varid = varid - 1
  c_ndims = size(start)
  ierr = nc_inq_varndims(c_ncid,c_varid,C_LOC(c_ndims))
  ndims = c_ndims
  do i = 1,ndims
    c_start(ndims-i+1) = start(i) - 1
    c_ncount(ndims-i+1) = ncount(i)
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
  integer (C_INT), target :: c_ndims
  integer (C_SIZE_T), dimension(size(start)) :: c_start, c_ncount
  real (C_DOUBLE), dimension(size(dp,1),size(dp,2),size(dp,3),size(dp,4)), target :: c_dp
  integer i, ndims
  c_ncid = ncid
  c_varid = varid - 1
  c_ndims = size(start)
  ierr = nc_inq_varndims(c_ncid,c_varid,C_LOC(c_ndims))
  ndims = c_ndims
  do i = 1,ndims
    c_start(ndims-i+1) = start(i) - 1
    c_ncount(ndims-i+1) = ncount(i)
  end do
  ierr = nc_get_vara_double(c_ncid,c_varid,c_start,c_ncount,C_LOC(c_dp))
  dp = c_dp
end function nf_get_vara_double_d4

integer function nf_get_vara_double_d5(ncid,varid,start,ncount,dp) result(ierr)
  implicit none
  integer, intent(in) :: ncid, varid
  integer, dimension(:), intent(in) :: start, ncount
  real(kind=8), dimension(:,:,:,:,:), intent(out) :: dp
  integer (C_INT) :: c_ncid, c_varid
  integer (C_INT), target :: c_ndims
  integer (C_SIZE_T), dimension(size(start)) :: c_start, c_ncount
  real (C_DOUBLE), dimension(size(dp,1),size(dp,2),size(dp,3),size(dp,4),size(dp,5)), &
      target :: c_dp
  integer i, ndims
  c_ncid = ncid
  c_varid = varid - 1
  c_ndims = size(start)
  ierr = nc_inq_varndims(c_ncid,c_varid,C_LOC(c_ndims))
  ndims = c_ndims
  do i = 1,ndims
    c_start(ndims-i+1) = start(i) - 1
    c_ncount(ndims-i+1) = ncount(i)
  end do
  ierr = nc_get_vara_double(c_ncid,c_varid,c_start,c_ncount,C_LOC(c_dp))
  dp = c_dp
end function nf_get_vara_double_d5

integer function nf_get_var_real_d1(ncid,varid,fp) result(ierr)
  implicit none
  integer, intent(in) :: ncid, varid
  real(kind=4), dimension(:), intent(out) :: fp
  integer (C_INT) :: c_ncid, c_varid
  real (C_FLOAT), dimension(size(fp)), target :: c_fp
  c_ncid = ncid
  c_varid = varid - 1
  ierr = nc_get_var_float(c_ncid,c_varid,C_LOC(c_fp))
  fp = c_fp
end function nf_get_var_real_d1

integer function nf_get_var_real_d2(ncid,varid,fp) result(ierr)
  implicit none
  integer, intent(in) :: ncid, varid
  real(kind=4), dimension(:,:), intent(out) :: fp
  integer (C_INT) :: c_ncid, c_varid
  real (C_FLOAT), dimension(size(fp,1),size(fp,2)), target :: c_fp
  c_ncid = ncid
  c_varid = varid - 1
  ierr = nc_get_var_float(c_ncid,c_varid,C_LOC(c_fp))
  fp = c_fp
end function nf_get_var_real_d2

integer function nf_get_var_real_d3(ncid,varid,fp) result(ierr)
  implicit none
  integer, intent(in) :: ncid, varid
  real(kind=4), dimension(:,:,:), intent(out) :: fp
  integer (C_INT) :: c_ncid, c_varid
  real (C_FLOAT), dimension(size(fp,1),size(fp,2),size(fp,3)), target :: c_fp
  c_ncid = ncid
  c_varid = varid - 1
  ierr = nc_get_var_float(c_ncid,c_varid,C_LOC(c_fp))
  fp = c_fp
end function nf_get_var_real_d3

integer function nf_get_var_real_d4(ncid,varid,fp) result(ierr)
  implicit none
  integer, intent(in) :: ncid, varid
  real(kind=4), dimension(:,:,:,:), intent(out) :: fp
  integer (C_INT) :: c_ncid, c_varid
  real (C_FLOAT), dimension(size(fp,1),size(fp,2),size(fp,3),size(fp,4)), target :: c_fp
  c_ncid = ncid
  c_varid = varid - 1
  ierr = nc_get_var_float(c_ncid,c_varid,C_LOC(c_fp))
  fp = c_fp
end function nf_get_var_real_d4

integer function nf_get_var_real_d5(ncid,varid,fp) result(ierr)
  implicit none
  integer, intent(in) :: ncid, varid
  real(kind=4), dimension(:,:,:,:,:), intent(out) :: fp
  integer (C_INT) :: c_ncid, c_varid
  real (C_FLOAT), dimension(size(fp,1),size(fp,2),size(fp,3),size(fp,4),size(fp,5)), &
      target :: c_fp
  c_ncid = ncid
  c_varid = varid - 1
  ierr = nc_get_var_float(c_ncid,c_varid,C_LOC(c_fp))
  fp = c_fp
end function nf_get_var_real_d5

integer function nf_get_var_double_d1(ncid,varid,dp) result(ierr)
  implicit none
  integer, intent(in) :: ncid, varid
  real(kind=8), dimension(:), intent(out) :: dp
  integer (C_INT) :: c_ncid, c_varid
  real (C_DOUBLE), dimension(size(dp)), target :: c_dp
  c_ncid = ncid
  c_varid = varid - 1
  ierr = nc_get_var_double(c_ncid,c_varid,C_LOC(c_dp))
  dp = c_dp
end function nf_get_var_double_d1

integer function nf_get_var_double_d2(ncid,varid,dp) result(ierr)
  implicit none
  integer, intent(in) :: ncid, varid
  real(kind=8), dimension(:,:), intent(out) :: dp
  integer (C_INT) :: c_ncid, c_varid
  real (C_DOUBLE), dimension(size(dp,1),size(dp,2)), target :: c_dp
  c_ncid = ncid
  c_varid = varid - 1
  ierr = nc_get_var_double(c_ncid,c_varid,C_LOC(c_dp))
  dp = c_dp
end function nf_get_var_double_d2

integer function nf_get_var_double_d3(ncid,varid,dp) result(ierr)
  implicit none
  integer, intent(in) :: ncid, varid
  real(kind=8), dimension(:,:,:), intent(out) :: dp
  integer (C_INT) :: c_ncid, c_varid
  real (C_DOUBLE), dimension(size(dp,1),size(dp,2),size(dp,3)), target :: c_dp
  c_ncid = ncid
  c_varid = varid - 1
  ierr = nc_get_var_double(c_ncid,c_varid,C_LOC(c_dp))
  dp = c_dp
end function nf_get_var_double_d3

integer function nf_get_var_double_d4(ncid,varid,dp) result(ierr)
  implicit none
  integer, intent(in) :: ncid, varid
  real(kind=8), dimension(:,:,:,:), intent(out) :: dp
  integer (C_INT) :: c_ncid, c_varid
  real (C_DOUBLE), dimension(size(dp,1),size(dp,2),size(dp,3),size(dp,4)), target :: c_dp
  c_ncid = ncid
  c_varid = varid - 1
  ierr = nc_get_var_double(c_ncid,c_varid,C_LOC(c_dp))
  dp = c_dp
end function nf_get_var_double_d4

integer function nf_get_var_double_d5(ncid,varid,dp) result(ierr)
  implicit none
  integer, intent(in) :: ncid, varid
  real(kind=8), dimension(:,:,:,:,:), intent(out) :: dp
  integer (C_INT) :: c_ncid, c_varid
  real (C_DOUBLE), dimension(size(dp,1),size(dp,2),size(dp,3),size(dp,4),size(dp,5)), &
      target :: c_dp
  c_ncid = ncid
  c_varid = varid - 1
  ierr = nc_get_var_double(c_ncid,c_varid,C_LOC(c_dp))
  dp = c_dp
end function nf_get_var_double_d5

integer function nf_get_var1_real_s(ncid,varid,start,fp) result(ierr)
  implicit none
  integer, intent(in) :: ncid, varid
  integer, intent(in) :: start
  real(kind=4), intent(out) :: fp
  integer (C_INT) :: c_ncid, c_varid
  integer (C_SIZE_T), dimension(1) :: c_start
  real (C_FLOAT), target :: c_fp
  c_ncid = ncid
  c_varid = varid - 1
  c_start = start - 1
  ierr = nc_get_var1_float(c_ncid,c_varid,c_start,C_LOC(c_fp))
  fp = c_fp
end function nf_get_var1_real_s

integer function nf_get_var1_real_v(ncid,varid,start,fp) result(ierr)
  implicit none
  integer, intent(in) :: ncid, varid
  integer, dimension(:), intent(in) :: start
  real(kind=4), intent(out) :: fp
  integer (C_INT) :: c_ncid, c_varid
  integer (C_INT), target :: c_ndims
  integer (C_SIZE_T), dimension(size(start)) :: c_start
  real (C_FLOAT), target :: c_fp
  integer i, ndims
  c_ncid = ncid
  c_varid = varid - 1
  c_ndims = size(start)
  ierr = nc_inq_varndims(c_ncid,c_varid,C_LOC(c_ndims))
  ndims = c_ndims
  do i = 1,ndims
    c_start(ndims-i+1) = start(i) - 1
  end do
  ierr = nc_get_var1_float(c_ncid,c_varid,c_start,C_LOC(c_fp))
  fp = c_fp
end function nf_get_var1_real_v

integer function nf_get_var1_int_s(ncid,varid,start,ip) result(ierr)
  implicit none
  integer, intent(in) :: ncid, varid
  integer, intent(in) :: start
  integer(kind=4), intent(out) :: ip
  integer (C_INT) :: c_ncid, c_varid
  integer (C_SIZE_T), dimension(1) :: c_start
  integer (C_INT), target :: c_ip
  c_ncid = ncid
  c_varid = varid - 1
  c_start = start - 1
  ierr = nc_get_var1_int(c_ncid,c_varid,c_start,C_LOC(c_ip))
  ip = c_ip
end function nf_get_var1_int_s

integer function nf_get_var1_int_v(ncid,varid,start,ip) result(ierr)
  implicit none
  integer, intent(in) :: ncid, varid
  integer, dimension(:), intent(in) :: start
  integer(kind=4), intent(out) :: ip
  integer (C_INT) :: c_ncid, c_varid
  integer (C_INT), target :: c_ndims
  integer (C_SIZE_T), dimension(size(start)) :: c_start
  integer (C_INT), target :: c_ip
  integer i, ndims
  c_ncid = ncid
  c_varid = varid - 1
  c_ndims = size(start)
  ierr = nc_inq_varndims(c_ncid,c_varid,C_LOC(c_ndims))
  ndims = c_ndims
  do i = 1,ndims
    c_start(ndims-i+1) = start(i) - 1
  end do
  ierr = nc_get_var1_int(c_ncid,c_varid,c_start,C_LOC(c_ip))
  ip = c_ip
end function nf_get_var1_int_v

integer function nf_get_var1_int1_s(ncid,varid,start,bp) result(ierr)
  implicit none
  integer, intent(in) :: ncid, varid
  integer, intent(in) :: start
  integer(kind=1), intent(out) :: bp
  integer (C_INT) :: c_ncid, c_varid
  integer (C_SIZE_T), dimension(1) :: c_start
  integer (C_INT), target :: c_bp
  c_ncid = ncid
  c_varid = varid - 1
  c_start = start - 1
  ierr = nc_get_var1_schar(c_ncid,c_varid,c_start,C_LOC(c_bp))
  bp = c_bp
end function nf_get_var1_int1_s

integer function nf_get_var1_int1_v(ncid,varid,start,bp) result(ierr)
  implicit none
  integer, intent(in) :: ncid, varid
  integer, dimension(:), intent(in) :: start
  integer(kind=1), intent(out) :: bp
  integer (C_INT) :: c_ncid, c_varid
  integer (C_INT), target :: c_ndims
  integer (C_SIZE_T), dimension(size(start)) :: c_start
  integer (C_INT), target :: c_bp
  integer i, ndims
  c_ncid = ncid
  c_varid = varid - 1
  c_ndims = size(start)
  ierr = nc_inq_varndims(c_ncid,c_varid,C_LOC(c_ndims))
  ndims = c_ndims
  do i = 1,ndims
    c_start(ndims-i+1) = start(i) - 1
  end do
  ierr = nc_get_var1_schar(c_ncid,c_varid,c_start,C_LOC(c_bp))
  bp = c_bp
end function nf_get_var1_int1_v

integer function nf_get_var1_int2_s(ncid,varid,start,sp) result(ierr)
  implicit none
  integer, intent(in) :: ncid, varid
  integer, intent(in) :: start
  integer(kind=2), intent(out) :: sp
  integer (C_INT) :: c_ncid, c_varid
  integer (C_SIZE_T), dimension(1) :: c_start
  integer (C_INT), target :: c_sp
  c_ncid = ncid
  c_varid = varid - 1
  c_start = start - 1
  ierr = nc_get_var1_short(c_ncid,c_varid,c_start,C_LOC(c_sp))
  sp = c_sp
end function nf_get_var1_int2_s

integer function nf_get_var1_int2_v(ncid,varid,start,sp) result(ierr)
  implicit none
  integer, intent(in) :: ncid, varid
  integer, dimension(:), intent(in) :: start
  integer(kind=2), intent(out) :: sp
  integer (C_INT) :: c_ncid, c_varid
  integer (C_INT), target :: c_ndims
  integer (C_SIZE_T), dimension(size(start)) :: c_start
  integer (C_INT), target :: c_sp
  integer i, ndims
  c_ncid = ncid
  c_varid = varid - 1
  c_ndims = size(start)
  ierr = nc_inq_varndims(c_ncid,c_varid,C_LOC(c_ndims))
  ndims = c_ndims
  do i = 1,ndims
    c_start(ndims-i+1) = start(i) - 1
  end do
  ierr = nc_get_var1_short(c_ncid,c_varid,c_start,C_LOC(c_sp))
  sp = c_sp
end function nf_get_var1_int2_v

integer function nf_get_var1_double_s(ncid,varid,start,dp) result(ierr)
  implicit none
  integer, intent(in) :: ncid, varid
  integer, intent(in) :: start
  real(kind=8), intent(out) :: dp
  integer (C_INT) :: c_ncid, c_varid
  integer (C_SIZE_T), dimension(1) :: c_start
  real (C_DOUBLE), target :: c_dp
  c_ncid = ncid
  c_varid = varid - 1
  c_start = start - 1
  ierr = nc_get_var1_double(c_ncid,c_varid,c_start,C_LOC(c_dp))
  dp = c_dp
end function nf_get_var1_double_s

integer function nf_get_var1_double_v(ncid,varid,start,dp) result(ierr)
  implicit none
  integer, intent(in) :: ncid, varid
  integer, dimension(:), intent(in) :: start
  real(kind=8), intent(out) :: dp
  integer (C_INT) :: c_ncid, c_varid
  integer (C_INT), target :: c_ndims
  integer (C_SIZE_T), dimension(size(start)) :: c_start
  real (C_DOUBLE), target :: c_dp
  integer i, ndims
  c_ncid = ncid
  c_varid = varid - 1
  c_ndims = size(start)
  ierr = nc_inq_varndims(c_ncid,c_varid,C_LOC(c_ndims))
  ndims = c_ndims
  do i = 1,ndims
    c_start(ndims-i+1) = start(i) - 1
  end do
  ierr = nc_get_var1_double(c_ncid,c_varid,c_start,C_LOC(c_dp))
  dp = c_dp
end function nf_get_var1_double_v

integer function nf_get_vars_text(ncid,varid,start,ncount,stride,tp) result(ierr)
  implicit none
  integer, intent(in) :: ncid, varid
  integer, dimension(:), intent(in) :: start, ncount, stride  
  character(len=*), intent(out) :: tp
  integer (C_INT) :: c_ncid, c_varid
  integer (C_INT), target :: c_ndims
  integer (C_SIZE_T), dimension(size(start)) :: c_start, c_ncount
  integer (C_INTPTR_T), dimension(size(stride)) :: c_stride
  character, dimension(len(tp)), target :: c_tp
  integer i, ndims
  c_ncid = ncid
  c_varid = varid - 1
  c_ndims = size(start)
  ierr = nc_inq_varndims(c_ncid,c_varid,C_LOC(c_ndims))
  ndims = c_ndims
  do i = 1,ndims
    c_start(ndims-i+1) = start(i) - 1
    c_ncount(ndims-i+1) = ncount(i)
    c_stride(ndims-i+1) = stride(i)
  end do
  ierr = nc_get_vars_text(c_ncid,c_varid,c_start,c_ncount,c_stride,C_LOC(c_tp))
  call fc_strcopy(c_tp,tp)
end function nf_get_vars_text

integer function nf_get_vars_int1_d1(ncid,varid,start,ncount,stride,bp) result(ierr)
  implicit none
  integer, intent(in) :: ncid, varid
  integer, dimension(:), intent(in) :: start, ncount, stride
  integer(kind=1), dimension(:), intent(out) :: bp
  integer (C_INT) :: c_ncid, c_varid
  integer (C_INT), target :: c_ndims
  integer (C_SIZE_T), dimension(size(start)) :: c_start
  integer (C_SIZE_T), dimension(size(ncount)) :: c_ncount
  integer (C_INTPTR_T), dimension(size(stride)) :: c_stride
  integer (C_SIGNED_CHAR), dimension(size(bp)), target :: c_bp
  integer i, ndims
  c_ncid = ncid
  c_varid = varid - 1
  c_ndims = size(start)
  ierr = nc_inq_varndims(c_ncid,c_varid,C_LOC(c_ndims))
  ndims = c_ndims
  do i = 1,ndims
    c_start(ndims-i+1) = start(i) - 1
    c_ncount(ndims-i+1) = ncount(i)
    c_stride(ndims-i+1) = stride(i)
  end do
  ierr = nc_get_vars_schar(c_ncid,c_varid,c_start,c_ncount,c_stride,C_LOC(c_bp))
  bp = c_bp
end function nf_get_vars_int1_d1

integer function nf_get_vars_int1_d2(ncid,varid,start,ncount,stride,bp) result(ierr)
  implicit none
  integer, intent(in) :: ncid, varid
  integer, dimension(:), intent(in) :: start, ncount, stride
  integer(kind=1), dimension(:,:), intent(out) :: bp
  integer (C_INT) :: c_ncid, c_varid
  integer (C_INT), target :: c_ndims
  integer (C_SIZE_T), dimension(size(start)) :: c_start
  integer (C_SIZE_T), dimension(size(ncount)) :: c_ncount
  integer (C_INTPTR_T), dimension(size(stride)) :: c_stride
  integer (C_SIGNED_CHAR), dimension(size(bp,1),size(bp,2)), target :: c_bp
  integer i, ndims
  c_ncid = ncid
  c_varid = varid - 1
  c_ndims = size(start)
  ierr = nc_inq_varndims(c_ncid,c_varid,C_LOC(c_ndims))
  ndims = c_ndims
  do i = 1,ndims
    c_start(ndims-i+1) = start(i) - 1
    c_ncount(ndims-i+1) = ncount(i)
    c_stride(ndims-i+1) = stride(i)
  end do
  ierr = nc_get_vars_schar(c_ncid,c_varid,c_start,c_ncount,c_stride,C_LOC(c_bp))
  bp = c_bp
end function nf_get_vars_int1_d2

integer function nf_get_vars_int1_d3(ncid,varid,start,ncount,stride,bp) result(ierr)
  implicit none
  integer, intent(in) :: ncid, varid
  integer, dimension(:), intent(in) :: start, ncount, stride
  integer(kind=1), dimension(:,:,:), intent(out) :: bp
  integer (C_INT) :: c_ncid, c_varid
  integer (C_INT), target :: c_ndims
  integer (C_SIZE_T), dimension(size(start)) :: c_start
  integer (C_SIZE_T), dimension(size(ncount)) :: c_ncount
  integer (C_INTPTR_T), dimension(size(stride)) :: c_stride
  integer (C_SIGNED_CHAR), dimension(size(bp,1),size(bp,2),size(bp,3)), target :: c_bp
  integer i, ndims
  c_ncid = ncid
  c_varid = varid - 1
  c_ndims = size(start)
  ierr = nc_inq_varndims(c_ncid,c_varid,C_LOC(c_ndims))
  ndims = c_ndims
  do i = 1,ndims
    c_start(ndims-i+1) = start(i) - 1
    c_ncount(ndims-i+1) = ncount(i)
    c_stride(ndims-i+1) = stride(i)
  end do
  ierr = nc_get_vars_schar(c_ncid,c_varid,c_start,c_ncount,c_stride,C_LOC(c_bp))
  bp = c_bp
end function nf_get_vars_int1_d3

integer function nf_get_vars_int1_d4(ncid,varid,start,ncount,stride,bp) result(ierr)
  implicit none
  integer, intent(in) :: ncid, varid
  integer, dimension(:), intent(in) :: start, ncount, stride
  integer(kind=1), dimension(:,:,:,:), intent(out) :: bp
  integer (C_INT) :: c_ncid, c_varid
  integer (C_INT), target :: c_ndims
  integer (C_SIZE_T), dimension(size(start)) :: c_start
  integer (C_SIZE_T), dimension(size(ncount)) :: c_ncount
  integer (C_INTPTR_T), dimension(size(stride)) :: c_stride
  integer (C_SIGNED_CHAR), dimension(size(bp,1),size(bp,2),size(bp,3),size(bp,4)), target :: c_bp
  integer i, ndims
  c_ncid = ncid
  c_varid = varid - 1
  c_ndims = size(start)
  ierr = nc_inq_varndims(c_ncid,c_varid,C_LOC(c_ndims))
  ndims = c_ndims
  do i = 1,ndims
    c_start(ndims-i+1) = start(i) - 1
    c_ncount(ndims-i+1) = ncount(i)
    c_stride(ndims-i+1) = stride(i)
  end do
  ierr = nc_get_vars_schar(c_ncid,c_varid,c_start,c_ncount,c_stride,C_LOC(c_bp))
  bp = c_bp
end function nf_get_vars_int1_d4

integer function nf_get_vars_int1_d5(ncid,varid,start,ncount,stride,bp) result(ierr)
  implicit none
  integer, intent(in) :: ncid, varid
  integer, dimension(:), intent(in) :: start, ncount, stride
  integer(kind=1), dimension(:,:,:,:,:), intent(out) :: bp
  integer (C_INT) :: c_ncid, c_varid
  integer (C_INT), target :: c_ndims
  integer (C_SIZE_T), dimension(size(start)) :: c_start
  integer (C_SIZE_T), dimension(size(ncount)) :: c_ncount
  integer (C_INTPTR_T), dimension(size(stride)) :: c_stride
  integer (C_SIGNED_CHAR), dimension(size(bp,1),size(bp,2),size(bp,3),size(bp,4),size(bp,5)), &
      target :: c_bp
  integer i, ndims
  c_ncid = ncid
  c_varid = varid - 1
  c_ndims = size(start)
  ierr = nc_inq_varndims(c_ncid,c_varid,C_LOC(c_ndims))
  ndims = c_ndims
  do i = 1,ndims
    c_start(ndims-i+1) = start(i) - 1
    c_ncount(ndims-i+1) = ncount(i)
    c_stride(ndims-i+1) = stride(i)
  end do
  ierr = nc_get_vars_schar(c_ncid,c_varid,c_start,c_ncount,c_stride,C_LOC(c_bp))
  bp = c_bp
end function nf_get_vars_int1_d5

integer function nf_get_vars_int1_d6(ncid,varid,start,ncount,stride,bp) result(ierr)
  implicit none
  integer, intent(in) :: ncid, varid
  integer, dimension(:), intent(in) :: start, ncount, stride
  integer(kind=1), dimension(:,:,:,:,:,:), intent(out) :: bp
  integer (C_INT) :: c_ncid, c_varid
  integer (C_INT), target :: c_ndims
  integer (C_SIZE_T), dimension(size(start)) :: c_start
  integer (C_SIZE_T), dimension(size(ncount)) :: c_ncount
  integer (C_INTPTR_T), dimension(size(stride)) :: c_stride
  integer (C_SIGNED_CHAR), dimension(size(bp,1),size(bp,2),size(bp,3),size(bp,4),size(bp,5), &
      size(bp,6)), target :: c_bp
  integer i, ndims
  c_ncid = ncid
  c_varid = varid - 1
  c_ndims = size(start)
  ierr = nc_inq_varndims(c_ncid,c_varid,C_LOC(c_ndims))
  ndims = c_ndims
  do i = 1,ndims
    c_start(ndims-i+1) = start(i) - 1
    c_ncount(ndims-i+1) = ncount(i)
    c_stride(ndims-i+1) = stride(i)
  end do
  ierr = nc_get_vars_schar(c_ncid,c_varid,c_start,c_ncount,c_stride,C_LOC(c_bp))
  bp = c_bp
end function nf_get_vars_int1_d6

integer function nf_get_vars_int1_d7(ncid,varid,start,ncount,stride,bp) result(ierr)
  implicit none
  integer, intent(in) :: ncid, varid
  integer, dimension(:), intent(in) :: start, ncount, stride
  integer(kind=1), dimension(:,:,:,:,:,:,:), intent(out) :: bp
  integer (C_INT) :: c_ncid, c_varid
  integer (C_INT), target :: c_ndims
  integer (C_SIZE_T), dimension(size(start)) :: c_start
  integer (C_SIZE_T), dimension(size(ncount)) :: c_ncount
  integer (C_INTPTR_T), dimension(size(stride)) :: c_stride
  integer (C_SIGNED_CHAR), dimension(size(bp,1),size(bp,2),size(bp,3),size(bp,4),size(bp,5), &
      size(bp,6),size(bp,7)), target :: c_bp
  integer i, ndims
  c_ncid = ncid
  c_varid = varid - 1
  c_ndims = size(start)
  ierr = nc_inq_varndims(c_ncid,c_varid,C_LOC(c_ndims))
  ndims = c_ndims
  do i = 1,ndims
    c_start(ndims-i+1) = start(i) - 1
    c_ncount(ndims-i+1) = ncount(i)
    c_stride(ndims-i+1) = stride(i)
  end do
  ierr = nc_get_vars_schar(c_ncid,c_varid,c_start,c_ncount,c_stride,C_LOC(c_bp))
  bp = c_bp
end function nf_get_vars_int1_d7

integer function nf_get_vars_int2_d1(ncid,varid,start,ncount,stride,sp) result(ierr)
  implicit none
  integer, intent(in) :: ncid, varid
  integer, dimension(:), intent(in) :: start, ncount, stride
  integer(kind=2), dimension(:), intent(out) :: sp
  integer (C_INT) :: c_ncid, c_varid
  integer (C_INT), target :: c_ndims
  integer (C_SIZE_T), dimension(size(start)) :: c_start
  integer (C_SIZE_T), dimension(size(ncount)) :: c_ncount
  integer (C_INTPTR_T), dimension(size(stride)) :: c_stride
  integer (C_SHORT), dimension(size(sp)), target :: c_sp
  integer i, ndims
  c_ncid = ncid
  c_varid = varid - 1
  c_ndims = size(start)
  ierr = nc_inq_varndims(c_ncid,c_varid,C_LOC(c_ndims))
  ndims = c_ndims
  do i = 1,ndims
    c_start(ndims-i+1) = start(i) - 1
    c_ncount(ndims-i+1) = ncount(i)
    c_stride(ndims-i+1) = stride(i)
  end do
  ierr = nc_get_vars_short(c_ncid,c_varid,c_start,c_ncount,c_stride,C_LOC(c_sp))
  sp = c_sp
end function nf_get_vars_int2_d1

integer function nf_get_vars_int2_d2(ncid,varid,start,ncount,stride,sp) result(ierr)
  implicit none
  integer, intent(in) :: ncid, varid
  integer, dimension(:), intent(in) :: start, ncount, stride
  integer(kind=2), dimension(:,:), intent(out) :: sp
  integer (C_INT) :: c_ncid, c_varid
  integer (C_INT), target :: c_ndims
  integer (C_SIZE_T), dimension(size(start)) :: c_start
  integer (C_SIZE_T), dimension(size(ncount)) :: c_ncount
  integer (C_INTPTR_T), dimension(size(stride)) :: c_stride
  integer (C_SHORT), dimension(size(sp,1),size(sp,2)), target :: c_sp
  integer i, ndims
  c_ncid = ncid
  c_varid = varid - 1
  c_ndims = size(start)
  ierr = nc_inq_varndims(c_ncid,c_varid,C_LOC(c_ndims))
  ndims = c_ndims
  do i = 1,ndims
    c_start(ndims-i+1) = start(i) - 1
    c_ncount(ndims-i+1) = ncount(i)
    c_stride(ndims-i+1) = stride(i)
  end do
  ierr = nc_get_vars_short(c_ncid,c_varid,c_start,c_ncount,c_stride,C_LOC(c_sp))
  sp = c_sp
end function nf_get_vars_int2_d2

integer function nf_get_vars_int2_d3(ncid,varid,start,ncount,stride,sp) result(ierr)
  implicit none
  integer, intent(in) :: ncid, varid
  integer, dimension(:), intent(in) :: start, ncount, stride
  integer(kind=2), dimension(:,:,:), intent(out) :: sp
  integer (C_INT) :: c_ncid, c_varid
  integer (C_INT), target :: c_ndims
  integer (C_SIZE_T), dimension(size(start)) :: c_start
  integer (C_SIZE_T), dimension(size(ncount)) :: c_ncount
  integer (C_INTPTR_T), dimension(size(stride)) :: c_stride
  integer (C_SHORT), dimension(size(sp,1),size(sp,2),size(sp,3)), target :: c_sp
  integer i, ndims
  c_ncid = ncid
  c_varid = varid - 1
  c_ndims = size(start)
  ierr = nc_inq_varndims(c_ncid,c_varid,C_LOC(c_ndims))
  ndims = c_ndims
  do i = 1,ndims
    c_start(ndims-i+1) = start(i) - 1
    c_ncount(ndims-i+1) = ncount(i)
    c_stride(ndims-i+1) = stride(i)
  end do
  ierr = nc_get_vars_short(c_ncid,c_varid,c_start,c_ncount,c_stride,C_LOC(c_sp))
  sp = c_sp
end function nf_get_vars_int2_d3

integer function nf_get_vars_int2_d4(ncid,varid,start,ncount,stride,sp) result(ierr)
  implicit none
  integer, intent(in) :: ncid, varid
  integer, dimension(:), intent(in) :: start, ncount, stride
  integer(kind=2), dimension(:,:,:,:), intent(out) :: sp
  integer (C_INT) :: c_ncid, c_varid
  integer (C_INT), target :: c_ndims
  integer (C_SIZE_T), dimension(size(start)) :: c_start
  integer (C_SIZE_T), dimension(size(ncount)) :: c_ncount
  integer (C_INTPTR_T), dimension(size(stride)) :: c_stride
  integer (C_SHORT), dimension(size(sp,1),size(sp,2),size(sp,3),size(sp,4)), target :: c_sp
  integer i, ndims
  c_ncid = ncid
  c_varid = varid - 1
  c_ndims = size(start)
  ierr = nc_inq_varndims(c_ncid,c_varid,C_LOC(c_ndims))
  ndims = c_ndims
  do i = 1,ndims
    c_start(ndims-i+1) = start(i) - 1
    c_ncount(ndims-i+1) = ncount(i)
    c_stride(ndims-i+1) = stride(i)
  end do
  ierr = nc_get_vars_short(c_ncid,c_varid,c_start,c_ncount,c_stride,C_LOC(c_sp))
  sp = c_sp
end function nf_get_vars_int2_d4

integer function nf_get_vars_int2_d5(ncid,varid,start,ncount,stride,sp) result(ierr)
  implicit none
  integer, intent(in) :: ncid, varid
  integer, dimension(:), intent(in) :: start, ncount, stride
  integer(kind=2), dimension(:,:,:,:,:), intent(out) :: sp
  integer (C_INT) :: c_ncid, c_varid
  integer (C_INT), target :: c_ndims
  integer (C_SIZE_T), dimension(size(start)) :: c_start
  integer (C_SIZE_T), dimension(size(ncount)) :: c_ncount
  integer (C_INTPTR_T), dimension(size(stride)) :: c_stride
  integer (C_SHORT), dimension(size(sp,1),size(sp,2),size(sp,3),size(sp,4),size(sp,5)), &
      target :: c_sp
  integer i, ndims
  c_ncid = ncid
  c_varid = varid - 1
  c_ndims = size(start)
  ierr = nc_inq_varndims(c_ncid,c_varid,C_LOC(c_ndims))
  ndims = c_ndims
  do i = 1,ndims
    c_start(ndims-i+1) = start(i) - 1
    c_ncount(ndims-i+1) = ncount(i)
    c_stride(ndims-i+1) = stride(i)
  end do
  ierr = nc_get_vars_short(c_ncid,c_varid,c_start,c_ncount,c_stride,C_LOC(c_sp))
  sp = c_sp
end function nf_get_vars_int2_d5

integer function nf_get_vars_int2_d6(ncid,varid,start,ncount,stride,sp) result(ierr)
  implicit none
  integer, intent(in) :: ncid, varid
  integer, dimension(:), intent(in) :: start, ncount, stride
  integer(kind=2), dimension(:,:,:,:,:,:), intent(out) :: sp
  integer (C_INT) :: c_ncid, c_varid
  integer (C_INT), target :: c_ndims
  integer (C_SIZE_T), dimension(size(start)) :: c_start
  integer (C_SIZE_T), dimension(size(ncount)) :: c_ncount
  integer (C_INTPTR_T), dimension(size(stride)) :: c_stride
  integer (C_SHORT), dimension(size(sp,1),size(sp,2),size(sp,3),size(sp,4),size(sp,5), &
      size(sp,6)), target :: c_sp
  integer i, ndims
  c_ncid = ncid
  c_varid = varid - 1
  c_ndims = size(start)
  ierr = nc_inq_varndims(c_ncid,c_varid,C_LOC(c_ndims))
  ndims = c_ndims
  do i = 1,ndims
    c_start(ndims-i+1) = start(i) - 1
    c_ncount(ndims-i+1) = ncount(i)
    c_stride(ndims-i+1) = stride(i)
  end do
  ierr = nc_get_vars_short(c_ncid,c_varid,c_start,c_ncount,c_stride,C_LOC(c_sp))
  sp = c_sp
end function nf_get_vars_int2_d6

integer function nf_get_vars_int2_d7(ncid,varid,start,ncount,stride,sp) result(ierr)
  implicit none
  integer, intent(in) :: ncid, varid
  integer, dimension(:), intent(in) :: start, ncount, stride
  integer(kind=2), dimension(:,:,:,:,:,:,:), intent(out) :: sp
  integer (C_INT) :: c_ncid, c_varid
  integer (C_INT), target :: c_ndims
  integer (C_SIZE_T), dimension(size(start)) :: c_start
  integer (C_SIZE_T), dimension(size(ncount)) :: c_ncount
  integer (C_INTPTR_T), dimension(size(stride)) :: c_stride
  integer (C_SHORT), dimension(size(sp,1),size(sp,2),size(sp,3),size(sp,4),size(sp,5), &
      size(sp,6),size(sp,7)), target :: c_sp
  integer i, ndims
  c_ncid = ncid
  c_varid = varid - 1
  c_ndims = size(start)
  ierr = nc_inq_varndims(c_ncid,c_varid,C_LOC(c_ndims))
  ndims = c_ndims
  do i = 1,ndims
    c_start(ndims-i+1) = start(i) - 1
    c_ncount(ndims-i+1) = ncount(i)
    c_stride(ndims-i+1) = stride(i)
  end do
  ierr = nc_get_vars_short(c_ncid,c_varid,c_start,c_ncount,c_stride,C_LOC(c_sp))
  sp = c_sp
end function nf_get_vars_int2_d7

integer function nf_get_vars_int_d1(ncid,varid,start,ncount,stride,ip) result(ierr)
  implicit none
  integer, intent(in) :: ncid, varid
  integer, dimension(:), intent(in) :: start, ncount, stride
  integer(kind=4), dimension(:), intent(out) :: ip
  integer (C_INT) :: c_ncid, c_varid
  integer (C_INT), target :: c_ndims
  integer (C_SIZE_T), dimension(size(start)) :: c_start
  integer (C_SIZE_T), dimension(size(ncount)) :: c_ncount
  integer (C_INTPTR_T), dimension(size(stride)) :: c_stride
  integer (C_INT), dimension(size(ip)), target :: c_ip
  integer i, ndims
  c_ncid = ncid
  c_varid = varid - 1
  c_ndims = size(start)
  ierr = nc_inq_varndims(c_ncid,c_varid,C_LOC(c_ndims))
  ndims = c_ndims
  do i = 1,ndims
    c_start(ndims-i+1) = start(i) - 1
    c_ncount(ndims-i+1) = ncount(i)
    c_stride(ndims-i+1) = stride(i)
  end do
  ierr = nc_get_vars_int(c_ncid,c_varid,c_start,c_ncount,c_stride,C_LOC(c_ip))
  ip = c_ip
end function nf_get_vars_int_d1

integer function nf_get_vars_int_d2(ncid,varid,start,ncount,stride,ip) result(ierr)
  implicit none
  integer, intent(in) :: ncid, varid
  integer, dimension(:), intent(in) :: start, ncount, stride
  integer(kind=4), dimension(:,:), intent(out) :: ip
  integer (C_INT) :: c_ncid, c_varid
  integer (C_INT), target :: c_ndims
  integer (C_SIZE_T), dimension(size(start)) :: c_start
  integer (C_SIZE_T), dimension(size(ncount)) :: c_ncount
  integer (C_INTPTR_T), dimension(size(stride)) :: c_stride
  integer (C_INT), dimension(size(ip,1),size(ip,2)), target :: c_ip
  integer i, ndims
  c_ncid = ncid
  c_varid = varid - 1
  c_ndims = size(start)
  ierr = nc_inq_varndims(c_ncid,c_varid,C_LOC(c_ndims))
  ndims = c_ndims
  do i = 1,ndims
    c_start(ndims-i+1) = start(i) - 1
    c_ncount(ndims-i+1) = ncount(i)
    c_stride(ndims-i+1) = stride(i)
  end do
  ierr = nc_get_vars_int(c_ncid,c_varid,c_start,c_ncount,c_stride,C_LOC(c_ip))
  ip = c_ip
end function nf_get_vars_int_d2

integer function nf_get_vars_int_d3(ncid,varid,start,ncount,stride,ip) result(ierr)
  implicit none
  integer, intent(in) :: ncid, varid
  integer, dimension(:), intent(in) :: start, ncount, stride
  integer(kind=4), dimension(:,:,:), intent(out) :: ip
  integer (C_INT) :: c_ncid, c_varid
  integer (C_INT), target :: c_ndims
  integer (C_SIZE_T), dimension(size(start)) :: c_start
  integer (C_SIZE_T), dimension(size(ncount)) :: c_ncount
  integer (C_INTPTR_T), dimension(size(stride)) :: c_stride
  integer (C_INT), dimension(size(ip,1),size(ip,2),size(ip,3)), target :: c_ip
  integer i, ndims
  c_ncid = ncid
  c_varid = varid - 1
  c_ndims = size(start)
  ierr = nc_inq_varndims(c_ncid,c_varid,C_LOC(c_ndims))
  ndims = c_ndims
  do i = 1,ndims
    c_start(ndims-i+1) = start(i) - 1
    c_ncount(ndims-i+1) = ncount(i)
    c_stride(ndims-i+1) = stride(i)
  end do
  ierr = nc_get_vars_int(c_ncid,c_varid,c_start,c_ncount,c_stride,C_LOC(c_ip))
  ip = c_ip
end function nf_get_vars_int_d3

integer function nf_get_vars_int_d4(ncid,varid,start,ncount,stride,ip) result(ierr)
  implicit none
  integer, intent(in) :: ncid, varid
  integer, dimension(:), intent(in) :: start, ncount, stride
  integer(kind=4), dimension(:,:,:,:), intent(out) :: ip
  integer (C_INT) :: c_ncid, c_varid
  integer (C_INT), target :: c_ndims
  integer (C_SIZE_T), dimension(size(start)) :: c_start
  integer (C_SIZE_T), dimension(size(ncount)) :: c_ncount
  integer (C_INTPTR_T), dimension(size(stride)) :: c_stride
  integer (C_INT), dimension(size(ip,1),size(ip,2),size(ip,3),size(ip,4)), target :: c_ip
  integer i, ndims
  c_ncid = ncid
  c_varid = varid - 1
  c_ndims = size(start)
  ierr = nc_inq_varndims(c_ncid,c_varid,C_LOC(c_ndims))
  ndims = c_ndims
  do i = 1,ndims
    c_start(ndims-i+1) = start(i) - 1
    c_ncount(ndims-i+1) = ncount(i)
    c_stride(ndims-i+1) = stride(i)
  end do
  ierr = nc_get_vars_int(c_ncid,c_varid,c_start,c_ncount,c_stride,C_LOC(c_ip))
  ip = c_ip
end function nf_get_vars_int_d4

integer function nf_get_vars_int_d5(ncid,varid,start,ncount,stride,ip) result(ierr)
  implicit none
  integer, intent(in) :: ncid, varid
  integer, dimension(:), intent(in) :: start, ncount, stride
  integer(kind=4), dimension(:,:,:,:,:), intent(out) :: ip
  integer (C_INT) :: c_ncid, c_varid
  integer (C_INT), target :: c_ndims
  integer (C_SIZE_T), dimension(size(start)) :: c_start
  integer (C_SIZE_T), dimension(size(ncount)) :: c_ncount
  integer (C_INTPTR_T), dimension(size(stride)) :: c_stride
  integer (C_INT), dimension(size(ip,1),size(ip,2),size(ip,3),size(ip,4),size(ip,5)), &
      target :: c_ip
  integer i, ndims
  c_ncid = ncid
  c_varid = varid - 1
  c_ndims = size(start)
  ierr = nc_inq_varndims(c_ncid,c_varid,C_LOC(c_ndims))
  ndims = c_ndims
  do i = 1,ndims
    c_start(ndims-i+1) = start(i) - 1
    c_ncount(ndims-i+1) = ncount(i)
    c_stride(ndims-i+1) = stride(i)
  end do
  ierr = nc_get_vars_int(c_ncid,c_varid,c_start,c_ncount,c_stride,C_LOC(c_ip))
  ip = c_ip
end function nf_get_vars_int_d5

integer function nf_get_vars_int_d6(ncid,varid,start,ncount,stride,ip) result(ierr)
  implicit none
  integer, intent(in) :: ncid, varid
  integer, dimension(:), intent(in) :: start, ncount, stride
  integer(kind=4), dimension(:,:,:,:,:,:), intent(out) :: ip
  integer (C_INT) :: c_ncid, c_varid
  integer (C_INT), target :: c_ndims
  integer (C_SIZE_T), dimension(size(start)) :: c_start
  integer (C_SIZE_T), dimension(size(ncount)) :: c_ncount
  integer (C_INTPTR_T), dimension(size(stride)) :: c_stride
  integer (C_INT), dimension(size(ip,1),size(ip,2),size(ip,3),size(ip,4),size(ip,5), &
      size(ip,6)), target :: c_ip
  integer i, ndims
  c_ncid = ncid
  c_varid = varid - 1
  c_ndims = size(start)
  ierr = nc_inq_varndims(c_ncid,c_varid,C_LOC(c_ndims))
  ndims = c_ndims
  do i = 1,ndims
    c_start(ndims-i+1) = start(i) - 1
    c_ncount(ndims-i+1) = ncount(i)
    c_stride(ndims-i+1) = stride(i)
  end do
  ierr = nc_get_vars_int(c_ncid,c_varid,c_start,c_ncount,c_stride,C_LOC(c_ip))
  ip = c_ip
end function nf_get_vars_int_d6

integer function nf_get_vars_int_d7(ncid,varid,start,ncount,stride,ip) result(ierr)
  implicit none
  integer, intent(in) :: ncid, varid
  integer, dimension(:), intent(in) :: start, ncount, stride
  integer(kind=4), dimension(:,:,:,:,:,:,:), intent(out) :: ip
  integer (C_INT) :: c_ncid, c_varid
  integer (C_INT), target :: c_ndims
  integer (C_SIZE_T), dimension(size(start)) :: c_start
  integer (C_SIZE_T), dimension(size(ncount)) :: c_ncount
  integer (C_INTPTR_T), dimension(size(stride)) :: c_stride
  integer (C_INT), dimension(size(ip,1),size(ip,2),size(ip,3),size(ip,4),size(ip,5), &
      size(ip,6),size(ip,7)), target :: c_ip
  integer i, ndims
  c_ncid = ncid
  c_varid = varid - 1
  c_ndims = size(start)
  ierr = nc_inq_varndims(c_ncid,c_varid,C_LOC(c_ndims))
  ndims = c_ndims
  do i = 1,ndims
    c_start(ndims-i+1) = start(i) - 1
    c_ncount(ndims-i+1) = ncount(i)
    c_stride(ndims-i+1) = stride(i)
  end do
  ierr = nc_get_vars_int(c_ncid,c_varid,c_start,c_ncount,c_stride,C_LOC(c_ip))
  ip = c_ip
end function nf_get_vars_int_d7

integer function nf_get_vars_real_d1(ncid,varid,start,ncount,stride,fp) result(ierr)
  implicit none
  integer, intent(in) :: ncid, varid
  integer, dimension(:), intent(in) :: start, ncount, stride
  real(kind=4), dimension(:), intent(out) :: fp
  integer (C_INT) :: c_ncid, c_varid
  integer (C_INT), target :: c_ndims
  integer (C_SIZE_T), dimension(size(start)) :: c_start
  integer (C_SIZE_T), dimension(size(ncount)) :: c_ncount
  integer (C_INTPTR_T), dimension(size(stride)) :: c_stride
  real (C_FLOAT), dimension(size(fp)), target :: c_fp
  integer i, ndims
  c_ncid = ncid
  c_varid = varid - 1
  c_ndims = size(start)
  ierr = nc_inq_varndims(c_ncid,c_varid,C_LOC(c_ndims))
  ndims = c_ndims
  do i = 1,ndims
    c_start(ndims-i+1) = start(i) - 1
    c_ncount(ndims-i+1) = ncount(i)
    c_stride(ndims-i+1) = stride(i)
  end do
  ierr = nc_get_vars_float(c_ncid,c_varid,c_start,c_ncount,c_stride,C_LOC(c_fp))
  fp = c_fp
end function nf_get_vars_real_d1

integer function nf_get_vars_real_d2(ncid,varid,start,ncount,stride,fp) result(ierr)
  implicit none
  integer, intent(in) :: ncid, varid
  integer, dimension(:), intent(in) :: start, ncount, stride
  real(kind=4), dimension(:,:), intent(out) :: fp
  integer (C_INT) :: c_ncid, c_varid
  integer (C_INT), target :: c_ndims
  integer (C_SIZE_T), dimension(size(start)) :: c_start
  integer (C_SIZE_T), dimension(size(ncount)) :: c_ncount
  integer (C_INTPTR_T), dimension(size(stride)) :: c_stride
  real (C_FLOAT), dimension(size(fp,1),size(fp,2)), target :: c_fp
  integer i, ndims
  c_ncid = ncid
  c_varid = varid - 1
  c_ndims = size(start)
  ierr = nc_inq_varndims(c_ncid,c_varid,C_LOC(c_ndims))
  ndims = c_ndims
  do i = 1,ndims
    c_start(ndims-i+1) = start(i) - 1
    c_ncount(ndims-i+1) = ncount(i)
    c_stride(ndims-i+1) = stride(i)
  end do
  ierr = nc_get_vars_float(c_ncid,c_varid,c_start,c_ncount,c_stride,C_LOC(c_fp))
  fp = c_fp
end function nf_get_vars_real_d2

integer function nf_get_vars_real_d3(ncid,varid,start,ncount,stride,fp) result(ierr)
  implicit none
  integer, intent(in) :: ncid, varid
  integer, dimension(:), intent(in) :: start, ncount, stride
  real(kind=4), dimension(:,:,:), intent(out) :: fp
  integer (C_INT) :: c_ncid, c_varid
  integer (C_INT), target :: c_ndims
  integer (C_SIZE_T), dimension(size(start)) :: c_start
  integer (C_SIZE_T), dimension(size(ncount)) :: c_ncount
  integer (C_INTPTR_T), dimension(size(stride)) :: c_stride
  real (C_FLOAT), dimension(size(fp,1),size(fp,2),size(fp,3)), target :: c_fp
  integer i, ndims
  c_ncid = ncid
  c_varid = varid - 1
  c_ndims = size(start)
  ierr = nc_inq_varndims(c_ncid,c_varid,C_LOC(c_ndims))
  ndims = c_ndims
  do i = 1,ndims
    c_start(ndims-i+1) = start(i) - 1
    c_ncount(ndims-i+1) = ncount(i)
    c_stride(ndims-i+1) = stride(i)
  end do
  ierr = nc_get_vars_float(c_ncid,c_varid,c_start,c_ncount,c_stride,C_LOC(c_fp))
  fp = c_fp
end function nf_get_vars_real_d3

integer function nf_get_vars_real_d4(ncid,varid,start,ncount,stride,fp) result(ierr)
  implicit none
  integer, intent(in) :: ncid, varid
  integer, dimension(:), intent(in) :: start, ncount, stride
  real(kind=4), dimension(:,:,:,:), intent(out) :: fp
  integer (C_INT) :: c_ncid, c_varid
  integer (C_INT), target :: c_ndims
  integer (C_SIZE_T), dimension(size(start)) :: c_start
  integer (C_SIZE_T), dimension(size(ncount)) :: c_ncount
  integer (C_INTPTR_T), dimension(size(stride)) :: c_stride
  real (C_FLOAT), dimension(size(fp,1),size(fp,2),size(fp,3),size(fp,4)), target :: c_fp
  integer i, ndims
  c_ncid = ncid
  c_varid = varid - 1
  c_ndims = size(start)
  ierr = nc_inq_varndims(c_ncid,c_varid,C_LOC(c_ndims))
  ndims = c_ndims
  do i = 1,ndims
    c_start(ndims-i+1) = start(i) - 1
    c_ncount(ndims-i+1) = ncount(i)
    c_stride(ndims-i+1) = stride(i)
  end do
  ierr = nc_get_vars_float(c_ncid,c_varid,c_start,c_ncount,c_stride,C_LOC(c_fp))
  fp = c_fp
end function nf_get_vars_real_d4

integer function nf_get_vars_real_d5(ncid,varid,start,ncount,stride,fp) result(ierr)
  implicit none
  integer, intent(in) :: ncid, varid
  integer, dimension(:), intent(in) :: start, ncount, stride
  real(kind=4), dimension(:,:,:,:,:), intent(out) :: fp
  integer (C_INT) :: c_ncid, c_varid
  integer (C_INT), target :: c_ndims
  integer (C_SIZE_T), dimension(size(start)) :: c_start
  integer (C_SIZE_T), dimension(size(ncount)) :: c_ncount
  integer (C_INTPTR_T), dimension(size(stride)) :: c_stride
  real (C_FLOAT), dimension(size(fp,1),size(fp,2),size(fp,3),size(fp,4),size(fp,5)), &
      target :: c_fp
  integer i, ndims
  c_ncid = ncid
  c_varid = varid - 1
  c_ndims = size(start)
  ierr = nc_inq_varndims(c_ncid,c_varid,C_LOC(c_ndims))
  ndims = c_ndims
  do i = 1,ndims
    c_start(ndims-i+1) = start(i) - 1
    c_ncount(ndims-i+1) = ncount(i)
    c_stride(ndims-i+1) = stride(i)
  end do
  ierr = nc_get_vars_float(c_ncid,c_varid,c_start,c_ncount,c_stride,C_LOC(c_fp))
  fp = c_fp
end function nf_get_vars_real_d5

integer function nf_get_vars_real_d6(ncid,varid,start,ncount,stride,fp) result(ierr)
  implicit none
  integer, intent(in) :: ncid, varid
  integer, dimension(:), intent(in) :: start, ncount, stride
  real(kind=4), dimension(:,:,:,:,:,:), intent(out) :: fp
  integer (C_INT) :: c_ncid, c_varid
  integer (C_INT), target :: c_ndims
  integer (C_SIZE_T), dimension(size(start)) :: c_start
  integer (C_SIZE_T), dimension(size(ncount)) :: c_ncount
  integer (C_INTPTR_T), dimension(size(stride)) :: c_stride
  real (C_FLOAT), dimension(size(fp,1),size(fp,2),size(fp,3),size(fp,4),size(fp,5), &
      size(fp,6)), target :: c_fp
  integer i, ndims
  c_ncid = ncid
  c_varid = varid - 1
  c_ndims = size(start)
  ierr = nc_inq_varndims(c_ncid,c_varid,C_LOC(c_ndims))
  ndims = c_ndims
  do i = 1,ndims
    c_start(ndims-i+1) = start(i) - 1
    c_ncount(ndims-i+1) = ncount(i)
    c_stride(ndims-i+1) = stride(i)
  end do
  ierr = nc_get_vars_float(c_ncid,c_varid,c_start,c_ncount,c_stride,C_LOC(c_fp))
  fp = c_fp
end function nf_get_vars_real_d6

integer function nf_get_vars_real_d7(ncid,varid,start,ncount,stride,fp) result(ierr)
  implicit none
  integer, intent(in) :: ncid, varid
  integer, dimension(:), intent(in) :: start, ncount, stride
  real(kind=4), dimension(:,:,:,:,:,:,:), intent(out) :: fp
  integer (C_INT) :: c_ncid, c_varid
  integer (C_INT), target :: c_ndims
  integer (C_SIZE_T), dimension(size(start)) :: c_start
  integer (C_SIZE_T), dimension(size(ncount)) :: c_ncount
  integer (C_INTPTR_T), dimension(size(stride)) :: c_stride
  real (C_FLOAT), dimension(size(fp,1),size(fp,2),size(fp,3),size(fp,4),size(fp,5), &
      size(fp,6),size(fp,7)), target :: c_fp
  integer i, ndims
  c_ncid = ncid
  c_varid = varid - 1
  c_ndims = size(start)
  ierr = nc_inq_varndims(c_ncid,c_varid,C_LOC(c_ndims))
  ndims = c_ndims
  do i = 1,ndims
    c_start(ndims-i+1) = start(i) - 1
    c_ncount(ndims-i+1) = ncount(i)
    c_stride(ndims-i+1) = stride(i)
  end do
  ierr = nc_get_vars_float(c_ncid,c_varid,c_start,c_ncount,c_stride,C_LOC(c_fp))
  fp = c_fp
end function nf_get_vars_real_d7

integer function nf_get_vars_double_d1(ncid,varid,start,ncount,stride,dp) result(ierr)
  implicit none
  integer, intent(in) :: ncid, varid
  integer, dimension(:), intent(in) :: start, ncount, stride
  real(kind=8), dimension(:), intent(out) :: dp
  integer (C_INT) :: c_ncid, c_varid
  integer (C_INT), target :: c_ndims
  integer (C_SIZE_T), dimension(size(start)) :: c_start
  integer (C_SIZE_T), dimension(size(ncount)) :: c_ncount
  integer (C_INTPTR_T), dimension(size(stride)) :: c_stride
  real (C_DOUBLE), dimension(size(dp)), target :: c_dp
  integer i, ndims
  c_ncid = ncid
  c_varid = varid - 1
  c_ndims = size(start)
  ierr = nc_inq_varndims(c_ncid,c_varid,C_LOC(c_ndims))
  ndims = c_ndims
  do i = 1,ndims
    c_start(ndims-i+1) = start(i) - 1
    c_ncount(ndims-i+1) = ncount(i)
    c_stride(ndims-i+1) = stride(i)
  end do
  ierr = nc_get_vars_double(c_ncid,c_varid,c_start,c_ncount,c_stride,C_LOC(c_dp))
  dp = c_dp
end function nf_get_vars_double_d1

integer function nf_get_vars_double_d2(ncid,varid,start,ncount,stride,dp) result(ierr)
  implicit none
  integer, intent(in) :: ncid, varid
  integer, dimension(:), intent(in) :: start, ncount, stride
  real(kind=8), dimension(:,:), intent(out) :: dp
  integer (C_INT) :: c_ncid, c_varid
  integer (C_INT), target :: c_ndims
  integer (C_SIZE_T), dimension(size(start)) :: c_start
  integer (C_SIZE_T), dimension(size(ncount)) :: c_ncount
  integer (C_INTPTR_T), dimension(size(stride)) :: c_stride
  real (C_DOUBLE), dimension(size(dp,1),size(dp,2)), target :: c_dp
  integer i, ndims
  c_ncid = ncid
  c_varid = varid - 1
  c_ndims = size(start)
  ierr = nc_inq_varndims(c_ncid,c_varid,C_LOC(c_ndims))
  ndims = c_ndims
  do i = 1,ndims
    c_start(ndims-i+1) = start(i) - 1
    c_ncount(ndims-i+1) = ncount(i)
    c_stride(ndims-i+1) = stride(i)
  end do
  ierr = nc_get_vars_double(c_ncid,c_varid,c_start,c_ncount,c_stride,C_LOC(c_dp))
  dp = c_dp
end function nf_get_vars_double_d2

integer function nf_get_vars_double_d3(ncid,varid,start,ncount,stride,dp) result(ierr)
  implicit none
  integer, intent(in) :: ncid, varid
  integer, dimension(:), intent(in) :: start, ncount, stride
  real(kind=8), dimension(:,:,:), intent(out) :: dp
  integer (C_INT) :: c_ncid, c_varid
  integer (C_INT), target :: c_ndims
  integer (C_SIZE_T), dimension(size(start)) :: c_start
  integer (C_SIZE_T), dimension(size(ncount)) :: c_ncount
  integer (C_INTPTR_T), dimension(size(stride)) :: c_stride
  real (C_DOUBLE), dimension(size(dp,1),size(dp,2),size(dp,3)), target :: c_dp
  integer i, ndims
  c_ncid = ncid
  c_varid = varid - 1
  c_ndims = size(start)
  ierr = nc_inq_varndims(c_ncid,c_varid,C_LOC(c_ndims))
  ndims = c_ndims
  do i = 1,ndims
    c_start(ndims-i+1) = start(i) - 1
    c_ncount(ndims-i+1) = ncount(i)
    c_stride(ndims-i+1) = stride(i)
  end do
  ierr = nc_get_vars_double(c_ncid,c_varid,c_start,c_ncount,c_stride,C_LOC(c_dp))
  dp = c_dp
end function nf_get_vars_double_d3

integer function nf_get_vars_double_d4(ncid,varid,start,ncount,stride,dp) result(ierr)
  implicit none
  integer, intent(in) :: ncid, varid
  integer, dimension(:), intent(in) :: start, ncount, stride
  real(kind=8), dimension(:,:,:,:), intent(out) :: dp
  integer (C_INT) :: c_ncid, c_varid
  integer (C_INT), target :: c_ndims
  integer (C_SIZE_T), dimension(size(start)) :: c_start
  integer (C_SIZE_T), dimension(size(ncount)) :: c_ncount
  integer (C_INTPTR_T), dimension(size(stride)) :: c_stride
  real (C_DOUBLE), dimension(size(dp,1),size(dp,2),size(dp,3),size(dp,4)), target :: c_dp
  integer i, ndims
  c_ncid = ncid
  c_varid = varid - 1
  c_ndims = size(start)
  ierr = nc_inq_varndims(c_ncid,c_varid,C_LOC(c_ndims))
  ndims = c_ndims
  do i = 1,ndims
    c_start(ndims-i+1) = start(i) - 1
    c_ncount(ndims-i+1) = ncount(i)
    c_stride(ndims-i+1) = stride(i)
  end do
  ierr = nc_get_vars_double(c_ncid,c_varid,c_start,c_ncount,c_stride,C_LOC(c_dp))
  dp = c_dp
end function nf_get_vars_double_d4

integer function nf_get_vars_double_d5(ncid,varid,start,ncount,stride,dp) result(ierr)
  implicit none
  integer, intent(in) :: ncid, varid
  integer, dimension(:), intent(in) :: start, ncount, stride
  real(kind=8), dimension(:,:,:,:,:), intent(out) :: dp
  integer (C_INT) :: c_ncid, c_varid
  integer (C_INT), target :: c_ndims
  integer (C_SIZE_T), dimension(size(start)) :: c_start
  integer (C_SIZE_T), dimension(size(ncount)) :: c_ncount
  integer (C_INTPTR_T), dimension(size(stride)) :: c_stride
  real (C_DOUBLE), dimension(size(dp,1),size(dp,2),size(dp,3),size(dp,4),size(dp,5)), &
      target :: c_dp
  integer i, ndims
  c_ncid = ncid
  c_varid = varid - 1
  c_ndims = size(start)
  ierr = nc_inq_varndims(c_ncid,c_varid,C_LOC(c_ndims))
  ndims = c_ndims
  do i = 1,ndims
    c_start(ndims-i+1) = start(i) - 1
    c_ncount(ndims-i+1) = ncount(i)
    c_stride(ndims-i+1) = stride(i)
  end do
  ierr = nc_get_vars_double(c_ncid,c_varid,c_start,c_ncount,c_stride,C_LOC(c_dp))
  dp = c_dp
end function nf_get_vars_double_d5

integer function nf_get_vars_double_d6(ncid,varid,start,ncount,stride,dp) result(ierr)
  implicit none
  integer, intent(in) :: ncid, varid
  integer, dimension(:), intent(in) :: start, ncount, stride
  real(kind=8), dimension(:,:,:,:,:,:), intent(out) :: dp
  integer (C_INT) :: c_ncid, c_varid
  integer (C_INT), target :: c_ndims
  integer (C_SIZE_T), dimension(size(start)) :: c_start
  integer (C_SIZE_T), dimension(size(ncount)) :: c_ncount
  integer (C_INTPTR_T), dimension(size(stride)) :: c_stride
  real (C_DOUBLE), dimension(size(dp,1),size(dp,2),size(dp,3),size(dp,4),size(dp,5), &
      size(dp,6)), target :: c_dp
  integer i, ndims
  c_ncid = ncid
  c_varid = varid - 1
  c_ndims = size(start)
  ierr = nc_inq_varndims(c_ncid,c_varid,C_LOC(c_ndims))
  ndims = c_ndims
  do i = 1,ndims
    c_start(ndims-i+1) = start(i) - 1
    c_ncount(ndims-i+1) = ncount(i)
    c_stride(ndims-i+1) = stride(i)
  end do
  ierr = nc_get_vars_double(c_ncid,c_varid,c_start,c_ncount,c_stride,C_LOC(c_dp))
  dp = c_dp
end function nf_get_vars_double_d6

integer function nf_get_vars_double_d7(ncid,varid,start,ncount,stride,dp) result(ierr)
  implicit none
  integer, intent(in) :: ncid, varid
  integer, dimension(:), intent(in) :: start, ncount, stride
  real(kind=8), dimension(:,:,:,:,:,:,:), intent(out) :: dp
  integer (C_INT) :: c_ncid, c_varid
  integer (C_INT), target :: c_ndims
  integer (C_SIZE_T), dimension(size(start)) :: c_start
  integer (C_SIZE_T), dimension(size(ncount)) :: c_ncount
  integer (C_INTPTR_T), dimension(size(stride)) :: c_stride
  real (C_DOUBLE), dimension(size(dp,1),size(dp,2),size(dp,3),size(dp,4),size(dp,5), &
      size(dp,6),size(dp,7)), target :: c_dp
  integer i, ndims
  c_ncid = ncid
  c_varid = varid - 1
  c_ndims = size(start)
  ierr = nc_inq_varndims(c_ncid,c_varid,C_LOC(c_ndims))
  ndims = c_ndims
  do i = 1,ndims
    c_start(ndims-i+1) = start(i) - 1
    c_ncount(ndims-i+1) = ncount(i)
    c_stride(ndims-i+1) = stride(i)
  end do
  ierr = nc_get_vars_double(c_ncid,c_varid,c_start,c_ncount,c_stride,C_LOC(c_dp))
  dp = c_dp
end function nf_get_vars_double_d7

integer function nf_get_varm_int1_d1(ncid,varid,start,ncount,stride,imap,bp) result(ierr)
  implicit none
  integer, intent(in) :: ncid, varid
  integer, dimension(:), intent(in) :: start, ncount, stride, imap
  integer(kind=1), dimension(:), intent(out) :: bp
  integer (C_INT) :: c_ncid, c_varid
  integer (C_INT), target :: c_ndims
  integer (C_SIZE_T), dimension(size(start)) :: c_start
  integer (C_SIZE_T), dimension(size(ncount)) :: c_ncount
  integer (C_INTPTR_T), dimension(size(stride)) :: c_stride
  integer (C_INTPTR_T), dimension(size(imap)) :: c_imap
  integer (C_SIGNED_CHAR), dimension(size(bp)), target :: c_bp
  integer i, ndims
  c_ncid = ncid
  c_varid = varid - 1
  c_ndims = size(start)
  ierr = nc_inq_varndims(c_ncid,c_varid,C_LOC(c_ndims))
  ndims = c_ndims
  do i = 1,ndims
    c_start(ndims-i+1) = start(i) - 1
    c_ncount(ndims-i+1) = ncount(i)
    c_stride(ndims-i+1) = stride(i)
    c_imap(ndims-i+1) = imap(i)
  end do
  ierr = nc_get_varm_schar(c_ncid,c_varid,c_start,c_ncount,c_stride,c_imap,C_LOC(c_bp))
  bp = c_bp
end function nf_get_varm_int1_d1

integer function nf_get_varm_int1_d2(ncid,varid,start,ncount,stride,imap,bp) result(ierr)
  implicit none
  integer, intent(in) :: ncid, varid
  integer, dimension(:), intent(in) :: start, ncount, stride, imap
  integer(kind=1), dimension(:,:), intent(out) :: bp
  integer (C_INT) :: c_ncid, c_varid
  integer (C_INT), target :: c_ndims
  integer (C_SIZE_T), dimension(size(start)) :: c_start
  integer (C_SIZE_T), dimension(size(ncount)) :: c_ncount
  integer (C_INTPTR_T), dimension(size(stride)) :: c_stride
  integer (C_INTPTR_T), dimension(size(imap)) :: c_imap
  integer (C_SIGNED_CHAR), dimension(size(bp,1),size(bp,2)), target :: c_bp
  integer i, ndims
  c_ncid = ncid
  c_varid = varid - 1
  c_ndims = size(start)
  ierr = nc_inq_varndims(c_ncid,c_varid,C_LOC(c_ndims))
  ndims = c_ndims
  do i = 1,ndims
    c_start(ndims-i+1) = start(i) - 1
    c_ncount(ndims-i+1) = ncount(i)
    c_stride(ndims-i+1) = stride(i)
    c_imap(ndims-i+1) = imap(i)
  end do
  ierr = nc_get_varm_schar(c_ncid,c_varid,c_start,c_ncount,c_stride,c_imap,C_LOC(c_bp))
  bp = c_bp
end function nf_get_varm_int1_d2

integer function nf_get_varm_int1_d3(ncid,varid,start,ncount,stride,imap,bp) result(ierr)
  implicit none
  integer, intent(in) :: ncid, varid
  integer, dimension(:), intent(in) :: start, ncount, stride, imap
  integer(kind=1), dimension(:,:,:), intent(out) :: bp
  integer (C_INT) :: c_ncid, c_varid
  integer (C_INT), target :: c_ndims
  integer (C_SIZE_T), dimension(size(start)) :: c_start
  integer (C_SIZE_T), dimension(size(ncount)) :: c_ncount
  integer (C_INTPTR_T), dimension(size(stride)) :: c_stride
  integer (C_INTPTR_T), dimension(size(imap)) :: c_imap
  integer (C_SIGNED_CHAR), dimension(size(bp,1),size(bp,2),size(bp,3)), target :: c_bp
  integer i, ndims
  c_ncid = ncid
  c_varid = varid - 1
  c_ndims = size(start)
  ierr = nc_inq_varndims(c_ncid,c_varid,C_LOC(c_ndims))
  ndims = c_ndims
  do i = 1,ndims
    c_start(ndims-i+1) = start(i) - 1
    c_ncount(ndims-i+1) = ncount(i)
    c_stride(ndims-i+1) = stride(i)
    c_imap(ndims-i+1) = imap(i)
  end do
  ierr = nc_get_varm_schar(c_ncid,c_varid,c_start,c_ncount,c_stride,c_imap,C_LOC(c_bp))
  bp = c_bp
end function nf_get_varm_int1_d3

integer function nf_get_varm_int1_d4(ncid,varid,start,ncount,stride,imap,bp) result(ierr)
  implicit none
  integer, intent(in) :: ncid, varid
  integer, dimension(:), intent(in) :: start, ncount, stride, imap
  integer(kind=1), dimension(:,:,:,:), intent(out) :: bp
  integer (C_INT) :: c_ncid, c_varid
  integer (C_INT), target :: c_ndims
  integer (C_SIZE_T), dimension(size(start)) :: c_start
  integer (C_SIZE_T), dimension(size(ncount)) :: c_ncount
  integer (C_INTPTR_T), dimension(size(stride)) :: c_stride
  integer (C_INTPTR_T), dimension(size(imap)) :: c_imap
  integer (C_SIGNED_CHAR), dimension(size(bp,1),size(bp,2),size(bp,3),size(bp,4)), target :: c_bp
  integer i, ndims
  c_ncid = ncid
  c_varid = varid - 1
  c_ndims = size(start)
  ierr = nc_inq_varndims(c_ncid,c_varid,C_LOC(c_ndims))
  ndims = c_ndims
  do i = 1,ndims
    c_start(ndims-i+1) = start(i) - 1
    c_ncount(ndims-i+1) = ncount(i)
    c_stride(ndims-i+1) = stride(i)
    c_imap(ndims-i+1) = imap(i)
  end do
  ierr = nc_get_varm_schar(c_ncid,c_varid,c_start,c_ncount,c_stride,c_imap,C_LOC(c_bp))
  bp = c_bp
end function nf_get_varm_int1_d4

integer function nf_get_varm_int1_d5(ncid,varid,start,ncount,stride,imap,bp) result(ierr)
  implicit none
  integer, intent(in) :: ncid, varid
  integer, dimension(:), intent(in) :: start, ncount, stride, imap
  integer(kind=1), dimension(:,:,:,:,:), intent(out) :: bp
  integer (C_INT) :: c_ncid, c_varid
  integer (C_INT), target :: c_ndims
  integer (C_SIZE_T), dimension(size(start)) :: c_start
  integer (C_SIZE_T), dimension(size(ncount)) :: c_ncount
  integer (C_INTPTR_T), dimension(size(stride)) :: c_stride
  integer (C_INTPTR_T), dimension(size(imap)) :: c_imap
  integer (C_SIGNED_CHAR), dimension(size(bp,1),size(bp,2),size(bp,3),size(bp,4),size(bp,5)), &
      target :: c_bp
  integer i, ndims
  c_ncid = ncid
  c_varid = varid - 1
  c_ndims = size(start)
  ierr = nc_inq_varndims(c_ncid,c_varid,C_LOC(c_ndims))
  ndims = c_ndims
  do i = 1,ndims
    c_start(ndims-i+1) = start(i) - 1
    c_ncount(ndims-i+1) = ncount(i)
    c_stride(ndims-i+1) = stride(i)
    c_imap(ndims-i+1) = imap(i)
  end do
  ierr = nc_get_varm_schar(c_ncid,c_varid,c_start,c_ncount,c_stride,c_imap,C_LOC(c_bp))
  bp = c_bp
end function nf_get_varm_int1_d5

integer function nf_get_varm_int1_d6(ncid,varid,start,ncount,stride,imap,bp) result(ierr)
  implicit none
  integer, intent(in) :: ncid, varid
  integer, dimension(:), intent(in) :: start, ncount, stride, imap
  integer(kind=1), dimension(:,:,:,:,:,:), intent(out) :: bp
  integer (C_INT) :: c_ncid, c_varid
  integer (C_INT), target :: c_ndims
  integer (C_SIZE_T), dimension(size(start)) :: c_start
  integer (C_SIZE_T), dimension(size(ncount)) :: c_ncount
  integer (C_INTPTR_T), dimension(size(stride)) :: c_stride
  integer (C_INTPTR_T), dimension(size(imap)) :: c_imap
  integer (C_SIGNED_CHAR), dimension(size(bp,1),size(bp,2),size(bp,3),size(bp,4),size(bp,5), &
      size(bp,6)), target :: c_bp
  integer i, ndims
  c_ncid = ncid
  c_varid = varid - 1
  c_ndims = size(start)
  ierr = nc_inq_varndims(c_ncid,c_varid,C_LOC(c_ndims))
  ndims = c_ndims
  do i = 1,ndims
    c_start(ndims-i+1) = start(i) - 1
    c_ncount(ndims-i+1) = ncount(i)
    c_stride(ndims-i+1) = stride(i)
    c_imap(ndims-i+1) = imap(i)
  end do
  ierr = nc_get_varm_schar(c_ncid,c_varid,c_start,c_ncount,c_stride,c_imap,C_LOC(c_bp))
  bp = c_bp
end function nf_get_varm_int1_d6

integer function nf_get_varm_int1_d7(ncid,varid,start,ncount,stride,imap,bp) result(ierr)
  implicit none
  integer, intent(in) :: ncid, varid
  integer, dimension(:), intent(in) :: start, ncount, stride, imap
  integer(kind=1), dimension(:,:,:,:,:,:,:), intent(out) :: bp
  integer (C_INT) :: c_ncid, c_varid
  integer (C_INT), target :: c_ndims
  integer (C_SIZE_T), dimension(size(start)) :: c_start
  integer (C_SIZE_T), dimension(size(ncount)) :: c_ncount
  integer (C_INTPTR_T), dimension(size(stride)) :: c_stride
  integer (C_INTPTR_T), dimension(size(imap)) :: c_imap
  integer (C_SIGNED_CHAR), dimension(size(bp,1),size(bp,2),size(bp,3),size(bp,4),size(bp,5), &
      size(bp,6),size(bp,7)), target :: c_bp
  integer i, ndims
  c_ncid = ncid
  c_varid = varid - 1
  c_ndims = size(start)
  ierr = nc_inq_varndims(c_ncid,c_varid,C_LOC(c_ndims))
  ndims = c_ndims
  do i = 1,ndims
    c_start(ndims-i+1) = start(i) - 1
    c_ncount(ndims-i+1) = ncount(i)
    c_stride(ndims-i+1) = stride(i)
    c_imap(ndims-i+1) = imap(i)
  end do
  ierr = nc_get_varm_schar(c_ncid,c_varid,c_start,c_ncount,c_stride,c_imap,C_LOC(c_bp))
  bp = c_bp
end function nf_get_varm_int1_d7

integer function nf_get_varm_int2_d1(ncid,varid,start,ncount,stride,imap,sp) result(ierr)
  implicit none
  integer, intent(in) :: ncid, varid
  integer, dimension(:), intent(in) :: start, ncount, stride, imap
  integer(kind=2), dimension(:), intent(out) :: sp
  integer (C_INT) :: c_ncid, c_varid
  integer (C_INT), target :: c_ndims
  integer (C_SIZE_T), dimension(size(start)) :: c_start
  integer (C_SIZE_T), dimension(size(ncount)) :: c_ncount
  integer (C_INTPTR_T), dimension(size(stride)) :: c_stride
  integer (C_INTPTR_T), dimension(size(imap)) :: c_imap
  integer (C_SHORT), dimension(size(sp)), target :: c_sp
  integer i, ndims
  c_ncid = ncid
  c_varid = varid - 1
  c_ndims = size(start)
  ierr = nc_inq_varndims(c_ncid,c_varid,C_LOC(c_ndims))
  ndims = c_ndims
  do i = 1,ndims
    c_start(ndims-i+1) = start(i) - 1
    c_ncount(ndims-i+1) = ncount(i)
    c_stride(ndims-i+1) = stride(i)
    c_imap(ndims-i+1) = imap(i)
  end do
  ierr = nc_get_varm_short(c_ncid,c_varid,c_start,c_ncount,c_stride,c_imap,C_LOC(c_sp))
  sp = c_sp
end function nf_get_varm_int2_d1

integer function nf_get_varm_int2_d2(ncid,varid,start,ncount,stride,imap,sp) result(ierr)
  implicit none
  integer, intent(in) :: ncid, varid
  integer, dimension(:), intent(in) :: start, ncount, stride, imap
  integer(kind=2), dimension(:,:), intent(out) :: sp
  integer (C_INT) :: c_ncid, c_varid
  integer (C_INT), target :: c_ndims
  integer (C_SIZE_T), dimension(size(start)) :: c_start
  integer (C_SIZE_T), dimension(size(ncount)) :: c_ncount
  integer (C_INTPTR_T), dimension(size(stride)) :: c_stride
  integer (C_INTPTR_T), dimension(size(imap)) :: c_imap
  integer (C_SHORT), dimension(size(sp,1),size(sp,2)), target :: c_sp
  integer i, ndims
  c_ncid = ncid
  c_varid = varid - 1
  c_ndims = size(start)
  ierr = nc_inq_varndims(c_ncid,c_varid,C_LOC(c_ndims))
  ndims = c_ndims
  do i = 1,ndims
    c_start(ndims-i+1) = start(i) - 1
    c_ncount(ndims-i+1) = ncount(i)
    c_stride(ndims-i+1) = stride(i)
    c_imap(ndims-i+1) = imap(i)
  end do
  ierr = nc_get_varm_short(c_ncid,c_varid,c_start,c_ncount,c_stride,c_imap,C_LOC(c_sp))
  sp = c_sp
end function nf_get_varm_int2_d2

integer function nf_get_varm_int2_d3(ncid,varid,start,ncount,stride,imap,sp) result(ierr)
  implicit none
  integer, intent(in) :: ncid, varid
  integer, dimension(:), intent(in) :: start, ncount, stride, imap
  integer(kind=2), dimension(:,:,:), intent(out) :: sp
  integer (C_INT) :: c_ncid, c_varid
  integer (C_INT), target :: c_ndims
  integer (C_SIZE_T), dimension(size(start)) :: c_start
  integer (C_SIZE_T), dimension(size(ncount)) :: c_ncount
  integer (C_INTPTR_T), dimension(size(stride)) :: c_stride
  integer (C_INTPTR_T), dimension(size(imap)) :: c_imap
  integer (C_SHORT), dimension(size(sp,1),size(sp,2),size(sp,3)), target :: c_sp
  integer i, ndims
  c_ncid = ncid
  c_varid = varid - 1
  c_ndims = size(start)
  ierr = nc_inq_varndims(c_ncid,c_varid,C_LOC(c_ndims))
  ndims = c_ndims
  do i = 1,ndims
    c_start(ndims-i+1) = start(i) - 1
    c_ncount(ndims-i+1) = ncount(i)
    c_stride(ndims-i+1) = stride(i)
    c_imap(ndims-i+1) = imap(i)
  end do
  ierr = nc_get_varm_short(c_ncid,c_varid,c_start,c_ncount,c_stride,c_imap,C_LOC(c_sp))
  sp = c_sp
end function nf_get_varm_int2_d3

integer function nf_get_varm_int2_d4(ncid,varid,start,ncount,stride,imap,sp) result(ierr)
  implicit none
  integer, intent(in) :: ncid, varid
  integer, dimension(:), intent(in) :: start, ncount, stride, imap
  integer(kind=2), dimension(:,:,:,:), intent(out) :: sp
  integer (C_INT) :: c_ncid, c_varid
  integer (C_INT), target :: c_ndims
  integer (C_SIZE_T), dimension(size(start)) :: c_start
  integer (C_SIZE_T), dimension(size(ncount)) :: c_ncount
  integer (C_INTPTR_T), dimension(size(stride)) :: c_stride
  integer (C_INTPTR_T), dimension(size(imap)) :: c_imap
  integer (C_SHORT), dimension(size(sp,1),size(sp,2),size(sp,3),size(sp,4)), target :: c_sp
  integer i, ndims
  c_ncid = ncid
  c_varid = varid - 1
  c_ndims = size(start)
  ierr = nc_inq_varndims(c_ncid,c_varid,C_LOC(c_ndims))
  ndims = c_ndims
  do i = 1,ndims
    c_start(ndims-i+1) = start(i) - 1
    c_ncount(ndims-i+1) = ncount(i)
    c_stride(ndims-i+1) = stride(i)
    c_imap(ndims-i+1) = imap(i)
  end do
  ierr = nc_get_varm_short(c_ncid,c_varid,c_start,c_ncount,c_stride,c_imap,C_LOC(c_sp))
  sp = c_sp
end function nf_get_varm_int2_d4

integer function nf_get_varm_int2_d5(ncid,varid,start,ncount,stride,imap,sp) result(ierr)
  implicit none
  integer, intent(in) :: ncid, varid
  integer, dimension(:), intent(in) :: start, ncount, stride, imap
  integer(kind=2), dimension(:,:,:,:,:), intent(out) :: sp
  integer (C_INT) :: c_ncid, c_varid
  integer (C_INT), target :: c_ndims
  integer (C_SIZE_T), dimension(size(start)) :: c_start
  integer (C_SIZE_T), dimension(size(ncount)) :: c_ncount
  integer (C_INTPTR_T), dimension(size(stride)) :: c_stride
  integer (C_INTPTR_T), dimension(size(imap)) :: c_imap
  integer (C_SHORT), dimension(size(sp,1),size(sp,2),size(sp,3),size(sp,4),size(sp,5)), &
      target :: c_sp
  integer i, ndims
  c_ncid = ncid
  c_varid = varid - 1
  c_ndims = size(start)
  ierr = nc_inq_varndims(c_ncid,c_varid,C_LOC(c_ndims))
  ndims = c_ndims
  do i = 1,ndims
    c_start(ndims-i+1) = start(i) - 1
    c_ncount(ndims-i+1) = ncount(i)
    c_stride(ndims-i+1) = stride(i)
    c_imap(ndims-i+1) = imap(i)
  end do
  ierr = nc_get_varm_short(c_ncid,c_varid,c_start,c_ncount,c_stride,c_imap,C_LOC(c_sp))
  sp = c_sp
end function nf_get_varm_int2_d5

integer function nf_get_varm_int2_d6(ncid,varid,start,ncount,stride,imap,sp) result(ierr)
  implicit none
  integer, intent(in) :: ncid, varid
  integer, dimension(:), intent(in) :: start, ncount, stride, imap
  integer(kind=2), dimension(:,:,:,:,:,:), intent(out) :: sp
  integer (C_INT) :: c_ncid, c_varid
  integer (C_INT), target :: c_ndims
  integer (C_SIZE_T), dimension(size(start)) :: c_start
  integer (C_SIZE_T), dimension(size(ncount)) :: c_ncount
  integer (C_INTPTR_T), dimension(size(stride)) :: c_stride
  integer (C_INTPTR_T), dimension(size(imap)) :: c_imap
  integer (C_SHORT), dimension(size(sp,1),size(sp,2),size(sp,3),size(sp,4),size(sp,5), &
      size(sp,6)), target :: c_sp
  integer i, ndims
  c_ncid = ncid
  c_varid = varid - 1
  c_ndims = size(start)
  ierr = nc_inq_varndims(c_ncid,c_varid,C_LOC(c_ndims))
  ndims = c_ndims
  do i = 1,ndims
    c_start(ndims-i+1) = start(i) - 1
    c_ncount(ndims-i+1) = ncount(i)
    c_stride(ndims-i+1) = stride(i)
    c_imap(ndims-i+1) = imap(i)
  end do
  ierr = nc_get_varm_short(c_ncid,c_varid,c_start,c_ncount,c_stride,c_imap,C_LOC(c_sp))
  sp = c_sp
end function nf_get_varm_int2_d6

integer function nf_get_varm_int2_d7(ncid,varid,start,ncount,stride,imap,sp) result(ierr)
  implicit none
  integer, intent(in) :: ncid, varid
  integer, dimension(:), intent(in) :: start, ncount, stride, imap
  integer(kind=2), dimension(:,:,:,:,:,:,:), intent(out) :: sp
  integer (C_INT) :: c_ncid, c_varid
  integer (C_INT), target :: c_ndims
  integer (C_SIZE_T), dimension(size(start)) :: c_start
  integer (C_SIZE_T), dimension(size(ncount)) :: c_ncount
  integer (C_INTPTR_T), dimension(size(stride)) :: c_stride
  integer (C_INTPTR_T), dimension(size(imap)) :: c_imap
  integer (C_SHORT), dimension(size(sp,1),size(sp,2),size(sp,3),size(sp,4),size(sp,5), &
      size(sp,6),size(sp,7)), target :: c_sp
  integer i, ndims
  c_ncid = ncid
  c_varid = varid - 1
  c_ndims = size(start)
  ierr = nc_inq_varndims(c_ncid,c_varid,C_LOC(c_ndims))
  ndims = c_ndims
  do i = 1,ndims
    c_start(ndims-i+1) = start(i) - 1
    c_ncount(ndims-i+1) = ncount(i)
    c_stride(ndims-i+1) = stride(i)
    c_imap(ndims-i+1) = imap(i)
  end do
  ierr = nc_get_varm_short(c_ncid,c_varid,c_start,c_ncount,c_stride,c_imap,C_LOC(c_sp))
  sp = c_sp
end function nf_get_varm_int2_d7

integer function nf_get_varm_int_d1(ncid,varid,start,ncount,stride,imap,ip) result(ierr)
  implicit none
  integer, intent(in) :: ncid, varid
  integer, dimension(:), intent(in) :: start, ncount, stride, imap
  integer(kind=4), dimension(:), intent(out) :: ip
  integer (C_INT) :: c_ncid, c_varid
  integer (C_INT), target :: c_ndims
  integer (C_SIZE_T), dimension(size(start)) :: c_start
  integer (C_SIZE_T), dimension(size(ncount)) :: c_ncount
  integer (C_INTPTR_T), dimension(size(stride)) :: c_stride
  integer (C_INTPTR_T), dimension(size(imap)) :: c_imap
  integer (C_INT), dimension(size(ip)), target :: c_ip
  integer i, ndims
  c_ncid = ncid
  c_varid = varid - 1
  c_ndims = size(start)
  ierr = nc_inq_varndims(c_ncid,c_varid,C_LOC(c_ndims))
  ndims = c_ndims
  do i = 1,ndims
    c_start(ndims-i+1) = start(i) - 1
    c_ncount(ndims-i+1) = ncount(i)
    c_stride(ndims-i+1) = stride(i)
    c_imap(ndims-i+1) = imap(i)
  end do
  ierr = nc_get_varm_int(c_ncid,c_varid,c_start,c_ncount,c_stride,c_imap,C_LOC(c_ip))
  ip = c_ip
end function nf_get_varm_int_d1

integer function nf_get_varm_int_d2(ncid,varid,start,ncount,stride,imap,ip) result(ierr)
  implicit none
  integer, intent(in) :: ncid, varid
  integer, dimension(:), intent(in) :: start, ncount, stride, imap
  integer(kind=4), dimension(:,:), intent(out) :: ip
  integer (C_INT) :: c_ncid, c_varid
  integer (C_INT), target :: c_ndims
  integer (C_SIZE_T), dimension(size(start)) :: c_start
  integer (C_SIZE_T), dimension(size(ncount)) :: c_ncount
  integer (C_INTPTR_T), dimension(size(stride)) :: c_stride
  integer (C_INTPTR_T), dimension(size(imap)) :: c_imap
  integer (C_INT), dimension(size(ip,1),size(ip,2)), target :: c_ip
  integer i, ndims
  c_ncid = ncid
  c_varid = varid - 1
  c_ndims = size(start)
  ierr = nc_inq_varndims(c_ncid,c_varid,C_LOC(c_ndims))
  ndims = c_ndims
  do i = 1,ndims
    c_start(ndims-i+1) = start(i) - 1
    c_ncount(ndims-i+1) = ncount(i)
    c_stride(ndims-i+1) = stride(i)
    c_imap(ndims-i+1) = imap(i)
  end do
  ierr = nc_get_varm_int(c_ncid,c_varid,c_start,c_ncount,c_stride,c_imap,C_LOC(c_ip))
  ip = c_ip
end function nf_get_varm_int_d2

integer function nf_get_varm_int_d3(ncid,varid,start,ncount,stride,imap,ip) result(ierr)
  implicit none
  integer, intent(in) :: ncid, varid
  integer, dimension(:), intent(in) :: start, ncount, stride, imap
  integer(kind=4), dimension(:,:,:), intent(out) :: ip
  integer (C_INT) :: c_ncid, c_varid
  integer (C_INT), target :: c_ndims
  integer (C_SIZE_T), dimension(size(start)) :: c_start
  integer (C_SIZE_T), dimension(size(ncount)) :: c_ncount
  integer (C_INTPTR_T), dimension(size(stride)) :: c_stride
  integer (C_INTPTR_T), dimension(size(imap)) :: c_imap
  integer (C_INT), dimension(size(ip,1),size(ip,2),size(ip,3)), target :: c_ip
  integer i, ndims
  c_ncid = ncid
  c_varid = varid - 1
  c_ndims = size(start)
  ierr = nc_inq_varndims(c_ncid,c_varid,C_LOC(c_ndims))
  ndims = c_ndims
  do i = 1,ndims
    c_start(ndims-i+1) = start(i) - 1
    c_ncount(ndims-i+1) = ncount(i)
    c_stride(ndims-i+1) = stride(i)
    c_imap(ndims-i+1) = imap(i)
  end do
  ierr = nc_get_varm_int(c_ncid,c_varid,c_start,c_ncount,c_stride,c_imap,C_LOC(c_ip))
  ip = c_ip
end function nf_get_varm_int_d3

integer function nf_get_varm_int_d4(ncid,varid,start,ncount,stride,imap,ip) result(ierr)
  implicit none
  integer, intent(in) :: ncid, varid
  integer, dimension(:), intent(in) :: start, ncount, stride, imap
  integer(kind=4), dimension(:,:,:,:), intent(out) :: ip
  integer (C_INT) :: c_ncid, c_varid
  integer (C_INT), target :: c_ndims
  integer (C_SIZE_T), dimension(size(start)) :: c_start
  integer (C_SIZE_T), dimension(size(ncount)) :: c_ncount
  integer (C_INTPTR_T), dimension(size(stride)) :: c_stride
  integer (C_INTPTR_T), dimension(size(imap)) :: c_imap
  integer (C_INT), dimension(size(ip,1),size(ip,2),size(ip,3),size(ip,4)), target :: c_ip
  integer i, ndims
  c_ncid = ncid
  c_varid = varid - 1
  c_ndims = size(start)
  ierr = nc_inq_varndims(c_ncid,c_varid,C_LOC(c_ndims))
  ndims = c_ndims
  do i = 1,ndims
    c_start(ndims-i+1) = start(i) - 1
    c_ncount(ndims-i+1) = ncount(i)
    c_stride(ndims-i+1) = stride(i)
    c_imap(ndims-i+1) = imap(i)
  end do
  ierr = nc_get_varm_int(c_ncid,c_varid,c_start,c_ncount,c_stride,c_imap,C_LOC(c_ip))
  ip = c_ip
end function nf_get_varm_int_d4

integer function nf_get_varm_int_d5(ncid,varid,start,ncount,stride,imap,ip) result(ierr)
  implicit none
  integer, intent(in) :: ncid, varid
  integer, dimension(:), intent(in) :: start, ncount, stride, imap
  integer(kind=4), dimension(:,:,:,:,:), intent(out) :: ip
  integer (C_INT) :: c_ncid, c_varid
  integer (C_INT), target :: c_ndims
  integer (C_SIZE_T), dimension(size(start)) :: c_start
  integer (C_SIZE_T), dimension(size(ncount)) :: c_ncount
  integer (C_INTPTR_T), dimension(size(stride)) :: c_stride
  integer (C_INTPTR_T), dimension(size(imap)) :: c_imap
  integer (C_INT), dimension(size(ip,1),size(ip,2),size(ip,3),size(ip,4),size(ip,5)), &
      target :: c_ip
  integer i, ndims
  c_ncid = ncid
  c_varid = varid - 1
  c_ndims = size(start)
  ierr = nc_inq_varndims(c_ncid,c_varid,C_LOC(c_ndims))
  ndims = c_ndims
  do i = 1,ndims
    c_start(ndims-i+1) = start(i) - 1
    c_ncount(ndims-i+1) = ncount(i)
    c_stride(ndims-i+1) = stride(i)
    c_imap(ndims-i+1) = imap(i)
  end do
  ierr = nc_get_varm_int(c_ncid,c_varid,c_start,c_ncount,c_stride,c_imap,C_LOC(c_ip))
  ip = c_ip
end function nf_get_varm_int_d5

integer function nf_get_varm_int_d6(ncid,varid,start,ncount,stride,imap,ip) result(ierr)
  implicit none
  integer, intent(in) :: ncid, varid
  integer, dimension(:), intent(in) :: start, ncount, stride, imap
  integer(kind=4), dimension(:,:,:,:,:,:), intent(out) :: ip
  integer (C_INT) :: c_ncid, c_varid
  integer (C_INT), target :: c_ndims
  integer (C_SIZE_T), dimension(size(start)) :: c_start
  integer (C_SIZE_T), dimension(size(ncount)) :: c_ncount
  integer (C_INTPTR_T), dimension(size(stride)) :: c_stride
  integer (C_INTPTR_T), dimension(size(imap)) :: c_imap
  integer (C_INT), dimension(size(ip,1),size(ip,2),size(ip,3),size(ip,4),size(ip,5), &
      size(ip,6)), target :: c_ip
  integer i, ndims
  c_ncid = ncid
  c_varid = varid - 1
  c_ndims = size(start)
  ierr = nc_inq_varndims(c_ncid,c_varid,C_LOC(c_ndims))
  ndims = c_ndims
  do i = 1,ndims
    c_start(ndims-i+1) = start(i) - 1
    c_ncount(ndims-i+1) = ncount(i)
    c_stride(ndims-i+1) = stride(i)
    c_imap(ndims-i+1) = imap(i)
  end do
  ierr = nc_get_varm_int(c_ncid,c_varid,c_start,c_ncount,c_stride,c_imap,C_LOC(c_ip))
  ip = c_ip
end function nf_get_varm_int_d6

integer function nf_get_varm_int_d7(ncid,varid,start,ncount,stride,imap,ip) result(ierr)
  implicit none
  integer, intent(in) :: ncid, varid
  integer, dimension(:), intent(in) :: start, ncount, stride, imap
  integer(kind=4), dimension(:,:,:,:,:,:,:), intent(out) :: ip
  integer (C_INT) :: c_ncid, c_varid
  integer (C_INT), target :: c_ndims
  integer (C_SIZE_T), dimension(size(start)) :: c_start
  integer (C_SIZE_T), dimension(size(ncount)) :: c_ncount
  integer (C_INTPTR_T), dimension(size(stride)) :: c_stride
  integer (C_INTPTR_T), dimension(size(imap)) :: c_imap
  integer (C_INT), dimension(size(ip,1),size(ip,2),size(ip,3),size(ip,4),size(ip,5), &
      size(ip,6),size(ip,7)), target :: c_ip
  integer i, ndims
  c_ncid = ncid
  c_varid = varid - 1
  c_ndims = size(start)
  ierr = nc_inq_varndims(c_ncid,c_varid,C_LOC(c_ndims))
  ndims = c_ndims
  do i = 1,ndims
    c_start(ndims-i+1) = start(i) - 1
    c_ncount(ndims-i+1) = ncount(i)
    c_stride(ndims-i+1) = stride(i)
    c_imap(ndims-i+1) = imap(i)
  end do
  ierr = nc_get_varm_int(c_ncid,c_varid,c_start,c_ncount,c_stride,c_imap,C_LOC(c_ip))
  ip = c_ip
end function nf_get_varm_int_d7

integer function nf_get_varm_real_d1(ncid,varid,start,ncount,stride,imap,fp) result(ierr)
  implicit none
  integer, intent(in) :: ncid, varid
  integer, dimension(:), intent(in) :: start, ncount, stride, imap
  real(kind=4), dimension(:), intent(out) :: fp
  integer (C_INT) :: c_ncid, c_varid
  integer (C_INT), target :: c_ndims
  integer (C_SIZE_T), dimension(size(start)) :: c_start
  integer (C_SIZE_T), dimension(size(ncount)) :: c_ncount
  integer (C_INTPTR_T), dimension(size(stride)) :: c_stride
  integer (C_INTPTR_T), dimension(size(imap)) :: c_imap
  real (C_FLOAT), dimension(size(fp)), target :: c_fp
  integer i, ndims
  c_ncid = ncid
  c_varid = varid - 1
  c_ndims = size(start)
  ierr = nc_inq_varndims(c_ncid,c_varid,C_LOC(c_ndims))
  ndims = c_ndims
  do i = 1,ndims
    c_start(ndims-i+1) = start(i) - 1
    c_ncount(ndims-i+1) = ncount(i)
    c_stride(ndims-i+1) = stride(i)
    c_imap(ndims-i+1) = imap(i)
  end do
  ierr = nc_get_varm_float(c_ncid,c_varid,c_start,c_ncount,c_stride,c_imap,C_LOC(c_fp))
  fp = c_fp
end function nf_get_varm_real_d1

integer function nf_get_varm_real_d2(ncid,varid,start,ncount,stride,imap,fp) result(ierr)
  implicit none
  integer, intent(in) :: ncid, varid
  integer, dimension(:), intent(in) :: start, ncount, stride, imap
  real(kind=4), dimension(:,:), intent(out) :: fp
  integer (C_INT) :: c_ncid, c_varid
  integer (C_INT), target :: c_ndims
  integer (C_SIZE_T), dimension(size(start)) :: c_start
  integer (C_SIZE_T), dimension(size(ncount)) :: c_ncount
  integer (C_INTPTR_T), dimension(size(stride)) :: c_stride
  integer (C_INTPTR_T), dimension(size(imap)) :: c_imap
  real (C_FLOAT), dimension(size(fp,1),size(fp,2)), target :: c_fp
  integer i, ndims
  c_ncid = ncid
  c_varid = varid - 1
  c_ndims = size(start)
  ierr = nc_inq_varndims(c_ncid,c_varid,C_LOC(c_ndims))
  ndims = c_ndims
  do i = 1,ndims
    c_start(ndims-i+1) = start(i) - 1
    c_ncount(ndims-i+1) = ncount(i)
    c_stride(ndims-i+1) = stride(i)
    c_imap(ndims-i+1) = imap(i)
  end do
  ierr = nc_get_varm_float(c_ncid,c_varid,c_start,c_ncount,c_stride,c_imap,C_LOC(c_fp))
  fp = c_fp
end function nf_get_varm_real_d2

integer function nf_get_varm_real_d3(ncid,varid,start,ncount,stride,imap,fp) result(ierr)
  implicit none
  integer, intent(in) :: ncid, varid
  integer, dimension(:), intent(in) :: start, ncount, stride, imap
  real(kind=4), dimension(:,:,:), intent(out) :: fp
  integer (C_INT) :: c_ncid, c_varid
  integer (C_INT), target :: c_ndims
  integer (C_SIZE_T), dimension(size(start)) :: c_start
  integer (C_SIZE_T), dimension(size(ncount)) :: c_ncount
  integer (C_INTPTR_T), dimension(size(stride)) :: c_stride
  integer (C_INTPTR_T), dimension(size(imap)) :: c_imap
  real (C_FLOAT), dimension(size(fp,1),size(fp,2),size(fp,3)), target :: c_fp
  integer i, ndims
  c_ncid = ncid
  c_varid = varid - 1
  c_ndims = size(start)
  ierr = nc_inq_varndims(c_ncid,c_varid,C_LOC(c_ndims))
  ndims = c_ndims
  do i = 1,ndims
    c_start(ndims-i+1) = start(i) - 1
    c_ncount(ndims-i+1) = ncount(i)
    c_stride(ndims-i+1) = stride(i)
    c_imap(ndims-i+1) = imap(i)
  end do
  ierr = nc_get_varm_float(c_ncid,c_varid,c_start,c_ncount,c_stride,c_imap,C_LOC(c_fp))
  fp = c_fp
end function nf_get_varm_real_d3

integer function nf_get_varm_real_d4(ncid,varid,start,ncount,stride,imap,fp) result(ierr)
  implicit none
  integer, intent(in) :: ncid, varid
  integer, dimension(:), intent(in) :: start, ncount, stride, imap
  real(kind=4), dimension(:,:,:,:), intent(out) :: fp
  integer (C_INT) :: c_ncid, c_varid
  integer (C_INT), target :: c_ndims
  integer (C_SIZE_T), dimension(size(start)) :: c_start
  integer (C_SIZE_T), dimension(size(ncount)) :: c_ncount
  integer (C_INTPTR_T), dimension(size(stride)) :: c_stride
  integer (C_INTPTR_T), dimension(size(imap)) :: c_imap
  real (C_FLOAT), dimension(size(fp,1),size(fp,2),size(fp,3),size(fp,4)), target :: c_fp
  integer i, ndims
  c_ncid = ncid
  c_varid = varid - 1
  c_ndims = size(start)
  ierr = nc_inq_varndims(c_ncid,c_varid,C_LOC(c_ndims))
  ndims = c_ndims
  do i = 1,ndims
    c_start(ndims-i+1) = start(i) - 1
    c_ncount(ndims-i+1) = ncount(i)
    c_stride(ndims-i+1) = stride(i)
    c_imap(ndims-i+1) = imap(i)
  end do
  ierr = nc_get_varm_float(c_ncid,c_varid,c_start,c_ncount,c_stride,c_imap,C_LOC(c_fp))
  fp = c_fp
end function nf_get_varm_real_d4

integer function nf_get_varm_real_d5(ncid,varid,start,ncount,stride,imap,fp) result(ierr)
  implicit none
  integer, intent(in) :: ncid, varid
  integer, dimension(:), intent(in) :: start, ncount, stride, imap
  real(kind=4), dimension(:,:,:,:,:), intent(out) :: fp
  integer (C_INT) :: c_ncid, c_varid
  integer (C_INT), target :: c_ndims
  integer (C_SIZE_T), dimension(size(start)) :: c_start
  integer (C_SIZE_T), dimension(size(ncount)) :: c_ncount
  integer (C_INTPTR_T), dimension(size(stride)) :: c_stride
  integer (C_INTPTR_T), dimension(size(imap)) :: c_imap
  real (C_FLOAT), dimension(size(fp,1),size(fp,2),size(fp,3),size(fp,4),size(fp,5)), &
      target :: c_fp
  integer i, ndims
  c_ncid = ncid
  c_varid = varid - 1
  c_ndims = size(start)
  ierr = nc_inq_varndims(c_ncid,c_varid,C_LOC(c_ndims))
  ndims = c_ndims
  do i = 1,ndims
    c_start(ndims-i+1) = start(i) - 1
    c_ncount(ndims-i+1) = ncount(i)
    c_stride(ndims-i+1) = stride(i)
    c_imap(ndims-i+1) = imap(i)
  end do
  ierr = nc_get_varm_float(c_ncid,c_varid,c_start,c_ncount,c_stride,c_imap,C_LOC(c_fp))
  fp = c_fp
end function nf_get_varm_real_d5

integer function nf_get_varm_real_d6(ncid,varid,start,ncount,stride,imap,fp) result(ierr)
  implicit none
  integer, intent(in) :: ncid, varid
  integer, dimension(:), intent(in) :: start, ncount, stride, imap
  real(kind=4), dimension(:,:,:,:,:,:), intent(out) :: fp
  integer (C_INT) :: c_ncid, c_varid
  integer (C_INT), target :: c_ndims
  integer (C_SIZE_T), dimension(size(start)) :: c_start
  integer (C_SIZE_T), dimension(size(ncount)) :: c_ncount
  integer (C_INTPTR_T), dimension(size(stride)) :: c_stride
  integer (C_INTPTR_T), dimension(size(imap)) :: c_imap
  real (C_FLOAT), dimension(size(fp,1),size(fp,2),size(fp,3),size(fp,4),size(fp,5), &
      size(fp,6)), target :: c_fp
  integer i, ndims
  c_ncid = ncid
  c_varid = varid - 1
  c_ndims = size(start)
  ierr = nc_inq_varndims(c_ncid,c_varid,C_LOC(c_ndims))
  ndims = c_ndims
  do i = 1,ndims
    c_start(ndims-i+1) = start(i) - 1
    c_ncount(ndims-i+1) = ncount(i)
    c_stride(ndims-i+1) = stride(i)
    c_imap(ndims-i+1) = imap(i)
  end do
  ierr = nc_get_varm_float(c_ncid,c_varid,c_start,c_ncount,c_stride,c_imap,C_LOC(c_fp))
  fp = c_fp
end function nf_get_varm_real_d6

integer function nf_get_varm_real_d7(ncid,varid,start,ncount,stride,imap,fp) result(ierr)
  implicit none
  integer, intent(in) :: ncid, varid
  integer, dimension(:), intent(in) :: start, ncount, stride, imap
  real(kind=4), dimension(:,:,:,:,:,:,:), intent(out) :: fp
  integer (C_INT) :: c_ncid, c_varid
  integer (C_INT), target :: c_ndims
  integer (C_SIZE_T), dimension(size(start)) :: c_start
  integer (C_SIZE_T), dimension(size(ncount)) :: c_ncount
  integer (C_INTPTR_T), dimension(size(stride)) :: c_stride
  integer (C_INTPTR_T), dimension(size(imap)) :: c_imap
  real (C_FLOAT), dimension(size(fp,1),size(fp,2),size(fp,3),size(fp,4),size(fp,5), &
      size(fp,6),size(fp,7)), target :: c_fp
  integer i, ndims
  c_ncid = ncid
  c_varid = varid - 1
  c_ndims = size(start)
  ierr = nc_inq_varndims(c_ncid,c_varid,C_LOC(c_ndims))
  ndims = c_ndims
  do i = 1,ndims
    c_start(ndims-i+1) = start(i) - 1
    c_ncount(ndims-i+1) = ncount(i)
    c_stride(ndims-i+1) = stride(i)
    c_imap(ndims-i+1) = imap(i)
  end do
  ierr = nc_get_varm_float(c_ncid,c_varid,c_start,c_ncount,c_stride,c_imap,C_LOC(c_fp))
  fp = c_fp
end function nf_get_varm_real_d7

integer function nf_get_varm_double_d1(ncid,varid,start,ncount,stride,imap,dp) result(ierr)
  implicit none
  integer, intent(in) :: ncid, varid
  integer, dimension(:), intent(in) :: start, ncount, stride, imap
  real(kind=8), dimension(:), intent(out) :: dp
  integer (C_INT) :: c_ncid, c_varid
  integer (C_INT), target :: c_ndims
  integer (C_SIZE_T), dimension(size(start)) :: c_start
  integer (C_SIZE_T), dimension(size(ncount)) :: c_ncount
  integer (C_INTPTR_T), dimension(size(stride)) :: c_stride
  integer (C_INTPTR_T), dimension(size(imap)) :: c_imap
  real (C_DOUBLE), dimension(size(dp)), target :: c_dp
  integer i, ndims
  c_ncid = ncid
  c_varid = varid - 1
  c_ndims = size(start)
  ierr = nc_inq_varndims(c_ncid,c_varid,C_LOC(c_ndims))
  ndims = c_ndims
  do i = 1,ndims
    c_start(ndims-i+1) = start(i) - 1
    c_ncount(ndims-i+1) = ncount(i)
    c_stride(ndims-i+1) = stride(i)
    c_imap(ndims-i+1) = imap(i)
  end do
  ierr = nc_get_varm_double(c_ncid,c_varid,c_start,c_ncount,c_stride,c_imap,C_LOC(c_dp))
  dp = c_dp
end function nf_get_varm_double_d1

integer function nf_get_varm_double_d2(ncid,varid,start,ncount,stride,imap,dp) result(ierr)
  implicit none
  integer, intent(in) :: ncid, varid
  integer, dimension(:), intent(in) :: start, ncount, stride, imap
  real(kind=8), dimension(:,:), intent(out) :: dp
  integer (C_INT) :: c_ncid, c_varid
  integer (C_INT), target :: c_ndims
  integer (C_SIZE_T), dimension(size(start)) :: c_start
  integer (C_SIZE_T), dimension(size(ncount)) :: c_ncount
  integer (C_INTPTR_T), dimension(size(stride)) :: c_stride
  integer (C_INTPTR_T), dimension(size(imap)) :: c_imap
  real (C_DOUBLE), dimension(size(dp,1),size(dp,2)), target :: c_dp
  integer i, ndims
  c_ncid = ncid
  c_varid = varid - 1
  c_ndims = size(start)
  ierr = nc_inq_varndims(c_ncid,c_varid,C_LOC(c_ndims))
  ndims = c_ndims
  do i = 1,ndims
    c_start(ndims-i+1) = start(i) - 1
    c_ncount(ndims-i+1) = ncount(i)
    c_stride(ndims-i+1) = stride(i)
    c_imap(ndims-i+1) = imap(i)
  end do
  ierr = nc_get_varm_double(c_ncid,c_varid,c_start,c_ncount,c_stride,c_imap,C_LOC(c_dp))
  dp = c_dp
end function nf_get_varm_double_d2

integer function nf_get_varm_double_d3(ncid,varid,start,ncount,stride,imap,dp) result(ierr)
  implicit none
  integer, intent(in) :: ncid, varid
  integer, dimension(:), intent(in) :: start, ncount, stride, imap
  real(kind=8), dimension(:,:,:), intent(out) :: dp
  integer (C_INT) :: c_ncid, c_varid
  integer (C_INT), target :: c_ndims
  integer (C_SIZE_T), dimension(size(start)) :: c_start
  integer (C_SIZE_T), dimension(size(ncount)) :: c_ncount
  integer (C_INTPTR_T), dimension(size(stride)) :: c_stride
  integer (C_INTPTR_T), dimension(size(imap)) :: c_imap
  real (C_DOUBLE), dimension(size(dp,1),size(dp,2),size(dp,3)), target :: c_dp
  integer i, ndims
  c_ncid = ncid
  c_varid = varid - 1
  c_ndims = size(start)
  ierr = nc_inq_varndims(c_ncid,c_varid,C_LOC(c_ndims))
  ndims = c_ndims
  do i = 1,ndims
    c_start(ndims-i+1) = start(i) - 1
    c_ncount(ndims-i+1) = ncount(i)
    c_stride(ndims-i+1) = stride(i)
    c_imap(ndims-i+1) = imap(i)
  end do
  ierr = nc_get_varm_double(c_ncid,c_varid,c_start,c_ncount,c_stride,c_imap,C_LOC(c_dp))
  dp = c_dp
end function nf_get_varm_double_d3

integer function nf_get_varm_double_d4(ncid,varid,start,ncount,stride,imap,dp) result(ierr)
  implicit none
  integer, intent(in) :: ncid, varid
  integer, dimension(:), intent(in) :: start, ncount, stride, imap
  real(kind=8), dimension(:,:,:,:), intent(out) :: dp
  integer (C_INT) :: c_ncid, c_varid
  integer (C_INT), target :: c_ndims
  integer (C_SIZE_T), dimension(size(start)) :: c_start
  integer (C_SIZE_T), dimension(size(ncount)) :: c_ncount
  integer (C_INTPTR_T), dimension(size(stride)) :: c_stride
  integer (C_INTPTR_T), dimension(size(imap)) :: c_imap
  real (C_DOUBLE), dimension(size(dp,1),size(dp,2),size(dp,3),size(dp,4)), target :: c_dp
  integer i, ndims
  c_ncid = ncid
  c_varid = varid - 1
  c_ndims = size(start)
  ierr = nc_inq_varndims(c_ncid,c_varid,C_LOC(c_ndims))
  ndims = c_ndims
  do i = 1,ndims
    c_start(ndims-i+1) = start(i) - 1
    c_ncount(ndims-i+1) = ncount(i)
    c_stride(ndims-i+1) = stride(i)
    c_imap(ndims-i+1) = imap(i)
  end do
  ierr = nc_get_varm_double(c_ncid,c_varid,c_start,c_ncount,c_stride,c_imap,C_LOC(c_dp))
  dp = c_dp
end function nf_get_varm_double_d4

integer function nf_get_varm_double_d5(ncid,varid,start,ncount,stride,imap,dp) result(ierr)
  implicit none
  integer, intent(in) :: ncid, varid
  integer, dimension(:), intent(in) :: start, ncount, stride, imap
  real(kind=8), dimension(:,:,:,:,:), intent(out) :: dp
  integer (C_INT) :: c_ncid, c_varid
  integer (C_INT), target :: c_ndims
  integer (C_SIZE_T), dimension(size(start)) :: c_start
  integer (C_SIZE_T), dimension(size(ncount)) :: c_ncount
  integer (C_INTPTR_T), dimension(size(stride)) :: c_stride
  integer (C_INTPTR_T), dimension(size(imap)) :: c_imap
  real (C_DOUBLE), dimension(size(dp,1),size(dp,2),size(dp,3),size(dp,4),size(dp,5)), &
      target :: c_dp
  integer i, ndims
  c_ncid = ncid
  c_varid = varid - 1
  c_ndims = size(start)
  ierr = nc_inq_varndims(c_ncid,c_varid,C_LOC(c_ndims))
  ndims = c_ndims
  do i = 1,ndims
    c_start(ndims-i+1) = start(i) - 1
    c_ncount(ndims-i+1) = ncount(i)
    c_stride(ndims-i+1) = stride(i)
    c_imap(ndims-i+1) = imap(i)
  end do
  ierr = nc_get_varm_double(c_ncid,c_varid,c_start,c_ncount,c_stride,c_imap,C_LOC(c_dp))
  dp = c_dp
end function nf_get_varm_double_d5

integer function nf_get_varm_double_d6(ncid,varid,start,ncount,stride,imap,dp) result(ierr)
  implicit none
  integer, intent(in) :: ncid, varid
  integer, dimension(:), intent(in) :: start, ncount, stride, imap
  real(kind=8), dimension(:,:,:,:,:,:), intent(out) :: dp
  integer (C_INT) :: c_ncid, c_varid
  integer (C_INT), target :: c_ndims
  integer (C_SIZE_T), dimension(size(start)) :: c_start
  integer (C_SIZE_T), dimension(size(ncount)) :: c_ncount
  integer (C_INTPTR_T), dimension(size(stride)) :: c_stride
  integer (C_INTPTR_T), dimension(size(imap)) :: c_imap
  real (C_DOUBLE), dimension(size(dp,1),size(dp,2),size(dp,3),size(dp,4),size(dp,5), &
      size(dp,6)), target :: c_dp
  integer i, ndims
  c_ncid = ncid
  c_varid = varid - 1
  c_ndims = size(start)
  ierr = nc_inq_varndims(c_ncid,c_varid,C_LOC(c_ndims))
  ndims = c_ndims
  do i = 1,ndims
    c_start(ndims-i+1) = start(i) - 1
    c_ncount(ndims-i+1) = ncount(i)
    c_stride(ndims-i+1) = stride(i)
    c_imap(ndims-i+1) = imap(i)
  end do
  ierr = nc_get_varm_double(c_ncid,c_varid,c_start,c_ncount,c_stride,c_imap,C_LOC(c_dp))
  dp = c_dp
end function nf_get_varm_double_d6

integer function nf_get_varm_double_d7(ncid,varid,start,ncount,stride,imap,dp) result(ierr)
  implicit none
  integer, intent(in) :: ncid, varid
  integer, dimension(:), intent(in) :: start, ncount, stride, imap
  real(kind=8), dimension(:,:,:,:,:,:,:), intent(out) :: dp
  integer (C_INT) :: c_ncid, c_varid
  integer (C_INT), target :: c_ndims
  integer (C_SIZE_T), dimension(size(start)) :: c_start
  integer (C_SIZE_T), dimension(size(ncount)) :: c_ncount
  integer (C_INTPTR_T), dimension(size(stride)) :: c_stride
  integer (C_INTPTR_T), dimension(size(imap)) :: c_imap
  real (C_DOUBLE), dimension(size(dp,1),size(dp,2),size(dp,3),size(dp,4),size(dp,5), &
      size(dp,6),size(dp,7)), target :: c_dp
  integer i, ndims
  c_ncid = ncid
  c_varid = varid - 1
  c_ndims = size(start)
  ierr = nc_inq_varndims(c_ncid,c_varid,C_LOC(c_ndims))
  ndims = c_ndims
  do i = 1,ndims
    c_start(ndims-i+1) = start(i) - 1
    c_ncount(ndims-i+1) = ncount(i)
    c_stride(ndims-i+1) = stride(i)
    c_imap(ndims-i+1) = imap(i)
  end do
  ierr = nc_get_varm_double(c_ncid,c_varid,c_start,c_ncount,c_stride,c_imap,C_LOC(c_dp))
  dp = c_dp
end function nf_get_varm_double_d7

integer function nf_def_dim(ncid,name,size,dimid) result(ierr)
  implicit none
  integer, intent(in) :: ncid, size
  integer, intent(out) :: dimid
  character(len=*), intent(in) :: name
  integer (C_INT) :: c_ncid
  integer (C_SIZE_T) :: c_size
  integer (C_INT), target :: c_dimid
  character, dimension(len(name)) :: c_name
  c_ncid = ncid
  call cf_strcopy(name,c_name)
  c_size = size
  ierr = nc_def_dim(c_ncid,c_name,c_size,C_LOC(c_dimid))
  dimid = c_dimid + 1
end function nf_def_dim

integer function nf_def_var_s(ncid,name,xtype,ndims,dimids,varid) result(ierr)
  implicit none
  integer, intent(in) :: ncid, xtype, ndims
  integer, intent(out) :: varid
  integer, intent(in) :: dimids
  character(len=*), intent(in) :: name
  integer (C_INT) :: c_ncid, c_xtype, c_ndims
  integer (C_INT), target :: c_varid
  integer (C_INT), dimension(1) :: c_dimids
  character, dimension(len(name)) :: c_name
  c_ncid = ncid
  call cf_strcopy(name,c_name)
  c_xtype = xtype
  c_ndims = ndims
  c_dimids = dimids - 1
  ierr = nc_def_var(c_ncid,c_name,c_xtype,c_ndims,c_dimids,C_LOC(c_varid))
  varid = c_varid + 1
end function nf_def_var_s

integer function nf_def_var_v(ncid,name,xtype,ndims,dimids,varid) result(ierr)
  implicit none
  integer, intent(in) :: ncid, xtype, ndims
  integer, intent(out) :: varid
  integer, dimension(:), intent(in) :: dimids
  character(len=*), intent(in) :: name
  integer (C_INT) :: c_ncid, c_xtype, c_ndims
  integer (C_INT), target :: c_varid
  integer (C_INT), dimension(size(dimids)) :: c_dimids
  character, dimension(len(name)) :: c_name
  integer i
  c_ncid = ncid
  call cf_strcopy(name,c_name)
  c_xtype = xtype
  c_ndims = ndims
  do i = 1,ndims
    c_dimids(ndims-i+1) = dimids(i) - 1
  end do
  ierr = nc_def_var(c_ncid,c_name,c_xtype,c_ndims,c_dimids,C_LOC(c_varid))
  varid = c_varid + 1
end function nf_def_var_v

integer function nf_rename_dim(ncid,dimid,name) result(ierr)
  implicit none
  integer, intent(in) :: ncid, dimid
  character(len=*), intent(in) :: name
  integer (C_INT) :: c_ncid, c_dimid
  character, dimension(len(name)) :: c_name
  c_ncid = ncid
  c_dimid = dimid - 1
  call cf_strcopy(name,c_name)
  ierr = nc_rename_dim(c_ncid,c_dimid,c_name)
end function nf_rename_dim

integer function nf_rename_att(ncid,varid,name,newname) result(ierr)
  implicit none
  integer, intent(in) :: ncid, varid
  character(len=*), intent(in) :: name
  character(len=*), intent(in) :: newname
  integer (C_INT) :: c_ncid, c_varid
  character, dimension(len(name)) :: c_name
  character, dimension(len(newname)) :: c_newname
  c_ncid = ncid
  c_varid = varid - 1
  call cf_strcopy(name,c_name)
  call cf_strcopy(newname,c_newname)
  ierr = nc_rename_att(c_ncid,c_varid,c_name,c_newname)
end function nf_rename_att

integer function nf_rename_var(ncid,varid,name) result(ierr)
  implicit none
  integer, intent(in) :: ncid, varid
  character(len=*), intent(in) :: name
  integer (C_INT) :: c_ncid, c_varid
  character, dimension(len(name)) :: c_name
  c_ncid = ncid
  c_varid = varid - 1
  call cf_strcopy(name,c_name)
  ierr = nc_rename_var(c_ncid,c_varid,c_name)
end function nf_rename_var

integer function nf_put_att_text(ncid,varid,name,clen,tp) result(ierr)
  implicit none
  integer, intent(in) :: ncid, varid, clen
  character(len=*), intent(in) :: name
  character(len=*), intent(in) :: tp
  integer (C_INT) :: c_ncid, c_varid
  integer (C_SIZE_T) :: c_clen
  character, dimension(len(name)) :: c_name
  character, dimension(len(tp)) :: c_tp
  c_ncid = ncid
  c_varid = varid - 1
  call cf_strcopy(name,c_name)
  c_clen = clen
  call cf_strcopy(tp,c_tp)
  ierr = nc_put_att_text(c_ncid,c_varid,c_name,c_clen,c_tp)
end function nf_put_att_text

integer function nf_put_att_int2_s(ncid,varid,name,xtype,slen,sp) result(ierr)
  implicit none
  integer, intent(in) :: ncid, varid, xtype, slen
  integer(kind=2), intent(in) :: sp
  character(len=*), intent(in) :: name
  integer (C_INT) :: c_ncid, c_varid, c_xtype
  integer (C_SIZE_T) :: c_slen
  integer (C_SHORT), target :: c_sp  
  character, dimension(len(name)) :: c_name
  c_ncid = ncid
  c_varid = varid - 1
  call cf_strcopy(name,c_name)
  c_xtype = xtype
  c_slen = slen
  c_sp = sp
  ierr = nc_put_att_short(c_ncid,c_varid,c_name,c_xtype,c_slen,C_LOC(c_sp))
end function nf_put_att_int2_s

integer function nf_put_att_int2_v(ncid,varid,name,xtype,slen,sp) result(ierr)
  implicit none
  integer, intent(in) :: ncid, varid, xtype, slen
  integer(kind=2), dimension(:), intent(in) :: sp
  character(len=*), intent(in) :: name
  integer (C_INT) :: c_ncid, c_varid, c_xtype
  integer (C_SIZE_T) :: c_slen
  integer (C_SHORT), dimension(slen), target :: c_sp  
  character, dimension(len(name)) :: c_name
  c_ncid = ncid
  c_varid = varid - 1
  call cf_strcopy(name,c_name)
  c_xtype = xtype
  c_slen = slen
  c_sp = sp(1:slen)
  ierr = nc_put_att_short(c_ncid,c_varid,c_name,c_xtype,c_slen,C_LOC(c_sp))
end function nf_put_att_int2_v

integer function nf_put_att_real_s(ncid,varid,name,xtype,rlen,fp) result(ierr)
  implicit none
  integer, intent(in) :: ncid, varid, xtype, rlen
  real(kind=4), intent(in) :: fp
  character(len=*), intent(in) :: name
  integer (C_INT) :: c_ncid, c_varid, c_xtype
  integer (C_SIZE_T) :: c_rlen
  real (C_FLOAT), target :: c_fp  
  character, dimension(len(name)) :: c_name
  c_ncid = ncid
  c_varid = varid - 1
  call cf_strcopy(name,c_name)
  c_xtype = xtype
  c_rlen = rlen
  c_fp = fp
  ierr = nc_put_att_float(c_ncid,c_varid,c_name,c_xtype,c_rlen,C_LOC(c_fp))
end function nf_put_att_real_s

integer function nf_put_att_real_v(ncid,varid,name,xtype,rlen,fp) result(ierr)
  implicit none
  integer, intent(in) :: ncid, varid, xtype, rlen
  real(kind=4), dimension(:), intent(in) :: fp
  character(len=*), intent(in) :: name
  integer (C_INT) :: c_ncid, c_varid, c_xtype
  integer (C_SIZE_T) :: c_rlen
  real (C_FLOAT), dimension(size(fp)), target :: c_fp  
  character, dimension(len(name)) :: c_name
  c_ncid = ncid
  c_varid = varid - 1
  call cf_strcopy(name,c_name)
  c_xtype = xtype
  c_rlen = rlen
  c_fp = fp
  ierr = nc_put_att_float(c_ncid,c_varid,c_name,c_xtype,c_rlen,C_LOC(c_fp))
end function nf_put_att_real_v

integer function nf_put_att_int_s(ncid,varid,name,xtype,slen,ip) result(ierr)
  implicit none
  integer, intent(in) :: ncid, varid, xtype, slen
  integer(kind=4), intent(in) :: ip
  character(len=*), intent(in) :: name
  integer (C_INT) :: c_ncid, c_varid, c_xtype
  integer (C_SIZE_T) :: c_slen
  integer (C_INT), target :: c_ip  
  character, dimension(len(name)) :: c_name
  c_ncid = ncid
  c_varid = varid - 1
  call cf_strcopy(name,c_name)
  c_xtype = xtype
  c_slen = slen
  c_ip = ip
  ierr = nc_put_att_int(c_ncid,c_varid,c_name,c_xtype,c_slen,C_LOC(c_ip))
end function nf_put_att_int_s

integer function nf_put_att_int_v(ncid,varid,name,xtype,slen,ip) result(ierr)
  implicit none
  integer, intent(in) :: ncid, varid, xtype, slen
  integer(kind=4), dimension(:), intent(in) :: ip
  character(len=*), intent(in) :: name
  integer (C_INT) :: c_ncid, c_varid, c_xtype
  integer (C_SIZE_T) :: c_slen
  integer (C_INT), dimension(size(ip)), target :: c_ip  
  character, dimension(len(name)) :: c_name
  c_ncid = ncid
  c_varid = varid - 1
  call cf_strcopy(name,c_name)
  c_xtype = xtype
  c_slen = slen
  c_ip = ip
  ierr = nc_put_att_int(c_ncid,c_varid,c_name,c_xtype,c_slen,C_LOC(c_ip))
end function nf_put_att_int_v

integer function nf_put_att_int1_s(ncid,varid,name,xtype,blen,bp) result(ierr)
  implicit none
  integer, intent(in) :: ncid, varid, xtype, blen
  integer(kind=1), intent(in) :: bp
  character(len=*), intent(in) :: name
  integer (C_INT) :: c_ncid, c_varid, c_xtype
  integer (C_SIZE_T) :: c_blen
  integer (C_SIGNED_CHAR), dimension(1), target :: c_bp  
  character, dimension(len(name)) :: c_name
  c_ncid = ncid
  c_varid = varid - 1
  call cf_strcopy(name,c_name)
  c_xtype = xtype
  c_blen = blen
  c_bp = bp
  ierr = nc_put_att_schar(c_ncid,c_varid,c_name,c_xtype,c_blen,C_LOC(c_bp))
end function nf_put_att_int1_s

integer function nf_put_att_int1_v(ncid,varid,name,xtype,blen,bp) result(ierr)
  implicit none
  integer, intent(in) :: ncid, varid, xtype, blen
  integer(kind=1), dimension(:), intent(in) :: bp
  character(len=*), intent(in) :: name
  integer (C_INT) :: c_ncid, c_varid, c_xtype
  integer (C_SIZE_T) :: c_blen
  integer (C_SIGNED_CHAR), dimension(size(bp)), target :: c_bp  
  character, dimension(len(name)) :: c_name
  c_ncid = ncid
  c_varid = varid - 1
  call cf_strcopy(name,c_name)
  c_xtype = xtype
  c_blen = blen
  c_bp = bp
  ierr = nc_put_att_schar(c_ncid,c_varid,c_name,c_xtype,c_blen,C_LOC(c_bp))
end function nf_put_att_int1_v

integer function nf_put_att_double_s(ncid,varid,name,xtype,dlen,dp) result(ierr)
  implicit none
  integer, intent(in) :: ncid, varid, xtype, dlen
  real(kind=8), intent(in) :: dp
  character(len=*), intent(in) :: name
  integer (C_INT) :: c_ncid, c_varid, c_xtype
  integer (C_SIZE_T) :: c_dlen
  real (C_DOUBLE), target :: c_dp  
  character, dimension(len(name)) :: c_name
  c_ncid = ncid
  c_varid = varid - 1
  call cf_strcopy(name,c_name)
  c_xtype = xtype
  c_dlen = dlen
  c_dp = dp
  ierr = nc_put_att_double(c_ncid,c_varid,c_name,c_xtype,c_dlen,C_LOC(c_dp))
end function nf_put_att_double_s

integer function nf_put_att_double_v(ncid,varid,name,xtype,dlen,dp) result(ierr)
  implicit none
  integer, intent(in) :: ncid, varid, xtype, dlen
  real(kind=8), dimension(:), intent(in) :: dp
  character(len=*), intent(in) :: name
  integer (C_INT) :: c_ncid, c_varid, c_xtype
  integer (C_SIZE_T) :: c_dlen
  real (C_DOUBLE), dimension(size(dp)), target :: c_dp  
  character, dimension(len(name)) :: c_name
  c_ncid = ncid
  c_varid = varid - 1
  call cf_strcopy(name,c_name)
  c_xtype = xtype
  c_dlen = dlen
  c_dp = dp
  ierr = nc_put_att_double(c_ncid,c_varid,c_name,c_xtype,c_dlen,C_LOC(c_dp))
end function nf_put_att_double_v

integer function nf_put_vara_real_d1(ncid,varid,start,ncount,fp) result(ierr)
  implicit none
  integer, intent(in) :: ncid, varid
  integer, dimension(:), intent(in) :: start, ncount
  real(kind=4), dimension(:), intent(in) :: fp
  integer (C_INT) :: c_ncid, c_varid
  integer (C_INT), target :: c_ndims
  integer (C_SIZE_T), dimension(size(start)) :: c_start, c_ncount
  real (C_FLOAT), dimension(size(fp)), target :: c_fp
  integer i, ndims
  c_ncid = ncid
  c_varid = varid - 1
  c_ndims = size(start)
  ierr = nc_inq_varndims(c_ncid,c_varid,C_LOC(c_ndims))
  ndims = c_ndims
  do i = 1,ndims
    c_start(ndims-i+1) = start(i) - 1
    c_ncount(ndims-i+1) = ncount(i)
  end do
  c_fp = fp
  ierr = nc_put_vara_float(c_ncid,c_varid,c_start,c_ncount,C_LOC(c_fp))
end function nf_put_vara_real_d1

integer function nf_put_vara_real_d2(ncid,varid,start,ncount,fp) result(ierr)
  implicit none
  integer, intent(in) :: ncid, varid
  integer, dimension(:), intent(in) :: start, ncount
  real(kind=4), dimension(:,:), intent(in) :: fp
  integer (C_INT) :: c_ncid, c_varid
  integer (C_INT), target :: c_ndims
  integer (C_SIZE_T), dimension(size(start)) :: c_start, c_ncount
  real (C_FLOAT), dimension(size(fp,1),size(fp,2)), target :: c_fp
  integer i, ndims
  c_ncid = ncid
  c_varid = varid - 1
  c_ndims = size(start)
  ierr = nc_inq_varndims(c_ncid,c_varid,C_LOC(c_ndims))
  ndims = c_ndims
  do i = 1,ndims
    c_start(ndims-i+1) = start(i) - 1
    c_ncount(ndims-i+1) = ncount(i)
  end do
  c_fp = fp
  ierr = nc_put_vara_float(c_ncid,c_varid,c_start,c_ncount,C_LOC(c_fp))
end function nf_put_vara_real_d2

integer function nf_put_vara_real_d3(ncid,varid,start,ncount,fp) result(ierr)
  implicit none
  integer, intent(in) :: ncid, varid
  integer, dimension(:), intent(in) :: start, ncount
  real(kind=4), dimension(:,:,:), intent(in) :: fp
  integer (C_INT) :: c_ncid, c_varid
  integer (C_INT), target :: c_ndims
  integer (C_SIZE_T), dimension(size(start)) :: c_start, c_ncount
  real (C_FLOAT), dimension(size(fp,1),size(fp,2),size(fp,3)), target :: c_fp
  integer i, ndims
  c_ncid = ncid
  c_varid = varid - 1
  c_ndims = size(start)
  ierr = nc_inq_varndims(c_ncid,c_varid,C_LOC(c_ndims))
  ndims = c_ndims
  do i = 1,ndims
    c_start(ndims-i+1) = start(i) - 1
    c_ncount(ndims-i+1) = ncount(i)
  end do
  c_fp = fp
  ierr = nc_put_vara_float(c_ncid,c_varid,c_start,c_ncount,C_LOC(c_fp))
end function nf_put_vara_real_d3

integer function nf_put_vara_real_d4(ncid,varid,start,ncount,fp) result(ierr)
  implicit none
  integer, intent(in) :: ncid, varid
  integer, dimension(:), intent(in) :: start, ncount
  real(kind=4), dimension(:,:,:,:), intent(in) :: fp
  integer (C_INT) :: c_ncid, c_varid
  integer (C_INT), target :: c_ndims
  integer (C_SIZE_T), dimension(size(start)) :: c_start, c_ncount
  real (C_FLOAT), dimension(size(fp,1),size(fp,2),size(fp,3),size(fp,4)), target :: c_fp
  integer i, ndims
  c_ncid = ncid
  c_varid = varid - 1
  c_ndims = size(start)
  ierr = nc_inq_varndims(c_ncid,c_varid,C_LOC(c_ndims))
  ndims = c_ndims
  do i = 1,ndims
    c_start(ndims-i+1) = start(i) - 1
    c_ncount(ndims-i+1) = ncount(i)
  end do
  c_fp = fp
  ierr = nc_put_vara_float(c_ncid,c_varid,c_start,c_ncount,C_LOC(c_fp))
end function nf_put_vara_real_d4

integer function nf_put_vara_real_d5(ncid,varid,start,ncount,fp) result(ierr)
  implicit none
  integer, intent(in) :: ncid, varid
  integer, dimension(:), intent(in) :: start, ncount
  real(kind=4), dimension(:,:,:,:,:), intent(in) :: fp
  integer (C_INT) :: c_ncid, c_varid
  integer (C_INT), target :: c_ndims
  integer (C_SIZE_T), dimension(size(start)) :: c_start, c_ncount
  real (C_FLOAT), dimension(size(fp,1),size(fp,2),size(fp,3),size(fp,4),size(fp,5)), &
      target :: c_fp
  integer i, ndims
  c_ncid = ncid
  c_varid = varid - 1
  c_ndims = size(start)
  ierr = nc_inq_varndims(c_ncid,c_varid,C_LOC(c_ndims))
  ndims = c_ndims
  do i = 1,ndims
    c_start(ndims-i+1) = start(i) - 1
    c_ncount(ndims-i+1) = ncount(i)
  end do
  c_fp = fp
  ierr = nc_put_vara_float(c_ncid,c_varid,c_start,c_ncount,C_LOC(c_fp))
end function nf_put_vara_real_d5

integer function nf_put_vara_int_d1(ncid,varid,start,ncount,ip) result(ierr)
  implicit none
  integer, intent(in) :: ncid, varid
  integer, dimension(:), intent(in) :: start, ncount
  integer(kind=4), dimension(:), intent(in) :: ip
  integer (C_INT) :: c_ncid, c_varid
  integer (C_INT), target :: c_ndims
  integer (C_SIZE_T), dimension(size(start)) :: c_start, c_ncount
  integer (C_INT), dimension(size(ip)), target :: c_ip
  integer i, ndims
  c_ncid = ncid
  c_varid = varid - 1
  c_ndims = size(start)
  ierr = nc_inq_varndims(c_ncid,c_varid,C_LOC(c_ndims))
  ndims = c_ndims
  do i = 1,ndims
    c_start(ndims-i+1) = start(i) - 1
    c_ncount(ndims-i+1) = ncount(i)
  end do
  c_ip = ip
  ierr = nc_put_vara_int(c_ncid,c_varid,c_start,c_ncount,C_LOC(c_ip))
end function nf_put_vara_int_d1

integer function nf_put_vara_int_d2(ncid,varid,start,ncount,ip) result(ierr)
  implicit none
  integer, intent(in) :: ncid, varid
  integer, dimension(:), intent(in) :: start, ncount
  integer(kind=4), dimension(:,:), intent(in) :: ip
  integer (C_INT) :: c_ncid, c_varid
  integer (C_INT), target :: c_ndims
  integer (C_SIZE_T), dimension(size(start)) :: c_start, c_ncount
  integer (C_INT), dimension(size(ip,1),size(ip,2)), target :: c_ip
  integer i, ndims
  c_ncid = ncid
  c_varid = varid - 1
  c_ndims = size(start)
  ierr = nc_inq_varndims(c_ncid,c_varid,C_LOC(c_ndims))
  ndims = c_ndims
  do i = 1,ndims
    c_start(ndims-i+1) = start(i) - 1
    c_ncount(ndims-i+1) = ncount(i)
  end do
  c_ip = ip
  ierr = nc_put_vara_int(c_ncid,c_varid,c_start,c_ncount,C_LOC(c_ip))
end function nf_put_vara_int_d2

integer function nf_put_vara_int_d3(ncid,varid,start,ncount,ip) result(ierr)
  implicit none
  integer, intent(in) :: ncid, varid
  integer, dimension(:), intent(in) :: start, ncount
  integer(kind=4), dimension(:,:,:), intent(in) :: ip
  integer (C_INT) :: c_ncid, c_varid
  integer (C_INT), target :: c_ndims
  integer (C_SIZE_T), dimension(size(start)) :: c_start, c_ncount
  integer (C_INT), dimension(size(ip,1),size(ip,2),size(ip,3)), target :: c_ip
  integer i, ndims
  c_ncid = ncid
  c_varid = varid - 1
  c_ndims = size(start)
  ierr = nc_inq_varndims(c_ncid,c_varid,C_LOC(c_ndims))
  ndims = c_ndims
  do i = 1,ndims
    c_start(ndims-i+1) = start(i) - 1
    c_ncount(ndims-i+1) = ncount(i)
  end do
  c_ip = ip
  ierr = nc_put_vara_int(c_ncid,c_varid,c_start,c_ncount,C_LOC(c_ip))
end function nf_put_vara_int_d3

integer function nf_put_vara_int_d4(ncid,varid,start,ncount,ip) result(ierr)
  implicit none
  integer, intent(in) :: ncid, varid
  integer, dimension(:), intent(in) :: start, ncount
  integer(kind=4), dimension(:,:,:,:), intent(in) :: ip
  integer (C_INT) :: c_ncid, c_varid
  integer (C_INT), target :: c_ndims
  integer (C_SIZE_T), dimension(size(start)) :: c_start, c_ncount
  integer (C_INT), dimension(size(ip,1),size(ip,2),size(ip,3),size(ip,4)), target :: c_ip
  integer i, ndims
  c_ncid = ncid
  c_varid = varid - 1
  c_ndims = size(start)
  ierr = nc_inq_varndims(c_ncid,c_varid,C_LOC(c_ndims))
  ndims = c_ndims
  do i = 1,ndims
    c_start(ndims-i+1) = start(i) - 1
    c_ncount(ndims-i+1) = ncount(i)
  end do
  c_ip = ip
  ierr = nc_put_vara_int(c_ncid,c_varid,c_start,c_ncount,C_LOC(c_ip))
end function nf_put_vara_int_d4

integer function nf_put_vara_int_d5(ncid,varid,start,ncount,ip) result(ierr)
  implicit none
  integer, intent(in) :: ncid, varid
  integer, dimension(:), intent(in) :: start, ncount
  integer(kind=4), dimension(:,:,:,:,:), intent(in) :: ip
  integer (C_INT) :: c_ncid, c_varid
  integer (C_INT), target :: c_ndims
  integer (C_SIZE_T), dimension(size(start)) :: c_start, c_ncount
  integer (C_INT), dimension(size(ip,1),size(ip,2),size(ip,3),size(ip,4),size(ip,5)), &
      target :: c_ip
  integer i, ndims
  c_ncid = ncid
  c_varid = varid - 1
  c_ndims = size(start)
  ierr = nc_inq_varndims(c_ncid,c_varid,C_LOC(c_ndims))
  ndims = c_ndims
  do i = 1,ndims
    c_start(ndims-i+1) = start(i) - 1
    c_ncount(ndims-i+1) = ncount(i)
  end do
  c_ip = ip
  ierr = nc_put_vara_int(c_ncid,c_varid,c_start,c_ncount,C_LOC(c_ip))
end function nf_put_vara_int_d5

integer function nf_put_vara_int2_d1(ncid,varid,start,ncount,sp) result(ierr)
  implicit none
  integer, intent(in) :: ncid, varid
  integer, dimension(:), intent(in) :: start, ncount
  integer(kind=2), dimension(:), intent(in) :: sp
  integer (C_INT) :: c_ncid, c_varid
  integer (C_INT), target :: c_ndims
  integer (C_SIZE_T), dimension(size(start)) :: c_start, c_ncount
  integer (C_SHORT), dimension(size(sp)), target :: c_sp
  integer i, ndims
  c_ncid = ncid
  c_varid = varid - 1
  c_ndims = size(start)
  ierr = nc_inq_varndims(c_ncid,c_varid,C_LOC(c_ndims))
  ndims = c_ndims
  do i = 1,ndims
    c_start(ndims-i+1) = start(i) - 1
    c_ncount(ndims-i+1) = ncount(i)
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
  integer (C_INT), target :: c_ndims
  integer (C_SIZE_T), dimension(size(start)) :: c_start, c_ncount
  integer (C_SHORT), dimension(size(sp,1),size(sp,2)), target :: c_sp
  integer i, ndims
  c_ncid = ncid
  c_varid = varid - 1
  c_ndims = size(start)
  ierr = nc_inq_varndims(c_ncid,c_varid,C_LOC(c_ndims))
  ndims = c_ndims
  do i = 1,ndims
    c_start(ndims-i+1) = start(i) - 1
    c_ncount(ndims-i+1) = ncount(i)
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
  integer (C_INT), target :: c_ndims
  integer (C_SIZE_T), dimension(size(start)) :: c_start, c_ncount
  integer (C_SHORT), dimension(size(sp,1),size(sp,2),size(sp,3)), target :: c_sp
  integer i, ndims
  c_ncid = ncid
  c_varid = varid - 1
  c_ndims = size(start)
  ierr = nc_inq_varndims(c_ncid,c_varid,C_LOC(c_ndims))
  ndims = c_ndims
  do i = 1,ndims
    c_start(ndims-i+1) = start(i) - 1
    c_ncount(ndims-i+1) = ncount(i)
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
  integer (C_INT), target :: c_ndims
  integer (C_SIZE_T), dimension(size(start)) :: c_start, c_ncount
  integer (C_SHORT), dimension(size(sp,1),size(sp,2),size(sp,3),size(sp,4)), target :: c_sp
  integer i, ndims
  c_ncid = ncid
  c_varid = varid - 1
  c_ndims = size(start)
  ierr = nc_inq_varndims(c_ncid,c_varid,C_LOC(c_ndims))
  ndims = c_ndims
  do i = 1,ndims
    c_start(ndims-i+1) = start(i) - 1
    c_ncount(ndims-i+1) = ncount(i)
  end do
  c_sp = sp
  ierr = nc_put_vara_short(c_ncid,c_varid,c_start,c_ncount,C_LOC(c_sp))
end function nf_put_vara_int2_d4

integer function nf_put_vara_int2_d5(ncid,varid,start,ncount,sp) result(ierr)
  implicit none
  integer, intent(in) :: ncid, varid
  integer, dimension(:), intent(in) :: start, ncount
  integer(kind=2), dimension(:,:,:,:,:), intent(in) :: sp
  integer (C_INT) :: c_ncid, c_varid
  integer (C_INT), target :: c_ndims
  integer (C_SIZE_T), dimension(size(start)) :: c_start, c_ncount
  integer (C_SHORT), dimension(size(sp,1),size(sp,2),size(sp,3),size(sp,4),size(sp,5)), &
      target :: c_sp
  integer i, ndims
  c_ncid = ncid
  c_varid = varid - 1
  c_ndims = size(start)
  ierr = nc_inq_varndims(c_ncid,c_varid,C_LOC(c_ndims))
  ndims = c_ndims
  do i = 1,ndims
    c_start(ndims-i+1) = start(i) - 1
    c_ncount(ndims-i+1) = ncount(i)
  end do
  c_sp = sp
  ierr = nc_put_vara_short(c_ncid,c_varid,c_start,c_ncount,C_LOC(c_sp))
end function nf_put_vara_int2_d5

integer function nf_put_vara_double_d1(ncid,varid,start,ncount,dp) result(ierr)
  implicit none
  integer, intent(in) :: ncid, varid
  integer, dimension(:), intent(in) :: start, ncount
  real(kind=8), dimension(:), intent(in) :: dp
  integer (C_INT) :: c_ncid, c_varid
  integer (C_INT), target :: c_ndims
  integer (C_SIZE_T), dimension(size(start)) :: c_start, c_ncount
  real (C_DOUBLE), dimension(size(dp)), target :: c_dp
  integer i, ndims
  c_ncid = ncid
  c_varid = varid - 1
  c_ndims = size(start)
  ierr = nc_inq_varndims(c_ncid,c_varid,C_LOC(c_ndims))
  ndims = c_ndims
  do i = 1,ndims
    c_start(ndims-i+1) = start(i) - 1
    c_ncount(ndims-i+1) = ncount(i)
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
  integer (C_INT), target :: c_ndims
  integer (C_SIZE_T), dimension(size(start)) :: c_start, c_ncount
  real (C_DOUBLE), dimension(size(dp,1),size(dp,2)), target :: c_dp
  integer i, ndims
  c_ncid = ncid
  c_varid = varid - 1
  c_ndims = size(start)
  ierr = nc_inq_varndims(c_ncid,c_varid,C_LOC(c_ndims))
  ndims = c_ndims
  do i = 1,ndims
    c_start(ndims-i+1) = start(i) - 1
    c_ncount(ndims-i+1) = ncount(i)
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
  integer (C_INT), target :: c_ndims
  integer (C_SIZE_T), dimension(size(start)) :: c_start, c_ncount
  real (C_DOUBLE), dimension(size(dp,1),size(dp,2),size(dp,3)), target :: c_dp
  integer i, ndims
  c_ncid = ncid
  c_varid = varid - 1
  c_ndims = size(start)
  ierr = nc_inq_varndims(c_ncid,c_varid,C_LOC(c_ndims))
  ndims = c_ndims
  do i = 1,ndims
    c_start(ndims-i+1) = start(i) - 1
    c_ncount(ndims-i+1) = ncount(i)
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
  integer (C_INT), target :: c_ndims
  integer (C_SIZE_T), dimension(size(start)) :: c_start, c_ncount
  real (C_DOUBLE), dimension(size(dp,1),size(dp,2),size(dp,3),size(dp,4)), target :: c_dp
  integer i, ndims
  c_ncid = ncid
  c_varid = varid - 1
  c_ndims = size(start)
  ierr = nc_inq_varndims(c_ncid,c_varid,C_LOC(c_ndims))
  ndims = c_ndims
  do i = 1,ndims
    c_start(ndims-i+1) = start(i) - 1
    c_ncount(ndims-i+1) = ncount(i)
  end do
  c_dp = dp
  ierr = nc_put_vara_double(c_ncid,c_varid,c_start,c_ncount,C_LOC(c_dp))
end function nf_put_vara_double_d4

integer function nf_put_vara_double_d5(ncid,varid,start,ncount,dp) result(ierr)
  implicit none
  integer, intent(in) :: ncid, varid
  integer, dimension(:), intent(in) :: start, ncount
  real(kind=8), dimension(:,:,:,:,:), intent(in) :: dp
  integer (C_INT) :: c_ncid, c_varid
  integer (C_INT), target :: c_ndims
  integer (C_SIZE_T), dimension(size(start)) :: c_start, c_ncount
  real (C_DOUBLE), dimension(size(dp,1),size(dp,2),size(dp,3),size(dp,4),size(dp,5)), target :: c_dp
  integer i, ndims
  c_ncid = ncid
  c_varid = varid - 1
  c_ndims = size(start)
  ierr = nc_inq_varndims(c_ncid,c_varid,C_LOC(c_ndims))
  ndims = c_ndims
  do i = 1,ndims
    c_start(ndims-i+1) = start(i) - 1
    c_ncount(ndims-i+1) = ncount(i)
  end do
  c_dp = dp
  ierr = nc_put_vara_double(c_ncid,c_varid,c_start,c_ncount,C_LOC(c_dp))
end function nf_put_vara_double_d5

integer function nf_put_var_real_d1(ncid,varid,fp) result(ierr)
  implicit none
  integer, intent(in) :: ncid, varid
  real(kind=4), dimension(:), intent(in) :: fp
  integer (C_INT) :: c_ncid, c_varid
  real (C_FLOAT), dimension(size(fp)), target :: c_fp
  c_ncid = ncid
  c_varid = varid - 1
  c_fp = fp
  ierr = nc_put_var_float(c_ncid,c_varid,C_LOC(c_fp))
end function nf_put_var_real_d1

integer function nf_put_var_real_d2(ncid,varid,fp) result(ierr)
  implicit none
  integer, intent(in) :: ncid, varid
  real(kind=4), dimension(:,:), intent(in) :: fp
  integer (C_INT) :: c_ncid, c_varid
  real (C_FLOAT), dimension(size(fp,1),size(fp,2)), target :: c_fp
  c_ncid = ncid
  c_varid = varid - 1
  c_fp = fp
  ierr = nc_put_var_float(c_ncid,c_varid,C_LOC(c_fp))
end function nf_put_var_real_d2

integer function nf_put_var_real_d3(ncid,varid,fp) result(ierr)
  implicit none
  integer, intent(in) :: ncid, varid
  real(kind=4), dimension(:,:,:), intent(in) :: fp
  integer (C_INT) :: c_ncid, c_varid
  real (C_FLOAT), dimension(size(fp,1),size(fp,2),size(fp,3)), target :: c_fp
  c_ncid = ncid
  c_varid = varid - 1
  c_fp = fp
  ierr = nc_put_var_float(c_ncid,c_varid,C_LOC(c_fp))
end function nf_put_var_real_d3

integer function nf_put_var_real_d4(ncid,varid,fp) result(ierr)
  implicit none
  integer, intent(in) :: ncid, varid
  real(kind=4), dimension(:,:,:,:), intent(in) :: fp
  integer (C_INT) :: c_ncid, c_varid
  real (C_FLOAT), dimension(size(fp,1),size(fp,2),size(fp,3),size(fp,4)), target :: c_fp
  c_ncid = ncid
  c_varid = varid - 1
  c_fp = fp
  ierr = nc_put_var_float(c_ncid,c_varid,C_LOC(c_fp))
end function nf_put_var_real_d4

integer function nf_put_var_real_d5(ncid,varid,fp) result(ierr)
  implicit none
  integer, intent(in) :: ncid, varid
  real(kind=4), dimension(:,:,:,:,:), intent(in) :: fp
  integer (C_INT) :: c_ncid, c_varid
  real (C_FLOAT), dimension(size(fp,1),size(fp,2),size(fp,3),size(fp,4),size(fp,5)), &
      target :: c_fp
  c_ncid = ncid
  c_varid = varid - 1
  c_fp = fp
  ierr = nc_put_var_float(c_ncid,c_varid,C_LOC(c_fp))
end function nf_put_var_real_d5

integer function nf_put_var1_int(ncid,varid,start,ip) result(ierr)
  implicit none
  integer, intent(in) :: ncid, varid
  integer, dimension(:), intent(in) :: start
  integer(kind=4), intent(in) :: ip
  integer (C_INT) :: c_ncid, c_varid
  integer (C_INT), target :: c_ndims
  integer (C_SIZE_T), dimension(size(start)) :: c_start
  integer (C_INT), target :: c_ip
  integer i, ndims
  c_ncid = ncid
  c_varid = varid - 1
  c_ndims = size(start)
  ierr = nc_inq_varndims(c_ncid,c_varid,C_LOC(c_ndims))
  ndims = c_ndims
  do i = 1,ndims
    c_start(ndims-i+1) = start(i) - 1
  end do
  c_ip = ip
  ierr = nc_put_var1_int(c_ncid,c_varid,c_start,C_LOC(c_ip))
end function nf_put_var1_int

integer function nf_put_var1_real(ncid,varid,start,rp) result(ierr)
  implicit none
  integer, intent(in) :: ncid, varid
  integer, dimension(:), intent(in) :: start
  real(kind=4), intent(in) :: rp
  integer (C_INT) :: c_ncid, c_varid
  integer (C_INT), target :: c_ndims
  integer (C_SIZE_T), dimension(size(start)) :: c_start
  real (C_FLOAT), target :: c_rp
  integer i, ndims
  c_ncid = ncid
  c_varid = varid - 1
  c_ndims = size(start)
  ierr = nc_inq_varndims(c_ncid,c_varid,C_LOC(c_ndims))
  ndims = c_ndims
  do i = 1,ndims
    c_start(ndims-i+1) = start(i) - 1
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
  integer (C_INT), target :: c_ndims
  integer (C_SIZE_T), dimension(size(start)) :: c_start
  real (C_DOUBLE), target :: c_dp
  integer i, ndims
  c_ncid = ncid
  c_varid = varid - 1
  c_ndims = size(start)
  ierr = nc_inq_varndims(c_ncid,c_varid,C_LOC(c_ndims))
  ndims = c_ndims
  do i = 1,ndims
    c_start(ndims-i+1) = start(i) - 1
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
  integer (C_INT), target :: c_ndims
  integer (C_SIZE_T), dimension(size(start)) :: c_start
  character, dimension(len(tp)) :: c_tp
  integer i, ndims
  c_ncid = ncid
  c_varid = varid - 1
  c_ndims = size(start)
  ierr = nc_inq_varndims(c_ncid,c_varid,C_LOC(c_ndims))
  ndims = c_ndims
  do i = 1,ndims
    c_start(ndims-i+1) = start(i) - 1
  end do
  call cf_strcopy(tp,c_tp)
  ierr = nc_put_var1_text(c_ncid,c_varid,c_start,c_tp)
end function nf_put_var1_text

integer function nf_put_var1_int1(ncid,varid,start,bp) result(ierr)
  implicit none
  integer, intent(in) :: ncid, varid
  integer, dimension(:), intent(in) :: start
  integer(kind=1), intent(in) :: bp
  integer (C_INT) :: c_ncid, c_varid
  integer (C_INT), target :: c_ndims
  integer (C_SIZE_T), dimension(size(start)) :: c_start
  integer (C_SIGNED_CHAR), target :: c_bp
  integer i, ndims
  c_ncid = ncid
  c_varid = varid - 1
  c_ndims = size(start)
  ierr = nc_inq_varndims(c_ncid,c_varid,C_LOC(c_ndims))
  ndims = c_ndims
  do i = 1,ndims
    c_start(ndims-i+1) = start(i) - 1
  end do
  c_bp = bp
  ierr = nc_put_var1_schar(c_ncid,c_varid,c_start,C_LOC(c_bp))
end function nf_put_var1_int1

integer function nf_put_var1_int2(ncid,varid,start,sp) result(ierr)
  implicit none
  integer, intent(in) :: ncid, varid
  integer, dimension(:), intent(in) :: start
  integer(kind=2), intent(in) :: sp
  integer (C_INT) :: c_ncid, c_varid
  integer (C_INT), target :: c_ndims
  integer (C_SIZE_T), dimension(size(start)) :: c_start
  integer (C_SHORT), target :: c_sp
  integer i, ndims
  c_ncid = ncid
  c_varid = varid - 1
  c_ndims = size(start)
  ierr = nc_inq_varndims(c_ncid,c_varid,C_LOC(c_ndims))
  ndims = c_ndims
  do i = 1,ndims
    c_start(ndims-i+1) = start(i) - 1
  end do
  c_sp = sp
  ierr = nc_put_var1_short(c_ncid,c_varid,c_start,C_LOC(c_sp))
end function nf_put_var1_int2

integer function nf_put_vars_text(ncid,varid,start,ncount,stride,tp) result(ierr)
  implicit none
  integer, intent(in) :: ncid, varid
  integer, dimension(:), intent(in) :: start, ncount, stride  
  character(len=*), intent(in) :: tp
  integer (C_INT) :: c_ncid, c_varid
  integer (C_INT), target :: c_ndims
  integer (C_SIZE_T), dimension(size(start)) :: c_start, c_ncount
  integer (C_INTPTR_T), dimension(size(stride)) :: c_stride
  character, dimension(len(tp)) :: c_tp
  integer i, ndims
  c_ncid = ncid
  c_varid = varid - 1
  c_ndims = size(start)
  ierr = nc_inq_varndims(c_ncid,c_varid,C_LOC(c_ndims))
  ndims = c_ndims
  do i = 1,ndims
    c_start(ndims-i+1) = start(i) - 1
    c_ncount(ndims-i+1) = ncount(i)
    c_stride(ndims-i+1) = stride(i)
  end do
  call cf_strcopy(tp,c_tp)
  ierr = nc_put_vars_text(c_ncid,c_varid,c_start,c_ncount,c_stride,c_tp)
end function nf_put_vars_text

integer function nf_put_vars_int1_d1(ncid,varid,start,ncount,stride,bp) result(ierr)
  implicit none
  integer, intent(in) :: ncid, varid
  integer, dimension(:), intent(in) :: start, ncount, stride
  integer(kind=1), dimension(:), intent(in) :: bp
  integer (C_INT) :: c_ncid, c_varid
  integer (C_INT), target :: c_ndims
  integer (C_SIZE_T), dimension(size(start)) :: c_start
  integer (C_SIZE_T), dimension(size(ncount)) :: c_ncount
  integer (C_INTPTR_T), dimension(size(stride)) :: c_stride
  integer (C_SIGNED_CHAR), dimension(size(bp)), target :: c_bp
  integer i, ndims
  c_ncid = ncid
  c_varid = varid - 1
  c_ndims = size(start)
  ierr = nc_inq_varndims(c_ncid,c_varid,C_LOC(c_ndims))
  ndims = c_ndims
  do i = 1,ndims
    c_start(ndims-i+1) = start(i) - 1
    c_ncount(ndims-i+1) = ncount(i)
    c_stride(ndims-i+1) = stride(i)
  end do
  c_bp = bp
  ierr = nc_put_vars_schar(c_ncid,c_varid,c_start,c_ncount,c_stride,C_LOC(c_bp))
end function nf_put_vars_int1_d1

integer function nf_put_vars_int1_d2(ncid,varid,start,ncount,stride,bp) result(ierr)
  implicit none
  integer, intent(in) :: ncid, varid
  integer, dimension(:), intent(in) :: start, ncount, stride
  integer(kind=1), dimension(:,:), intent(in) :: bp
  integer (C_INT) :: c_ncid, c_varid
  integer (C_INT), target :: c_ndims
  integer (C_SIZE_T), dimension(size(start)) :: c_start
  integer (C_SIZE_T), dimension(size(ncount)) :: c_ncount
  integer (C_INTPTR_T), dimension(size(stride)) :: c_stride
  integer (C_SIGNED_CHAR), dimension(size(bp,1),size(bp,2)), target :: c_bp
  integer i, ndims
  c_ncid = ncid
  c_varid = varid - 1
  c_ndims = size(start)
  ierr = nc_inq_varndims(c_ncid,c_varid,C_LOC(c_ndims))
  ndims = c_ndims
  do i = 1,ndims
    c_start(ndims-i+1) = start(i) - 1
    c_ncount(ndims-i+1) = ncount(i)
    c_stride(ndims-i+1) = stride(i)
  end do
  c_bp = bp
  ierr = nc_put_vars_schar(c_ncid,c_varid,c_start,c_ncount,c_stride,C_LOC(c_bp))
end function nf_put_vars_int1_d2

integer function nf_put_vars_int1_d3(ncid,varid,start,ncount,stride,bp) result(ierr)
  implicit none
  integer, intent(in) :: ncid, varid
  integer, dimension(:), intent(in) :: start, ncount, stride
  integer(kind=1), dimension(:,:,:), intent(in) :: bp
  integer (C_INT) :: c_ncid, c_varid
  integer (C_INT), target :: c_ndims
  integer (C_SIZE_T), dimension(size(start)) :: c_start
  integer (C_SIZE_T), dimension(size(ncount)) :: c_ncount
  integer (C_INTPTR_T), dimension(size(stride)) :: c_stride
  integer (C_SIGNED_CHAR), dimension(size(bp,1),size(bp,2),size(bp,3)), target :: c_bp
  integer i, ndims
  c_ncid = ncid
  c_varid = varid - 1
  c_ndims = size(start)
  ierr = nc_inq_varndims(c_ncid,c_varid,C_LOC(c_ndims))
  ndims = c_ndims
  do i = 1,ndims
    c_start(ndims-i+1) = start(i) - 1
    c_ncount(ndims-i+1) = ncount(i)
    c_stride(ndims-i+1) = stride(i)
  end do
  c_bp = bp
  ierr = nc_put_vars_schar(c_ncid,c_varid,c_start,c_ncount,c_stride,C_LOC(c_bp))
end function nf_put_vars_int1_d3

integer function nf_put_vars_int1_d4(ncid,varid,start,ncount,stride,bp) result(ierr)
  implicit none
  integer, intent(in) :: ncid, varid
  integer, dimension(:), intent(in) :: start, ncount, stride
  integer(kind=1), dimension(:,:,:,:), intent(in) :: bp
  integer (C_INT) :: c_ncid, c_varid
  integer (C_INT), target :: c_ndims
  integer (C_SIZE_T), dimension(size(start)) :: c_start
  integer (C_SIZE_T), dimension(size(ncount)) :: c_ncount
  integer (C_INTPTR_T), dimension(size(stride)) :: c_stride
  integer (C_SIGNED_CHAR), dimension(size(bp,1),size(bp,2),size(bp,3),size(bp,4)), target :: c_bp
  integer i, ndims
  c_ncid = ncid
  c_varid = varid - 1
  c_ndims = size(start)
  ierr = nc_inq_varndims(c_ncid,c_varid,C_LOC(c_ndims))
  ndims = c_ndims
  do i = 1,ndims
    c_start(ndims-i+1) = start(i) - 1
    c_ncount(ndims-i+1) = ncount(i)
    c_stride(ndims-i+1) = stride(i)
  end do
  c_bp = bp
  ierr = nc_put_vars_schar(c_ncid,c_varid,c_start,c_ncount,c_stride,C_LOC(c_bp))
end function nf_put_vars_int1_d4

integer function nf_put_vars_int1_d5(ncid,varid,start,ncount,stride,bp) result(ierr)
  implicit none
  integer, intent(in) :: ncid, varid
  integer, dimension(:), intent(in) :: start, ncount, stride
  integer(kind=1), dimension(:,:,:,:,:), intent(in) :: bp
  integer (C_INT) :: c_ncid, c_varid
  integer (C_INT), target :: c_ndims
  integer (C_SIZE_T), dimension(size(start)) :: c_start
  integer (C_SIZE_T), dimension(size(ncount)) :: c_ncount
  integer (C_INTPTR_T), dimension(size(stride)) :: c_stride
  integer (C_SIGNED_CHAR), dimension(size(bp,1),size(bp,2),size(bp,3),size(bp,4),size(bp,5)), &
      target :: c_bp
  integer i, ndims
  c_ncid = ncid
  c_varid = varid - 1
  c_ndims = size(start)
  ierr = nc_inq_varndims(c_ncid,c_varid,C_LOC(c_ndims))
  ndims = c_ndims
  do i = 1,ndims
    c_start(ndims-i+1) = start(i) - 1
    c_ncount(ndims-i+1) = ncount(i)
    c_stride(ndims-i+1) = stride(i)
  end do
  c_bp = bp
  ierr = nc_put_vars_schar(c_ncid,c_varid,c_start,c_ncount,c_stride,C_LOC(c_bp))
end function nf_put_vars_int1_d5

integer function nf_put_vars_int1_d6(ncid,varid,start,ncount,stride,bp) result(ierr)
  implicit none
  integer, intent(in) :: ncid, varid
  integer, dimension(:), intent(in) :: start, ncount, stride
  integer(kind=1), dimension(:,:,:,:,:,:), intent(in) :: bp
  integer (C_INT) :: c_ncid, c_varid
  integer (C_INT), target :: c_ndims
  integer (C_SIZE_T), dimension(size(start)) :: c_start
  integer (C_SIZE_T), dimension(size(ncount)) :: c_ncount
  integer (C_INTPTR_T), dimension(size(stride)) :: c_stride
  integer (C_SIGNED_CHAR), dimension(size(bp,1),size(bp,2),size(bp,3),size(bp,4),size(bp,5), &
      size(bp,6)), target :: c_bp
  integer i, ndims
  c_ncid = ncid
  c_varid = varid - 1
  c_ndims = size(start)
  ierr = nc_inq_varndims(c_ncid,c_varid,C_LOC(c_ndims))
  ndims = c_ndims
  do i = 1,ndims
    c_start(ndims-i+1) = start(i) - 1
    c_ncount(ndims-i+1) = ncount(i)
    c_stride(ndims-i+1) = stride(i)
  end do
  c_bp = bp
  ierr = nc_put_vars_schar(c_ncid,c_varid,c_start,c_ncount,c_stride,C_LOC(c_bp))
end function nf_put_vars_int1_d6

integer function nf_put_vars_int1_d7(ncid,varid,start,ncount,stride,bp) result(ierr)
  implicit none
  integer, intent(in) :: ncid, varid
  integer, dimension(:), intent(in) :: start, ncount, stride
  integer(kind=1), dimension(:,:,:,:,:,:,:), intent(in) :: bp
  integer (C_INT) :: c_ncid, c_varid
  integer (C_INT), target :: c_ndims
  integer (C_SIZE_T), dimension(size(start)) :: c_start
  integer (C_SIZE_T), dimension(size(ncount)) :: c_ncount
  integer (C_INTPTR_T), dimension(size(stride)) :: c_stride
  integer (C_SIGNED_CHAR), dimension(size(bp,1),size(bp,2),size(bp,3),size(bp,4),size(bp,5), &
      size(bp,6),size(bp,6)), target :: c_bp
  integer i, ndims
  c_ncid = ncid
  c_varid = varid - 1
  c_ndims = size(start)
  ierr = nc_inq_varndims(c_ncid,c_varid,C_LOC(c_ndims))
  ndims = c_ndims
  do i = 1,ndims
    c_start(ndims-i+1) = start(i) - 1
    c_ncount(ndims-i+1) = ncount(i)
    c_stride(ndims-i+1) = stride(i)
  end do
  c_bp = bp
  ierr = nc_put_vars_schar(c_ncid,c_varid,c_start,c_ncount,c_stride,C_LOC(c_bp))
end function nf_put_vars_int1_d7

integer function nf_put_vars_int2_d1(ncid,varid,start,ncount,stride,sp) result(ierr)
  implicit none
  integer, intent(in) :: ncid, varid
  integer, dimension(:), intent(in) :: start, ncount, stride
  integer(kind=2), dimension(:), intent(in) :: sp
  integer (C_INT) :: c_ncid, c_varid
  integer (C_INT), target :: c_ndims
  integer (C_SIZE_T), dimension(size(start)) :: c_start
  integer (C_SIZE_T), dimension(size(ncount)) :: c_ncount
  integer (C_INTPTR_T), dimension(size(stride)) :: c_stride
  integer (C_SHORT), dimension(size(sp)), target :: c_sp
  integer i, ndims
  c_ncid = ncid
  c_varid = varid - 1
  c_ndims = size(start)
  ierr = nc_inq_varndims(c_ncid,c_varid,C_LOC(c_ndims))
  ndims = c_ndims
  do i = 1,ndims
    c_start(ndims-i+1) = start(i) - 1
    c_ncount(ndims-i+1) = ncount(i)
    c_stride(ndims-i+1) = stride(i)
  end do
  c_sp = sp
  ierr = nc_put_vars_short(c_ncid,c_varid,c_start,c_ncount,c_stride,C_LOC(c_sp))
end function nf_put_vars_int2_d1

integer function nf_put_vars_int2_d2(ncid,varid,start,ncount,stride,sp) result(ierr)
  implicit none
  integer, intent(in) :: ncid, varid
  integer, dimension(:), intent(in) :: start, ncount, stride
  integer(kind=2), dimension(:,:), intent(in) :: sp
  integer (C_INT) :: c_ncid, c_varid
  integer (C_INT), target :: c_ndims
  integer (C_SIZE_T), dimension(size(start)) :: c_start
  integer (C_SIZE_T), dimension(size(ncount)) :: c_ncount
  integer (C_INTPTR_T), dimension(size(stride)) :: c_stride
  integer (C_SHORT), dimension(size(sp,1),size(sp,2)), target :: c_sp
  integer i, ndims
  c_ncid = ncid
  c_varid = varid - 1
  c_ndims = size(start)
  ierr = nc_inq_varndims(c_ncid,c_varid,C_LOC(c_ndims))
  ndims = c_ndims
  do i = 1,ndims
    c_start(ndims-i+1) = start(i) - 1
    c_ncount(ndims-i+1) = ncount(i)
    c_stride(ndims-i+1) = stride(i)
  end do
  c_sp = sp
  ierr = nc_put_vars_short(c_ncid,c_varid,c_start,c_ncount,c_stride,C_LOC(c_sp))
end function nf_put_vars_int2_d2

integer function nf_put_vars_int2_d3(ncid,varid,start,ncount,stride,sp) result(ierr)
  implicit none
  integer, intent(in) :: ncid, varid
  integer, dimension(:), intent(in) :: start, ncount, stride
  integer(kind=2), dimension(:,:,:), intent(in) :: sp
  integer (C_INT) :: c_ncid, c_varid
  integer (C_INT), target :: c_ndims
  integer (C_SIZE_T), dimension(size(start)) :: c_start
  integer (C_SIZE_T), dimension(size(ncount)) :: c_ncount
  integer (C_INTPTR_T), dimension(size(stride)) :: c_stride
  integer (C_SHORT), dimension(size(sp,1),size(sp,2),size(sp,3)), target :: c_sp
  integer i, ndims
  c_ncid = ncid
  c_varid = varid - 1
  c_ndims = size(start)
  ierr = nc_inq_varndims(c_ncid,c_varid,C_LOC(c_ndims))
  ndims = c_ndims
  do i = 1,ndims
    c_start(ndims-i+1) = start(i) - 1
    c_ncount(ndims-i+1) = ncount(i)
    c_stride(ndims-i+1) = stride(i)
  end do
  c_sp = sp
  ierr = nc_put_vars_short(c_ncid,c_varid,c_start,c_ncount,c_stride,C_LOC(c_sp))
end function nf_put_vars_int2_d3

integer function nf_put_vars_int2_d4(ncid,varid,start,ncount,stride,sp) result(ierr)
  implicit none
  integer, intent(in) :: ncid, varid
  integer, dimension(:), intent(in) :: start, ncount, stride
  integer(kind=2), dimension(:,:,:,:), intent(in) :: sp
  integer (C_INT) :: c_ncid, c_varid
  integer (C_INT), target :: c_ndims
  integer (C_SIZE_T), dimension(size(start)) :: c_start
  integer (C_SIZE_T), dimension(size(ncount)) :: c_ncount
  integer (C_INTPTR_T), dimension(size(stride)) :: c_stride
  integer (C_SHORT), dimension(size(sp,1),size(sp,2),size(sp,3),size(sp,4)), target :: c_sp
  integer i, ndims
  c_ncid = ncid
  c_varid = varid - 1
  c_ndims = size(start)
  ierr = nc_inq_varndims(c_ncid,c_varid,C_LOC(c_ndims))
  ndims = c_ndims
  do i = 1,ndims
    c_start(ndims-i+1) = start(i) - 1
    c_ncount(ndims-i+1) = ncount(i)
    c_stride(ndims-i+1) = stride(i)
  end do
  c_sp = sp
  ierr = nc_put_vars_short(c_ncid,c_varid,c_start,c_ncount,c_stride,C_LOC(c_sp))
end function nf_put_vars_int2_d4

integer function nf_put_vars_int2_d5(ncid,varid,start,ncount,stride,sp) result(ierr)
  implicit none
  integer, intent(in) :: ncid, varid
  integer, dimension(:), intent(in) :: start, ncount, stride
  integer(kind=2), dimension(:,:,:,:,:), intent(in) :: sp
  integer (C_INT) :: c_ncid, c_varid
  integer (C_INT), target :: c_ndims
  integer (C_SIZE_T), dimension(size(start)) :: c_start
  integer (C_SIZE_T), dimension(size(ncount)) :: c_ncount
  integer (C_INTPTR_T), dimension(size(stride)) :: c_stride
  integer (C_SHORT), dimension(size(sp,1),size(sp,2),size(sp,3),size(sp,4),size(sp,5)), &
      target :: c_sp
  integer i, ndims
  c_ncid = ncid
  c_varid = varid - 1
  c_ndims = size(start)
  ierr = nc_inq_varndims(c_ncid,c_varid,C_LOC(c_ndims))
  ndims = c_ndims
  do i = 1,ndims
    c_start(ndims-i+1) = start(i) - 1
    c_ncount(ndims-i+1) = ncount(i)
    c_stride(ndims-i+1) = stride(i)
  end do
  c_sp = sp
  ierr = nc_put_vars_short(c_ncid,c_varid,c_start,c_ncount,c_stride,C_LOC(c_sp))
end function nf_put_vars_int2_d5

integer function nf_put_vars_int2_d6(ncid,varid,start,ncount,stride,sp) result(ierr)
  implicit none
  integer, intent(in) :: ncid, varid
  integer, dimension(:), intent(in) :: start, ncount, stride
  integer(kind=2), dimension(:,:,:,:,:,:), intent(in) :: sp
  integer (C_INT) :: c_ncid, c_varid
  integer (C_INT), target :: c_ndims
  integer (C_SIZE_T), dimension(size(start)) :: c_start
  integer (C_SIZE_T), dimension(size(ncount)) :: c_ncount
  integer (C_INTPTR_T), dimension(size(stride)) :: c_stride
  integer (C_SHORT), dimension(size(sp,1),size(sp,2),size(sp,3),size(sp,4),size(sp,5), &
      size(sp,6)), target :: c_sp
  integer i, ndims
  c_ncid = ncid
  c_varid = varid - 1
  c_ndims = size(start)
  ierr = nc_inq_varndims(c_ncid,c_varid,C_LOC(c_ndims))
  ndims = c_ndims
  do i = 1,ndims
    c_start(ndims-i+1) = start(i) - 1
    c_ncount(ndims-i+1) = ncount(i)
    c_stride(ndims-i+1) = stride(i)
  end do
  c_sp = sp
  ierr = nc_put_vars_short(c_ncid,c_varid,c_start,c_ncount,c_stride,C_LOC(c_sp))
end function nf_put_vars_int2_d6

integer function nf_put_vars_int2_d7(ncid,varid,start,ncount,stride,sp) result(ierr)
  implicit none
  integer, intent(in) :: ncid, varid
  integer, dimension(:), intent(in) :: start, ncount, stride
  integer(kind=2), dimension(:,:,:,:,:,:,:), intent(in) :: sp
  integer (C_INT) :: c_ncid, c_varid
  integer (C_INT), target :: c_ndims
  integer (C_SIZE_T), dimension(size(start)) :: c_start
  integer (C_SIZE_T), dimension(size(ncount)) :: c_ncount
  integer (C_INTPTR_T), dimension(size(stride)) :: c_stride
  integer (C_SHORT), dimension(size(sp,1),size(sp,2),size(sp,3),size(sp,4),size(sp,5), &
      size(sp,6),size(sp,7)), target :: c_sp
  integer i, ndims
  c_ncid = ncid
  c_varid = varid - 1
  c_ndims = size(start)
  ierr = nc_inq_varndims(c_ncid,c_varid,C_LOC(c_ndims))
  ndims = c_ndims
  do i = 1,ndims
    c_start(ndims-i+1) = start(i) - 1
    c_ncount(ndims-i+1) = ncount(i)
    c_stride(ndims-i+1) = stride(i)
  end do
  c_sp = sp
  ierr = nc_put_vars_short(c_ncid,c_varid,c_start,c_ncount,c_stride,C_LOC(c_sp))
end function nf_put_vars_int2_d7

integer function nf_put_vars_int_d1(ncid,varid,start,ncount,stride,ip) result(ierr)
  implicit none
  integer, intent(in) :: ncid, varid
  integer, dimension(:), intent(in) :: start, ncount, stride
  integer(kind=4), dimension(:), intent(in) :: ip
  integer (C_INT) :: c_ncid, c_varid
  integer (C_INT), target :: c_ndims
  integer (C_SIZE_T), dimension(size(start)) :: c_start
  integer (C_SIZE_T), dimension(size(ncount)) :: c_ncount
  integer (C_INTPTR_T), dimension(size(stride)) :: c_stride
  integer (C_INT), dimension(size(ip)), target :: c_ip
  integer i, ndims
  c_ncid = ncid
  c_varid = varid - 1
  c_ndims = size(start)
  ierr = nc_inq_varndims(c_ncid,c_varid,C_LOC(c_ndims))
  ndims = c_ndims
  do i = 1,ndims
    c_start(ndims-i+1) = start(i) - 1
    c_ncount(ndims-i+1) = ncount(i)
    c_stride(ndims-i+1) = stride(i)
  end do
  c_ip = ip
  ierr = nc_put_vars_int(c_ncid,c_varid,c_start,c_ncount,c_stride,C_LOC(c_ip))
end function nf_put_vars_int_d1

integer function nf_put_vars_int_d2(ncid,varid,start,ncount,stride,ip) result(ierr)
  implicit none
  integer, intent(in) :: ncid, varid
  integer, dimension(:), intent(in) :: start, ncount, stride
  integer(kind=4), dimension(:,:), intent(in) :: ip
  integer (C_INT) :: c_ncid, c_varid
  integer (C_INT), target :: c_ndims
  integer (C_SIZE_T), dimension(size(start)) :: c_start
  integer (C_SIZE_T), dimension(size(ncount)) :: c_ncount
  integer (C_INTPTR_T), dimension(size(stride)) :: c_stride
  integer (C_INT), dimension(size(ip,1),size(ip,2)), target :: c_ip
  integer i, ndims
  c_ncid = ncid
  c_varid = varid - 1
  c_ndims = size(start)
  ierr = nc_inq_varndims(c_ncid,c_varid,C_LOC(c_ndims))
  ndims = c_ndims
  do i = 1,ndims
    c_start(ndims-i+1) = start(i) - 1
    c_ncount(ndims-i+1) = ncount(i)
    c_stride(ndims-i+1) = stride(i)
  end do
  c_ip = ip
  ierr = nc_put_vars_int(c_ncid,c_varid,c_start,c_ncount,c_stride,C_LOC(c_ip))
end function nf_put_vars_int_d2

integer function nf_put_vars_int_d3(ncid,varid,start,ncount,stride,ip) result(ierr)
  implicit none
  integer, intent(in) :: ncid, varid
  integer, dimension(:), intent(in) :: start, ncount, stride
  integer(kind=4), dimension(:,:,:), intent(in) :: ip
  integer (C_INT) :: c_ncid, c_varid
  integer (C_INT), target :: c_ndims
  integer (C_SIZE_T), dimension(size(start)) :: c_start
  integer (C_SIZE_T), dimension(size(ncount)) :: c_ncount
  integer (C_INTPTR_T), dimension(size(stride)) :: c_stride
  integer (C_INT), dimension(size(ip,1),size(ip,2),size(ip,3)), target :: c_ip
  integer i, ndims
  c_ncid = ncid
  c_varid = varid - 1
  c_ndims = size(start)
  ierr = nc_inq_varndims(c_ncid,c_varid,C_LOC(c_ndims))
  ndims = c_ndims
  do i = 1,ndims
    c_start(ndims-i+1) = start(i) - 1
    c_ncount(ndims-i+1) = ncount(i)
    c_stride(ndims-i+1) = stride(i)
  end do
  c_ip = ip
  ierr = nc_put_vars_int(c_ncid,c_varid,c_start,c_ncount,c_stride,C_LOC(c_ip))
end function nf_put_vars_int_d3

integer function nf_put_vars_int_d4(ncid,varid,start,ncount,stride,ip) result(ierr)
  implicit none
  integer, intent(in) :: ncid, varid
  integer, dimension(:), intent(in) :: start, ncount, stride
  integer(kind=4), dimension(:,:,:,:), intent(in) :: ip
  integer (C_INT) :: c_ncid, c_varid
  integer (C_INT), target :: c_ndims
  integer (C_SIZE_T), dimension(size(start)) :: c_start
  integer (C_SIZE_T), dimension(size(ncount)) :: c_ncount
  integer (C_INTPTR_T), dimension(size(stride)) :: c_stride
  integer (C_INT), dimension(size(ip,1),size(ip,2),size(ip,3),size(ip,4)), target :: c_ip
  integer i, ndims
  c_ncid = ncid
  c_varid = varid - 1
  c_ndims = size(start)
  ierr = nc_inq_varndims(c_ncid,c_varid,C_LOC(c_ndims))
  ndims = c_ndims
  do i = 1,ndims
    c_start(ndims-i+1) = start(i) - 1
    c_ncount(ndims-i+1) = ncount(i)
    c_stride(ndims-i+1) = stride(i)
  end do
  c_ip = ip
  ierr = nc_put_vars_int(c_ncid,c_varid,c_start,c_ncount,c_stride,C_LOC(c_ip))
end function nf_put_vars_int_d4

integer function nf_put_vars_int_d5(ncid,varid,start,ncount,stride,ip) result(ierr)
  implicit none
  integer, intent(in) :: ncid, varid
  integer, dimension(:), intent(in) :: start, ncount, stride
  integer(kind=4), dimension(:,:,:,:,:), intent(in) :: ip
  integer (C_INT) :: c_ncid, c_varid
  integer (C_INT), target :: c_ndims
  integer (C_SIZE_T), dimension(size(start)) :: c_start
  integer (C_SIZE_T), dimension(size(ncount)) :: c_ncount
  integer (C_INTPTR_T), dimension(size(stride)) :: c_stride
  integer (C_INT), dimension(size(ip,1),size(ip,2),size(ip,3),size(ip,4),size(ip,5)), &
      target :: c_ip
  integer i, ndims
  c_ncid = ncid
  c_varid = varid - 1
  c_ndims = size(start)
  ierr = nc_inq_varndims(c_ncid,c_varid,C_LOC(c_ndims))
  ndims = c_ndims
  do i = 1,ndims
    c_start(ndims-i+1) = start(i) - 1
    c_ncount(ndims-i+1) = ncount(i)
    c_stride(ndims-i+1) = stride(i)
  end do
  c_ip = ip
  ierr = nc_put_vars_int(c_ncid,c_varid,c_start,c_ncount,c_stride,C_LOC(c_ip))
end function nf_put_vars_int_d5

integer function nf_put_vars_int_d6(ncid,varid,start,ncount,stride,ip) result(ierr)
  implicit none
  integer, intent(in) :: ncid, varid
  integer, dimension(:), intent(in) :: start, ncount, stride
  integer(kind=4), dimension(:,:,:,:,:,:), intent(in) :: ip
  integer (C_INT) :: c_ncid, c_varid
  integer (C_INT), target :: c_ndims
  integer (C_SIZE_T), dimension(size(start)) :: c_start
  integer (C_SIZE_T), dimension(size(ncount)) :: c_ncount
  integer (C_INTPTR_T), dimension(size(stride)) :: c_stride
  integer (C_INT), dimension(size(ip,1),size(ip,2),size(ip,3),size(ip,4),size(ip,5), &
      size(ip,6)), target :: c_ip
  integer i, ndims
  c_ncid = ncid
  c_varid = varid - 1
  c_ndims = size(start)
  ierr = nc_inq_varndims(c_ncid,c_varid,C_LOC(c_ndims))
  ndims = c_ndims
  do i = 1,ndims
    c_start(ndims-i+1) = start(i) - 1
    c_ncount(ndims-i+1) = ncount(i)
    c_stride(ndims-i+1) = stride(i)
  end do
  c_ip = ip
  ierr = nc_put_vars_int(c_ncid,c_varid,c_start,c_ncount,c_stride,C_LOC(c_ip))
end function nf_put_vars_int_d6

integer function nf_put_vars_int_d7(ncid,varid,start,ncount,stride,ip) result(ierr)
  implicit none
  integer, intent(in) :: ncid, varid
  integer, dimension(:), intent(in) :: start, ncount, stride
  integer(kind=4), dimension(:,:,:,:,:,:,:), intent(in) :: ip
  integer (C_INT) :: c_ncid, c_varid
  integer (C_INT), target :: c_ndims
  integer (C_SIZE_T), dimension(size(start)) :: c_start
  integer (C_SIZE_T), dimension(size(ncount)) :: c_ncount
  integer (C_INTPTR_T), dimension(size(stride)) :: c_stride
  integer (C_INT), dimension(size(ip,1),size(ip,2),size(ip,3),size(ip,4),size(ip,5), &
      size(ip,6),size(ip,7)), target :: c_ip
  integer i, ndims
  c_ncid = ncid
  c_varid = varid - 1
  c_ndims = size(start)
  ierr = nc_inq_varndims(c_ncid,c_varid,C_LOC(c_ndims))
  ndims = c_ndims
  do i = 1,ndims
    c_start(ndims-i+1) = start(i) - 1
    c_ncount(ndims-i+1) = ncount(i)
    c_stride(ndims-i+1) = stride(i)
  end do
  c_ip = ip
  ierr = nc_put_vars_int(c_ncid,c_varid,c_start,c_ncount,c_stride,C_LOC(c_ip))
end function nf_put_vars_int_d7

integer function nf_put_vars_real_d1(ncid,varid,start,ncount,stride,fp) result(ierr)
  implicit none
  integer, intent(in) :: ncid, varid
  integer, dimension(:), intent(in) :: start, ncount, stride
  real(kind=4), dimension(:), intent(in) :: fp
  integer (C_INT) :: c_ncid, c_varid
  integer (C_INT), target :: c_ndims
  integer (C_SIZE_T), dimension(size(start)) :: c_start
  integer (C_SIZE_T), dimension(size(ncount)) :: c_ncount
  integer (C_INTPTR_T), dimension(size(stride)) :: c_stride
  real (C_FLOAT), dimension(size(fp)), target :: c_fp
  integer i, ndims
  c_ncid = ncid
  c_varid = varid - 1
  c_ndims = size(start)
  ierr = nc_inq_varndims(c_ncid,c_varid,C_LOC(c_ndims))
  ndims = c_ndims
  do i = 1,ndims
    c_start(ndims-i+1) = start(i) - 1
    c_ncount(ndims-i+1) = ncount(i)
    c_stride(ndims-i+1) = stride(i)
  end do
  c_fp = fp
  ierr = nc_put_vars_float(c_ncid,c_varid,c_start,c_ncount,c_stride,C_LOC(c_fp))
end function nf_put_vars_real_d1

integer function nf_put_vars_real_d2(ncid,varid,start,ncount,stride,fp) result(ierr)
  implicit none
  integer, intent(in) :: ncid, varid
  integer, dimension(:), intent(in) :: start, ncount, stride
  real(kind=4), dimension(:,:), intent(in) :: fp
  integer (C_INT) :: c_ncid, c_varid
  integer (C_INT), target :: c_ndims
  integer (C_SIZE_T), dimension(size(start)) :: c_start
  integer (C_SIZE_T), dimension(size(ncount)) :: c_ncount
  integer (C_INTPTR_T), dimension(size(stride)) :: c_stride
  real (C_FLOAT), dimension(size(fp,1),size(fp,2)), target :: c_fp
  integer i, ndims
  c_ncid = ncid
  c_varid = varid - 1
  c_ndims = size(start)
  ierr = nc_inq_varndims(c_ncid,c_varid,C_LOC(c_ndims))
  ndims = c_ndims
  do i = 1,ndims
    c_start(ndims-i+1) = start(i) - 1
    c_ncount(ndims-i+1) = ncount(i)
    c_stride(ndims-i+1) = stride(i)
  end do
  c_fp = fp
  ierr = nc_put_vars_float(c_ncid,c_varid,c_start,c_ncount,c_stride,C_LOC(c_fp))
end function nf_put_vars_real_d2

integer function nf_put_vars_real_d3(ncid,varid,start,ncount,stride,fp) result(ierr)
  implicit none
  integer, intent(in) :: ncid, varid
  integer, dimension(:), intent(in) :: start, ncount, stride
  real(kind=4), dimension(:,:,:), intent(in) :: fp
  integer (C_INT) :: c_ncid, c_varid
  integer (C_INT), target :: c_ndims
  integer (C_SIZE_T), dimension(size(start)) :: c_start
  integer (C_SIZE_T), dimension(size(ncount)) :: c_ncount
  integer (C_INTPTR_T), dimension(size(stride)) :: c_stride
  real (C_FLOAT), dimension(size(fp,1),size(fp,2),size(fp,3)), target :: c_fp
  integer i, ndims
  c_ncid = ncid
  c_varid = varid - 1
  c_ndims = size(start)
  ierr = nc_inq_varndims(c_ncid,c_varid,C_LOC(c_ndims))
  ndims = c_ndims
  do i = 1,ndims
    c_start(ndims-i+1) = start(i) - 1
    c_ncount(ndims-i+1) = ncount(i)
    c_stride(ndims-i+1) = stride(i)
  end do
  c_fp = fp
  ierr = nc_put_vars_float(c_ncid,c_varid,c_start,c_ncount,c_stride,C_LOC(c_fp))
end function nf_put_vars_real_d3

integer function nf_put_vars_real_d4(ncid,varid,start,ncount,stride,fp) result(ierr)
  implicit none
  integer, intent(in) :: ncid, varid
  integer, dimension(:), intent(in) :: start, ncount, stride
  real(kind=4), dimension(:,:,:,:), intent(in) :: fp
  integer (C_INT) :: c_ncid, c_varid
  integer (C_INT), target :: c_ndims
  integer (C_SIZE_T), dimension(size(start)) :: c_start
  integer (C_SIZE_T), dimension(size(ncount)) :: c_ncount
  integer (C_INTPTR_T), dimension(size(stride)) :: c_stride
  real (C_FLOAT), dimension(size(fp,1),size(fp,2),size(fp,3),size(fp,4)), target :: c_fp
  integer i, ndims
  c_ncid = ncid
  c_varid = varid - 1
  c_ndims = size(start)
  ierr = nc_inq_varndims(c_ncid,c_varid,C_LOC(c_ndims))
  ndims = c_ndims
  do i = 1,ndims
    c_start(ndims-i+1) = start(i) - 1
    c_ncount(ndims-i+1) = ncount(i)
    c_stride(ndims-i+1) = stride(i)
  end do
  c_fp = fp
  ierr = nc_put_vars_float(c_ncid,c_varid,c_start,c_ncount,c_stride,C_LOC(c_fp))
end function nf_put_vars_real_d4

integer function nf_put_vars_real_d5(ncid,varid,start,ncount,stride,fp) result(ierr)
  implicit none
  integer, intent(in) :: ncid, varid
  integer, dimension(:), intent(in) :: start, ncount, stride
  real(kind=4), dimension(:,:,:,:,:), intent(in) :: fp
  integer (C_INT) :: c_ncid, c_varid
  integer (C_INT), target :: c_ndims
  integer (C_SIZE_T), dimension(size(start)) :: c_start
  integer (C_SIZE_T), dimension(size(ncount)) :: c_ncount
  integer (C_INTPTR_T), dimension(size(stride)) :: c_stride
  real (C_FLOAT), dimension(size(fp,1),size(fp,2),size(fp,3),size(fp,4),size(fp,5)), &
      target :: c_fp
  integer i, ndims
  c_ncid = ncid
  c_varid = varid - 1
  c_ndims = size(start)
  ierr = nc_inq_varndims(c_ncid,c_varid,C_LOC(c_ndims))
  ndims = c_ndims
  do i = 1,ndims
    c_start(ndims-i+1) = start(i) - 1
    c_ncount(ndims-i+1) = ncount(i)
    c_stride(ndims-i+1) = stride(i)
  end do
  c_fp = fp
  ierr = nc_put_vars_float(c_ncid,c_varid,c_start,c_ncount,c_stride,C_LOC(c_fp))
end function nf_put_vars_real_d5

integer function nf_put_vars_real_d6(ncid,varid,start,ncount,stride,fp) result(ierr)
  implicit none
  integer, intent(in) :: ncid, varid
  integer, dimension(:), intent(in) :: start, ncount, stride
  real(kind=4), dimension(:,:,:,:,:,:), intent(in) :: fp
  integer (C_INT) :: c_ncid, c_varid
  integer (C_INT), target :: c_ndims
  integer (C_SIZE_T), dimension(size(start)) :: c_start
  integer (C_SIZE_T), dimension(size(ncount)) :: c_ncount
  integer (C_INTPTR_T), dimension(size(stride)) :: c_stride
  real (C_FLOAT), dimension(size(fp,1),size(fp,2),size(fp,3),size(fp,4),size(fp,5), &
      size(fp,6)), target :: c_fp
  integer i, ndims
  c_ncid = ncid
  c_varid = varid - 1
  c_ndims = size(start)
  ierr = nc_inq_varndims(c_ncid,c_varid,C_LOC(c_ndims))
  ndims = c_ndims
  do i = 1,ndims
    c_start(ndims-i+1) = start(i) - 1
    c_ncount(ndims-i+1) = ncount(i)
    c_stride(ndims-i+1) = stride(i)
  end do
  c_fp = fp
  ierr = nc_put_vars_float(c_ncid,c_varid,c_start,c_ncount,c_stride,C_LOC(c_fp))
end function nf_put_vars_real_d6

integer function nf_put_vars_real_d7(ncid,varid,start,ncount,stride,fp) result(ierr)
  implicit none
  integer, intent(in) :: ncid, varid
  integer, dimension(:), intent(in) :: start, ncount, stride
  real(kind=4), dimension(:,:,:,:,:,:,:), intent(in) :: fp
  integer (C_INT) :: c_ncid, c_varid
  integer (C_INT), target :: c_ndims
  integer (C_SIZE_T), dimension(size(start)) :: c_start
  integer (C_SIZE_T), dimension(size(ncount)) :: c_ncount
  integer (C_INTPTR_T), dimension(size(stride)) :: c_stride
  real (C_FLOAT), dimension(size(fp,1),size(fp,2),size(fp,3),size(fp,4),size(fp,5), &
      size(fp,6),size(fp,7)), target :: c_fp
  integer i, ndims
  c_ncid = ncid
  c_varid = varid - 1
  c_ndims = size(start)
  ierr = nc_inq_varndims(c_ncid,c_varid,C_LOC(c_ndims))
  ndims = c_ndims
  do i = 1,ndims
    c_start(ndims-i+1) = start(i) - 1
    c_ncount(ndims-i+1) = ncount(i)
    c_stride(ndims-i+1) = stride(i)
  end do
  c_fp = fp
  ierr = nc_put_vars_float(c_ncid,c_varid,c_start,c_ncount,c_stride,C_LOC(c_fp))
end function nf_put_vars_real_d7

integer function nf_put_vars_double_d1(ncid,varid,start,ncount,stride,dp) result(ierr)
  implicit none
  integer, intent(in) :: ncid, varid
  integer, dimension(:), intent(in) :: start, ncount, stride
  real(kind=8), dimension(:), intent(in) :: dp
  integer (C_INT) :: c_ncid, c_varid
  integer (C_INT), target :: c_ndims
  integer (C_SIZE_T), dimension(size(start)) :: c_start
  integer (C_SIZE_T), dimension(size(ncount)) :: c_ncount
  integer (C_INTPTR_T), dimension(size(stride)) :: c_stride
  real (C_DOUBLE), dimension(size(dp)), target :: c_dp
  integer i, ndims
  c_ncid = ncid
  c_varid = varid - 1
  c_ndims = size(start)
  ierr = nc_inq_varndims(c_ncid,c_varid,C_LOC(c_ndims))
  ndims = c_ndims
  do i = 1,ndims
    c_start(ndims-i+1) = start(i) - 1
    c_ncount(ndims-i+1) = ncount(i)
    c_stride(ndims-i+1) = stride(i)
  end do
  c_dp = dp
  ierr = nc_put_vars_double(c_ncid,c_varid,c_start,c_ncount,c_stride,C_LOC(c_dp))
end function nf_put_vars_double_d1

integer function nf_put_vars_double_d2(ncid,varid,start,ncount,stride,dp) result(ierr)
  implicit none
  integer, intent(in) :: ncid, varid
  integer, dimension(:), intent(in) :: start, ncount, stride
  real(kind=8), dimension(:,:), intent(in) :: dp
  integer (C_INT) :: c_ncid, c_varid
  integer (C_INT), target :: c_ndims
  integer (C_SIZE_T), dimension(size(start)) :: c_start
  integer (C_SIZE_T), dimension(size(ncount)) :: c_ncount
  integer (C_INTPTR_T), dimension(size(stride)) :: c_stride
  real (C_DOUBLE), dimension(size(dp,1),size(dp,2)), target :: c_dp
  integer i, ndims
  c_ncid = ncid
  c_varid = varid - 1
  c_ndims = size(start)
  ierr = nc_inq_varndims(c_ncid,c_varid,C_LOC(c_ndims))
  ndims = c_ndims
  do i = 1,ndims
    c_start(ndims-i+1) = start(i) - 1
    c_ncount(ndims-i+1) = ncount(i)
    c_stride(ndims-i+1) = stride(i)
  end do
  c_dp = dp
  ierr = nc_put_vars_double(c_ncid,c_varid,c_start,c_ncount,c_stride,C_LOC(c_dp))
end function nf_put_vars_double_d2

integer function nf_put_vars_double_d3(ncid,varid,start,ncount,stride,dp) result(ierr)
  implicit none
  integer, intent(in) :: ncid, varid
  integer, dimension(:), intent(in) :: start, ncount, stride
  real(kind=8), dimension(:,:,:), intent(in) :: dp
  integer (C_INT) :: c_ncid, c_varid
  integer (C_INT), target :: c_ndims
  integer (C_SIZE_T), dimension(size(start)) :: c_start
  integer (C_SIZE_T), dimension(size(ncount)) :: c_ncount
  integer (C_INTPTR_T), dimension(size(stride)) :: c_stride
  real (C_DOUBLE), dimension(size(dp,1),size(dp,2),size(dp,3)), target :: c_dp
  integer i, ndims
  c_ncid = ncid
  c_varid = varid - 1
  c_ndims = size(start)
  ierr = nc_inq_varndims(c_ncid,c_varid,C_LOC(c_ndims))
  ndims = c_ndims
  do i = 1,ndims
    c_start(ndims-i+1) = start(i) - 1
    c_ncount(ndims-i+1) = ncount(i)
    c_stride(ndims-i+1) = stride(i)
  end do
  c_dp = dp
  ierr = nc_put_vars_double(c_ncid,c_varid,c_start,c_ncount,c_stride,C_LOC(c_dp))
end function nf_put_vars_double_d3

integer function nf_put_vars_double_d4(ncid,varid,start,ncount,stride,dp) result(ierr)
  implicit none
  integer, intent(in) :: ncid, varid
  integer, dimension(:), intent(in) :: start, ncount, stride
  real(kind=8), dimension(:,:,:,:), intent(in) :: dp
  integer (C_INT) :: c_ncid, c_varid
  integer (C_INT), target :: c_ndims
  integer (C_SIZE_T), dimension(size(start)) :: c_start
  integer (C_SIZE_T), dimension(size(ncount)) :: c_ncount
  integer (C_INTPTR_T), dimension(size(stride)) :: c_stride
  real (C_DOUBLE), dimension(size(dp,1),size(dp,2),size(dp,3),size(dp,4)), target :: c_dp
  integer i, ndims
  c_ncid = ncid
  c_varid = varid - 1
  c_ndims = size(start)
  ierr = nc_inq_varndims(c_ncid,c_varid,C_LOC(c_ndims))
  ndims = c_ndims
  do i = 1,ndims
    c_start(ndims-i+1) = start(i) - 1
    c_ncount(ndims-i+1) = ncount(i)
    c_stride(ndims-i+1) = stride(i)
  end do
  c_dp = dp
  ierr = nc_put_vars_double(c_ncid,c_varid,c_start,c_ncount,c_stride,C_LOC(c_dp))
end function nf_put_vars_double_d4

integer function nf_put_vars_double_d5(ncid,varid,start,ncount,stride,dp) result(ierr)
  implicit none
  integer, intent(in) :: ncid, varid
  integer, dimension(:), intent(in) :: start, ncount, stride
  real(kind=8), dimension(:,:,:,:,:), intent(in) :: dp
  integer (C_INT) :: c_ncid, c_varid
  integer (C_INT), target :: c_ndims
  integer (C_SIZE_T), dimension(size(start)) :: c_start
  integer (C_SIZE_T), dimension(size(ncount)) :: c_ncount
  integer (C_INTPTR_T), dimension(size(stride)) :: c_stride
  real (C_DOUBLE), dimension(size(dp,1),size(dp,2),size(dp,3),size(dp,4),size(dp,5)), &
      target :: c_dp
  integer i, ndims
  c_ncid = ncid
  c_varid = varid - 1
  c_ndims = size(start)
  ierr = nc_inq_varndims(c_ncid,c_varid,C_LOC(c_ndims))
  ndims = c_ndims
  do i = 1,ndims
    c_start(ndims-i+1) = start(i) - 1
    c_ncount(ndims-i+1) = ncount(i)
    c_stride(ndims-i+1) = stride(i)
  end do
  c_dp = dp
  ierr = nc_put_vars_double(c_ncid,c_varid,c_start,c_ncount,c_stride,C_LOC(c_dp))
end function nf_put_vars_double_d5

integer function nf_put_vars_double_d6(ncid,varid,start,ncount,stride,dp) result(ierr)
  implicit none
  integer, intent(in) :: ncid, varid
  integer, dimension(:), intent(in) :: start, ncount, stride
  real(kind=8), dimension(:,:,:,:,:,:), intent(in) :: dp
  integer (C_INT) :: c_ncid, c_varid
  integer (C_INT), target :: c_ndims
  integer (C_SIZE_T), dimension(size(start)) :: c_start
  integer (C_SIZE_T), dimension(size(ncount)) :: c_ncount
  integer (C_INTPTR_T), dimension(size(stride)) :: c_stride
  real (C_DOUBLE), dimension(size(dp,1),size(dp,2),size(dp,3),size(dp,4),size(dp,5), &
      size(dp,6)), target :: c_dp
  integer i, ndims
  c_ncid = ncid
  c_varid = varid - 1
  c_ndims = size(start)
  ierr = nc_inq_varndims(c_ncid,c_varid,C_LOC(c_ndims))
  ndims = c_ndims
  do i = 1,ndims
    c_start(ndims-i+1) = start(i) - 1
    c_ncount(ndims-i+1) = ncount(i)
    c_stride(ndims-i+1) = stride(i)
  end do
  c_dp = dp
  ierr = nc_put_vars_double(c_ncid,c_varid,c_start,c_ncount,c_stride,C_LOC(c_dp))
end function nf_put_vars_double_d6

integer function nf_put_vars_double_d7(ncid,varid,start,ncount,stride,dp) result(ierr)
  implicit none
  integer, intent(in) :: ncid, varid
  integer, dimension(:), intent(in) :: start, ncount, stride
  real(kind=8), dimension(:,:,:,:,:,:,:), intent(in) :: dp
  integer (C_INT) :: c_ncid, c_varid
  integer (C_INT), target :: c_ndims
  integer (C_SIZE_T), dimension(size(start)) :: c_start
  integer (C_SIZE_T), dimension(size(ncount)) :: c_ncount
  integer (C_INTPTR_T), dimension(size(stride)) :: c_stride
  real (C_DOUBLE), dimension(size(dp,1),size(dp,2),size(dp,3),size(dp,4),size(dp,5), &
      size(dp,6),size(dp,7)), target :: c_dp
  integer i, ndims
  c_ncid = ncid
  c_varid = varid - 1
  c_ndims = size(start)
  ierr = nc_inq_varndims(c_ncid,c_varid,C_LOC(c_ndims))
  ndims = c_ndims
  do i = 1,ndims
    c_start(ndims-i+1) = start(i) - 1
    c_ncount(ndims-i+1) = ncount(i)
    c_stride(ndims-i+1) = stride(i)
  end do
  c_dp = dp
  ierr = nc_put_vars_double(c_ncid,c_varid,c_start,c_ncount,c_stride,C_LOC(c_dp))
end function nf_put_vars_double_d7

integer function nf_put_varm_int1_d1(ncid,varid,start,ncount,stride,imap,bp) result(ierr)
  implicit none
  integer, intent(in) :: ncid, varid
  integer, dimension(:), intent(in) :: start, ncount, stride, imap
  integer(kind=1), dimension(:), intent(in) :: bp
  integer (C_INT) :: c_ncid, c_varid
  integer (C_INT), target :: c_ndims
  integer (C_SIZE_T), dimension(size(start)) :: c_start
  integer (C_SIZE_T), dimension(size(ncount)) :: c_ncount
  integer (C_INTPTR_T), dimension(size(stride)) :: c_stride
  integer (C_INTPTR_T), dimension(size(imap)) :: c_imap
  integer (C_SIGNED_CHAR), dimension(size(bp)), target :: c_bp
  integer i, ndims
  c_ncid = ncid
  c_varid = varid - 1
  c_ndims = size(start)
  ierr = nc_inq_varndims(c_ncid,c_varid,C_LOC(c_ndims))
  ndims = c_ndims
  do i = 1,ndims
    c_start(ndims-i+1) = start(i) - 1
    c_ncount(ndims-i+1) = ncount(i)
    c_stride(ndims-i+1) = stride(i)
    c_imap(ndims-i+1) = imap(i)
  end do
  c_bp = bp
  ierr = nc_put_varm_schar(c_ncid,c_varid,c_start,c_ncount,c_stride,c_imap,C_LOC(c_bp))
end function nf_put_varm_int1_d1

integer function nf_put_varm_int1_d2(ncid,varid,start,ncount,stride,imap,bp) result(ierr)
  implicit none
  integer, intent(in) :: ncid, varid
  integer, dimension(:), intent(in) :: start, ncount, stride, imap
  integer(kind=1), dimension(:,:), intent(in) :: bp
  integer (C_INT) :: c_ncid, c_varid
  integer (C_INT), target :: c_ndims
  integer (C_SIZE_T), dimension(size(start)) :: c_start
  integer (C_SIZE_T), dimension(size(ncount)) :: c_ncount
  integer (C_INTPTR_T), dimension(size(stride)) :: c_stride
  integer (C_INTPTR_T), dimension(size(imap)) :: c_imap
  integer (C_SIGNED_CHAR), dimension(size(bp,1),size(bp,2)), target :: c_bp
  integer i, ndims
  c_ncid = ncid
  c_varid = varid - 1
  c_ndims = size(start)
  ierr = nc_inq_varndims(c_ncid,c_varid,C_LOC(c_ndims))
  ndims = c_ndims
  do i = 1,ndims
    c_start(ndims-i+1) = start(i) - 1
    c_ncount(ndims-i+1) = ncount(i)
    c_stride(ndims-i+1) = stride(i)
    c_imap(ndims-i+1) = imap(i)
  end do
  c_bp = bp
  ierr = nc_put_varm_schar(c_ncid,c_varid,c_start,c_ncount,c_stride,c_imap,C_LOC(c_bp))
end function nf_put_varm_int1_d2

integer function nf_put_varm_int1_d3(ncid,varid,start,ncount,stride,imap,bp) result(ierr)
  implicit none
  integer, intent(in) :: ncid, varid
  integer, dimension(:), intent(in) :: start, ncount, stride, imap
  integer(kind=1), dimension(:,:,:), intent(in) :: bp
  integer (C_INT) :: c_ncid, c_varid
  integer (C_INT), target :: c_ndims
  integer (C_SIZE_T), dimension(size(start)) :: c_start
  integer (C_SIZE_T), dimension(size(ncount)) :: c_ncount
  integer (C_INTPTR_T), dimension(size(stride)) :: c_stride
  integer (C_INTPTR_T), dimension(size(imap)) :: c_imap
  integer (C_SIGNED_CHAR), dimension(size(bp,1),size(bp,2),size(bp,3)), target :: c_bp
  integer i, ndims
  c_ncid = ncid
  c_varid = varid - 1
  c_ndims = size(start)
  ierr = nc_inq_varndims(c_ncid,c_varid,C_LOC(c_ndims))
  ndims = c_ndims
  do i = 1,ndims
    c_start(ndims-i+1) = start(i) - 1
    c_ncount(ndims-i+1) = ncount(i)
    c_stride(ndims-i+1) = stride(i)
    c_imap(ndims-i+1) = imap(i)
  end do
  c_bp = bp
  ierr = nc_put_varm_schar(c_ncid,c_varid,c_start,c_ncount,c_stride,c_imap,C_LOC(c_bp))
end function nf_put_varm_int1_d3

integer function nf_put_varm_int1_d4(ncid,varid,start,ncount,stride,imap,bp) result(ierr)
  implicit none
  integer, intent(in) :: ncid, varid
  integer, dimension(:), intent(in) :: start, ncount, stride, imap
  integer(kind=1), dimension(:,:,:,:), intent(in) :: bp
  integer (C_INT) :: c_ncid, c_varid
  integer (C_INT), target :: c_ndims
  integer (C_SIZE_T), dimension(size(start)) :: c_start
  integer (C_SIZE_T), dimension(size(ncount)) :: c_ncount
  integer (C_INTPTR_T), dimension(size(stride)) :: c_stride
  integer (C_INTPTR_T), dimension(size(imap)) :: c_imap
  integer (C_SIGNED_CHAR), dimension(size(bp,1),size(bp,2),size(bp,3),size(bp,4)), target :: c_bp
  integer i, ndims
  c_ncid = ncid
  c_varid = varid - 1
  c_ndims = size(start)
  ierr = nc_inq_varndims(c_ncid,c_varid,C_LOC(c_ndims))
  ndims = c_ndims
  do i = 1,ndims
    c_start(ndims-i+1) = start(i) - 1
    c_ncount(ndims-i+1) = ncount(i)
    c_stride(ndims-i+1) = stride(i)
    c_imap(ndims-i+1) = imap(i)
  end do
  c_bp = bp
  ierr = nc_put_varm_schar(c_ncid,c_varid,c_start,c_ncount,c_stride,c_imap,C_LOC(c_bp))
end function nf_put_varm_int1_d4

integer function nf_put_varm_int1_d5(ncid,varid,start,ncount,stride,imap,bp) result(ierr)
  implicit none
  integer, intent(in) :: ncid, varid
  integer, dimension(:), intent(in) :: start, ncount, stride, imap
  integer(kind=1), dimension(:,:,:,:,:), intent(in) :: bp
  integer (C_INT) :: c_ncid, c_varid
  integer (C_INT), target :: c_ndims
  integer (C_SIZE_T), dimension(size(start)) :: c_start
  integer (C_SIZE_T), dimension(size(ncount)) :: c_ncount
  integer (C_INTPTR_T), dimension(size(stride)) :: c_stride
  integer (C_INTPTR_T), dimension(size(imap)) :: c_imap
  integer (C_SIGNED_CHAR), dimension(size(bp,1),size(bp,2),size(bp,3),size(bp,4),size(bp,5)), &
    target :: c_bp
  integer i, ndims
  c_ncid = ncid
  c_varid = varid - 1
  c_ndims = size(start)
  ierr = nc_inq_varndims(c_ncid,c_varid,C_LOC(c_ndims))
  ndims = c_ndims
  do i = 1,ndims
    c_start(ndims-i+1) = start(i) - 1
    c_ncount(ndims-i+1) = ncount(i)
    c_stride(ndims-i+1) = stride(i)
    c_imap(ndims-i+1) = imap(i)
  end do
  c_bp = bp
  ierr = nc_put_varm_schar(c_ncid,c_varid,c_start,c_ncount,c_stride,c_imap,C_LOC(c_bp))
end function nf_put_varm_int1_d5

integer function nf_put_varm_int1_d6(ncid,varid,start,ncount,stride,imap,bp) result(ierr)
  implicit none
  integer, intent(in) :: ncid, varid
  integer, dimension(:), intent(in) :: start, ncount, stride, imap
  integer(kind=1), dimension(:,:,:,:,:,:), intent(in) :: bp
  integer (C_INT) :: c_ncid, c_varid
  integer (C_INT), target :: c_ndims
  integer (C_SIZE_T), dimension(size(start)) :: c_start
  integer (C_SIZE_T), dimension(size(ncount)) :: c_ncount
  integer (C_INTPTR_T), dimension(size(stride)) :: c_stride
  integer (C_INTPTR_T), dimension(size(imap)) :: c_imap
  integer (C_SIGNED_CHAR), dimension(size(bp,1),size(bp,2),size(bp,3),size(bp,4),size(bp,5), &
      size(bp,6)), target :: c_bp
  integer i, ndims
  c_ncid = ncid
  c_varid = varid - 1
  c_ndims = size(start)
  ierr = nc_inq_varndims(c_ncid,c_varid,C_LOC(c_ndims))
  ndims = c_ndims
  do i = 1,ndims
    c_start(ndims-i+1) = start(i) - 1
    c_ncount(ndims-i+1) = ncount(i)
    c_stride(ndims-i+1) = stride(i)
    c_imap(ndims-i+1) = imap(i)
  end do
  c_bp = bp
  ierr = nc_put_varm_schar(c_ncid,c_varid,c_start,c_ncount,c_stride,c_imap,C_LOC(c_bp))
end function nf_put_varm_int1_d6

integer function nf_put_varm_int1_d7(ncid,varid,start,ncount,stride,imap,bp) result(ierr)
  implicit none
  integer, intent(in) :: ncid, varid
  integer, dimension(:), intent(in) :: start, ncount, stride, imap
  integer(kind=1), dimension(:,:,:,:,:,:,:), intent(in) :: bp
  integer (C_INT) :: c_ncid, c_varid
  integer (C_INT), target :: c_ndims
  integer (C_SIZE_T), dimension(size(start)) :: c_start
  integer (C_SIZE_T), dimension(size(ncount)) :: c_ncount
  integer (C_INTPTR_T), dimension(size(stride)) :: c_stride
  integer (C_INTPTR_T), dimension(size(imap)) :: c_imap
  integer (C_SIGNED_CHAR), dimension(size(bp,1),size(bp,2),size(bp,3),size(bp,4),size(bp,5), &
      size(bp,6),size(bp,7)), target :: c_bp
  integer i, ndims
  c_ncid = ncid
  c_varid = varid - 1
  c_ndims = size(start)
  ierr = nc_inq_varndims(c_ncid,c_varid,C_LOC(c_ndims))
  ndims = c_ndims
  do i = 1,ndims
    c_start(ndims-i+1) = start(i) - 1
    c_ncount(ndims-i+1) = ncount(i)
    c_stride(ndims-i+1) = stride(i)
    c_imap(ndims-i+1) = imap(i)
  end do
  c_bp = bp
  ierr = nc_put_varm_schar(c_ncid,c_varid,c_start,c_ncount,c_stride,c_imap,C_LOC(c_bp))
end function nf_put_varm_int1_d7

integer function nf_put_varm_int2_d1(ncid,varid,start,ncount,stride,imap,sp) result(ierr)
  implicit none
  integer, intent(in) :: ncid, varid
  integer, dimension(:), intent(in) :: start, ncount, stride, imap
  integer(kind=2), dimension(:), intent(in) :: sp
  integer (C_INT) :: c_ncid, c_varid
  integer (C_INT), target :: c_ndims
  integer (C_SIZE_T), dimension(size(start)) :: c_start
  integer (C_SIZE_T), dimension(size(ncount)) :: c_ncount
  integer (C_INTPTR_T), dimension(size(stride)) :: c_stride
  integer (C_INTPTR_T), dimension(size(imap)) :: c_imap
  integer (C_SHORT), dimension(size(sp)), target :: c_sp
  integer i, ndims
  c_ncid = ncid
  c_varid = varid - 1
  c_ndims = size(start)
  ierr = nc_inq_varndims(c_ncid,c_varid,C_LOC(c_ndims))
  ndims = c_ndims
  do i = 1,ndims
    c_start(ndims-i+1) = start(i) - 1
    c_ncount(ndims-i+1) = ncount(i)
    c_stride(ndims-i+1) = stride(i)
    c_imap(ndims-i+1) = imap(i)
  end do
  c_sp = sp
  ierr = nc_put_varm_short(c_ncid,c_varid,c_start,c_ncount,c_stride,c_imap,C_LOC(c_sp))
end function nf_put_varm_int2_d1

integer function nf_put_varm_int2_d2(ncid,varid,start,ncount,stride,imap,sp) result(ierr)
  implicit none
  integer, intent(in) :: ncid, varid
  integer, dimension(:), intent(in) :: start, ncount, stride, imap
  integer(kind=2), dimension(:,:), intent(in) :: sp
  integer (C_INT) :: c_ncid, c_varid
  integer (C_INT), target :: c_ndims
  integer (C_SIZE_T), dimension(size(start)) :: c_start
  integer (C_SIZE_T), dimension(size(ncount)) :: c_ncount
  integer (C_INTPTR_T), dimension(size(stride)) :: c_stride
  integer (C_INTPTR_T), dimension(size(imap)) :: c_imap
  integer (C_SHORT), dimension(size(sp,1),size(sp,2)), target :: c_sp
  integer i, ndims
  c_ncid = ncid
  c_varid = varid - 1
  c_ndims = size(start)
  ierr = nc_inq_varndims(c_ncid,c_varid,C_LOC(c_ndims))
  ndims = c_ndims
  do i = 1,ndims
    c_start(ndims-i+1) = start(i) - 1
    c_ncount(ndims-i+1) = ncount(i)
    c_stride(ndims-i+1) = stride(i)
    c_imap(ndims-i+1) = imap(i)
  end do
  c_sp = sp
  ierr = nc_put_varm_short(c_ncid,c_varid,c_start,c_ncount,c_stride,c_imap,C_LOC(c_sp))
end function nf_put_varm_int2_d2

integer function nf_put_varm_int2_d3(ncid,varid,start,ncount,stride,imap,sp) result(ierr)
  implicit none
  integer, intent(in) :: ncid, varid
  integer, dimension(:), intent(in) :: start, ncount, stride, imap
  integer(kind=2), dimension(:,:,:), intent(in) :: sp
  integer (C_INT) :: c_ncid, c_varid
  integer (C_INT), target :: c_ndims
  integer (C_SIZE_T), dimension(size(start)) :: c_start
  integer (C_SIZE_T), dimension(size(ncount)) :: c_ncount
  integer (C_INTPTR_T), dimension(size(stride)) :: c_stride
  integer (C_INTPTR_T), dimension(size(imap)) :: c_imap
  integer (C_SHORT), dimension(size(sp,1),size(sp,2),size(sp,3)), target :: c_sp
  integer i, ndims
  c_ncid = ncid
  c_varid = varid - 1
  c_ndims = size(start)
  ierr = nc_inq_varndims(c_ncid,c_varid,C_LOC(c_ndims))
  ndims = c_ndims
  do i = 1,ndims
    c_start(ndims-i+1) = start(i) - 1
    c_ncount(ndims-i+1) = ncount(i)
    c_stride(ndims-i+1) = stride(i)
    c_imap(ndims-i+1) = imap(i)
  end do
  c_sp = sp
  ierr = nc_put_varm_short(c_ncid,c_varid,c_start,c_ncount,c_stride,c_imap,C_LOC(c_sp))
end function nf_put_varm_int2_d3

integer function nf_put_varm_int2_d4(ncid,varid,start,ncount,stride,imap,sp) result(ierr)
  implicit none
  integer, intent(in) :: ncid, varid
  integer, dimension(:), intent(in) :: start, ncount, stride, imap
  integer(kind=2), dimension(:,:,:,:), intent(in) :: sp
  integer (C_INT) :: c_ncid, c_varid
  integer (C_INT), target :: c_ndims
  integer (C_SIZE_T), dimension(size(start)) :: c_start
  integer (C_SIZE_T), dimension(size(ncount)) :: c_ncount
  integer (C_INTPTR_T), dimension(size(stride)) :: c_stride
  integer (C_INTPTR_T), dimension(size(imap)) :: c_imap
  integer (C_SHORT), dimension(size(sp,1),size(sp,2),size(sp,3),size(sp,4)), target :: c_sp
  integer i, ndims
  c_ncid = ncid
  c_varid = varid - 1
  c_ndims = size(start)
  ierr = nc_inq_varndims(c_ncid,c_varid,C_LOC(c_ndims))
  ndims = c_ndims
  do i = 1,ndims
    c_start(ndims-i+1) = start(i) - 1
    c_ncount(ndims-i+1) = ncount(i)
    c_stride(ndims-i+1) = stride(i)
    c_imap(ndims-i+1) = imap(i)
  end do
  c_sp = sp
  ierr = nc_put_varm_short(c_ncid,c_varid,c_start,c_ncount,c_stride,c_imap,C_LOC(c_sp))
end function nf_put_varm_int2_d4

integer function nf_put_varm_int2_d5(ncid,varid,start,ncount,stride,imap,sp) result(ierr)
  implicit none
  integer, intent(in) :: ncid, varid
  integer, dimension(:), intent(in) :: start, ncount, stride, imap
  integer(kind=2), dimension(:,:,:,:,:), intent(in) :: sp
  integer (C_INT) :: c_ncid, c_varid
  integer (C_INT), target :: c_ndims
  integer (C_SIZE_T), dimension(size(start)) :: c_start
  integer (C_SIZE_T), dimension(size(ncount)) :: c_ncount
  integer (C_INTPTR_T), dimension(size(stride)) :: c_stride
  integer (C_INTPTR_T), dimension(size(imap)) :: c_imap
  integer (C_SHORT), dimension(size(sp,1),size(sp,2),size(sp,3),size(sp,4),size(sp,5)), &
      target :: c_sp
  integer i, ndims
  c_ncid = ncid
  c_varid = varid - 1
  c_ndims = size(start)
  ierr = nc_inq_varndims(c_ncid,c_varid,C_LOC(c_ndims))
  ndims = c_ndims
  do i = 1,ndims
    c_start(ndims-i+1) = start(i) - 1
    c_ncount(ndims-i+1) = ncount(i)
    c_stride(ndims-i+1) = stride(i)
    c_imap(ndims-i+1) = imap(i)
  end do
  c_sp = sp
  ierr = nc_put_varm_short(c_ncid,c_varid,c_start,c_ncount,c_stride,c_imap,C_LOC(c_sp))
end function nf_put_varm_int2_d5

integer function nf_put_varm_int2_d6(ncid,varid,start,ncount,stride,imap,sp) result(ierr)
  implicit none
  integer, intent(in) :: ncid, varid
  integer, dimension(:), intent(in) :: start, ncount, stride, imap
  integer(kind=2), dimension(:,:,:,:,:,:), intent(in) :: sp
  integer (C_INT) :: c_ncid, c_varid
  integer (C_INT), target :: c_ndims
  integer (C_SIZE_T), dimension(size(start)) :: c_start
  integer (C_SIZE_T), dimension(size(ncount)) :: c_ncount
  integer (C_INTPTR_T), dimension(size(stride)) :: c_stride
  integer (C_INTPTR_T), dimension(size(imap)) :: c_imap
  integer (C_SHORT), dimension(size(sp,1),size(sp,2),size(sp,3),size(sp,4),size(sp,5), &
      size(sp,6)), target :: c_sp
  integer i, ndims
  c_ncid = ncid
  c_varid = varid - 1
  c_ndims = size(start)
  ierr = nc_inq_varndims(c_ncid,c_varid,C_LOC(c_ndims))
  ndims = c_ndims
  do i = 1,ndims
    c_start(ndims-i+1) = start(i) - 1
    c_ncount(ndims-i+1) = ncount(i)
    c_stride(ndims-i+1) = stride(i)
    c_imap(ndims-i+1) = imap(i)
  end do
  c_sp = sp
  ierr = nc_put_varm_short(c_ncid,c_varid,c_start,c_ncount,c_stride,c_imap,C_LOC(c_sp))
end function nf_put_varm_int2_d6

integer function nf_put_varm_int2_d7(ncid,varid,start,ncount,stride,imap,sp) result(ierr)
  implicit none
  integer, intent(in) :: ncid, varid
  integer, dimension(:), intent(in) :: start, ncount, stride, imap
  integer(kind=2), dimension(:,:,:,:,:,:,:), intent(in) :: sp
  integer (C_INT) :: c_ncid, c_varid
  integer (C_INT), target :: c_ndims
  integer (C_SIZE_T), dimension(size(start)) :: c_start
  integer (C_SIZE_T), dimension(size(ncount)) :: c_ncount
  integer (C_INTPTR_T), dimension(size(stride)) :: c_stride
  integer (C_INTPTR_T), dimension(size(imap)) :: c_imap
  integer (C_SHORT), dimension(size(sp,1),size(sp,2),size(sp,3),size(sp,4),size(sp,5), &
      size(sp,6),size(sp,7)), target :: c_sp
  integer i, ndims
  c_ncid = ncid
  c_varid = varid - 1
  c_ndims = size(start)
  ierr = nc_inq_varndims(c_ncid,c_varid,C_LOC(c_ndims))
  ndims = c_ndims
  do i = 1,ndims
    c_start(ndims-i+1) = start(i) - 1
    c_ncount(ndims-i+1) = ncount(i)
    c_stride(ndims-i+1) = stride(i)
    c_imap(ndims-i+1) = imap(i)
  end do
  c_sp = sp
  ierr = nc_put_varm_short(c_ncid,c_varid,c_start,c_ncount,c_stride,c_imap,C_LOC(c_sp))
end function nf_put_varm_int2_d7

integer function nf_put_varm_int_d1(ncid,varid,start,ncount,stride,imap,ip) result(ierr)
  implicit none
  integer, intent(in) :: ncid, varid
  integer, dimension(:), intent(in) :: start, ncount, stride, imap
  integer(kind=4), dimension(:), intent(in) :: ip
  integer (C_INT) :: c_ncid, c_varid
  integer (C_INT), target :: c_ndims
  integer (C_SIZE_T), dimension(size(start)) :: c_start
  integer (C_SIZE_T), dimension(size(ncount)) :: c_ncount
  integer (C_INTPTR_T), dimension(size(stride)) :: c_stride
  integer (C_INTPTR_T), dimension(size(imap)) :: c_imap
  integer (C_INT), dimension(size(ip)), target :: c_ip
  integer i, ndims
  c_ncid = ncid
  c_varid = varid - 1
  c_ndims = size(start)
  ierr = nc_inq_varndims(c_ncid,c_varid,C_LOC(c_ndims))
  ndims = c_ndims
  do i = 1,ndims
    c_start(ndims-i+1) = start(i) - 1
    c_ncount(ndims-i+1) = ncount(i)
    c_stride(ndims-i+1) = stride(i)
    c_imap(ndims-i+1) = imap(i)
  end do
  c_ip = ip
  ierr = nc_put_varm_int(c_ncid,c_varid,c_start,c_ncount,c_stride,c_imap,C_LOC(c_ip))
end function nf_put_varm_int_d1

integer function nf_put_varm_int_d2(ncid,varid,start,ncount,stride,imap,ip) result(ierr)
  implicit none
  integer, intent(in) :: ncid, varid
  integer, dimension(:), intent(in) :: start, ncount, stride, imap
  integer(kind=4), dimension(:,:), intent(in) :: ip
  integer (C_INT) :: c_ncid, c_varid
  integer (C_INT), target :: c_ndims
  integer (C_SIZE_T), dimension(size(start)) :: c_start
  integer (C_SIZE_T), dimension(size(ncount)) :: c_ncount
  integer (C_INTPTR_T), dimension(size(stride)) :: c_stride
  integer (C_INTPTR_T), dimension(size(imap)) :: c_imap
  integer (C_INT), dimension(size(ip,1),size(ip,2)), target :: c_ip
  integer i, ndims
  c_ncid = ncid
  c_varid = varid - 1
  c_ndims = size(start)
  ierr = nc_inq_varndims(c_ncid,c_varid,C_LOC(c_ndims))
  ndims = c_ndims
  do i = 1,ndims
    c_start(ndims-i+1) = start(i) - 1
    c_ncount(ndims-i+1) = ncount(i)
    c_stride(ndims-i+1) = stride(i)
    c_imap(ndims-i+1) = imap(i)
  end do
  c_ip = ip
  ierr = nc_put_varm_int(c_ncid,c_varid,c_start,c_ncount,c_stride,c_imap,C_LOC(c_ip))
end function nf_put_varm_int_d2

integer function nf_put_varm_int_d3(ncid,varid,start,ncount,stride,imap,ip) result(ierr)
  implicit none
  integer, intent(in) :: ncid, varid
  integer, dimension(:), intent(in) :: start, ncount, stride, imap
  integer(kind=4), dimension(:,:,:), intent(in) :: ip
  integer (C_INT) :: c_ncid, c_varid
  integer (C_INT), target :: c_ndims
  integer (C_SIZE_T), dimension(size(start)) :: c_start
  integer (C_SIZE_T), dimension(size(ncount)) :: c_ncount
  integer (C_INTPTR_T), dimension(size(stride)) :: c_stride
  integer (C_INTPTR_T), dimension(size(imap)) :: c_imap
  integer (C_INT), dimension(size(ip,1),size(ip,2),size(ip,3)), target :: c_ip
  integer i, ndims
  c_ncid = ncid
  c_varid = varid - 1
  c_ndims = size(start)
  ierr = nc_inq_varndims(c_ncid,c_varid,C_LOC(c_ndims))
  ndims = c_ndims
  do i = 1,ndims
    c_start(ndims-i+1) = start(i) - 1
    c_ncount(ndims-i+1) = ncount(i)
    c_stride(ndims-i+1) = stride(i)
    c_imap(ndims-i+1) = imap(i)
  end do
  c_ip = ip
  ierr = nc_put_varm_int(c_ncid,c_varid,c_start,c_ncount,c_stride,c_imap,C_LOC(c_ip))
end function nf_put_varm_int_d3

integer function nf_put_varm_int_d4(ncid,varid,start,ncount,stride,imap,ip) result(ierr)
  implicit none
  integer, intent(in) :: ncid, varid
  integer, dimension(:), intent(in) :: start, ncount, stride, imap
  integer(kind=4), dimension(:,:,:,:), intent(in) :: ip
  integer (C_INT) :: c_ncid, c_varid
  integer (C_INT), target :: c_ndims
  integer (C_SIZE_T), dimension(size(start)) :: c_start
  integer (C_SIZE_T), dimension(size(ncount)) :: c_ncount
  integer (C_INTPTR_T), dimension(size(stride)) :: c_stride
  integer (C_INTPTR_T), dimension(size(imap)) :: c_imap
  integer (C_INT), dimension(size(ip,1),size(ip,2),size(ip,3),size(ip,4)), target :: c_ip
  integer i, ndims
  c_ncid = ncid
  c_varid = varid - 1
  c_ndims = size(start)
  ierr = nc_inq_varndims(c_ncid,c_varid,C_LOC(c_ndims))
  ndims = c_ndims
  do i = 1,ndims
    c_start(ndims-i+1) = start(i) - 1
    c_ncount(ndims-i+1) = ncount(i)
    c_stride(ndims-i+1) = stride(i)
    c_imap(ndims-i+1) = imap(i)
  end do
  c_ip = ip
  ierr = nc_put_varm_int(c_ncid,c_varid,c_start,c_ncount,c_stride,c_imap,C_LOC(c_ip))
end function nf_put_varm_int_d4

integer function nf_put_varm_int_d5(ncid,varid,start,ncount,stride,imap,ip) result(ierr)
  implicit none
  integer, intent(in) :: ncid, varid
  integer, dimension(:), intent(in) :: start, ncount, stride, imap
  integer(kind=4), dimension(:,:,:,:,:), intent(in) :: ip
  integer (C_INT) :: c_ncid, c_varid
  integer (C_INT), target :: c_ndims
  integer (C_SIZE_T), dimension(size(start)) :: c_start
  integer (C_SIZE_T), dimension(size(ncount)) :: c_ncount
  integer (C_INTPTR_T), dimension(size(stride)) :: c_stride
  integer (C_INTPTR_T), dimension(size(imap)) :: c_imap
  integer (C_INT), dimension(size(ip,1),size(ip,2),size(ip,3),size(ip,4),size(ip,5)), &
      target :: c_ip
  integer i, ndims
  c_ncid = ncid
  c_varid = varid - 1
  c_ndims = size(start)
  ierr = nc_inq_varndims(c_ncid,c_varid,C_LOC(c_ndims))
  ndims = c_ndims
  do i = 1,ndims
    c_start(ndims-i+1) = start(i) - 1
    c_ncount(ndims-i+1) = ncount(i)
    c_stride(ndims-i+1) = stride(i)
    c_imap(ndims-i+1) = imap(i)
  end do
  c_ip = ip
  ierr = nc_put_varm_int(c_ncid,c_varid,c_start,c_ncount,c_stride,c_imap,C_LOC(c_ip))
end function nf_put_varm_int_d5

integer function nf_put_varm_int_d6(ncid,varid,start,ncount,stride,imap,ip) result(ierr)
  implicit none
  integer, intent(in) :: ncid, varid
  integer, dimension(:), intent(in) :: start, ncount, stride, imap
  integer(kind=4), dimension(:,:,:,:,:,:), intent(in) :: ip
  integer (C_INT) :: c_ncid, c_varid
  integer (C_INT), target :: c_ndims
  integer (C_SIZE_T), dimension(size(start)) :: c_start
  integer (C_SIZE_T), dimension(size(ncount)) :: c_ncount
  integer (C_INTPTR_T), dimension(size(stride)) :: c_stride
  integer (C_INTPTR_T), dimension(size(imap)) :: c_imap
  integer (C_INT), dimension(size(ip,1),size(ip,2),size(ip,3),size(ip,4),size(ip,5), &
      size(ip,6)), target :: c_ip
  integer i, ndims
  c_ncid = ncid
  c_varid = varid - 1
  c_ndims = size(start)
  ierr = nc_inq_varndims(c_ncid,c_varid,C_LOC(c_ndims))
  ndims = c_ndims
  do i = 1,ndims
    c_start(ndims-i+1) = start(i) - 1
    c_ncount(ndims-i+1) = ncount(i)
    c_stride(ndims-i+1) = stride(i)
    c_imap(ndims-i+1) = imap(i)
  end do
  c_ip = ip
  ierr = nc_put_varm_int(c_ncid,c_varid,c_start,c_ncount,c_stride,c_imap,C_LOC(c_ip))
end function nf_put_varm_int_d6

integer function nf_put_varm_int_d7(ncid,varid,start,ncount,stride,imap,ip) result(ierr)
  implicit none
  integer, intent(in) :: ncid, varid
  integer, dimension(:), intent(in) :: start, ncount, stride, imap
  integer(kind=4), dimension(:,:,:,:,:,:,:), intent(in) :: ip
  integer (C_INT) :: c_ncid, c_varid
  integer (C_INT), target :: c_ndims
  integer (C_SIZE_T), dimension(size(start)) :: c_start
  integer (C_SIZE_T), dimension(size(ncount)) :: c_ncount
  integer (C_INTPTR_T), dimension(size(stride)) :: c_stride
  integer (C_INTPTR_T), dimension(size(imap)) :: c_imap
  integer (C_INT), dimension(size(ip,1),size(ip,2),size(ip,3),size(ip,4),size(ip,5), &
      size(ip,6),size(ip,7)), target :: c_ip
  integer i, ndims
  c_ncid = ncid
  c_varid = varid - 1
  c_ndims = size(start)
  ierr = nc_inq_varndims(c_ncid,c_varid,C_LOC(c_ndims))
  ndims = c_ndims
  do i = 1,ndims
    c_start(ndims-i+1) = start(i) - 1
    c_ncount(ndims-i+1) = ncount(i)
    c_stride(ndims-i+1) = stride(i)
    c_imap(ndims-i+1) = imap(i)
  end do
  c_ip = ip
  ierr = nc_put_varm_int(c_ncid,c_varid,c_start,c_ncount,c_stride,c_imap,C_LOC(c_ip))
end function nf_put_varm_int_d7

integer function nf_put_varm_real_d1(ncid,varid,start,ncount,stride,imap,fp) result(ierr)
  implicit none
  integer, intent(in) :: ncid, varid
  integer, dimension(:), intent(in) :: start, ncount, stride, imap
  real(kind=4), dimension(:), intent(in) :: fp
  integer (C_INT) :: c_ncid, c_varid
  integer (C_INT), target :: c_ndims
  integer (C_SIZE_T), dimension(size(start)) :: c_start
  integer (C_SIZE_T), dimension(size(ncount)) :: c_ncount
  integer (C_INTPTR_T), dimension(size(stride)) :: c_stride
  integer (C_INTPTR_T), dimension(size(imap)) :: c_imap
  real (C_FLOAT), dimension(size(fp)), target :: c_fp
  integer i, ndims
  c_ncid = ncid
  c_varid = varid - 1
  c_ndims = size(start)
  ierr = nc_inq_varndims(c_ncid,c_varid,C_LOC(c_ndims))
  ndims = c_ndims
  do i = 1,ndims
    c_start(ndims-i+1) = start(i) - 1
    c_ncount(ndims-i+1) = ncount(i)
    c_stride(ndims-i+1) = stride(i)
    c_imap(ndims-i+1) = imap(i)
  end do
  c_fp = fp
  ierr = nc_put_varm_float(c_ncid,c_varid,c_start,c_ncount,c_stride,c_imap,C_LOC(c_fp))
end function nf_put_varm_real_d1

integer function nf_put_varm_real_d2(ncid,varid,start,ncount,stride,imap,fp) result(ierr)
  implicit none
  integer, intent(in) :: ncid, varid
  integer, dimension(:), intent(in) :: start, ncount, stride, imap
  real(kind=4), dimension(:,:), intent(in) :: fp
  integer (C_INT) :: c_ncid, c_varid
  integer (C_INT), target :: c_ndims
  integer (C_SIZE_T), dimension(size(start)) :: c_start
  integer (C_SIZE_T), dimension(size(ncount)) :: c_ncount
  integer (C_INTPTR_T), dimension(size(stride)) :: c_stride
  integer (C_INTPTR_T), dimension(size(imap)) :: c_imap
  real (C_FLOAT), dimension(size(fp,1),size(fp,2)), target :: c_fp
  integer i, ndims
  c_ncid = ncid
  c_varid = varid - 1
  c_ndims = size(start)
  ierr = nc_inq_varndims(c_ncid,c_varid,C_LOC(c_ndims))
  ndims = c_ndims
  do i = 1,ndims
    c_start(ndims-i+1) = start(i) - 1
    c_ncount(ndims-i+1) = ncount(i)
    c_stride(ndims-i+1) = stride(i)
    c_imap(ndims-i+1) = imap(i)
  end do
  c_fp = fp
  ierr = nc_put_varm_float(c_ncid,c_varid,c_start,c_ncount,c_stride,c_imap,C_LOC(c_fp))
end function nf_put_varm_real_d2

integer function nf_put_varm_real_d3(ncid,varid,start,ncount,stride,imap,fp) result(ierr)
  implicit none
  integer, intent(in) :: ncid, varid
  integer, dimension(:), intent(in) :: start, ncount, stride, imap
  real(kind=4), dimension(:,:,:), intent(in) :: fp
  integer (C_INT) :: c_ncid, c_varid
  integer (C_INT), target :: c_ndims
  integer (C_SIZE_T), dimension(size(start)) :: c_start
  integer (C_SIZE_T), dimension(size(ncount)) :: c_ncount
  integer (C_INTPTR_T), dimension(size(stride)) :: c_stride
  integer (C_INTPTR_T), dimension(size(imap)) :: c_imap
  real (C_FLOAT), dimension(size(fp,1),size(fp,2),size(fp,3)), target :: c_fp
  integer i, ndims
  c_ncid = ncid
  c_varid = varid - 1
  c_ndims = size(start)
  ierr = nc_inq_varndims(c_ncid,c_varid,C_LOC(c_ndims))
  ndims = c_ndims
  do i = 1,ndims
    c_start(ndims-i+1) = start(i) - 1
    c_ncount(ndims-i+1) = ncount(i)
    c_stride(ndims-i+1) = stride(i)
    c_imap(ndims-i+1) = imap(i)
  end do
  c_fp = fp
  ierr = nc_put_varm_float(c_ncid,c_varid,c_start,c_ncount,c_stride,c_imap,C_LOC(c_fp))
end function nf_put_varm_real_d3

integer function nf_put_varm_real_d4(ncid,varid,start,ncount,stride,imap,fp) result(ierr)
  implicit none
  integer, intent(in) :: ncid, varid
  integer, dimension(:), intent(in) :: start, ncount, stride, imap
  real(kind=4), dimension(:,:,:,:), intent(in) :: fp
  integer (C_INT) :: c_ncid, c_varid
  integer (C_INT), target :: c_ndims
  integer (C_SIZE_T), dimension(size(start)) :: c_start
  integer (C_SIZE_T), dimension(size(ncount)) :: c_ncount
  integer (C_INTPTR_T), dimension(size(stride)) :: c_stride
  integer (C_INTPTR_T), dimension(size(imap)) :: c_imap
  real (C_FLOAT), dimension(size(fp,1),size(fp,2),size(fp,3),size(fp,4)), target :: c_fp
  integer i, ndims
  c_ncid = ncid
  c_varid = varid - 1
  c_ndims = size(start)
  ierr = nc_inq_varndims(c_ncid,c_varid,C_LOC(c_ndims))
  ndims = c_ndims
  do i = 1,ndims
    c_start(ndims-i+1) = start(i) - 1
    c_ncount(ndims-i+1) = ncount(i)
    c_stride(ndims-i+1) = stride(i)
    c_imap(ndims-i+1) = imap(i)
  end do
  c_fp = fp
  ierr = nc_put_varm_float(c_ncid,c_varid,c_start,c_ncount,c_stride,c_imap,C_LOC(c_fp))
end function nf_put_varm_real_d4

integer function nf_put_varm_real_d5(ncid,varid,start,ncount,stride,imap,fp) result(ierr)
  implicit none
  integer, intent(in) :: ncid, varid
  integer, dimension(:), intent(in) :: start, ncount, stride, imap
  real(kind=4), dimension(:,:,:,:,:), intent(in) :: fp
  integer (C_INT) :: c_ncid, c_varid
  integer (C_INT), target :: c_ndims
  integer (C_SIZE_T), dimension(size(start)) :: c_start
  integer (C_SIZE_T), dimension(size(ncount)) :: c_ncount
  integer (C_INTPTR_T), dimension(size(stride)) :: c_stride
  integer (C_INTPTR_T), dimension(size(imap)) :: c_imap
  real (C_FLOAT), dimension(size(fp,1),size(fp,2),size(fp,3),size(fp,4),size(fp,5)), &
      target :: c_fp
  integer i, ndims
  c_ncid = ncid
  c_varid = varid - 1
  c_ndims = size(start)
  ierr = nc_inq_varndims(c_ncid,c_varid,C_LOC(c_ndims))
  ndims = c_ndims
  do i = 1,ndims
    c_start(ndims-i+1) = start(i) - 1
    c_ncount(ndims-i+1) = ncount(i)
    c_stride(ndims-i+1) = stride(i)
    c_imap(ndims-i+1) = imap(i)
  end do
  c_fp = fp
  ierr = nc_put_varm_float(c_ncid,c_varid,c_start,c_ncount,c_stride,c_imap,C_LOC(c_fp))
end function nf_put_varm_real_d5

integer function nf_put_varm_real_d6(ncid,varid,start,ncount,stride,imap,fp) result(ierr)
  implicit none
  integer, intent(in) :: ncid, varid
  integer, dimension(:), intent(in) :: start, ncount, stride, imap
  real(kind=4), dimension(:,:,:,:,:,:), intent(in) :: fp
  integer (C_INT) :: c_ncid, c_varid
  integer (C_INT), target :: c_ndims
  integer (C_SIZE_T), dimension(size(start)) :: c_start
  integer (C_SIZE_T), dimension(size(ncount)) :: c_ncount
  integer (C_INTPTR_T), dimension(size(stride)) :: c_stride
  integer (C_INTPTR_T), dimension(size(imap)) :: c_imap
  real (C_FLOAT), dimension(size(fp,1),size(fp,2),size(fp,3),size(fp,4),size(fp,5), &
      size(fp,6)), target :: c_fp
  integer i, ndims
  c_ncid = ncid
  c_varid = varid - 1
  c_ndims = size(start)
  ierr = nc_inq_varndims(c_ncid,c_varid,C_LOC(c_ndims))
  ndims = c_ndims
  do i = 1,ndims
    c_start(ndims-i+1) = start(i) - 1
    c_ncount(ndims-i+1) = ncount(i)
    c_stride(ndims-i+1) = stride(i)
    c_imap(ndims-i+1) = imap(i)
  end do
  c_fp = fp
  ierr = nc_put_varm_float(c_ncid,c_varid,c_start,c_ncount,c_stride,c_imap,C_LOC(c_fp))
end function nf_put_varm_real_d6

integer function nf_put_varm_real_d7(ncid,varid,start,ncount,stride,imap,fp) result(ierr)
  implicit none
  integer, intent(in) :: ncid, varid
  integer, dimension(:), intent(in) :: start, ncount, stride, imap
  real(kind=4), dimension(:,:,:,:,:,:,:), intent(in) :: fp
  integer (C_INT) :: c_ncid, c_varid
  integer (C_INT), target :: c_ndims
  integer (C_SIZE_T), dimension(size(start)) :: c_start
  integer (C_SIZE_T), dimension(size(ncount)) :: c_ncount
  integer (C_INTPTR_T), dimension(size(stride)) :: c_stride
  integer (C_INTPTR_T), dimension(size(imap)) :: c_imap
  real (C_FLOAT), dimension(size(fp,1),size(fp,2),size(fp,3),size(fp,4),size(fp,5), &
      size(fp,6),size(fp,7)), target :: c_fp
  integer i, ndims
  c_ncid = ncid
  c_varid = varid - 1
  c_ndims = size(start)
  ierr = nc_inq_varndims(c_ncid,c_varid,C_LOC(c_ndims))
  ndims = c_ndims
  do i = 1,ndims
    c_start(ndims-i+1) = start(i) - 1
    c_ncount(ndims-i+1) = ncount(i)
    c_stride(ndims-i+1) = stride(i)
    c_imap(ndims-i+1) = imap(i)
  end do
  c_fp = fp
  ierr = nc_put_varm_float(c_ncid,c_varid,c_start,c_ncount,c_stride,c_imap,C_LOC(c_fp))
end function nf_put_varm_real_d7

integer function nf_put_varm_double_d1(ncid,varid,start,ncount,stride,imap,dp) result(ierr)
  implicit none
  integer, intent(in) :: ncid, varid
  integer, dimension(:), intent(in) :: start, ncount, stride, imap
  real(kind=8), dimension(:), intent(in) :: dp
  integer (C_INT) :: c_ncid, c_varid
  integer (C_INT), target :: c_ndims
  integer (C_SIZE_T), dimension(size(start)) :: c_start
  integer (C_SIZE_T), dimension(size(ncount)) :: c_ncount
  integer (C_INTPTR_T), dimension(size(stride)) :: c_stride
  integer (C_INTPTR_T), dimension(size(imap)) :: c_imap
  real (C_DOUBLE), dimension(size(dp)), target :: c_dp
  integer i, ndims
  c_ncid = ncid
  c_varid = varid - 1
  c_ndims = size(start)
  ierr = nc_inq_varndims(c_ncid,c_varid,C_LOC(c_ndims))
  ndims = c_ndims
  do i = 1,ndims
    c_start(ndims-i+1) = start(i) - 1
    c_ncount(ndims-i+1) = ncount(i)
    c_stride(ndims-i+1) = stride(i)
    c_imap(ndims-i+1) = imap(i)
  end do
  c_dp = dp
  ierr = nc_put_varm_double(c_ncid,c_varid,c_start,c_ncount,c_stride,c_imap,C_LOC(c_dp))
end function nf_put_varm_double_d1

integer function nf_put_varm_double_d2(ncid,varid,start,ncount,stride,imap,dp) result(ierr)
  implicit none
  integer, intent(in) :: ncid, varid
  integer, dimension(:), intent(in) :: start, ncount, stride, imap
  real(kind=8), dimension(:,:), intent(in) :: dp
  integer (C_INT) :: c_ncid, c_varid
  integer (C_INT), target :: c_ndims
  integer (C_SIZE_T), dimension(size(start)) :: c_start
  integer (C_SIZE_T), dimension(size(ncount)) :: c_ncount
  integer (C_INTPTR_T), dimension(size(stride)) :: c_stride
  integer (C_INTPTR_T), dimension(size(imap)) :: c_imap
  real (C_DOUBLE), dimension(size(dp,1),size(dp,2)), target :: c_dp
  integer i, ndims
  c_ncid = ncid
  c_varid = varid - 1
  c_ndims = size(start)
  ierr = nc_inq_varndims(c_ncid,c_varid,C_LOC(c_ndims))
  ndims = c_ndims
  do i = 1,ndims
    c_start(ndims-i+1) = start(i) - 1
    c_ncount(ndims-i+1) = ncount(i)
    c_stride(ndims-i+1) = stride(i)
    c_imap(ndims-i+1) = imap(i)
  end do
  c_dp = dp
  ierr = nc_put_varm_double(c_ncid,c_varid,c_start,c_ncount,c_stride,c_imap,C_LOC(c_dp))
end function nf_put_varm_double_d2

integer function nf_put_varm_double_d3(ncid,varid,start,ncount,stride,imap,dp) result(ierr)
  implicit none
  integer, intent(in) :: ncid, varid
  integer, dimension(:), intent(in) :: start, ncount, stride, imap
  real(kind=8), dimension(:,:,:), intent(in) :: dp
  integer (C_INT) :: c_ncid, c_varid
  integer (C_INT), target :: c_ndims
  integer (C_SIZE_T), dimension(size(start)) :: c_start
  integer (C_SIZE_T), dimension(size(ncount)) :: c_ncount
  integer (C_INTPTR_T), dimension(size(stride)) :: c_stride
  integer (C_INTPTR_T), dimension(size(imap)) :: c_imap
  real (C_DOUBLE), dimension(size(dp,1),size(dp,2),size(dp,3)), target :: c_dp
  integer i, ndims
  c_ncid = ncid
  c_varid = varid - 1
  c_ndims = size(start)
  ierr = nc_inq_varndims(c_ncid,c_varid,C_LOC(c_ndims))
  ndims = c_ndims
  do i = 1,ndims
    c_start(ndims-i+1) = start(i) - 1
    c_ncount(ndims-i+1) = ncount(i)
    c_stride(ndims-i+1) = stride(i)
    c_imap(ndims-i+1) = imap(i)
  end do
  c_dp = dp
  ierr = nc_put_varm_double(c_ncid,c_varid,c_start,c_ncount,c_stride,c_imap,C_LOC(c_dp))
end function nf_put_varm_double_d3

integer function nf_put_varm_double_d4(ncid,varid,start,ncount,stride,imap,dp) result(ierr)
  implicit none
  integer, intent(in) :: ncid, varid
  integer, dimension(:), intent(in) :: start, ncount, stride, imap
  real(kind=8), dimension(:,:,:,:), intent(in) :: dp
  integer (C_INT) :: c_ncid, c_varid
  integer (C_INT), target :: c_ndims
  integer (C_SIZE_T), dimension(size(start)) :: c_start
  integer (C_SIZE_T), dimension(size(ncount)) :: c_ncount
  integer (C_INTPTR_T), dimension(size(stride)) :: c_stride
  integer (C_INTPTR_T), dimension(size(imap)) :: c_imap
  real (C_DOUBLE), dimension(size(dp,1),size(dp,2),size(dp,3),size(dp,4)), target :: c_dp
  integer i, ndims
  c_ncid = ncid
  c_varid = varid - 1
  c_ndims = size(start)
  ierr = nc_inq_varndims(c_ncid,c_varid,C_LOC(c_ndims))
  ndims = c_ndims
  do i = 1,ndims
    c_start(ndims-i+1) = start(i) - 1
    c_ncount(ndims-i+1) = ncount(i)
    c_stride(ndims-i+1) = stride(i)
    c_imap(ndims-i+1) = imap(i)
  end do
  c_dp = dp
  ierr = nc_put_varm_double(c_ncid,c_varid,c_start,c_ncount,c_stride,c_imap,C_LOC(c_dp))
end function nf_put_varm_double_d4

integer function nf_put_varm_double_d5(ncid,varid,start,ncount,stride,imap,dp) result(ierr)
  implicit none
  integer, intent(in) :: ncid, varid
  integer, dimension(:), intent(in) :: start, ncount, stride, imap
  real(kind=8), dimension(:,:,:,:,:), intent(in) :: dp
  integer (C_INT) :: c_ncid, c_varid
  integer (C_INT), target :: c_ndims
  integer (C_SIZE_T), dimension(size(start)) :: c_start
  integer (C_SIZE_T), dimension(size(ncount)) :: c_ncount
  integer (C_INTPTR_T), dimension(size(stride)) :: c_stride
  integer (C_INTPTR_T), dimension(size(imap)) :: c_imap
  real (C_DOUBLE), dimension(size(dp,1),size(dp,2),size(dp,3),size(dp,4),size(dp,5)), &
      target :: c_dp
  integer i, ndims
  c_ncid = ncid
  c_varid = varid - 1
  c_ndims = size(start)
  ierr = nc_inq_varndims(c_ncid,c_varid,C_LOC(c_ndims))
  ndims = c_ndims
  do i = 1,ndims
    c_start(ndims-i+1) = start(i) - 1
    c_ncount(ndims-i+1) = ncount(i)
    c_stride(ndims-i+1) = stride(i)
    c_imap(ndims-i+1) = imap(i)
  end do
  c_dp = dp
  ierr = nc_put_varm_double(c_ncid,c_varid,c_start,c_ncount,c_stride,c_imap,C_LOC(c_dp))
end function nf_put_varm_double_d5

integer function nf_put_varm_double_d6(ncid,varid,start,ncount,stride,imap,dp) result(ierr)
  implicit none
  integer, intent(in) :: ncid, varid
  integer, dimension(:), intent(in) :: start, ncount, stride, imap
  real(kind=8), dimension(:,:,:,:,:,:), intent(in) :: dp
  integer (C_INT) :: c_ncid, c_varid
  integer (C_INT), target :: c_ndims
  integer (C_SIZE_T), dimension(size(start)) :: c_start
  integer (C_SIZE_T), dimension(size(ncount)) :: c_ncount
  integer (C_INTPTR_T), dimension(size(stride)) :: c_stride
  integer (C_INTPTR_T), dimension(size(imap)) :: c_imap
  real (C_DOUBLE), dimension(size(dp,1),size(dp,2),size(dp,3),size(dp,4),size(dp,5), &
      size(dp,6)), target :: c_dp
  integer i, ndims
  c_ncid = ncid
  c_varid = varid - 1
  c_ndims = size(start)
  ierr = nc_inq_varndims(c_ncid,c_varid,C_LOC(c_ndims))
  ndims = c_ndims
  do i = 1,ndims
    c_start(ndims-i+1) = start(i) - 1
    c_ncount(ndims-i+1) = ncount(i)
    c_stride(ndims-i+1) = stride(i)
    c_imap(ndims-i+1) = imap(i)
  end do
  c_dp = dp
  ierr = nc_put_varm_double(c_ncid,c_varid,c_start,c_ncount,c_stride,c_imap,C_LOC(c_dp))
end function nf_put_varm_double_d6

integer function nf_put_varm_double_d7(ncid,varid,start,ncount,stride,imap,dp) result(ierr)
  implicit none
  integer, intent(in) :: ncid, varid
  integer, dimension(:), intent(in) :: start, ncount, stride, imap
  real(kind=8), dimension(:,:,:,:,:,:,:), intent(in) :: dp
  integer (C_INT) :: c_ncid, c_varid
  integer (C_INT), target :: c_ndims
  integer (C_SIZE_T), dimension(size(start)) :: c_start
  integer (C_SIZE_T), dimension(size(ncount)) :: c_ncount
  integer (C_INTPTR_T), dimension(size(stride)) :: c_stride
  integer (C_INTPTR_T), dimension(size(imap)) :: c_imap
  real (C_DOUBLE), dimension(size(dp,1),size(dp,2),size(dp,3),size(dp,4),size(dp,5), &
      size(dp,6),size(dp,7)), target :: c_dp
  integer i, ndims
  c_ncid = ncid
  c_varid = varid - 1
  c_ndims = size(start)
  ierr = nc_inq_varndims(c_ncid,c_varid,C_LOC(c_ndims))
  ndims = c_ndims
  do i = 1,ndims
    c_start(ndims-i+1) = start(i) - 1
    c_ncount(ndims-i+1) = ncount(i)
    c_stride(ndims-i+1) = stride(i)
    c_imap(ndims-i+1) = imap(i)
  end do
  c_dp = dp
  ierr = nc_put_varm_double(c_ncid,c_varid,c_start,c_ncount,c_stride,c_imap,C_LOC(c_dp))
end function nf_put_varm_double_d7

integer function nf_copy_att(ncidin,varidin,name,ncidout,varidout) result(ierr)
  implicit none
  integer, intent(in) :: ncidin, varidin, ncidout, varidout
  character(len=*), intent(in) :: name
  integer (C_INT) :: c_ncidin, c_ncidout, c_varidin, c_varidout
  character, dimension(len(name)) :: c_name
  c_ncidin = ncidin
  c_varidin = varidin - 1
  c_ncidout = ncidout
  c_varidout = varidout - 1
  call cf_strcopy(name,c_name)
  ierr = nc_copy_att(c_ncidin,c_varidin,c_name,c_ncidout,c_varidout)
end function nf_copy_att

integer function nf_del_att(ncid,varid,name) result(ierr)
  implicit none
  integer, intent(in) :: ncid, varid
  character(len=*), intent(in) :: name
  integer (C_INT) :: c_ncid, c_varid
  character, dimension(len(name)) :: c_name
  c_ncid = ncid
  c_varid = varid - 1
  call cf_strcopy(name,c_name)
  ierr = nc_del_att(c_ncid,c_varid,c_name)
end function nf_del_att

subroutine cf_strcopy(fname,cname)
  implicit none
  character(len=*), intent(in) :: fname
  character, dimension(:), intent(out) :: cname
  integer i, ix
  ix = len_trim(fname)
  do i = 1,len_trim(fname)
    cname(i) = fname(i:i)
  end do
  cname(ix+1:ix+1) = C_NULL_CHAR
end subroutine cf_strcopy

subroutine fc_strcopy(cname,fname)
  implicit none
  character(len=*), intent(out) :: fname
  character, dimension(:), intent(in) :: cname
  character(len=len(fname)) temp
  integer ix
  temp = ''
  temp = transfer(cname,temp)
  ix = index(temp,C_NULL_CHAR) - 1
  fname = ''    
  if ( ix>0 ) then
    fname = temp(1:ix)
  else
    fname = temp
  end if
end subroutine fc_strcopy

end module netcdf_m
#else
! Fortran 77 interface
module netcdf_m
public
include 'netcdf.inc'
end module netcdf_m
#endif
