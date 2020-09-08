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
    
module indices_m

implicit none

private
public iw_g,is_g,ise_g,ie_g,in_g,iwn_g,inw_g,isw_g,ies_g,iws_g
public ine_g,ien_g,inn_g,iss_g,iww_g,iee_g,iwu2_g,isv2_g
public lwws_g,lwss_g,lees_g,less_g,lwwn_g,lwnn_g,leen_g,lenn_g,lsww_g
public lssw_g,lsee_g,lsse_g,lnww_g,lnnw_g,lnee_g,lnne_g
public npann_g,npane_g,npanw_g,npans_g
public iw,is,ise,ie,ine,in,iwn,inw,isw,ies,iws,ien,inn,iss,iww,iee
public iwu,isv,ieu,inv,iwwu,issv,ieeu,innv
public iev,iwv,inu,isu,ieev,innu
public lwws,lwss,lees,less,lwwn,lwnn,leen,lenn,lsww
public lssw,lsee,lsse,lnww,lnnw,lnee,lnne
public indices_init,indices_end
public jn_g, je_g, js_g, jw_g
public jne_g, jen_g, jse_g, jes_g, jnw_g, jwn_g, jsw_g, jws_g
public unpack_nsew, unpack_ne
public unpack_nveu, unpack_svwu

integer, dimension(:), allocatable, save :: in,is,ie,iw                     ! default bounds
integer, dimension(:), allocatable, save :: ine,ien,ise,ies,isw,iws,inw,iwn ! corner=.true.
integer, dimension(:), allocatable, save :: inn,iss,iee,iww                 ! nrows=2
integer, dimension(:), allocatable, save :: inv,ieu,iwu,isv                 ! default boundsuv
integer, dimension(:), allocatable, save :: innv,ieeu,iwwu,issv             ! stag or allvec.and.nrows=2
integer, dimension(:), allocatable, save :: inu,iev,iwv,isu                 ! allvec
integer, dimension(:), allocatable, save :: innu,ieev                       ! allvec.and.nrows=2
integer, dimension(:), allocatable, save :: lwws,lwss,lees,less,lwwn,lwnn,leen,lenn,lsww
integer, dimension(:), allocatable, save :: lssw,lsee,lsse,lnww,lnnw,lnee,lnne
!$acc declare create(in,is,ie,iw)
!$acc declare create(iws,ies,iwn,ien,isw,inw,ise,ine)
!$acc declare create(inn,iss,iee,iww)
!$acc declare create(inv,isv,ieu,iwu)
!$acc declare create(innv,ieeu,issv,iwwu)
!$acc declare create(inu,isu,iev,iwv)
!$acc declare create(lwws,lwss,lees,less,lwwn,lwnn,leen,lenn,lsww)
!$acc declare create(lssw,lsee,lsse,lnww,lnnw,lnee,lnne)
integer, parameter, dimension(0:5) :: npann_g = (/ 1, 103, 3, 105, 5, 101 /)
integer, parameter, dimension(0:5) :: npans_g = (/ 104, 0, 100, 2, 102, 4 /)
integer, parameter, dimension(0:5) :: npane_g = (/ 102, 2, 104, 4, 100, 0 /)
integer, parameter, dimension(0:5) :: npanw_g = (/ 5, 105, 1, 101, 3, 103 /)

contains

subroutine indices_init(ifull,npan)

implicit none

integer, intent(in) :: ifull,npan

allocate(iw(ifull),is(ifull),ise(ifull))
allocate(ie(ifull),ine(ifull),in(ifull),iwn(ifull))
allocate(inw(ifull),isw(ifull),ies(ifull),iws(ifull))
allocate(ien(ifull),inn(ifull),iss(ifull),iww(ifull))
allocate(iee(ifull),iwu(ifull),isv(ifull),ieu(ifull))
allocate(inv(ifull),iwwu(ifull),issv(ifull),ieeu(ifull))
allocate(innv(ifull))
allocate(iev(ifull),iwv(ifull),inu(ifull),isu(ifull))
allocate(ieev(ifull),innu(ifull))
allocate(lwws(npan),lwss(npan))
allocate(lees(npan),less(npan),lwwn(npan),lwnn(npan))
allocate(leen(npan),lenn(npan),lsww(npan))
allocate(lssw(npan),lsee(npan),lsse(npan),lnww(npan))
allocate(lnnw(npan),lnee(npan),lnne(npan))

return
end subroutine indices_init

subroutine indices_end

implicit none

deallocate(iw,is,ise,ie,ine,in,iwn,inw,isw,ies,iws,ien,inn,iss,iww,iee,iwu,isv)
deallocate(ieu,inv,iwwu,issv,ieeu,innv)
deallocate(iev,iwv,inu,isu)
deallocate(ieev,innu)
deallocate(lwws,lwss,lees,less,lwwn,lwnn,leen,lenn,lsww)
deallocate(lssw,lsee,lsse,lnww,lnnw,lnee,lnne)

return
end subroutine indices_end

subroutine unpack_nsew(data_in,data_n,data_s,data_e,data_w)

use newmpar_m

implicit none

integer i, j, n, iq, ipan, jpan
real, dimension(ifull+iextra), intent(in) :: data_in
real, dimension(ifull), intent(out) :: data_n, data_s, data_e, data_w

ipan = il
jpan = jl/npan

data_e(1:ifull-1)    = data_in(2:ifull)
data_w(2:ifull)      = data_in(1:ifull-1)
data_n(1:ifull-ipan) = data_in(ipan+1:ifull)
data_s(ipan+1:ifull) = data_in(1:ifull-ipan)
do n = 1,npan
  do j = 1,jpan
    iq = 1 + (j-1)*ipan + (n-1)*ipan*jpan
    data_w(iq) = data_in(iw(iq))
    iq = j*ipan + (n-1)*ipan*jpan
    data_e(iq) = data_in(ie(iq))
  end do
  do i = 1,ipan
    iq = i + (n-1)*ipan*jpan
    data_s(iq) = data_in(is(iq))
    iq = i - ipan + n*ipan*jpan
    data_n(iq) = data_in(in(iq))
  end do
end do
    
return
end subroutine unpack_nsew

subroutine unpack_ne(data_in,data_n,data_e)

use newmpar_m

implicit none

integer i, j, n, iq, ipan, jpan
real, dimension(ifull+iextra), intent(in) :: data_in
real, dimension(ifull), intent(out) :: data_n, data_e

ipan = il
jpan = jl/npan

data_e(1:ifull-1)    = data_in(2:ifull)
data_n(1:ifull-ipan) = data_in(ipan+1:ifull)
do n = 1,npan
  do j = 1,jpan
    iq = j*ipan + (n-1)*ipan*jpan
    data_e(iq) = data_in(ie(iq))
  end do
  do i = 1,ipan
    iq = i - ipan + n*ipan*jpan
    data_n(iq) = data_in(in(iq))
  end do
end do
    
return
end subroutine unpack_ne

subroutine unpack_nveu(data_in_u,data_in_v,data_nv,data_eu)

use newmpar_m

implicit none

integer i, j, n, iq, ipan, jpan
real, dimension(ifull+iextra), intent(in) :: data_in_u, data_in_v
real, dimension(ifull), intent(out) :: data_nv, data_eu

ipan = il
jpan = jl/npan

data_eu(1:ifull-1) = data_in_u(2:ifull)
data_nv(1:ifull-ipan) = data_in_v(ipan+1:ifull)
do n = 1,npan
  do j = 1,jpan
    iq = j*ipan + (n-1)*ipan*jpan
    data_eu(iq) = data_in_u(ieu(iq))
  end do
  do i = 1,ipan
    iq = i - ipan + n*ipan*jpan
    data_nv(iq) = data_in_v(inv(iq))
  end do
end do

return
end subroutine unpack_nveu

subroutine unpack_svwu(data_in_u,data_in_v,data_sv,data_wu)

use newmpar_m

implicit none

integer i, j, n, iq, ipan, jpan
real, dimension(ifull+iextra), intent(in) :: data_in_u, data_in_v
real, dimension(ifull), intent(out) :: data_sv, data_wu

ipan = il
jpan = jl/npan

data_wu(2:ifull)      = data_in_u(1:ifull-1)
data_sv(ipan+1:ifull) = data_in_v(1:ifull-ipan)
do n = 1,npan
  do j = 1,jpan
    iq = 1 + (j-1)*ipan + (n-1)*ipan*jpan
    data_wu(iq) = data_in_u(iwu(iq))
  end do
  do i = 1,ipan
    iq = i + (n-1)*ipan*jpan
    data_sv(iq) = data_in_v(isv(iq))
  end do
end do

return
end subroutine unpack_svwu

function in_g(iq) result(iqq)
use newmpar_m
implicit none
integer, intent(in) :: iq
integer iqq
integer n,i,j
n=(iq-1)/(il_g*il_g)
j=(iq-1-n*il_g*il_g)/il_g+1
i=iq-(j-1)*il_g-n*il_g*il_g
if (j==il_g) then
  if (npann_g(n)<100) then
    iqq=i+npann_g(n)*il_g*il_g
  else
    iqq=1+(il_g-i)*il_g+(npann_g(n)-100)*il_g*il_g
  end if
else
  iqq=iq+il_g    
end if
end function in_g

function is_g(iq) result(iqq)
use newmpar_m
implicit none
integer, intent(in) :: iq
integer iqq
integer n,i,j
n=(iq-1)/(il_g*il_g)
j=(iq-1-n*il_g*il_g)/il_g+1
i=iq-(j-1)*il_g-n*il_g*il_g
if (j==1) then
  if (npans_g(n)<100) then
    iqq=i+(il_g-1)*il_g+npans_g(n)*il_g*il_g
  else
    iqq=il_g+(il_g-i)*il_g+(npans_g(n)-100)*il_g*il_g
  end if
else
  iqq=iq-il_g    
end if
end function is_g

function ie_g(iq) result(iqq)
use newmpar_m
implicit none
integer, intent(in) :: iq
integer iqq
integer n,i,j
n=(iq-1)/(il_g*il_g)
j=(iq-1-n*il_g*il_g)/il_g+1
i=iq-(j-1)*il_g-n*il_g*il_g
if (i==il_g) then
  if (npane_g(n)<100) then
    iqq=1+(j-1)*il_g+npane_g(n)*il_g*il_g
  else
    iqq=il_g+1-j+(npane_g(n)-100)*il_g*il_g
  end if
else
  iqq=iq+1
end if
end function ie_g

function iw_g(iq) result(iqq)
use newmpar_m
implicit none
integer, intent(in) :: iq
integer iqq
integer n,i,j
n=(iq-1)/(il_g*il_g)
j=(iq-1-n*il_g*il_g)/il_g+1
i=iq-(j-1)*il_g-n*il_g*il_g
if (i==1) then
  if (npanw_g(n)<100) then
    iqq=il_g+(j-1)*il_g+npanw_g(n)*il_g*il_g
  else
    iqq=il_g+1-j+(il_g-1)*il_g+(npanw_g(n)-100)*il_g*il_g
  end if
else
  iqq=iq-1
end if
end function iw_g

function inn_g(iq) result(iqq)
use newmpar_m
implicit none
integer, intent(in) :: iq
integer iqq
integer n,i,j
n=(iq-1)/(il_g*il_g)
j=(iq-1-n*il_g*il_g)/il_g+1
i=iq-(j-1)*il_g-n*il_g*il_g
if (npann_g(n)>=100.and.j==il_g) then
  iqq=ie_g(in_g(iq))
else
  iqq=in_g(in_g(iq))
end if
end function inn_g

function iss_g(iq) result(iqq)
use newmpar_m
implicit none
integer, intent(in) :: iq
integer iqq
integer n,i,j
n=(iq-1)/(il_g*il_g)
j=(iq-1-n*il_g*il_g)/il_g+1
i=iq-(j-1)*il_g-n*il_g*il_g
if (npans_g(n)>=100.and.j==1) then
  iqq=iw_g(is_g(iq))
else
  iqq=is_g(is_g(iq))
end if
end function iss_g

function iee_g(iq) result(iqq)
use newmpar_m
implicit none
integer, intent(in) :: iq
integer iqq
integer n,i,j
n=(iq-1)/(il_g*il_g)
j=(iq-1-n*il_g*il_g)/il_g+1
i=iq-(j-1)*il_g-n*il_g*il_g
if (npane_g(n)>=100.and.i==il_g) then
  iqq=in_g(ie_g(iq))
else
  iqq=ie_g(ie_g(iq))
end if
end function iee_g

function iww_g(iq) result(iqq)
use newmpar_m
implicit none
integer, intent(in) :: iq
integer iqq
integer n,i,j
n=(iq-1)/(il_g*il_g)
j=(iq-1-n*il_g*il_g)/il_g+1
i=iq-(j-1)*il_g-n*il_g*il_g
if (npanw_g(n)>=100.and.i==1) then
  iqq=is_g(iw_g(iq))
else
  iqq=iw_g(iw_g(iq))
end if
end function iww_g

function ine_g(iq) result(iqq)
use newmpar_m
implicit none
integer, intent(in) :: iq
integer iqq
integer n,i,j
n=(iq-1)/(il_g*il_g)
j=(iq-1-n*il_g*il_g)/il_g+1
i=iq-(j-1)*il_g-n*il_g*il_g
if (npane_g(n)>=100.and.i==il_g) then
  iqq=iw_g(ie_g(iq))
else
  iqq=in_g(ie_g(iq))
end if
end function ine_g

function ise_g(iq) result(iqq)
use newmpar_m
implicit none
integer, intent(in) :: iq
integer iqq
integer n,i,j
n=(iq-1)/(il_g*il_g)
j=(iq-1-n*il_g*il_g)/il_g+1
i=iq-(j-1)*il_g-n*il_g*il_g
if (npane_g(n)>=100.and.i==il_g) then
  iqq=ie_g(ie_g(iq))
else
  iqq=is_g(ie_g(iq))
end if
end function ise_g

function ien_g(iq) result(iqq)
use newmpar_m
implicit none
integer, intent(in) :: iq
integer iqq
integer n,i,j
n=(iq-1)/(il_g*il_g)
j=(iq-1-n*il_g*il_g)/il_g+1
i=iq-(j-1)*il_g-n*il_g*il_g
if (npann_g(n)>=100.and.j==il_g) then
  iqq=is_g(in_g(iq))
else
  iqq=ie_g(in_g(iq))
end if
end function ien_g

function iwn_g(iq) result(iqq)
use newmpar_m
implicit none
integer, intent(in) :: iq
integer iqq
integer n,i,j
n=(iq-1)/(il_g*il_g)
j=(iq-1-n*il_g*il_g)/il_g+1
i=iq-(j-1)*il_g-n*il_g*il_g
if (npann_g(n)>=100.and.j==il_g) then
  iqq=in_g(in_g(iq))
else
  iqq=iw_g(in_g(iq))
end if
end function iwn_g

function inw_g(iq) result(iqq)
use newmpar_m
implicit none
integer, intent(in) :: iq
integer iqq
integer n,i,j
n=(iq-1)/(il_g*il_g)
j=(iq-1-n*il_g*il_g)/il_g+1
i=iq-(j-1)*il_g-n*il_g*il_g
if (npanw_g(n)>=100.and.i==1) then
  iqq=iw_g(iw_g(iq))
else
  iqq=in_g(iw_g(iq))
end if
end function inw_g

function isw_g(iq) result(iqq)
use newmpar_m
implicit none
integer, intent(in) :: iq
integer iqq
integer n,i,j
n=(iq-1)/(il_g*il_g)
j=(iq-1-n*il_g*il_g)/il_g+1
i=iq-(j-1)*il_g-n*il_g*il_g
if (npanw_g(n)>=100.and.i==1) then
  iqq=ie_g(iw_g(iq))
else
  iqq=is_g(iw_g(iq))
end if
end function isw_g

function ies_g(iq) result(iqq)
use newmpar_m
implicit none
integer, intent(in) :: iq
integer iqq
integer n,i,j
n=(iq-1)/(il_g*il_g)
j=(iq-1-n*il_g*il_g)/il_g+1
i=iq-(j-1)*il_g-n*il_g*il_g
if (npans_g(n)>=100.and.j==1) then
  iqq=is_g(is_g(iq))
else
  iqq=ie_g(is_g(iq))
end if
end function ies_g

function iws_g(iq) result(iqq)
use newmpar_m
implicit none
integer, intent(in) :: iq
integer iqq
integer n,i,j
n=(iq-1)/(il_g*il_g)
j=(iq-1-n*il_g*il_g)/il_g+1
i=iq-(j-1)*il_g-n*il_g*il_g
if (npans_g(n)>=100.and.j==1) then
  iqq=in_g(is_g(iq))
else
  iqq=iw_g(is_g(iq))
end if
end function iws_g

function iwu2_g(iq) result(iqq)
use newmpar_m
implicit none
integer, intent(in) :: iq
integer iqq
integer n,i,j
n=(iq-1)/(il_g*il_g)
j=(iq-1-n*il_g*il_g)/il_g+1
i=iq-(j-1)*il_g-n*il_g*il_g
if (npanw_g(n)>=100.and.i==1) then
  iqq=iw_g(iq)+ifull_g
else
  iqq=iw_g(iq)
end if
end function iwu2_g

function isv2_g(iq) result(iqq)
use newmpar_m
implicit none
integer, intent(in) :: iq
integer iqq
integer n,i,j
n=(iq-1)/(il_g*il_g)
j=(iq-1-n*il_g*il_g)/il_g+1
i=iq-(j-1)*il_g-n*il_g*il_g
if (npans_g(n)>=100.and.j==1) then
  iqq=is_g(iq)-ifull_g
else
  iqq=is_g(iq)
end if
end function isv2_g

function leen_g(n) result(iqq)
use newmpar_m
implicit none
integer, intent(in) :: n
integer iqq
if (npann_g(n)>=100) then
  iqq=iss_g(in_g(il_g+(il_g-1)*il_g+n*il_g*il_g))
else
  iqq=iee_g(in_g(il_g+(il_g-1)*il_g+n*il_g*il_g))
end if
end function leen_g

function lenn_g(n) result(iqq)
use newmpar_m
implicit none
integer, intent(in) :: n
integer iqq
if (npann_g(n)>=100) then
  iqq=ise_g(in_g(il_g+(il_g-1)*il_g+n*il_g*il_g))
else
  iqq=ien_g(in_g(il_g+(il_g-1)*il_g+n*il_g*il_g))
end if
end function lenn_g

function lwnn_g(n) result(iqq)
use newmpar_m
implicit none
integer, intent(in) :: n
integer iqq
if (npann_g(n)>=100) then
  iqq=ine_g(in_g(1+(il_g-1)*il_g+n*il_g*il_g))
else
  iqq=iwn_g(in_g(1+(il_g-1)*il_g+n*il_g*il_g))
end if
end function lwnn_g

function lsee_g(n) result(iqq)
use newmpar_m
implicit none
integer, intent(in) :: n
integer iqq
if (npane_g(n)>=100) then
  iqq=iwn_g(ie_g(il_g+n*il_g*il_g))
else
  iqq=ise_g(ie_g(il_g+n*il_g*il_g))
end if
end function lsee_g

function lnee_g(n) result(iqq)
use newmpar_m
implicit none
integer, intent(in) :: n
integer iqq
if (npane_g(n)>=100) then
  iqq=iwn_g(ie_g(il_g+(il_g-1)*il_g+n*il_g*il_g))
else
  iqq=ine_g(ie_g(il_g+(il_g-1)*il_g+n*il_g*il_g))
end if
end function lnee_g

function lnne_g(n) result(iqq)
use newmpar_m
implicit none
integer, intent(in) :: n
integer iqq
if (npane_g(n)>=100) then
  iqq=iww_g(ie_g(il_g+(il_g-1)*il_g+n*il_g*il_g))
else
  iqq=inn_g(ie_g(il_g+(il_g-1)*il_g+n*il_g*il_g))
end if
end function lnne_g

function lsww_g(n) result(iqq)
use newmpar_m
implicit none
integer, intent(in) :: n
integer iqq
if (npanw_g(n)>=100) then
  iqq=ies_g(iw_g(1+n*il_g*il_g))
else
  iqq=isw_g(iw_g(1+n*il_g*il_g))
end if
end function lsww_g

function lssw_g(n) result(iqq)
use newmpar_m
implicit none
integer, intent(in) :: n
integer iqq
if (npanw_g(n)>=100) then
  iqq=iee_g(iw_g(1+n*il_g*il_g))
else
  iqq=iss_g(iw_g(1+n*il_g*il_g))
end if
end function lssw_g

function lnww_g(n) result(iqq)
use newmpar_m
implicit none
integer, intent(in) :: n
integer iqq
if (npanw_g(n)>=100) then
  iqq=iws_g(iw_g(1+(il_g-1)*il_g+n*il_g*il_g))
else
  iqq=inw_g(iw_g(1+(il_g-1)*il_g+n*il_g*il_g))
end if
end function lnww_g

function lwws_g(n) result(iqq)
use newmpar_m
implicit none
integer, intent(in) :: n
integer iqq
if (npans_g(n)>=100) then
  iqq=inn_g(is_g(1+n*il_g*il_g))
else
  iqq=iww_g(is_g(1+n*il_g*il_g))
end if
end function lwws_g

function lwss_g(n) result(iqq)
use newmpar_m
implicit none
integer, intent(in) :: n
integer iqq
if (npans_g(n)>=100) then
  iqq=inw_g(is_g(1+n*il_g*il_g))
else
  iqq=iws_g(is_g(1+n*il_g*il_g))
end if
end function lwss_g

function less_g(n) result(iqq)
use newmpar_m
implicit none
integer, intent(in) :: n
integer iqq
if (npans_g(n)>=100) then
  iqq=isw_g(is_g(il_g+n*il_g*il_g))
else
  iqq=ies_g(is_g(il_g+n*il_g*il_g))
end if
end function less_g

function lwwn_g(n) result(iqq)
use newmpar_m
implicit none
integer, intent(in) :: n
integer iqq
if (npann_g(n)>=100) then
  iqq=inn_g(in_g(1+(il_g-1)*il_g+n*il_g*il_g))
else
  iqq=iww_g(in_g(1+(il_g-1)*il_g+n*il_g*il_g))
end if
end function lwwn_g

function lsse_g(n) result(iqq)
use newmpar_m
implicit none
integer, intent(in) :: n
integer iqq
if (npane_g(n)>=100) then
  iqq=iee_g(ie_g(il_g+n*il_g*il_g))
else
  iqq=iss_g(ie_g(il_g+n*il_g*il_g))
end if
end function lsse_g

function lnnw_g(n) result(iqq)
use newmpar_m
implicit none
integer, intent(in) :: n
integer iqq
if (npanw_g(n)>=100) then
  iqq=iww_g(iw_g(1+(il_g-1)*il_g+n*il_g*il_g))
else
  iqq=inn_g(iw_g(1+(il_g-1)*il_g+n*il_g*il_g))
end if
end function lnnw_g

function lees_g(n) result(iqq)
use newmpar_m
implicit none
integer, intent(in) :: n
integer iqq
if (npans_g(n)>=100) then
  iqq=iss_g(is_g(il_g+n*il_g*il_g))
else
  iqq=iee_g(is_g(il_g+n*il_g*il_g))
end if
end function lees_g

function jn_g(iq,mil_g) result(iqq)
implicit none
integer, intent(in) :: iq, mil_g
integer iqq
integer n, i, j
n = (iq-1)/(mil_g*mil_g)
j = (iq-1-n*mil_g*mil_g)/mil_g + 1
i = iq - (j-1)*mil_g - n*mil_g*mil_g
if ( j==mil_g ) then
  if ( npann_g(n)<100 ) then
    iqq = i + npann_g(n)*mil_g*mil_g
  else
    iqq = 1 + (mil_g-i)*mil_g + (npann_g(n)-100)*mil_g*mil_g
  end if
else
  iqq = iq + mil_g    
end if
end function jn_g

function je_g(iq,mil_g) result(iqq)
implicit none
integer, intent(in) :: iq, mil_g
integer iqq
integer n, i, j
n = (iq-1)/(mil_g*mil_g)
j = (iq-1-n*mil_g*mil_g)/mil_g + 1
i = iq - (j-1)*mil_g - n*mil_g*mil_g
if ( i==mil_g ) then
  if ( npane_g(n)<100 ) then
    iqq = 1 + (j-1)*mil_g + npane_g(n)*mil_g*mil_g
  else
    iqq = mil_g + 1 - j + (npane_g(n)-100)*mil_g*mil_g
  end if
else
  iqq = iq + 1
end if
end function je_g

function js_g(iq,mil_g) result(iqq)
implicit none
integer, intent(in) :: iq, mil_g
integer iqq
integer n, i, j
n = (iq-1)/(mil_g*mil_g)
j = (iq-1-n*mil_g*mil_g)/mil_g + 1
i = iq - (j-1)*mil_g - n*mil_g*mil_g
if ( j==1 ) then
  if ( npans_g(n)<100 ) then
    iqq = i + (mil_g-1)*mil_g + npans_g(n)*mil_g*mil_g
  else
    iqq = mil_g + (mil_g-i)*mil_g + (npans_g(n)-100)*mil_g*mil_g
  end if
else
  iqq = iq - mil_g    
end if
end function js_g

function jw_g(iq,mil_g) result(iqq)
implicit none
integer, intent(in) :: iq, mil_g
integer iqq
integer n, i, j
n = (iq-1)/(mil_g*mil_g)
j = (iq-1-n*mil_g*mil_g)/mil_g + 1
i = iq - (j-1)*mil_g - n*mil_g*mil_g
if ( i==1 ) then
  if ( npanw_g(n)<100 ) then
    iqq = mil_g + (j-1)*mil_g + npanw_g(n)*mil_g*mil_g
  else
    iqq = mil_g + 1 - j + (mil_g-1)*mil_g + (npanw_g(n)-100)*mil_g*mil_g
  end if
else
  iqq = iq - 1
end if
end function jw_g

function jne_g(iq,mil_g) result(iqq)
implicit none
integer, intent(in) :: iq, mil_g
integer iqq
integer n, i, j
n = (iq-1)/(mil_g*mil_g)
j = (iq-1-n*mil_g*mil_g)/mil_g + 1
i = iq - (j-1)*mil_g - n*mil_g*mil_g
if ( npane_g(n)>=100 .and. i==mil_g ) then
  iqq = jw_g(je_g(iq,mil_g),mil_g)
else
  iqq = jn_g(je_g(iq,mil_g),mil_g)
end if
end function jne_g

function jen_g(iq,mil_g) result(iqq)
implicit none
integer, intent(in) :: iq, mil_g
integer iqq
integer n, i, j
n = (iq-1)/(mil_g*mil_g)
j = (iq-1-n*mil_g*mil_g)/mil_g + 1
i = iq - (j-1)*mil_g - n*mil_g*mil_g
if ( npann_g(n)>=100 .and. j==mil_g ) then
  iqq = js_g(jn_g(iq,mil_g),mil_g)
else
  iqq = je_g(jn_g(iq,mil_g),mil_g)
end if
end function jen_g

function jse_g(iq,mil_g) result(iqq)
implicit none
integer, intent(in) :: iq, mil_g
integer iqq
integer n, i, j
n = (iq-1)/(mil_g*mil_g)
j = (iq-1-n*mil_g*mil_g)/mil_g + 1
i = iq - (j-1)*mil_g - n*mil_g*mil_g
if ( npane_g(n)>=100 .and. i==mil_g ) then
  iqq = je_g(je_g(iq,mil_g),mil_g)
else
  iqq = js_g(je_g(iq,mil_g),mil_g)
end if
end function jse_g

function jes_g(iq,mil_g) result(iqq)
implicit none
integer, intent(in) :: iq, mil_g
integer iqq
integer n, i, j
n = (iq-1)/(mil_g*mil_g)
j = (iq-1-n*mil_g*mil_g)/mil_g + 1
i = iq - (j-1)*mil_g - n*mil_g*mil_g
if ( npans_g(n)>=100 .and. j==1 ) then
  iqq = js_g(js_g(iq,mil_g),mil_g)
else
  iqq = je_g(js_g(iq,mil_g),mil_g)
end if
end function jes_g

function jnw_g(iq,mil_g) result(iqq)
implicit none
integer, intent(in) :: iq, mil_g
integer iqq
integer n, i, j
n = (iq-1)/(mil_g*mil_g)
j = (iq-1-n*mil_g*mil_g)/mil_g + 1
i = iq - (j-1)*mil_g - n*mil_g*mil_g
if ( npanw_g(n)>=100 .and. i==1 ) then
  iqq = jw_g(jw_g(iq,mil_g),mil_g)
else
  iqq = jn_g(jw_g(iq,mil_g),mil_g)
end if
end function jnw_g

function jwn_g(iq,mil_g) result(iqq)
implicit none
integer, intent(in) :: iq, mil_g
integer iqq
integer n, i, j
n = (iq-1)/(mil_g*mil_g)
j = (iq-1-n*mil_g*mil_g)/mil_g + 1
i = iq - (j-1)*mil_g - n*mil_g*mil_g
if ( npann_g(n)>=100 .and. j==mil_g ) then
  iqq = jn_g(jn_g(iq,mil_g),mil_g)
else
  iqq = jw_g(jn_g(iq,mil_g),mil_g)
end if
end function jwn_g

function jsw_g(iq,mil_g) result(iqq)
implicit none
integer, intent(in) :: iq, mil_g
integer iqq
integer n, i, j
n = (iq-1)/(mil_g*mil_g)
j = (iq-1-n*mil_g*mil_g)/mil_g + 1
i = iq - (j-1)*mil_g - n*mil_g*mil_g
if ( npanw_g(n)>=100 .and. i==1 ) then
  iqq = je_g(jw_g(iq,mil_g),mil_g)
else
  iqq = js_g(jw_g(iq,mil_g),mil_g)
end if
end function jsw_g

function jws_g(iq,mil_g) result(iqq)
implicit none
integer, intent(in) :: iq, mil_g
integer iqq
integer n, i, j
n = (iq-1)/(mil_g*mil_g)
j = (iq-1-n*mil_g*mil_g)/mil_g + 1
i = iq - (j-1)*mil_g - n*mil_g*mil_g
if ( npans_g(n)>=100 .and. j==1 ) then
  iqq = jn_g(js_g(iq,mil_g),mil_g)
else
  iqq = jw_g(js_g(iq,mil_g),mil_g)
end if
end function jws_g

end module indices_m
