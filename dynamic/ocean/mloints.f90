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
    
module mloints

implicit none

private
public mlob2ints_bs, mlofill, mloclean
public mlo_miss, mlo_sal, mlo_zero 

integer, parameter :: mlo_miss = 0
integer, parameter :: mlo_sal = 1
integer, parameter :: mlo_zero = 2
real, parameter :: cxx = -9999. ! missing value flag

interface mlob2ints_bs
  module procedure mlob2ints_bs_2, mlob2ints_bs_3
end interface

interface mlofill
  module procedure mlofill_2, mlofill_3
end interface

interface mloclean
  module procedure mloclean_2, mloclean_3
end interface

contains


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Calculate depature points for MLO semi-Lagrangian advection

subroutine mlob2ints_bs_3(s,nface,xg,yg,wtr,bs_test,mlointschf)

use cc_acc
use cc_mpi
use newmpar_m
use parm_m
use sflux_m

implicit none

integer, intent(in) :: mlointschf
integer intsch, nn
integer ntr, nstart, nend, nlen
integer async_counter
integer, dimension(ifull,ol), intent(in) :: nface
real, dimension(ifull,ol), intent(in) :: xg, yg
real, dimension(:,:,:), intent(inout) :: s
real, dimension(-1:ipan+2,-1:jpan+2,1:npan,ol,nagg) :: sx
logical, dimension(ifull+iextra,ol), intent(in) :: wtr
logical, dimension(-1:ipan+2,-1:jpan+2,1:npan,ol) :: wx
logical, intent(in) :: bs_test

call START_LOG(waterints_begin)

ntr = size(s,3)

if ( mlointschf==0 ) then
  intsch = 0
else if ( mlointschf>0 ) then
  intsch = mod(ktau,mlointschf)
else
  write(6,*) "ERROR: Unknown value for mlointschf = ",mlointschf  
  call ccmpi_abort(-1)
end if

! bounds(nrows=2) must be called prior

!$acc enter data create(wx,sx)


call setwx(wx,wtr,intsch)
!$acc update device(wx)
  
do nstart = 1,ntr,nagg
  nend = min(nstart + nagg - 1, ntr)
  nlen = nend - nstart + 1

  do nn = 1,nlen
    async_counter = mod(nn-1, async_length)
    call setsx(sx(:,:,:,:,nn),s(:,:,nn-1+nstart),intsch)      
    !$acc update device(sx(:,:,:,:,nn)) async(async_counter)
  end do

  ! Loop over points that need to be calculated for other processes
  do nn = 1,nlen
    call calcdepproc(sx(:,:,:,:,nn),wx,nface,xg,yg,nn,bs_test,intsch)
  end do  

  ! Send messages to other processes.  We then start the calculation for this processor
  ! while waiting for the message to return, thereby overlapping computation with
  ! communication
  call intssync_send(nlen)
  
  ! Loop over points that need to be calculated for this process
  !$omp parallel
  do nn = 1,nlen
    async_counter = mod(nn-1, async_length)
    call calcdephost(sx(:,:,:,:,nn),wx,s(:,:,nn-1+nstart),nface,xg,yg,nn,bs_test,intsch,async_counter)
  end do
  !$acc wait
  !$omp end parallel

  call intssync_recv(s(:,:,nstart:nend))  
    
end do ! nstart


!$acc exit data delete(wx,sx)

call END_LOG(waterints_end)

return
end subroutine mlob2ints_bs_3

subroutine mlob2ints_bs_2(sa,sb,sc,sd,nface,xg,yg,wtr,bs_test,mlointschf)

use cc_acc
use cc_mpi
use newmpar_m
use parm_m
use sflux_m

implicit none

integer, intent(in) :: mlointschf
integer intsch, nn, nlen
integer async_counter
integer, dimension(ifull,ol), intent(in) :: nface
real, dimension(ifull,ol), intent(in) :: xg, yg
real, dimension(:,:), intent(inout) :: sa
real, dimension(:,:), intent(inout), optional :: sb, sc, sd
real, dimension(-1:ipan+2,-1:jpan+2,1:npan,ol,nagg) :: sx
logical, dimension(ifull+iextra,ol), intent(in) :: wtr
logical, dimension(-1:ipan+2,-1:jpan+2,1:npan,ol) :: wx
logical, intent(in) :: bs_test

call START_LOG(waterints_begin)

if ( mlointschf==0 ) then
  intsch = 0
else if ( mlointschf>0 ) then
  intsch = mod(ktau,mlointschf)
else
  write(6,*) "ERROR: Unknown value for mlointschf = ",mlointschf  
  call ccmpi_abort(-1)
end if

! bounds(nrows=2) must be called prior

!$acc enter data create(wx,sx)


call setwx(wx,wtr,intsch)
!$acc update device(wx)
  
async_counter = 0
nn = 1
call setsx(sx(:,:,:,:,nn),sa,intsch)      
!$acc update device(sx(:,:,:,:,nn)) async(async_counter)

if ( present(sb) ) then
  async_counter = mod(async_counter + 1, async_length)
  nn = nn + 1
  call setsx(sx(:,:,:,:,nn),sb,intsch)      
  !$acc update device(sx(:,:,:,:,nn)) async(async_counter)    
end if

if ( present(sc) ) then
  async_counter = mod(async_counter + 1, async_length)
  nn = nn + 1
  call setsx(sx(:,:,:,:,nn),sc,intsch)      
  !$acc update device(sx(:,:,:,:,nn)) async(async_counter)    
end if

if ( present(sd) ) then
  async_counter = mod(async_counter + 1, async_length)
  nn = nn + 1
  call setsx(sx(:,:,:,:,nn),sd,intsch)      
  !$acc update device(sx(:,:,:,:,nn)) async(async_counter)    
end if

! Loop over points that need to be calculated for other processes
nn = 1
call calcdepproc(sx(:,:,:,:,nn),wx,nface,xg,yg,nn,bs_test,intsch)

if ( present(sb) ) then
  nn = nn + 1
  call calcdepproc(sx(:,:,:,:,nn),wx,nface,xg,yg,nn,bs_test,intsch)   
end if

if ( present(sc) ) then
  nn = nn + 1
  call calcdepproc(sx(:,:,:,:,nn),wx,nface,xg,yg,nn,bs_test,intsch)   
end if

if ( present(sd) ) then
  nn = nn + 1
  call calcdepproc(sx(:,:,:,:,nn),wx,nface,xg,yg,nn,bs_test,intsch)   
end if

! Send messages to other processes.  We then start the calculation for this processor
! while waiting for the message to return, thereby overlapping computation with
! communication
nlen = nn
call intssync_send(nlen)
  
!$omp parallel

! Loop over points that need to be calculated for this process
async_counter = 0
nn = 1
call calcdephost(sx(:,:,:,:,nn),wx,sa,nface,xg,yg,nn,bs_test,intsch,async_counter)

if ( present(sb) ) then
  async_counter = mod(async_counter + 1, async_length)
  nn = nn + 1
  call calcdephost(sx(:,:,:,:,nn),wx,sb,nface,xg,yg,nn,bs_test,intsch,async_counter)    
end if

if ( present(sc) ) then
  async_counter = mod(async_counter + 1, async_length)
  nn = nn + 1
  call calcdephost(sx(:,:,:,:,nn),wx,sc,nface,xg,yg,nn,bs_test,intsch,async_counter)    
end if

if ( present(sd) ) then
  async_counter = mod(async_counter + 1, async_length)
  nn = nn + 1
  call calcdephost(sx(:,:,:,:,nn),wx,sd,nface,xg,yg,nn,bs_test,intsch,async_counter)    
end if

!$acc wait
!$omp end parallel

if ( present(sd) .and. present(sc) .and. present(sb) ) then
  call intssync_recv(ta=sa,tb=sb,tc=sc,td=sd)
else if ( present(sc) .and. present(sb) ) then
  call intssync_recv(ta=sa,tb=sb,tc=sc)    
else if ( present(sb) ) then
  call intssync_recv(ta=sa,tb=sb)
else
  call intssync_recv(ta=sa)
end if


!$acc exit data delete(wx,sx)

call END_LOG(waterints_end)


end subroutine mlob2ints_bs_2

subroutine mlofill_3(s,wtr,bc_test)

use cc_mpi
use indices_m
use newmpar_m

implicit none

integer, intent(in) :: bc_test
integer ntr, nn, k, ii, iq, nlen, nstart, nend, np
real, dimension(:,:,:), intent(inout) :: s
real, dimension(ifull+iextra,ol,nagg) :: s_old
real s_tot, s_count
logical, dimension(ifull+iextra,ol), intent(in) :: wtr
logical need_fill

call START_LOG(waterints_begin)

ntr = size(s,3)
need_fill = .false.

select case(bc_test)
  case(mlo_miss) ! fill
    do nn = 1,ntr
      do k = 1,ol
        where (.not.wtr(1:ifull,k))
          s(1:ifull,k,nn) = cxx - 1. ! missing value flag
        end where
      end do
    end do
    need_fill = .false.
  case(mlo_sal) ! salinity + fill
    do nn = 1,ntr
      do k = 1,ol  
        where ( s(1:ifull,k,nn)<2. .or. .not.wtr(1:ifull,k) )
          s(1:ifull,k,nn) = cxx - 1. ! missing value flag
        end where
      end do
    end do
    need_fill = .true.
  case(mlo_zero) ! zero
    do nn = 1,ntr
      do k = 1,ol
        where (.not.wtr(1:ifull,k))
          s(1:ifull,k,nn) = 0.
        end where
      end do
    end do
    need_fill = .false.
end select

! fill
if ( need_fill ) then
  do nstart = 1,ntr,nagg
    nend = min(nstart + nagg - 1, ntr)
    nlen = nend - nstart + 1
    do ii = 1,6 ! 6 iterations of fill should be enough
      s_old(1:ifull,:,1:nlen) = s(1:ifull,:,nstart:nend)
      call bounds(s_old(:,:,1:nlen))
      !$omp parallel do collapse(2) private(nn,k,iq,s_tot,s_count)
      do nn = 1,nlen
        do k = 1,ol
          np = nn - 1 + nstart  
          do iq = 1,ifull
            if ( s(iq,k,np)<cxx ) then
              s_tot = 0.
              s_count = 0.
              if ( s_old(is(iq),k,nn)>=cxx ) then
                s_tot = s_tot + s_old(is(iq),k,nn)
                s_count = s_count + 1.
              end if
              if ( s_old(in(iq),k,nn)>=cxx ) then
                s_tot = s_tot + s_old(in(iq),k,nn)
                s_count = s_count + 1.
              end if
              if ( s_old(iw(iq),k,nn)>=cxx ) then
                s_tot = s_tot + s_old(iw(iq),k,nn)
                s_count = s_count + 1.
              end if
              if ( s_old(ie(iq),k,nn)>=cxx ) then
                s_tot = s_tot + s_old(ie(iq),k,nn)
                s_count = s_count + 1.
              end if
              if ( s_count>0. ) then
                s(iq,k,np) = s_tot/s_count
              end if  
            end if
          end do  ! iq loop
        end do    ! k loop
      end do      ! nn loop
      !$omp end parallel do
    end do ! ii loop
  end do   ! nstart loop
end if   ! need_fill

call bounds(s,nrows=2)

call END_LOG(waterints_end)

end subroutine mlofill_3

subroutine mlofill_2(s,wtr,bc_test)

use newmpar_m

implicit none

integer, intent(in) :: bc_test
real, dimension(:,:), intent(inout) :: s
real, dimension(size(s,1),size(s,2),1) :: s_work
logical, dimension(ifull+iextra,ol), intent(in) :: wtr

s_work(:,:,1) = s(:,:)
call mlofill_3(s_work,wtr,bc_test)
s(:,:) = s_work(:,:,1)

end subroutine mlofill_2

subroutine mloclean_3(s,s_store,wtr,bc_test)

use cc_mpi
use newmpar_m

implicit none

integer, intent(in) :: bc_test
integer nn, ntr
real, dimension(:,:,:), intent(inout) :: s
real, dimension(:,:,:), intent(in) :: s_store
logical, dimension(ifull+iextra,ol), intent(in) :: wtr

call START_LOG(waterints_begin)

ntr = size(s,3)

do nn = 1,ntr
  where ( .not.wtr(1:ifull,1:ol) .or. s(1:ifull,1:ol,nn)<cxx )
    s(1:ifull,1:ol,nn) = s_store(1:ifull,1:ol,nn)
  end where
end do
if ( bc_test==mlo_sal ) then
  do nn = 1,ntr
    where ( s_store(1:ifull,1:ol,nn)<2. )
      s(1:ifull,1:ol,nn) = s_store(1:ifull,1:ol,nn)  
    end where
  end do
end if

call END_LOG(waterints_end)

end subroutine mloclean_3

subroutine mloclean_2(s,s_store,wtr,bc_test)

use newmpar_m

implicit none

integer, intent(in) :: bc_test
real, dimension(:,:), intent(inout) :: s
real, dimension(:,:), intent(in) :: s_store
real, dimension(size(s,1),size(s,2),1) :: s_work
real, dimension(size(s_store,1),size(s_store,2),1) :: s_store_work
logical, dimension(ifull+iextra,ol), intent(in) :: wtr

s_work(:,:,1) = s(:,:)
s_store_work(:,:,1) = s_store(:,:)
call mloclean_3(s_work,s_store_work,wtr,bc_test)
s(:,:) = s_work(:,:,1)

end subroutine mloclean_2

subroutine setwx(wx,wtr,intsch)

use cc_mpi             ! CC MPI routines
use indices_m          ! Grid index arrays
use newmpar_m          ! Grid parameters

implicit none

integer, intent(in) :: intsch
integer k, n, i, j, iq
logical, dimension(ifull+iextra,ol), intent(in) :: wtr
logical, dimension(-1:ipan+2,-1:jpan+2,1:npan,1:ol), intent(out) :: wx

! unpack wtr to wx for interpolation
wx = .false.

!======================== start of intsch=1 section ====================
if ( intsch==1 ) then
    
  do k = 1,ol
    wx(1:ipan,1:jpan,1:npan,k) = &
      reshape( wtr(1:ipan*jpan*npan,k), (/ ipan, jpan, npan /) )
    do n = 1,npan
      do j = 1,jpan
        iq = 1+(j-1)*ipan+(n-1)*ipan*jpan
        wx(0,j,n,k)      = wtr( iw(iq),k)
        wx(-1,j,n,k)     = wtr(iww(iq),k)
        iq = j*ipan+(n-1)*ipan*jpan
        wx(ipan+1,j,n,k) = wtr( ie(iq),k)
        wx(ipan+2,j,n,k) = wtr(iee(iq),k)
      end do            ! j loop
      do i = 1,ipan
        iq = i+(n-1)*ipan*jpan
        wx(i,0,n,k)      = wtr( is(iq),k)
        wx(i,-1,n,k)     = wtr(iss(iq),k)
        iq = i-ipan+n*ipan*jpan
        wx(i,jpan+1,n,k) = wtr( in(iq),k)
        wx(i,jpan+2,n,k) = wtr(inn(iq),k)
      end do            ! i loop
    end do
!   for ew interpolation, sometimes need (different from ns):
!       (-1,0),   (0,0),   (0,-1)   (-1,il+1),   (0,il+1),   (0,il+2)
!     (il+1,0),(il+2,0),(il+1,-1) (il+1,il+1),(il+2,il+1),(il+1,il+2)
    do n = 1,npan
      wx(-1,0,n,k)          = wtr(lwws(n),                  k)
      wx(0,0,n,k)           = wtr(iws(1+(n-1)*ipan*jpan),   k)
      wx(0,-1,n,k)          = wtr(lwss(n),                  k)
      wx(ipan+1,0,n,k)      = wtr(ies(ipan+(n-1)*ipan*jpan),k)
      wx(ipan+2,0,n,k)      = wtr(lees(n),                  k)
      wx(ipan+1,-1,n,k)     = wtr(less(n),                  k)
      wx(-1,jpan+1,n,k)     = wtr(lwwn(n),                  k)
      wx(0,jpan+2,n,k)      = wtr(lwnn(n),                  k)
      wx(ipan+2,jpan+1,n,k) = wtr(leen(n),                  k)
      wx(ipan+1,jpan+2,n,k) = wtr(lenn(n),                  k)
      wx(0,jpan+1,n,k)      = wtr(iwn(1-ipan+n*ipan*jpan),  k)
      wx(ipan+1,jpan+1,n,k) = wtr(ien(n*ipan*jpan),         k)
    end do           ! n loop
  end do             ! k loop

!========================   end of intsch=1 section ====================
else     ! if(intsch==1)then
!======================== start of intsch=2 section ====================
!       this is intsc           NS interps done first
!       first extend s arrays into sx - this one -1:il+2 & -1:il+2
    
  do k = 1,ol
    wx(1:ipan,1:jpan,1:npan,k) = &
      reshape( wtr(1:ipan*jpan*npan,k), (/ ipan, jpan, npan /) )
    do n = 1,npan
      do j = 1,jpan
        iq = 1+(j-1)*ipan+(n-1)*ipan*jpan
        wx(0,j,n,k)      = wtr( iw(iq),k)
        wx(-1,j,n,k)     = wtr(iww(iq),k)
        iq = j*ipan+(n-1)*ipan*jpan
        wx(ipan+1,j,n,k) = wtr( ie(iq),k)
        wx(ipan+2,j,n,k) = wtr(iee(iq),k)
      end do            ! j loop
      do i = 1,ipan
        iq = i+(n-1)*ipan*jpan
        wx(i,0,n,k)      = wtr( is(iq),k)
        wx(i,-1,n,k)     = wtr(iss(iq),k)
        iq = i-ipan+n*ipan*jpan
        wx(i,jpan+1,n,k) = wtr( in(iq),k)
        wx(i,jpan+2,n,k) = wtr(inn(iq),k)
      end do            ! i loop
    end do
!   for ns interpolation, sometimes need (different from ew):
!        (-1,0),   (0,0),   (0,-1)   (-1,il+1),   (0,il+1),   (0,il+2)
!      (il+1,0),(il+2,0),(il+1,-1) (il+1,il+1),(il+2,il+1),(il+1,il+2)
    do n = 1,npan
      wx(-1,0,n,k)          = wtr(lsww(n),k)
      wx(0,0,n,k)           = wtr(isw(1+(n-1)*ipan*jpan),   k)
      wx(0,-1,n,k)          = wtr(lssw(n),k)
      wx(ipan+2,0,n,k)      = wtr(lsee(n),k)
      wx(ipan+1,-1,n,k)     = wtr(lsse(n),k)
      wx(-1,jpan+1,n,k)     = wtr(lnww(n),k)
      wx(0,jpan+1,n,k)      = wtr(inw(1-ipan+n*ipan*jpan),  k)
      wx(0,jpan+2,n,k)      = wtr(lnnw(n),k)
      wx(ipan+2,jpan+1,n,k) = wtr(lnee(n),k)
      wx(ipan+1,jpan+2,n,k) = wtr(lnne(n),k)
      wx(ipan+1,0,n,k)      = wtr(ise(ipan+(n-1)*ipan*jpan),k)
      wx(ipan+1,jpan+1,n,k) = wtr(ine(n*ipan*jpan),         k)
    end do           ! n loop
  end do             ! k loop 
    
end if                     ! (intsch==1) .. else ..
!========================   end of intsch=1 section ====================

return
end subroutine setwx

subroutine setsx(sx,s,intsch)

use cc_mpi             ! CC MPI routines
use indices_m          ! Grid index arrays
use newmpar_m          ! Grid parameters

implicit none

integer, intent(in) :: intsch
integer k, n, i, j, iq
real, dimension(:,:), intent(in) :: s
real, dimension(-1:ipan+2,-1:jpan+2,1:npan,1:ol), intent(out) :: sx

! unpack s to sx for interpolation
sx = 0.

!======================== start of intsch=1 section ====================
if ( intsch==1 ) then

  do k = 1,ol
    sx(1:ipan,1:jpan,1:npan,k) = &
      reshape( s(1:ipan*jpan*npan,k), (/ ipan, jpan, npan /) )
    do n = 1,npan
      do j = 1,jpan
        iq = 1+(j-1)*ipan+(n-1)*ipan*jpan
        sx(0,j,n,k)      = s( iw(iq),k)
        sx(-1,j,n,k)     = s(iww(iq),k)
        iq = j*ipan+(n-1)*ipan*jpan
        sx(ipan+1,j,n,k) = s( ie(iq),k)
        sx(ipan+2,j,n,k) = s(iee(iq),k)
      end do            ! j loop
      do i = 1,ipan
        iq = i+(n-1)*ipan*jpan
        sx(i,0,n,k)      = s( is(iq),k)
        sx(i,-1,n,k)     = s(iss(iq),k)
        iq = i-ipan+n*ipan*jpan
        sx(i,jpan+1,n,k) = s( in(iq),k)
        sx(i,jpan+2,n,k) = s(inn(iq),k)
      end do            ! i loop
    end do
!   for ew interpolation, sometimes need (different from ns):
!       (-1,0),   (0,0),   (0,-1)   (-1,il+1),   (0,il+1),   (0,il+2)
!     (il+1,0),(il+2,0),(il+1,-1) (il+1,il+1),(il+2,il+1),(il+1,il+2)
    do n = 1,npan
      sx(-1,0,n,k)          = s(lwws(n),                  k)
      sx(0,0,n,k)           = s(iws(1+(n-1)*ipan*jpan),   k)
      sx(0,-1,n,k)          = s(lwss(n),                  k)
      sx(ipan+1,0,n,k)      = s(ies(ipan+(n-1)*ipan*jpan),k)
      sx(ipan+2,0,n,k)      = s(lees(n),                  k)
      sx(ipan+1,-1,n,k)     = s(less(n),                  k)
      sx(-1,jpan+1,n,k)     = s(lwwn(n),                  k)
      sx(0,jpan+2,n,k)      = s(lwnn(n),                  k)
      sx(ipan+2,jpan+1,n,k) = s(leen(n),                  k)
      sx(ipan+1,jpan+2,n,k) = s(lenn(n),                  k)
      sx(0,jpan+1,n,k)      = s(iwn(1-ipan+n*ipan*jpan),  k)
      sx(ipan+1,jpan+1,n,k) = s(ien(n*ipan*jpan),         k)
    end do           ! n loop
  end do             ! k loop
    
!========================   end of intsch=1 section ====================
else     ! if(intsch==1)then
!======================== start of intsch=2 section ====================
!       this is intsc           NS interps done first
!       first extend s arrays into sx - this one -1:il+2 & -1:il+2
    
  do k = 1,ol
    sx(1:ipan,1:jpan,1:npan,k) = &
      reshape( s(1:ipan*jpan*npan,k), (/ ipan, jpan, npan /) )
    do n = 1,npan
      do j = 1,jpan
        iq = 1+(j-1)*ipan+(n-1)*ipan*jpan
        sx(0,j,n,k)      = s( iw(iq),k)
        sx(-1,j,n,k)     = s(iww(iq),k)
        iq = j*ipan+(n-1)*ipan*jpan
        sx(ipan+1,j,n,k) = s( ie(iq),k)
        sx(ipan+2,j,n,k) = s(iee(iq),k)
      end do            ! j loop
      do i = 1,ipan
        iq = i+(n-1)*ipan*jpan
        sx(i,0,n,k)      = s( is(iq),k)
        sx(i,-1,n,k)     = s(iss(iq),k)
        iq = i-ipan+n*ipan*jpan
        sx(i,jpan+1,n,k) = s( in(iq),k)
        sx(i,jpan+2,n,k) = s(inn(iq),k)
      end do            ! i loop
    end do
!   for ns interpolation, sometimes need (different from ew):
!        (-1,0),   (0,0),   (0,-1)   (-1,il+1),   (0,il+1),   (0,il+2)
!      (il+1,0),(il+2,0),(il+1,-1) (il+1,il+1),(il+2,il+1),(il+1,il+2)
    do n = 1,npan
      sx(-1,0,n,k)          = s(lsww(n),k)
      sx(0,0,n,k)           = s(isw(1+(n-1)*ipan*jpan),   k)
      sx(0,-1,n,k)          = s(lssw(n),k)
      sx(ipan+2,0,n,k)      = s(lsee(n),k)
      sx(ipan+1,-1,n,k)     = s(lsse(n),k)
      sx(-1,jpan+1,n,k)     = s(lnww(n),k)
      sx(0,jpan+1,n,k)      = s(inw(1-ipan+n*ipan*jpan),  k)
      sx(0,jpan+2,n,k)      = s(lnnw(n),k)
      sx(ipan+2,jpan+1,n,k) = s(lnee(n),k)
      sx(ipan+1,jpan+2,n,k) = s(lnne(n),k)
      sx(ipan+1,0,n,k)      = s(ise(ipan+(n-1)*ipan*jpan),k)
      sx(ipan+1,jpan+1,n,k) = s(ine(n*ipan*jpan),         k)
    end do           ! n loop
  end do             ! k loop
    
end if                     ! (intsch==1) .. else ..
!========================   end of intsch=1 section ====================
    
return
end subroutine setsx

subroutine calcdepproc(sx,wx,nface,xg,yg,nn,bs_test,intsch)

use cc_mpi             ! CC MPI routines
use newmpar_m          ! Grid parameters
use parm_m             ! Model configuration

implicit none

integer, intent(in) :: nn, intsch
integer ii, iq, n, k, idel, jdel
integer, dimension(ifull,kl), intent(in) :: nface         ! interpolation coordinates
real, dimension(-1:ipan+2,-1:jpan+2,1:npan,1:ol), intent(in) :: sx
real, dimension(ifull,kl), intent(in) :: xg, yg           ! interpolation coordinates
real xxg, yyg, cmin, cmax, sx_ans
real dmul_2, dmul_3, cmul_1, cmul_2, cmul_3, cmul_4
real emul_1, emul_2, emul_3, emul_4, rmul_1, rmul_2, rmul_3, rmul_4
real sx_0m, sx_1m, sx_m0, sx_00, sx_10, sx_20, sx_m1, sx_01, sx_11, sx_21, sx_02, sx_12
logical, intent(in) :: bs_test
logical, dimension(-1:ipan+2,-1:jpan+2,1:npan,1:ol), intent(in) :: wx
logical bcub_water, blin_test


!======================== start of intsch=1 section ====================
if ( intsch==1 ) then
    
  ! Loop over points that need to be calculated for other processes
  if ( bs_test ) then

    do ii = 1,neighnum
      do iq = 1,drlen(ii)
        n = nint(dpoints(ii)%a(iq,1)) + noff ! Local index
        idel = int(dpoints(ii)%a(iq,2))
        xxg = dpoints(ii)%a(iq,2) - real(idel)
        jdel = int(dpoints(ii)%a(iq,3))
        yyg = dpoints(ii)%a(iq,3) - real(jdel)
        k = nint(dpoints(ii)%a(iq,4))
        idel = idel - ioff
        jdel = jdel - joff
            
        sx_0m = sx(idel,  jdel-1,n,k)
        sx_1m = sx(idel+1,jdel-1,n,k)
        sx_m0 = sx(idel-1,jdel  ,n,k)
        sx_00 = sx(idel,  jdel  ,n,k)
        sx_10 = sx(idel+1,jdel  ,n,k)
        sx_20 = sx(idel+2,jdel  ,n,k)
        sx_m1 = sx(idel-1,jdel+1,n,k)
        sx_01 = sx(idel,  jdel+1,n,k)
        sx_11 = sx(idel+1,jdel+1,n,k)
        sx_21 = sx(idel+2,jdel+1,n,k)
        sx_02 = sx(idel,  jdel+2,n,k)
        sx_12 = sx(idel+1,jdel+2,n,k)

        bcub_water = wx(idel,jdel-1,n,k) .and. wx(idel+1,jdel-1,n,k) .and.   &
                     wx(idel-1,jdel,n,k) .and. wx(idel,jdel,n,k) .and.       &
                     wx(idel+1,jdel,n,k) .and. wx(idel+2,jdel,n,k) .and.     &
                     wx(idel-1,jdel+1,n,k) .and. wx(idel,jdel+1,n,k) .and.   &
                     wx(idel+1,jdel+1,n,k) .and. wx(idel+2,jdel+1,n,k) .and. &
                     wx(idel,jdel+2,n,k) .and. wx(idel+1,jdel+2,n,k)

        blin_test = sx_00>=cxx .and. sx_10>=cxx .and. sx_01>=cxx .and. sx_11>=cxx

        if ( bcub_water ) then
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
          rmul_1 = sx_0m*dmul_2 + sx_1m*dmul_3
          rmul_2 = sx_m0*cmul_1 + sx_00*cmul_2 + &
                   sx_10*cmul_3 + sx_20*cmul_4
          rmul_3 = sx_m1*cmul_1 + sx_01*cmul_2 + &
                   sx_11*cmul_3 + sx_21*cmul_4
          rmul_4 = sx_02*dmul_2 + sx_12*dmul_3
          sx_ans = rmul_1*emul_1 + rmul_2*emul_2 + rmul_3*emul_3 + rmul_4*emul_4
          cmin = min(sx_00,sx_10,sx_01,sx_11)
          cmax = max(sx_00,sx_10,sx_01,sx_11)
          sextra(ii)%a(iq+(nn-1)*drlen(ii)) = min( max( cmin, sx_ans ), cmax ) ! Bermejo & Staniforth              
        else if ( blin_test ) then
          sextra(ii)%a(iq+(nn-1)*drlen(ii)) = (1.-xxg)*(1.-yyg)*sx_00 + xxg*(1.-yyg)*sx_10 &
                                            + (1.-xxg)*yyg*sx_01 + xxg*yyg*sx_11
        else
          sextra(ii)%a(iq+(nn-1)*drlen(ii)) = cxx - 1.
        end if
            
      end do          ! iq loop
    end do            ! ii loop
    
  else

    do ii = 1,neighnum
      do iq = 1,drlen(ii)
        n = nint(dpoints(ii)%a(iq,1)) + noff ! Local index
        idel = int(dpoints(ii)%a(iq,2))
        xxg = dpoints(ii)%a(iq,2) - real(idel)
        jdel = int(dpoints(ii)%a(iq,3))
        yyg = dpoints(ii)%a(iq,3) - real(jdel)
        k = nint(dpoints(ii)%a(iq,4))
        idel = idel - ioff
        jdel = jdel - joff

        sx_0m = sx(idel,  jdel-1,n,k)
        sx_1m = sx(idel+1,jdel-1,n,k)
        sx_m0 = sx(idel-1,jdel  ,n,k)
        sx_00 = sx(idel,  jdel  ,n,k)
        sx_10 = sx(idel+1,jdel  ,n,k)
        sx_20 = sx(idel+2,jdel  ,n,k)
        sx_m1 = sx(idel-1,jdel+1,n,k)
        sx_01 = sx(idel,  jdel+1,n,k)
        sx_11 = sx(idel+1,jdel+1,n,k)
        sx_21 = sx(idel+2,jdel+1,n,k)
        sx_02 = sx(idel,  jdel+2,n,k)
        sx_12 = sx(idel+1,jdel+2,n,k)

        bcub_water = wx(idel,jdel-1,n,k) .and. wx(idel+1,jdel-1,n,k) .and.   &
                     wx(idel-1,jdel,n,k) .and. wx(idel,jdel,n,k) .and.       &
                     wx(idel+1,jdel,n,k) .and. wx(idel+2,jdel,n,k) .and.     &
                     wx(idel-1,jdel+1,n,k) .and. wx(idel,jdel+1,n,k) .and.   &
                     wx(idel+1,jdel+1,n,k) .and. wx(idel+2,jdel+1,n,k) .and. &
                     wx(idel,jdel+2,n,k) .and. wx(idel+1,jdel+2,n,k)

        blin_test = sx_00>=cxx .and. sx_10>=cxx .and. sx_01>=cxx .and. sx_11>=cxx

        if ( bcub_water ) then
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
          rmul_1 = sx_0m*dmul_2 + sx_1m*dmul_3
          rmul_2 = sx_m0*cmul_1 + sx_00*cmul_2 + &
                   sx_10*cmul_3 + sx_20*cmul_4
          rmul_3 = sx_m1*cmul_1 + sx_01*cmul_2 + &
                   sx_11*cmul_3 + sx_21*cmul_4
          rmul_4 = sx_02*dmul_2 + sx_12*dmul_3
          sextra(ii)%a(iq+(nn-1)*drlen(ii)) = rmul_1*emul_1 + rmul_2*emul_2 + rmul_3*emul_3 + rmul_4*emul_4
        else if ( blin_test ) then
          sextra(ii)%a(iq+(nn-1)*drlen(ii)) = (1.-xxg)*(1.-yyg)*sx_00 + xxg*(1.-yyg)*sx_10 &
                                            + (1.-xxg)*yyg*sx_01 + xxg*yyg*sx_11
        else
          sextra(ii)%a(iq+(nn-1)*drlen(ii)) = cxx - 1.
        end if
        
      end do          ! iq loop
    end do            ! ii loop
      
  end if           ! bs_test ..else..  
        
!========================   end of intsch=1 section ====================
else     ! if(intsch==1)then
!======================== start of intsch=2 section ====================
!       this is intsc           NS interps done first
!       first extend s arrays into sx - this one -1:il+2 & -1:il+2
    
  if ( bs_test ) then
      
    do ii = neighnum,1,-1
      do iq = 1,drlen(ii)
        n = nint(dpoints(ii)%a(iq,1)) + noff ! Local index
        idel = int(dpoints(ii)%a(iq,2))
        xxg = dpoints(ii)%a(iq,2) - real(idel)
        jdel = int(dpoints(ii)%a(iq,3))
        yyg = dpoints(ii)%a(iq,3) - real(jdel)
        k = nint(dpoints(ii)%a(iq,4))
        idel = idel - ioff
        jdel = jdel - joff
            
        sx_0m = sx(idel,  jdel-1,n,k)
        sx_1m = sx(idel+1,jdel-1,n,k)
        sx_m0 = sx(idel-1,jdel  ,n,k)
        sx_00 = sx(idel,  jdel  ,n,k)
        sx_10 = sx(idel+1,jdel  ,n,k)
        sx_20 = sx(idel+2,jdel  ,n,k)
        sx_m1 = sx(idel-1,jdel+1,n,k)
        sx_01 = sx(idel,  jdel+1,n,k)
        sx_11 = sx(idel+1,jdel+1,n,k)
        sx_21 = sx(idel+2,jdel+1,n,k)
        sx_02 = sx(idel,  jdel+2,n,k)
        sx_12 = sx(idel+1,jdel+2,n,k)

        bcub_water = wx(idel,jdel-1,n,k) .and. wx(idel+1,jdel-1,n,k) .and.   &
                     wx(idel-1,jdel,n,k) .and. wx(idel,jdel,n,k) .and.       &
                     wx(idel+1,jdel,n,k) .and. wx(idel+2,jdel,n,k) .and.     &
                     wx(idel-1,jdel+1,n,k) .and. wx(idel,jdel+1,n,k) .and.   &
                     wx(idel+1,jdel+1,n,k) .and. wx(idel+2,jdel+1,n,k) .and. &
                     wx(idel,jdel+2,n,k) .and. wx(idel+1,jdel+2,n,k)

        blin_test = sx_00>=cxx .and. sx_10>=cxx .and. sx_01>=cxx .and. sx_11>=cxx

        if ( bcub_water ) then
          cmul_1 = (1.-yyg)*(2.-yyg)*(-yyg)/6.
          cmul_2 = (1.-yyg)*(2.-yyg)*(1.+yyg)/2.
          cmul_3 = yyg*(1.+yyg)*(2.-yyg)/2.
          cmul_4 = (1.-yyg)*(-yyg)*(1.+yyg)/6.
          dmul_2 = (1.-yyg)
          dmul_3 = yyg
          emul_1 = (1.-xxg)*(2.-xxg)*(-xxg)/6.
          emul_2 = (1.-xxg)*(2.-xxg)*(1.+xxg)/2.
          emul_3 = xxg*(1.+xxg)*(2.-xxg)/2.
          emul_4 = (1.-xxg)*(-xxg)*(1.+xxg)/6.
          rmul_1 = sx_m0*dmul_2 + sx_m1*dmul_3
          rmul_2 = sx_0m*cmul_1 + sx_00*cmul_2 + &
                   sx_01*cmul_3 + sx_02*cmul_4
          rmul_3 = sx_1m*cmul_1 + sx_10*cmul_2 + &
                   sx_11*cmul_3 + sx_12*cmul_4
          rmul_4 = sx_20*dmul_2 + sx_21*dmul_3
          sx_ans = rmul_1*emul_1 + rmul_2*emul_2 + rmul_3*emul_3 + rmul_4*emul_4
          cmin = min(sx_00,sx_10,sx_01,sx_11)
          cmax = max(sx_00,sx_10,sx_01,sx_11)
          sextra(ii)%a(iq+(nn-1)*drlen(ii)) = min( max( cmin, sx_ans ), cmax ) ! Bermejo & Staniforth              
        else if ( blin_test ) then
          sextra(ii)%a(iq+(nn-1)*drlen(ii)) = (1.-xxg)*(1.-yyg)*sx_00 + xxg*(1.-yyg)*sx_10 &
                                            + (1.-xxg)*yyg*sx_01 + xxg*yyg*sx_11
        else
          sextra(ii)%a(iq+(nn-1)*drlen(ii)) = cxx - 1.
        end if
            
      end do          ! iq loop
    end do            ! ii loop
      
  else
      
    do ii = neighnum,1,-1
      do iq = 1,drlen(ii)
        n = nint(dpoints(ii)%a(iq,1)) + noff ! Local index
        idel = int(dpoints(ii)%a(iq,2))
        xxg = dpoints(ii)%a(iq,2) - real(idel)
        jdel = int(dpoints(ii)%a(iq,3))
        yyg = dpoints(ii)%a(iq,3) - real(jdel)
        k = nint(dpoints(ii)%a(iq,4))
        idel = idel - ioff
        jdel = jdel - joff
            
        sx_0m = sx(idel,  jdel-1,n,k)
        sx_1m = sx(idel+1,jdel-1,n,k)
        sx_m0 = sx(idel-1,jdel  ,n,k)
        sx_00 = sx(idel,  jdel  ,n,k)
        sx_10 = sx(idel+1,jdel  ,n,k)
        sx_20 = sx(idel+2,jdel  ,n,k)
        sx_m1 = sx(idel-1,jdel+1,n,k)
        sx_01 = sx(idel,  jdel+1,n,k)
        sx_11 = sx(idel+1,jdel+1,n,k)
        sx_21 = sx(idel+2,jdel+1,n,k)
        sx_02 = sx(idel,  jdel+2,n,k)
        sx_12 = sx(idel+1,jdel+2,n,k)

        bcub_water = wx(idel,jdel-1,n,k) .and. wx(idel+1,jdel-1,n,k) .and.   &
                     wx(idel-1,jdel,n,k) .and. wx(idel,jdel,n,k) .and.       &
                     wx(idel+1,jdel,n,k) .and. wx(idel+2,jdel,n,k) .and.     &
                     wx(idel-1,jdel+1,n,k) .and. wx(idel,jdel+1,n,k) .and.   &
                     wx(idel+1,jdel+1,n,k) .and. wx(idel+2,jdel+1,n,k) .and. &
                     wx(idel,jdel+2,n,k) .and. wx(idel+1,jdel+2,n,k)

        blin_test = sx_00>=cxx .and. sx_10>=cxx .and. sx_01>=cxx .and. sx_11>=cxx

        if ( bcub_water ) then
          cmul_1 = (1.-yyg)*(2.-yyg)*(-yyg)/6.
          cmul_2 = (1.-yyg)*(2.-yyg)*(1.+yyg)/2.
          cmul_3 = yyg*(1.+yyg)*(2.-yyg)/2.
          cmul_4 = (1.-yyg)*(-yyg)*(1.+yyg)/6.
          dmul_2 = (1.-yyg)
          dmul_3 = yyg
          emul_1 = (1.-xxg)*(2.-xxg)*(-xxg)/6.
          emul_2 = (1.-xxg)*(2.-xxg)*(1.+xxg)/2.
          emul_3 = xxg*(1.+xxg)*(2.-xxg)/2.
          emul_4 = (1.-xxg)*(-xxg)*(1.+xxg)/6.
          rmul_1 = sx_m0*dmul_2 + sx_m1*dmul_3
          rmul_2 = sx_0m*cmul_1 + sx_00*cmul_2 + &
                   sx_01*cmul_3 + sx_02*cmul_4
          rmul_3 = sx_1m*cmul_1 + sx_10*cmul_2 + &
                   sx_11*cmul_3 + sx_12*cmul_4
          rmul_4 = sx_20*dmul_2 + sx_21*dmul_3
          sextra(ii)%a(iq+(nn-1)*drlen(ii)) = rmul_1*emul_1 + rmul_2*emul_2 + rmul_3*emul_3 + rmul_4*emul_4
        else if ( blin_test ) then
          sextra(ii)%a(iq+(nn-1)*drlen(ii)) = (1.-xxg)*(1.-yyg)*sx_00 + xxg*(1.-yyg)*sx_10 &
                                            + (1.-xxg)*yyg*sx_01 + xxg*yyg*sx_11
        else
          sextra(ii)%a(iq+(nn-1)*drlen(ii)) = cxx - 1.
        end if

      end do          ! iq loop
    end do            ! ii loop
      
  end if           ! bs_test ..else..  
    
end if                     ! (intsch==1) .. else ..
!========================   end of intsch=1 section ====================
        
return
end subroutine calcdepproc

subroutine calcdephost(sx,wx,s,nface,xg,yg,nn,bs_test,intsch,async_counter)

use cc_mpi             ! CC MPI routines
use newmpar_m          ! Grid parameters
use parm_m             ! Model configuration

implicit none

integer, dimension(ifull,kl), intent(in) :: nface         ! interpolation coordinates
integer, intent(in) :: nn, intsch, async_counter
integer iq, n, idel, jdel, k
real, dimension(:,:), intent(inout) :: s
real, dimension(-1:ipan+2,-1:jpan+2,1:npan,1:ol), intent(in) :: sx
real, dimension(ifull,kl), intent(in) :: xg, yg           ! interpolation coordinates
real xxg, yyg, cmin, cmax, sx_ans
real dmul_2, dmul_3, cmul_1, cmul_2, cmul_3, cmul_4
real emul_1, emul_2, emul_3, emul_4, rmul_1, rmul_2, rmul_3, rmul_4
real sx_0m, sx_1m, sx_m0, sx_00, sx_10, sx_20, sx_m1, sx_01, sx_11, sx_21, sx_02, sx_12
logical, intent(in) :: bs_test
logical, dimension(-1:ipan+2,-1:jpan+2,1:npan,1:ol), intent(in) :: wx
logical bcub_water, blin_test

!======================== start of intsch=1 section ====================
if ( intsch==1 ) then

  if ( bs_test ) then

#ifdef _OPENACC
    !$acc parallel loop collapse(2) copyout(s) present(sx,xg,yg,nface,wx) async(async_counter)
#else
    !$omp do schedule(static) private(k,iq,idel,jdel,n,xxg,yyg)                            &
    !$omp private(sx_0m,sx_1m,sx_m0,sx_00,sx_10,sx_20,sx_m1,sx_01,sx_11,sx_21,sx_02,sx_12) &
    !$omp private(cmul_1,cmul_2,cmul_3,cmul_4,dmul_2,dmul_3,emul_1,emul_2,emul_3,emul_4)   &
    !$omp private(rmul_1,rmul_2,rmul_3,rmul_4,bcub_water,blin_test,cmin,cmax)
#endif
    do k = 1,ol      
      do iq = 1,ifull
        idel = int(xg(iq,k))
        xxg = xg(iq,k) - real(idel)
        jdel = int(yg(iq,k))
        yyg = yg(iq,k) - real(jdel)
        idel = min( max(idel - ioff, 0), ipan)
        jdel = min( max(jdel - joff, 0), jpan)
        n = min( max(nface(iq,k) + noff, 1), npan)

        sx_0m = sx(idel,  jdel-1,n,k)
        sx_1m = sx(idel+1,jdel-1,n,k)
        sx_m0 = sx(idel-1,jdel  ,n,k)
        sx_00 = sx(idel,  jdel  ,n,k)
        sx_10 = sx(idel+1,jdel  ,n,k)
        sx_20 = sx(idel+2,jdel  ,n,k)
        sx_m1 = sx(idel-1,jdel+1,n,k)
        sx_01 = sx(idel,  jdel+1,n,k)
        sx_11 = sx(idel+1,jdel+1,n,k)
        sx_21 = sx(idel+2,jdel+1,n,k)
        sx_02 = sx(idel,  jdel+2,n,k)
        sx_12 = sx(idel+1,jdel+2,n,k)

        bcub_water = wx(idel,jdel-1,n,k) .and. wx(idel+1,jdel-1,n,k) .and.   &
                     wx(idel-1,jdel,n,k) .and. wx(idel,jdel,n,k) .and.       &
                     wx(idel+1,jdel,n,k) .and. wx(idel+2,jdel,n,k) .and.     &
                     wx(idel-1,jdel+1,n,k) .and. wx(idel,jdel+1,n,k) .and.   &
                     wx(idel+1,jdel+1,n,k) .and. wx(idel+2,jdel+1,n,k) .and. &
                     wx(idel,jdel+2,n,k) .and. wx(idel+1,jdel+2,n,k)

        blin_test = sx_00>=cxx .and. sx_10>=cxx .and. sx_01>=cxx .and. sx_11>=cxx

        if ( bcub_water ) then
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
          rmul_1 = sx_0m*dmul_2 + sx_1m*dmul_3
          rmul_2 = sx_m0*cmul_1 + sx_00*cmul_2 + &
                   sx_10*cmul_3 + sx_20*cmul_4
          rmul_3 = sx_m1*cmul_1 + sx_01*cmul_2 + &
                   sx_11*cmul_3 + sx_21*cmul_4
          rmul_4 = sx_02*dmul_2 + sx_12*dmul_3
          sx_ans = rmul_1*emul_1 + rmul_2*emul_2 + rmul_3*emul_3 + rmul_4*emul_4
          cmin = min(sx_00,sx_10,sx_01,sx_11)
          cmax = max(sx_00,sx_10,sx_01,sx_11)
          s(iq,k) = min( max( cmin, sx_ans ), cmax ) ! Bermejo & Staniforth              
        else if ( blin_test ) then
          s(iq,k) = (1.-xxg)*(1.-yyg)*sx_00 + xxg*(1.-yyg)*sx_10 &
                  + (1.-xxg)*yyg*sx_01 + xxg*yyg*sx_11
        else
          s(iq,k) = cxx - 1.
        end if            
      end do       ! iq loop
    end do         ! k loop
#ifdef _OPENACC
    !$acc end parallel loop
#else
    !$omp end do nowait
#endif

  else
      
#ifdef _OPENACC
    !$acc parallel loop collapse(2) copyout(s) present(sx,xg,yg,nface,wx) async(async_counter)
#else
    !$omp do schedule(static) private(k,iq,idel,jdel,n,xxg,yyg)                            &
    !$omp private(sx_0m,sx_1m,sx_m0,sx_00,sx_10,sx_20,sx_m1,sx_01,sx_11,sx_21,sx_02,sx_12) &
    !$omp private(cmul_1,cmul_2,cmul_3,cmul_4,dmul_2,dmul_3,emul_1,emul_2,emul_3,emul_4)   &
    !$omp private(rmul_1,rmul_2,rmul_3,rmul_4,bcub_water,blin_test)
#endif
    do k = 1,ol      
      do iq = 1,ifull
        idel = int(xg(iq,k))
        xxg = xg(iq,k) - real(idel)
        jdel = int(yg(iq,k))
        yyg = yg(iq,k) - real(jdel)
        idel = min( max(idel - ioff, 0), ipan)
        jdel = min( max(jdel - joff, 0), jpan)
        n = min( max(nface(iq,k) + noff, 1), npan)
            
        sx_0m = sx(idel,  jdel-1,n,k)
        sx_1m = sx(idel+1,jdel-1,n,k)
        sx_m0 = sx(idel-1,jdel  ,n,k)
        sx_00 = sx(idel,  jdel  ,n,k)
        sx_10 = sx(idel+1,jdel  ,n,k)
        sx_20 = sx(idel+2,jdel  ,n,k)
        sx_m1 = sx(idel-1,jdel+1,n,k)
        sx_01 = sx(idel,  jdel+1,n,k)
        sx_11 = sx(idel+1,jdel+1,n,k)
        sx_21 = sx(idel+2,jdel+1,n,k)
        sx_02 = sx(idel,  jdel+2,n,k)
        sx_12 = sx(idel+1,jdel+2,n,k)

        bcub_water = wx(idel,jdel-1,n,k) .and. wx(idel+1,jdel-1,n,k) .and.   &
                     wx(idel-1,jdel,n,k) .and. wx(idel,jdel,n,k) .and.       &
                     wx(idel+1,jdel,n,k) .and. wx(idel+2,jdel,n,k) .and.     &
                     wx(idel-1,jdel+1,n,k) .and. wx(idel,jdel+1,n,k) .and.   &
                     wx(idel+1,jdel+1,n,k) .and. wx(idel+2,jdel+1,n,k) .and. &
                     wx(idel,jdel+2,n,k) .and. wx(idel+1,jdel+2,n,k)

        blin_test = sx_00>=cxx .and. sx_10>=cxx .and. sx_01>=cxx .and. sx_11>=cxx

        if ( bcub_water ) then
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
          rmul_1 = sx_0m*dmul_2 + sx_1m*dmul_3
          rmul_2 = sx_m0*cmul_1 + sx_00*cmul_2 + &
                   sx_10*cmul_3 + sx_20*cmul_4
          rmul_3 = sx_m1*cmul_1 + sx_01*cmul_2 + &
                   sx_11*cmul_3 + sx_21*cmul_4
          rmul_4 = sx_02*dmul_2 + sx_12*dmul_3
          s(iq,k) = rmul_1*emul_1 + rmul_2*emul_2 + rmul_3*emul_3 + rmul_4*emul_4
        else if ( blin_test ) then
          s(iq,k) = (1.-xxg)*(1.-yyg)*sx_00 + xxg*(1.-yyg)*sx_10 &
                  + (1.-xxg)*yyg*sx_01 + xxg*yyg*sx_11
        else
          s(iq,k) = cxx - 1.
        end if

      end do       ! iq loop
    end do         ! k loop
#ifdef _OPENACC
    !$acc end parallel loop
#else
    !$omp end do nowait
#endif
          
  end if           ! bs_test ..else..    

!========================   end of intsch=1 section ====================
else     ! if(intsch==1)then
!======================== start of intsch=2 section ====================
!       this is intsc           NS interps done first
!       first extend s arrays into sx - this one -1:il+2 & -1:il+2
    
  if ( bs_test ) then
      
#ifdef _OPENACC
    !$acc parallel loop collapse(2) copyout(s) present(sx,xg,yg,nface,wx) async(async_counter)
#else
    !$omp do schedule(static) private(k,iq,idel,jdel,n,xxg,yyg)                            &
    !$omp private(sx_0m,sx_1m,sx_m0,sx_00,sx_10,sx_20,sx_m1,sx_01,sx_11,sx_21,sx_02,sx_12) &
    !$omp private(cmul_1,cmul_2,cmul_3,cmul_4,dmul_2,dmul_3,emul_1,emul_2,emul_3,emul_4)   &
    !$omp private(rmul_1,rmul_2,rmul_3,rmul_4,bcub_water,blin_test)
#endif
    do k = 1,ol
      do iq = 1,ifull
        idel = int(xg(iq,k))
        xxg = xg(iq,k) - real(idel)
        jdel = int(yg(iq,k))
        yyg = yg(iq,k) - real(jdel)
        idel = min( max(idel - ioff, 0), ipan)
        jdel = min( max(jdel - joff, 0), jpan)
        n = min( max(nface(iq,k) + noff, 1), npan)
            
        sx_0m = sx(idel,  jdel-1,n,k)
        sx_1m = sx(idel+1,jdel-1,n,k)
        sx_m0 = sx(idel-1,jdel  ,n,k)
        sx_00 = sx(idel,  jdel  ,n,k)
        sx_10 = sx(idel+1,jdel  ,n,k)
        sx_20 = sx(idel+2,jdel  ,n,k)
        sx_m1 = sx(idel-1,jdel+1,n,k)
        sx_01 = sx(idel,  jdel+1,n,k)
        sx_11 = sx(idel+1,jdel+1,n,k)
        sx_21 = sx(idel+2,jdel+1,n,k)
        sx_02 = sx(idel,  jdel+2,n,k)
        sx_12 = sx(idel+1,jdel+2,n,k)

        bcub_water = wx(idel,jdel-1,n,k) .and. wx(idel+1,jdel-1,n,k) .and.   &
                     wx(idel-1,jdel,n,k) .and. wx(idel,jdel,n,k) .and.       &
                     wx(idel+1,jdel,n,k) .and. wx(idel+2,jdel,n,k) .and.     &
                     wx(idel-1,jdel+1,n,k) .and. wx(idel,jdel+1,n,k) .and.   &
                     wx(idel+1,jdel+1,n,k) .and. wx(idel+2,jdel+1,n,k) .and. &
                     wx(idel,jdel+2,n,k) .and. wx(idel+1,jdel+2,n,k)

        blin_test = sx_00>=cxx .and. sx_10>=cxx .and. sx_01>=cxx .and. sx_11>=cxx

        if ( bcub_water ) then
          cmul_1 = (1.-yyg)*(2.-yyg)*(-yyg)/6.
          cmul_2 = (1.-yyg)*(2.-yyg)*(1.+yyg)/2.
          cmul_3 = yyg*(1.+yyg)*(2.-yyg)/2.
          cmul_4 = (1.-yyg)*(-yyg)*(1.+yyg)/6.
          dmul_2 = (1.-yyg)
          dmul_3 = yyg
          emul_1 = (1.-xxg)*(2.-xxg)*(-xxg)/6.
          emul_2 = (1.-xxg)*(2.-xxg)*(1.+xxg)/2.
          emul_3 = xxg*(1.+xxg)*(2.-xxg)/2.
          emul_4 = (1.-xxg)*(-xxg)*(1.+xxg)/6.
          rmul_1 = sx_m0*dmul_2 + sx_m1*dmul_3
          rmul_2 = sx_0m*cmul_1 + sx_00*cmul_2 + &
                   sx_01*cmul_3 + sx_02*cmul_4
          rmul_3 = sx_1m*cmul_1 + sx_10*cmul_2 + &
                   sx_11*cmul_3 + sx_12*cmul_4
          rmul_4 = sx_20*dmul_2 + sx_21*dmul_3
          sx_ans = rmul_1*emul_1 + rmul_2*emul_2 + rmul_3*emul_3 + rmul_4*emul_4
          cmin = min(sx_00,sx_10,sx_01,sx_11)
          cmax = max(sx_00,sx_10,sx_01,sx_11)
          s(iq,k) = min( max( cmin, sx_ans ), cmax ) ! Bermejo & Staniforth              
        else if ( blin_test ) then
          s(iq,k) = (1.-xxg)*(1.-yyg)*sx_00 + xxg*(1.-yyg)*sx_10 &
                  + (1.-xxg)*yyg*sx_01 + xxg*yyg*sx_11
        else
          s(iq,k) = cxx - 1.
        end if
            
      end do
    end do
#ifdef _OPENACC
    !$acc end parallel loop
#else
    !$omp end do nowait        
#endif
      
  else
      
#ifdef _OPENACC
    !$acc parallel loop collapse(2) copyout(s) present(sx,xg,yg,nface,wx) async(async_counter)
#else
    !$omp do schedule(static) private(k,iq,idel,jdel,n,xxg,yyg)                            &
    !$omp private(sx_0m,sx_1m,sx_m0,sx_00,sx_10,sx_20,sx_m1,sx_01,sx_11,sx_21,sx_02,sx_12) &
    !$omp private(cmul_1,cmul_2,cmul_3,cmul_4,dmul_2,dmul_3,emul_1,emul_2,emul_3,emul_4)   &
    !$omp private(rmul_1,rmul_2,rmul_3,rmul_4,bcub_water,blin_test)
#endif
    do k = 1,ol
      do iq = 1,ifull
        idel = int(xg(iq,k))
        xxg = xg(iq,k) - real(idel)
        jdel = int(yg(iq,k))
        yyg = yg(iq,k) - real(jdel)
        idel = min( max(idel - ioff, 0), ipan)
        jdel = min( max(jdel - joff, 0), jpan)
        n = min( max(nface(iq,k) + noff, 1), npan)

        sx_0m = sx(idel,  jdel-1,n,k)
        sx_1m = sx(idel+1,jdel-1,n,k)
        sx_m0 = sx(idel-1,jdel  ,n,k)
        sx_00 = sx(idel,  jdel  ,n,k)
        sx_10 = sx(idel+1,jdel  ,n,k)
        sx_20 = sx(idel+2,jdel  ,n,k)
        sx_m1 = sx(idel-1,jdel+1,n,k)
        sx_01 = sx(idel,  jdel+1,n,k)
        sx_11 = sx(idel+1,jdel+1,n,k)
        sx_21 = sx(idel+2,jdel+1,n,k)
        sx_02 = sx(idel,  jdel+2,n,k)
        sx_12 = sx(idel+1,jdel+2,n,k)

        bcub_water = wx(idel,jdel-1,n,k) .and. wx(idel+1,jdel-1,n,k) .and.   &
                     wx(idel-1,jdel,n,k) .and. wx(idel,jdel,n,k) .and.       &
                     wx(idel+1,jdel,n,k) .and. wx(idel+2,jdel,n,k) .and.     &
                     wx(idel-1,jdel+1,n,k) .and. wx(idel,jdel+1,n,k) .and.   &
                     wx(idel+1,jdel+1,n,k) .and. wx(idel+2,jdel+1,n,k) .and. &
                     wx(idel,jdel+2,n,k) .and. wx(idel+1,jdel+2,n,k)

        blin_test = sx_00>=cxx .and. sx_10>=cxx .and. sx_01>=cxx .and. sx_11>=cxx

        if ( bcub_water ) then
          cmul_1 = (1.-yyg)*(2.-yyg)*(-yyg)/6.
          cmul_2 = (1.-yyg)*(2.-yyg)*(1.+yyg)/2.
          cmul_3 = yyg*(1.+yyg)*(2.-yyg)/2.
          cmul_4 = (1.-yyg)*(-yyg)*(1.+yyg)/6.
          dmul_2 = (1.-yyg)
          dmul_3 = yyg
          emul_1 = (1.-xxg)*(2.-xxg)*(-xxg)/6.
          emul_2 = (1.-xxg)*(2.-xxg)*(1.+xxg)/2.
          emul_3 = xxg*(1.+xxg)*(2.-xxg)/2.
          emul_4 = (1.-xxg)*(-xxg)*(1.+xxg)/6.
          rmul_1 = sx_m0*dmul_2 + sx_m1*dmul_3
          rmul_2 = sx_0m*cmul_1 + sx_00*cmul_2 + &
                   sx_01*cmul_3 + sx_02*cmul_4
          rmul_3 = sx_1m*cmul_1 + sx_10*cmul_2 + &
                   sx_11*cmul_3 + sx_12*cmul_4
          rmul_4 = sx_20*dmul_2 + sx_21*dmul_3
          s(iq,k) = rmul_1*emul_1 + rmul_2*emul_2 + rmul_3*emul_3 + rmul_4*emul_4
        else if ( blin_test ) then
          s(iq,k) = (1.-xxg)*(1.-yyg)*sx_00 + xxg*(1.-yyg)*sx_10 &
                  + (1.-xxg)*yyg*sx_01 + xxg*yyg*sx_11
        else
          s(iq,k) = cxx - 1.
        end if
           
      end do
    end do
#ifdef _OPENACC
    !$acc end parallel loop
#else
    !$omp end do nowait
#endif
      
  end if           ! bs_test ..else..
    
end if                     ! (intsch==1) .. else ..
!========================   end of intsch=1 section ====================

return
end subroutine calcdephost
  
end module mloints
