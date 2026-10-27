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

module ints_m

private
public ints, ints_bl

interface ints
  module procedure ints3x, ints4
end interface

contains
    
! These routines interpolate depature points for the atmosphere semi-Lagrangian
! dynamics.  It is difficult to vectorise and to load-balance the message passing
! in these routines, although the timings do scale with increasing numbers of
! processors.

! The ints routines support Bermejo and Staniforth option to prevent overshooting
! for discontinuous fields (e.g., clouds, aerosols and tracers).
    
!     this one includes Bermejo & Staniforth option
!     parameter (mh_bs): B&S on/off depending on value of nfield
!     can put sx into work array (only 2d)
!     later may wish to save idel etc between array calls
!     this one does linear interp in x on outer y sides
!     doing x-interpolation before y-interpolation
!     nfield: 1 (psl), 2 (u, v), 3 (T), 4 (gases)
    
subroutine ints4(s,intsch,nface,xg,yg,nfield)

use cc_acc             ! CC OpenACC routines
use cc_mpi             ! CC MPI routines
use newmpar_m          ! Grid parameters
use parm_m             ! Model configuration

implicit none

integer, intent(in) :: intsch  ! method to interpolate panel corners
integer, intent(in) :: nfield  ! use B&S if nfield>=mh_bs
integer nstart, nend, nlen, np, nn, ntr
integer async_counter
integer, dimension(ifull,kl), intent(in) :: nface         ! interpolation coordinates
real, dimension(ifull,kl), intent(in) :: xg, yg           ! interpolation coordinates
real, dimension(:,:,:), intent(inout) :: s  ! array of tracers
real, dimension(-1:ipan+2,-1:jpan+2,1:npan,kl,nagg) :: sx ! unpacked tracer array

call START_LOG(ints_begin)

ntr = size(s,3)

! call bounds before calling ints
!call bounds(s,nrows=2)


!$acc enter data create(sx)

do nstart = 1,ntr,nagg
  nend = min(nstart + nagg - 1, ntr )
  nlen = nend - nstart + 1
    
  do nn = 1,nlen
    np = nn - 1 + nstart
    async_counter = mod(nn-1, async_length)
    call setsx(sx(:,:,:,:,nn),s(:,:,np),intsch)
    !$acc update device(sx(:,:,:,:,nn)) async(async_counter)      
  end do              ! nn loop  
    
  ! Loop over points that need to be calculated for other processes
  do nn = 1,nlen
    call calcdepproc(sx(:,:,:,:,nn),nface,xg,yg,nn,nfield,intsch)
  end do  
    
  ! Send messages to other processors.  We then start the calculation for this processor while waiting for
  ! the messages to return, thereby overlapping computation with communication.
  call intssync_send(nlen)

  ! Loop over points that need to be calculated for this process
  !$omp parallel
  do nn = 1,nlen
    np = nn - 1 + nstart
    async_counter = mod(nn-1, async_length)        
    call calcdephost(sx(:,:,:,:,nn),s(:,:,np),nface,xg,yg,nn,nfield,intsch,async_counter)
  end do
  !$acc wait
  !$omp end parallel
    
  call intssync_recv(s(:,:,nstart:nend))  

end do ! nstart

!$acc exit data delete(sx)


call END_LOG(ints_end)

return
end subroutine ints4

subroutine ints3x(sa,sb,sc,sd,se,intsch,nface,xg,yg,nfield)

use cc_acc             ! CC OpenACC routines
use cc_mpi             ! CC MPI routines
use newmpar_m          ! Grid parameters
use parm_m             ! Model configuration

implicit none

integer, intent(in) :: intsch  ! method to interpolate panel corners
integer, intent(in) :: nfield  ! use B&S if nfield>=mh_bs
integer nn, nlen
integer async_counter
integer, dimension(ifull,kl), intent(in) :: nface         ! interpolation coordinates
real, dimension(ifull,kl), intent(in) :: xg, yg           ! interpolation coordinates
real, dimension(:,:), intent(inout) :: sa  ! array of tracers
real, dimension(:,:), intent(inout), optional :: sb, sc, sd, se
real, dimension(-1:ipan+2,-1:jpan+2,1:npan,kl,nagg) :: sx ! unpacked tracer array

call START_LOG(ints_begin)

! call bounds before calling ints
!call bounds(s,nrows=2)


!$acc enter data create(sx)


! pack data into sx
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
call calcdepproc(sx(:,:,:,:,nn),nface,xg,yg,nn,nfield,intsch)

if ( present(sb) ) then
  nn = nn + 1
  call calcdepproc(sx(:,:,:,:,nn),nface,xg,yg,nn,nfield,intsch)    
end if

if ( present(sc) ) then
  nn = nn + 1
  call calcdepproc(sx(:,:,:,:,nn),nface,xg,yg,nn,nfield,intsch)    
end if

if ( present(sd) ) then
  nn = nn + 1
  call calcdepproc(sx(:,:,:,:,nn),nface,xg,yg,nn,nfield,intsch)    
end if

! Send messages to other processors.  We then start the calculation for this processor while waiting for
! the messages to return, thereby overlapping computation with communication.
nlen = nn
call intssync_send(nlen)

!$omp parallel

! Loop over points that need to be calculated for this process
async_counter = 0
nn = 1
call calcdephost(sx(:,:,:,:,nn),sa,nface,xg,yg,nn,nfield,intsch,async_counter)

if ( present(sb) ) then
  async_counter = mod(async_counter + 1, async_length)
  nn = nn + 1
  call calcdephost(sx(:,:,:,:,nn),sb,nface,xg,yg,nn,nfield,intsch,async_counter)    
end if

if ( present(sc) ) then
  async_counter = mod(async_counter + 1, async_length)
  nn = nn + 1
  call calcdephost(sx(:,:,:,:,nn),sc,nface,xg,yg,nn,nfield,intsch,async_counter)    
end if

if ( present(sd) ) then
  async_counter = mod(async_counter + 1, async_length)
  nn = nn + 1
  call calcdephost(sx(:,:,:,:,nn),sd,nface,xg,yg,nn,nfield,intsch,async_counter)    
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


if ( present(se) ) then

  ! pack data into sx
  async_counter = 0
  nn = 1
  call setsx(sx(:,:,:,:,nn),se,intsch)
  !$acc update device(sx(:,:,:,:,nn)) async(async_counter)

  ! Loop over points that need to be calculated for other processes
  nn = 1
  call calcdepproc(sx(:,:,:,:,nn),nface,xg,yg,nn,nfield,intsch)

  ! Send messages to other processors.  We then start the calculation for this processor while waiting for
  ! the messages to return, thereby overlapping computation with communication.
  nlen = nn
  call intssync_send(nlen)

  !$omp parallel

  ! Loop over points that need to be calculated for this process
  async_counter = 0
  nn = 1
  call calcdephost(sx(:,:,:,:,nn),se,nface,xg,yg,nn,nfield,intsch,async_counter)

  !$acc wait
  !$omp end parallel

  call intssync_recv(ta=se)
    
end if ! present(se)


!$acc exit data delete(sx)


call END_LOG(ints_end)

return
end subroutine ints3x

subroutine ints_bl(s,intsch,nface,xg,yg)  ! not usually called

use atmosphere_arrays_m ! Collection of atmosphere data
use cc_mpi              ! CC MPI routines
use indices_m           ! Grid index arrays
use newmpar_m           ! Grid parameters
use parm_m              ! Model configuration

implicit none

integer idel, iq, jdel
integer i, j, k, n
integer ii
integer, intent(in) :: intsch
integer, dimension(ifull,kl), intent(in) :: nface
real xxg, yyg
real, dimension(ifull,kl), intent(in) :: xg,yg      ! now passed through call
real, dimension(0:ipan+1,0:jpan+1,1:npan,kl) :: sx
real, dimension(ifull+iextra,kl), intent(inout) :: s

!     this one does bi-linear interpolation only
!     first extend s arrays into sx - this one -1:il+2 & -1:il+2
!                    but for bi-linear only need 0:il+1 &  0:il+1

call START_LOG(ints_begin)

! now call bounds before calling ints_bl
!call bounds(s,corner=.true.)

!$acc enter data create(sx)

sx(1:ipan,1:jpan,1:npan,1:kl) = reshape(s(1:ipan*jpan*npan,1:kl), (/ipan,jpan,npan,kl/))
do k = 1,kl
  do n = 1,npan
    do j = 1,jpan
      sx(0,j,n,k)      = s(iw(1+(j-1)*ipan+(n-1)*ipan*jpan),k)
      sx(ipan+1,j,n,k) = s(ie(j*ipan+(n-1)*ipan*jpan),k)
    end do               ! j loop
    do i = 1,ipan
      sx(i,0,n,k)      = s(is(i+(n-1)*ipan*jpan),k)
      sx(i,jpan+1,n,k) = s(in(i-ipan+n*ipan*jpan),k)
    end do               ! i loop
  end do
  do n = 1,npan
    sx(0,0,n,k)           = s(iws(1+(n-1)*ipan*jpan),k)
    sx(ipan+1,0,n,k)      = s(ies(ipan+(n-1)*ipan*jpan),k)
    sx(0,jpan+1,n,k)      = s(iwn(1-ipan+n*ipan*jpan),k)
    sx(ipan+1,jpan+1,n,k) = s(ien(n*ipan*jpan),k)
  end do                 ! n loop
end do                   ! k loop
!$acc update device(sx) async(0)

! Loop over points that need to be calculated for other processes
do ii = 1,neighnum
  do iq = 1,drlen(ii)
    n = nint(dpoints(ii)%a(iq,1)) + noff ! Local index
    !  Need global face index in fproc call
    idel = int(dpoints(ii)%a(iq,2))
    xxg = dpoints(ii)%a(iq,2) - real(idel)
    jdel = int(dpoints(ii)%a(iq,3))
    yyg = dpoints(ii)%a(iq,3) - real(jdel)
    k = nint(dpoints(ii)%a(iq,4))
    idel = idel - ioff
    jdel = jdel - joff
    sextra(ii)%a(iq) =      yyg*(xxg*sx(idel+1,jdel+1,n,k)+(1.-xxg)*sx(idel,jdel+1,n,k)) &
                     + (1.-yyg)*(xxg*sx(idel+1,  jdel,n,k)+(1.-xxg)*sx(idel,  jdel,n,k))
  end do
end do

call intssync_send

#ifdef _OPENACC
!$acc wait(0)
!$acc parallel loop collapse(2) copyout(s) present(sx,xg,yg,nface)
#else
!$omp parallel do schedule(static) private(k,iq,idel,xxg,jdel,yyg,n)
#endif
do k = 1,kl
  do iq = 1,ifull
    ! Convert face index from 0:npanels to array indices
    idel = int(xg(iq,k))
    xxg = xg(iq,k) - idel
    jdel = int(yg(iq,k))
    yyg = yg(iq,k) - jdel
    idel = min( max(idel - ioff, 0), ipan)
    jdel = min( max(jdel - joff, 0), jpan)
    n = min( max(nface(iq,k) + noff, 1), npan)
    s(iq,k) =      yyg*(xxg*sx(idel+1,jdel+1,n,k)+(1.-xxg)*sx(idel,jdel+1,n,k)) &
            + (1.-yyg)*(xxg*sx(idel+1,  jdel,n,k)+(1.-xxg)*sx(idel,  jdel,n,k))
  end do                  ! iq loop
end do                    ! k
#ifdef _OPENACC
!$acc end parallel loop
#else
!$omp end parallel do
#endif

call intssync_recv(s)

!$acc exit data delete(sx)

call END_LOG(ints_end)

return
end subroutine ints_bl
    
subroutine setsx(sx,s,intsch)

use atmosphere_arrays_m ! Collection of atmosphere data
use cc_mpi              ! CC MPI routines
use indices_m           ! Grid index arrays
use newmpar_m           ! Grid parameters

implicit none

integer, intent(in) :: intsch
integer k, n, i, j, iq
real, dimension(:,:), intent(in) :: s
real, dimension(-1:ipan+2,-1:jpan+2,1:npan,1:kl), intent(out) :: sx

! unpack s to sx for interpolation
sx = 0.

!======================== start of intsch=1 section ====================
if ( intsch==1 ) then

  ! this is intsb           EW interps done first
  ! first extend s arrays into sx - this one -1:il+2 & -1:il+2
  do k = 1,kl
    sx(1:ipan,1:jpan,1:npan,k) = reshape(s(1:ipan*jpan*npan,k), (/ipan,jpan,npan/))
    do n = 1,npan
      do j = 1,jpan
        iq = 1 + (j-1)*ipan + (n-1)*ipan*jpan
        sx(0,j,n,k)      = s(iw(iq),k)
        sx(-1,j,n,k)     = s(iww(iq),k)
        iq = j*ipan + (n-1)*ipan*jpan
        sx(ipan+1,j,n,k) = s(ie(iq),k)
        sx(ipan+2,j,n,k) = s(iee(iq),k)
      end do            ! j loop
      do i = 1,ipan
        iq = i + (n-1)*ipan*jpan
        sx(i,0,n,k)      = s(is(iq),k)
        sx(i,-1,n,k)     = s(iss(iq),k)
        iq = i - ipan + n*ipan*jpan
        sx(i,jpan+1,n,k) = s(in(iq),k)
        sx(i,jpan+2,n,k) = s(inn(iq),k)
      end do            ! i loop
      sx(-1,0,n,k)          = s(lwws(n),k)
      sx(0,0,n,k)           = s(iws(1+(n-1)*ipan*jpan),k)
      sx(0,-1,n,k)          = s(lwss(n),k)
      sx(ipan+1,0,n,k)      = s(ies(ipan+(n-1)*ipan*jpan),k)
      sx(ipan+2,0,n,k)      = s(lees(n),k)
      sx(ipan+1,-1,n,k)     = s(less(n),k)
      sx(-1,jpan+1,n,k)     = s(lwwn(n),k)
      sx(0,jpan+2,n,k)      = s(lwnn(n),k)
      sx(ipan+2,jpan+1,n,k) = s(leen(n),k)
      sx(ipan+1,jpan+2,n,k) = s(lenn(n),k)
      sx(0,jpan+1,n,k)      = s(iwn(1-ipan+n*ipan*jpan),k)
      sx(ipan+1,jpan+1,n,k) = s(ien(n*ipan*jpan),k)
    end do          ! n loop
  end do            ! k loop
 
!========================   end of intsch=1 section ====================
else     ! if(intsch==1)then
!======================== start of intsch=2 section ====================

  !       this is intsc           NS interps done first
  !       first extend s arrays into sx - this one -1:il+2 & -1:il+2
  do k = 1,kl
    sx(1:ipan,1:jpan,1:npan,k) = reshape(s(1:ipan*jpan*npan,k), (/ipan,jpan,npan/))  
    do n = 1,npan
      do j = 1,jpan
        iq = 1+(j-1)*ipan+(n-1)*ipan*jpan
        sx(0,j,n,k)      = s(iw(iq),k)
        sx(-1,j,n,k)     = s(iww(iq),k)
        iq = j*ipan+(n-1)*ipan*jpan
        sx(ipan+1,j,n,k) = s(ie(iq),k)
        sx(ipan+2,j,n,k) = s(iee(iq),k)
      end do            ! j loop
      do i = 1,ipan
        iq = i+(n-1)*ipan*jpan
        sx(i,0,n,k)      = s(is(iq),k)
        sx(i,-1,n,k)     = s(iss(iq),k)
        iq = i-ipan+n*ipan*jpan
        sx(i,jpan+1,n,k) = s(in(iq),k)
        sx(i,jpan+2,n,k) = s(inn(iq),k)
      end do            ! i loop
    end do
    do n = 1,npan
      sx(-1,0,n,k)          = s(lsww(n),k)
      sx(0,0,n,k)           = s(isw(1+(n-1)*ipan*jpan),k)
      sx(0,-1,n,k)          = s(lssw(n),k)
      sx(ipan+2,0,n,k)      = s(lsee(n),k)
      sx(ipan+1,-1,n,k)     = s(lsse(n),k)
      sx(-1,jpan+1,n,k)     = s(lnww(n),k)
      sx(0,jpan+1,n,k)      = s(inw(1-ipan+n*ipan*jpan),k)
      sx(0,jpan+2,n,k)      = s(lnnw(n),k)
      sx(ipan+2,jpan+1,n,k) = s(lnee(n),k)
      sx(ipan+1,jpan+2,n,k) = s(lnne(n),k)
      sx(ipan+1,0,n,k)      = s(ise(ipan+(n-1)*ipan*jpan),k)
      sx(ipan+1,jpan+1,n,k) = s(ine(n*ipan*jpan),k)
    end do              ! n loop
  end do                ! k loop

end if               ! (intsch==1) .. else ..
!========================   end of intsch=1 section ====================
    
return
end subroutine setsx
    
subroutine calcdepproc(sx,nface,xg,yg,nn,nfield,intsch)

use cc_mpi             ! CC MPI routines
use newmpar_m          ! Grid parameters
use parm_m             ! Model configuration
use parmhor_m          ! Horizontal advection parameters

implicit none

integer, intent(in) :: nn, nfield, intsch
integer ii, iq, n, idel, jdel, k
integer, dimension(ifull,kl), intent(in) :: nface         ! interpolation coordinates
real, dimension(-1:ipan+2,-1:jpan+2,1:npan,1:kl), intent(in) :: sx
real, dimension(ifull,kl), intent(in) :: xg, yg           ! interpolation coordinates
real xxg, yyg, cmin, cmax
real dmul_2, dmul_3, cmul_1, cmul_2, cmul_3, cmul_4
real emul_1, emul_2, emul_3, emul_4, rmul_1, rmul_2, rmul_3, rmul_4
real sx_0m, sx_1m, sx_m0, sx_00, sx_10, sx_20, sx_m1, sx_01, sx_11, sx_21, sx_02, sx_12

!======================== start of intsch=1 section ====================
if ( intsch==1 ) then

  ! Loop over points that need to be calculated for other processes
  if ( nfield<mh_bs ) then

    do ii = 1,neighnum
      do iq = 1,drlen(ii)
        ! depature point coordinates
        n = nint(dpoints(ii)%a(iq,1)) + noff ! Local index
        idel = int(dpoints(ii)%a(iq,2))
        xxg = dpoints(ii)%a(iq,2) - real(idel)
        jdel = int(dpoints(ii)%a(iq,3))
        yyg = dpoints(ii)%a(iq,3) - real(jdel)
        k = nint(dpoints(ii)%a(iq,4))
        idel = idel - ioff
        jdel = jdel - joff
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
        sx_0m = sx(idel,  jdel-1,n,k)
        sx_1m = sx(idel+1,jdel-1,n,k)
        sx_m0 = sx(idel-1,jdel,  n,k)
        sx_00 = sx(idel,  jdel,  n,k)
        sx_10 = sx(idel+1,jdel,  n,k)
        sx_20 = sx(idel+2,jdel,  n,k)
        sx_m1 = sx(idel-1,jdel+1,n,k)
        sx_01 = sx(idel,  jdel+1,n,k)
        sx_11 = sx(idel+1,jdel+1,n,k)
        sx_21 = sx(idel+2,jdel+1,n,k)
        sx_02 = sx(idel,  jdel+2,n,k)
        sx_12 = sx(idel+1,jdel+2,n,k)
        rmul_1 = sx_0m*dmul_2 + sx_1m*dmul_3
        rmul_2 = sx_m0*cmul_1 + sx_00*cmul_2 + &
                 sx_10*cmul_3 + sx_20*cmul_4
        rmul_3 = sx_m1*cmul_1 + sx_01*cmul_2 + &
                 sx_11*cmul_3 + sx_21*cmul_4
        rmul_4 = sx_02*dmul_2 + sx_12*dmul_3
        sextra(ii)%a(iq+(nn-1)*drlen(ii)) = rmul_1*emul_1 + rmul_2*emul_2 + rmul_3*emul_3 + rmul_4*emul_4
      end do      ! iq loop
    end do        ! ii loop

  else              ! (nfield<mh_bs)
      
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
        sx_0m = sx(idel,  jdel-1,n,k)
        sx_1m = sx(idel+1,jdel-1,n,k)
        sx_m0 = sx(idel-1,jdel,  n,k)
        sx_00 = sx(idel,  jdel,  n,k)
        sx_10 = sx(idel+1,jdel,  n,k)
        sx_20 = sx(idel+2,jdel,  n,k)
        sx_m1 = sx(idel-1,jdel+1,n,k)
        sx_01 = sx(idel,  jdel+1,n,k)
        sx_11 = sx(idel+1,jdel+1,n,k)
        sx_21 = sx(idel+2,jdel+1,n,k)
        sx_02 = sx(idel,  jdel+2,n,k)
        sx_12 = sx(idel+1,jdel+2,n,k)
        cmin = min(sx_00,sx_01,sx_10,sx_11)
        cmax = max(sx_00,sx_01,sx_10,sx_11)
        rmul_1 = sx_0m*dmul_2 + sx_1m*dmul_3
        rmul_2 = sx_m0*cmul_1 + sx_00*cmul_2 + &
                 sx_10*cmul_3 + sx_20*cmul_4
        rmul_3 = sx_m1*cmul_1 + sx_01*cmul_2 + &
                 sx_11*cmul_3 + sx_21*cmul_4
        rmul_4 = sx_02*dmul_2 + sx_12*dmul_3
        sextra(ii)%a(iq+(nn-1)*drlen(ii)) = min( max( cmin, &
            rmul_1*emul_1 + rmul_2*emul_2 + rmul_3*emul_3 + rmul_4*emul_4 ), cmax ) ! Bermejo & Staniforth
      end do      ! iq loop
    end do        ! ii loop
    
  end if

!========================   end of intsch=1 section ====================
else     ! if(intsch==1)then
!======================== start of intsch=2 section ====================

  if ( nfield<mh_bs ) then
      
    do ii = 1,neighnum
      do iq = 1,drlen(ii)
        n = nint(dpoints(ii)%a(iq,1)) + noff ! Local index
        !  Need global face index in fproc call
        idel = int(dpoints(ii)%a(iq,2))
        xxg = dpoints(ii)%a(iq,2) - real(idel)
        jdel = int(dpoints(ii)%a(iq,3))
        yyg = dpoints(ii)%a(iq,3) - real(jdel)
        k = nint(dpoints(ii)%a(iq,4))
        idel = idel - ioff
        jdel = jdel - joff
        ! bi-cubic
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
        sx_0m = sx(idel,  jdel-1,n,k)
        sx_1m = sx(idel+1,jdel-1,n,k)
        sx_m0 = sx(idel-1,jdel,  n,k)
        sx_00 = sx(idel,  jdel,  n,k)
        sx_10 = sx(idel+1,jdel,  n,k)
        sx_20 = sx(idel+2,jdel,  n,k)
        sx_m1 = sx(idel-1,jdel+1,n,k)
        sx_01 = sx(idel,  jdel+1,n,k)
        sx_11 = sx(idel+1,jdel+1,n,k)
        sx_21 = sx(idel+2,jdel+1,n,k)
        sx_02 = sx(idel,  jdel+2,n,k)
        sx_12 = sx(idel+1,jdel+2,n,k)
        rmul_1 = sx_m0*dmul_2 + sx_m1*dmul_3
        rmul_2 = sx_0m*cmul_1 + sx_00*cmul_2 + &
                 sx_01*cmul_3 + sx_02*cmul_4
        rmul_3 = sx_1m*cmul_1 + sx_10*cmul_2 + &
                 sx_11*cmul_3 + sx_12*cmul_4
        rmul_4 = sx_20*dmul_2 + sx_21*dmul_3
        sextra(ii)%a(iq+(nn-1)*drlen(ii)) = rmul_1*emul_1 + rmul_2*emul_2 + rmul_3*emul_3 + rmul_4*emul_4
      end do         ! iq loop
    end do           ! ii loop
      
  else
      
    do ii = 1,neighnum
      do iq = 1,drlen(ii)
        n = nint(dpoints(ii)%a(iq,1)) + noff ! Local index
        !  Need global face index in fproc call
        idel = int(dpoints(ii)%a(iq,2))
        xxg = dpoints(ii)%a(iq,2) - real(idel)
        jdel = int(dpoints(ii)%a(iq,3))
        yyg = dpoints(ii)%a(iq,3) - real(jdel)
        k = nint(dpoints(ii)%a(iq,4))
        idel = idel - ioff
        jdel = jdel - joff
        ! bi-cubic
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
        sx_0m = sx(idel,  jdel-1,n,k)
        sx_1m = sx(idel+1,jdel-1,n,k)
        sx_m0 = sx(idel-1,jdel,  n,k)
        sx_00 = sx(idel,  jdel,  n,k)
        sx_10 = sx(idel+1,jdel,  n,k)
        sx_20 = sx(idel+2,jdel,  n,k)
        sx_m1 = sx(idel-1,jdel+1,n,k)
        sx_01 = sx(idel,  jdel+1,n,k)
        sx_11 = sx(idel+1,jdel+1,n,k)
        sx_21 = sx(idel+2,jdel+1,n,k)
        sx_02 = sx(idel,  jdel+2,n,k)
        sx_12 = sx(idel+1,jdel+2,n,k)
        cmin = min(sx_00,sx_01,sx_10,sx_11)
        cmax = max(sx_00,sx_01,sx_10,sx_11)
        rmul_1 = sx_m0*dmul_2 + sx_m1*dmul_3
        rmul_2 = sx_0m*cmul_1 + sx_00*cmul_2 + &
                 sx_01*cmul_3 + sx_02*cmul_4
        rmul_3 = sx_1m*cmul_1 + sx_10*cmul_2 + &
                 sx_11*cmul_3 + sx_12*cmul_4
        rmul_4 = sx_20*dmul_2 + sx_21*dmul_3
        sextra(ii)%a(iq+(nn-1)*drlen(ii)) = min( max( cmin, &
            rmul_1*emul_1 + rmul_2*emul_2 + rmul_3*emul_3 + rmul_4*emul_4 ), cmax ) ! Bermejo & Staniforth
      end do      ! iq loop
    end do        ! ii loop
      
  end if
  
end if               ! (intsch==1) .. else ..
!========================   end of intsch=1 section ====================

return
end subroutine calcdepproc
  
subroutine calcdephost(sx,s,nface,xg,yg,nn,nfield,intsch,async_counter)

use atmosphere_arrays_m ! Collection of atmosphere data
use cc_mpi              ! CC MPI routines
use newmpar_m           ! Grid parameters
use parm_m              ! Model configuration

implicit none

integer, dimension(ifull,kl), intent(in) :: nface         ! interpolation coordinates
integer, intent(in) :: nn, nfield, intsch,async_counter
integer iq, n, idel, jdel, k
real, dimension(:,:), intent(inout) :: s
real, dimension(-1:ipan+2,-1:jpan+2,1:npan,1:kl), intent(in) :: sx
real, dimension(ifull,kl), intent(in) :: xg, yg           ! interpolation coordinates
real xxg, yyg, cmin, cmax
real dmul_2, dmul_3, cmul_1, cmul_2, cmul_3, cmul_4
real emul_1, emul_2, emul_3, emul_4, rmul_1, rmul_2, rmul_3, rmul_4
real sx_0m, sx_1m, sx_m0, sx_00, sx_10, sx_20, sx_m1, sx_01, sx_11, sx_21, sx_02, sx_12

!======================== start of intsch=1 section ====================
if ( intsch==1 ) then

  ! Loop over points that need to be calculated for this process
  if ( nfield<mh_bs ) then

#ifdef _OPENACC
    !$acc parallel loop collapse(2) copyout(s) present(sx,xg,yg,nface) async(async_counter)
#else
    !$omp do schedule(static) private(k,iq,idel,xxg,jdel,yyg)                        &
    !$omp private(n,cmul_1,cmul_2,cmul_3,cmul_4,dmul_2,dmul_3,emul_1,emul_2,emul_3)  &
    !$omp private(emul_4,rmul_1,rmul_2,rmul_3,rmul_4)                                &
    !$omp private(sx_0m,sx_1m,sx_m0,sx_00,sx_10,sx_20,sx_m1,sx_01)                   &
    !$omp private(sx_11,sx_21,sx_02,sx_12)
#endif
    do k = 1,kl
      do iq = 1,ifull    ! non Berm-Stan option
        idel = int(xg(iq,k))
        xxg = xg(iq,k) - real(idel)
        jdel = int(yg(iq,k))
        yyg = yg(iq,k) - real(jdel)
        idel = min( max(idel - ioff, 0), ipan)
        jdel = min( max(jdel - joff, 0), jpan)
        n = min( max(nface(iq,k) + noff, 1), npan)
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
        sx_0m = sx(idel,  jdel-1,n,k)
        sx_1m = sx(idel+1,jdel-1,n,k)
        sx_m0 = sx(idel-1,jdel,  n,k)
        sx_00 = sx(idel,  jdel,  n,k)
        sx_10 = sx(idel+1,jdel,  n,k)
        sx_20 = sx(idel+2,jdel,  n,k)
        sx_m1 = sx(idel-1,jdel+1,n,k)
        sx_01 = sx(idel,  jdel+1,n,k)
        sx_11 = sx(idel+1,jdel+1,n,k)
        sx_21 = sx(idel+2,jdel+1,n,k)
        sx_02 = sx(idel,  jdel+2,n,k)
        sx_12 = sx(idel+1,jdel+2,n,k)
        rmul_1 = sx_0m*dmul_2 + sx_1m*dmul_3
        rmul_2 = sx_m0*cmul_1 + sx_00*cmul_2 + &
                 sx_10*cmul_3 + sx_20*cmul_4
        rmul_3 = sx_m1*cmul_1 + sx_01*cmul_2 + &
                 sx_11*cmul_3 + sx_21*cmul_4
        rmul_4 = sx_02*dmul_2 + sx_12*dmul_3
        s(iq,k) = rmul_1*emul_1 + rmul_2*emul_2 + rmul_3*emul_3 + rmul_4*emul_4
      end do       ! iq loop
    end do         ! k loop
#ifdef _OPENACC
    !$acc end parallel loop
#else
    !$omp end do nowait
#endif

  else              ! (nfield<mh_bs)

#ifdef _OPENACC      
    !$acc parallel loop collapse(2) copyout(s) present(sx,xg,yg,nface) async(async_counter)
#else
    !$omp do schedule(static) private(k,iq,idel,xxg,jdel,yyg)                        &
    !$omp private(n,cmul_1,cmul_2,cmul_3,cmul_4,dmul_2,dmul_3,emul_1,emul_2,emul_3)  &
    !$omp private(emul_4,rmul_1,rmul_2,rmul_3,rmul_4,cmax,cmin)                      &
    !$omp private(sx_0m,sx_1m,sx_m0,sx_00,sx_10,sx_20,sx_m1,sx_01)                   &
    !$omp private(sx_11,sx_21,sx_02,sx_12)
#endif
    do k = 1,kl
      do iq = 1,ifull    ! Berm-Stan option here e.g. qg & gases
        idel = int(xg(iq,k))
        xxg = xg(iq,k) - real(idel)
        jdel = int(yg(iq,k))
        yyg = yg(iq,k) - real(jdel)
        idel = min( max(idel - ioff, 0), ipan)
        jdel = min( max(jdel - joff, 0), jpan)
        n = min( max(nface(iq,k) + noff, 1), npan)
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
        sx_0m = sx(idel,  jdel-1,n,k)
        sx_1m = sx(idel+1,jdel-1,n,k)
        sx_m0 = sx(idel-1,jdel,  n,k)
        sx_00 = sx(idel,  jdel,  n,k)
        sx_10 = sx(idel+1,jdel,  n,k)
        sx_20 = sx(idel+2,jdel,  n,k)
        sx_m1 = sx(idel-1,jdel+1,n,k)
        sx_01 = sx(idel,  jdel+1,n,k)
        sx_11 = sx(idel+1,jdel+1,n,k)
        sx_21 = sx(idel+2,jdel+1,n,k)
        sx_02 = sx(idel,  jdel+2,n,k)
        sx_12 = sx(idel+1,jdel+2,n,k)
        cmin = min(sx_00,sx_01,sx_10,sx_11)
        cmax = max(sx_00,sx_01,sx_10,sx_11)
        rmul_1 = sx_0m*dmul_2 + sx_1m*dmul_3
        rmul_2 = sx_m0*cmul_1 + sx_00*cmul_2 + &
                 sx_10*cmul_3 + sx_20*cmul_4
        rmul_3 = sx_m1*cmul_1 + sx_01*cmul_2 + &
                 sx_11*cmul_3 + sx_21*cmul_4
        rmul_4 = sx_02*dmul_2 + sx_12*dmul_3
        s(iq,k) = min( max( cmin, &
            rmul_1*emul_1 + rmul_2*emul_2 + rmul_3*emul_3 + rmul_4*emul_4 ), cmax ) ! Bermejo & Staniforth
      end do      ! iq loop
    end do        ! k loop
#ifdef _OPENACC
    !$acc end parallel loop
#else
    !$omp end do nowait
#endif

  end if ! nfield<mh_bs

!========================   end of intsch=1 section ====================
else     ! if(intsch==1)then
!======================== start of intsch=2 section ====================

  if ( nfield<mh_bs ) then
      
#ifdef _OPENACC
        !$acc parallel loop collapse(2) copyout(s) present(sx,xg,yg,nface) async(async_counter)
#else
        !$omp do schedule(static) private(k,iq,idel,xxg,jdel,yyg)                        &
        !$omp private(n,cmul_1,cmul_2,cmul_3,cmul_4,dmul_2,dmul_3,emul_1,emul_2,emul_3)  &
        !$omp private(emul_4,rmul_1,rmul_2,rmul_3,rmul_4)                                &
        !$omp private(sx_0m,sx_1m,sx_m0,sx_00,sx_10,sx_20,sx_m1,sx_01)                   &
        !$omp private(sx_11,sx_21,sx_02,sx_12)
#endif
        do k = 1,kl
          do iq = 1,ifull    ! non Berm-Stan option
            ! Convert face index from 0:npanels to array indices
            idel = int(xg(iq,k))
            xxg = xg(iq,k) - real(idel)
            jdel = int(yg(iq,k))
            yyg = yg(iq,k) - real(jdel)
            idel = min( max(idel - ioff, 0), ipan)
            jdel = min( max(jdel - joff, 0), jpan)
            n = min( max(nface(iq,k) + noff, 1), npan)
            ! bi-cubic
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
            sx_0m = sx(idel,  jdel-1,n,k)
            sx_1m = sx(idel+1,jdel-1,n,k)
            sx_m0 = sx(idel-1,jdel,  n,k)
            sx_00 = sx(idel,  jdel,  n,k)
            sx_10 = sx(idel+1,jdel,  n,k)
            sx_20 = sx(idel+2,jdel,  n,k)
            sx_m1 = sx(idel-1,jdel+1,n,k)
            sx_01 = sx(idel,  jdel+1,n,k)
            sx_11 = sx(idel+1,jdel+1,n,k)
            sx_21 = sx(idel+2,jdel+1,n,k)
            sx_02 = sx(idel,  jdel+2,n,k)
            sx_12 = sx(idel+1,jdel+2,n,k)
            rmul_1 = sx_m0*dmul_2 + sx_m1*dmul_3
            rmul_2 = sx_0m*cmul_1 + sx_00*cmul_2 + &
                     sx_01*cmul_3 + sx_02*cmul_4
            rmul_3 = sx_1m*cmul_1 + sx_10*cmul_2 + &
                     sx_11*cmul_3 + sx_12*cmul_4
            rmul_4 = sx_20*dmul_2 + sx_21*dmul_3
            s(iq,k) = rmul_1*emul_1 + rmul_2*emul_2 + rmul_3*emul_3 + rmul_4*emul_4
          end do       ! iq loop
        end do         ! k loop
#ifdef _OPENACC
        !$acc end parallel loop
#else
        !$omp end do nowait
#endif
      
  else

#ifdef _OPENACC
    !$acc parallel loop collapse(2) copyout(s) present(sx,xg,yg,nface) async(async_counter)
#else
    !$omp do schedule(static) private(k,iq,idel,xxg,jdel,yyg)                        &
    !$omp private(n,cmul_1,cmul_2,cmul_3,cmul_4,dmul_2,dmul_3,emul_1,emul_2,emul_3)  &
    !$omp private(emul_4,rmul_1,rmul_2,rmul_3,rmul_4,cmax,cmin)                      &
    !$omp private(sx_0m,sx_1m,sx_m0,sx_00,sx_10,sx_20,sx_m1,sx_01)                   &
    !$omp private(sx_11,sx_21,sx_02,sx_12)
#endif
    do k = 1,kl
      do iq = 1,ifull    ! Berm-Stan option here e.g. qg & gases
        idel = int(xg(iq,k))
        xxg = xg(iq,k) - real(idel)
        jdel = int(yg(iq,k))
        yyg = yg(iq,k) - real(jdel)
        idel = min( max(idel - ioff, 0), ipan)
        jdel = min( max(jdel - joff, 0), jpan)
        n = min( max(nface(iq,k) + noff, 1), npan)
        ! bi-cubic
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
        sx_0m = sx(idel,  jdel-1,n,k)
        sx_1m = sx(idel+1,jdel-1,n,k)
        sx_m0 = sx(idel-1,jdel,  n,k)
        sx_00 = sx(idel,  jdel,  n,k)
        sx_10 = sx(idel+1,jdel,  n,k)
        sx_20 = sx(idel+2,jdel,  n,k)
        sx_m1 = sx(idel-1,jdel+1,n,k)
        sx_01 = sx(idel,  jdel+1,n,k)
        sx_11 = sx(idel+1,jdel+1,n,k)
        sx_21 = sx(idel+2,jdel+1,n,k)
        sx_02 = sx(idel,  jdel+2,n,k)
        sx_12 = sx(idel+1,jdel+2,n,k)
        cmin = min(sx_00,sx_01,sx_10,sx_11)
        cmax = max(sx_00,sx_01,sx_10,sx_11)
        rmul_1 = sx_m0*dmul_2 + sx_m1*dmul_3
        rmul_2 = sx_0m*cmul_1 + sx_00*cmul_2 + &
                 sx_01*cmul_3 + sx_02*cmul_4
        rmul_3 = sx_1m*cmul_1 + sx_10*cmul_2 + &
                 sx_11*cmul_3 + sx_12*cmul_4
        rmul_4 = sx_20*dmul_2 + sx_21*dmul_3
        s(iq,k) = min( max( cmin, &
            rmul_1*emul_1 + rmul_2*emul_2 + rmul_3*emul_3 + rmul_4*emul_4 ), cmax ) ! Bermejo & Staniforth
      end do       ! iq loop
    end do         ! k loop
#ifdef _OPENACC
    !$acc end parallel loop
#else
    !$omp end do nowait
#endif
      
  end if
  
end if               ! (intsch==1) .. else ..
!========================   end of intsch=1 section ====================

return
end subroutine calcdephost
    
end module ints_m