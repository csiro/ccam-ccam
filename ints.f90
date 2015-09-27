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
    
subroutine ints(ntr,s,intsch,nface,xg,yg,nfield)

use cc_mpi             ! CC MPI routines
use indices_m          ! Grid index arrays

implicit none

include 'newmpar.h'    ! Grid parameters
include 'parm.h'       ! Model configuration
include 'parmhor.h'    ! Horizontal advection parameters

integer, intent(in) :: ntr     ! number of tracers to be interpolated
integer, intent(in) :: intsch  ! method to interpolate panel corners
integer, intent(in) :: nfield  ! use B&S if nfield>=mh_bs
integer idel, iq, jdel, nn
integer i, j, k, n, ii
integer, dimension(ifull,kl), intent(in) :: nface        ! interpolation coordinates
real xxg, yyg, cmin, cmax
real, dimension(ifull,kl), intent(in) :: xg, yg          ! interpolation coordinates
real, dimension(ifull+iextra,kl,ntr), intent(inout) :: s ! array of tracers
real, dimension(-1:ipan+2,-1:jpan+2,1:npan,kl,ntr) :: sx ! unpacked tracer array
real, dimension(4) :: cmul, emul, rmul
real, dimension(2:3) :: dmul

call START_LOG(ints_begin)

call bounds(s,nrows=2)

!======================== start of intsch=1 section ====================
if ( intsch == 1 ) then

  ! MJT notes - here we use JLM's unpacking of the indirect addressed array to a direct addressed
  ! array.
  sx(1:ipan,1:jpan,1:npan,1:kl,1:ntr) = reshape(s(1:ipan*jpan*npan,1:kl,1:ntr), (/ipan,jpan,npan,kl,ntr/))
  ! this is intsb           EW interps done first
  ! first extend s arrays into sx - this one -1:il+2 & -1:il+2
  do n = 1,npan
    do j = 1,jpan
      sx(0,j,n,:,:)      = s(iw(1+(j-1)*ipan+(n-1)*ipan*jpan),:,:)
      sx(-1,j,n,:,:)     = s(iww(1+(j-1)*ipan+(n-1)*ipan*jpan),:,:)
      sx(ipan+1,j,n,:,:) = s(ie(j*ipan+(n-1)*ipan*jpan),:,:)
      sx(ipan+2,j,n,:,:) = s(iee(j*ipan+(n-1)*ipan*jpan),:,:)
    end do            ! j loop
    do i = 1,ipan
      sx(i,0,n,:,:)      = s(is(i+(n-1)*ipan*jpan),:,:)
      sx(i,-1,n,:,:)     = s(iss(i+(n-1)*ipan*jpan),:,:)
      sx(i,jpan+1,n,:,:) = s(in(i-ipan+n*ipan*jpan),:,:)
      sx(i,jpan+2,n,:,:) = s(inn(i-ipan+n*ipan*jpan),:,:)
    end do            ! i loop
    ! for ew interpolation, sometimes need (different from ns):
    ! (-1,0),   (0,0),   (0,-1)   (-1,il+1),   (0,il+1),   (0,il+2)
    ! (il+1,0),(il+2,0),(il+1,-1) (il+1,il+1),(il+2,il+1),(il+1,il+2)
    sx(-1,0,n,:,:)          = s(lwws(n),:,:)
    sx(0,0,n,:,:)           = s(iws(1+(n-1)*ipan*jpan),:,:)
    sx(0,-1,n,:,:)          = s(lwss(n),:,:)
    sx(ipan+1,0,n,:,:)      = s(ies(ipan+(n-1)*ipan*jpan),:,:)
    sx(ipan+2,0,n,:,:)      = s(lees(n),:,:)
    sx(ipan+1,-1,n,:,:)     = s(less(n),:,:)
    sx(-1,jpan+1,n,:,:)     = s(lwwn(n),:,:)
    sx(0,jpan+2,n,:,:)      = s(lwnn(n),:,:)
    sx(ipan+2,jpan+1,n,:,:) = s(leen(n),:,:)
    sx(ipan+1,jpan+2,n,:,:) = s(lenn(n),:,:)
    sx(0,jpan+1,n,:,:)      = s(iwn(1-ipan+n*ipan*jpan),:,:)
    sx(ipan+1,jpan+1,n,:,:) = s(ien(n*ipan*jpan),:,:)
  end do              ! n loop

! Loop over points that need to be calculated for other processes
  if ( nfield < mh_bs ) then
    do ii = neighnum,1,-1
      do iq = 1,drlen(ii)
        
        ! depature point coordinates
        n = nint(dpoints(ii)%a(1,iq)) + noff ! Local index
        idel = int(dpoints(ii)%a(2,iq))
        xxg = dpoints(ii)%a(2,iq) - idel
        jdel = int(dpoints(ii)%a(3,iq))
        yyg = dpoints(ii)%a(3,iq) - jdel
        k = nint(dpoints(ii)%a(4,iq))
        idel = idel - ioff
        jdel = jdel - joff
              
        ! bi-cubic
        cmul(1) = (1.-xxg)*(2.-xxg)*(-xxg)/6.
        cmul(2) = (1.-xxg)*(2.-xxg)*(1.+xxg)/2.
        cmul(3) = xxg*(1.+xxg)*(2.-xxg)/2.
        cmul(4) = (1.-xxg)*(-xxg)*(1.+xxg)/6.
        dmul(2) = (1.-xxg)
        dmul(3) = xxg
        emul(1) = (1.-yyg)*(2.-yyg)*(-yyg)/6.
        emul(2) = (1.-yyg)*(2.-yyg)*(1.+yyg)/2.
        emul(3) = yyg*(1.+yyg)*(2.-yyg)/2.
        emul(4) = (1.-yyg)*(-yyg)*(1.+yyg)/6.
        do nn = 1,ntr
          rmul(1) = sum(sx(idel:idel+1,  jdel-1,n,k,nn)*dmul(2:3))
          rmul(2) = sum(sx(idel-1:idel+2,jdel,  n,k,nn)*cmul(1:4))
          rmul(3) = sum(sx(idel-1:idel+2,jdel+1,n,k,nn)*cmul(1:4))
          rmul(4) = sum(sx(idel:idel+1,  jdel+2,n,k,nn)*dmul(2:3))
          sextra(ii)%a(nn+(iq-1)*ntr) = sum(rmul(1:4)*emul(1:4))
        end do
        
      end do        ! iq loop
    end do          ! ii loop
  else              ! (nfield<mh_bs)
    do ii = neighnum,1,-1
      do iq = 1,drlen(ii)
        n = nint(dpoints(ii)%a(1,iq)) + noff ! Local index
        idel = int(dpoints(ii)%a(2,iq))
        xxg = dpoints(ii)%a(2,iq) - idel
        jdel = int(dpoints(ii)%a(3,iq))
        yyg = dpoints(ii)%a(3,iq) - jdel
        k = nint(dpoints(ii)%a(4,iq))
        idel = idel - ioff
        jdel = jdel - joff

        ! bi-cubic
        cmul(1) = (1.-xxg)*(2.-xxg)*(-xxg)/6.
        cmul(2) = (1.-xxg)*(2.-xxg)*(1.+xxg)/2.
        cmul(3) = xxg*(1.+xxg)*(2.-xxg)/2.
        cmul(4) = (1.-xxg)*(-xxg)*(1.+xxg)/6.
        dmul(2) = (1.-xxg)
        dmul(3) = xxg
        emul(1) = (1.-yyg)*(2.-yyg)*(-yyg)/6.
        emul(2) = (1.-yyg)*(2.-yyg)*(1.+yyg)/2.
        emul(3) = yyg*(1.+yyg)*(2.-yyg)/2.
        emul(4) = (1.-yyg)*(-yyg)*(1.+yyg)/6.
        do nn = 1,ntr
          cmin = minval(sx(idel:idel+1,jdel:jdel+1,n,k,nn))
          cmax = maxval(sx(idel:idel+1,jdel:jdel+1,n,k,nn))
          rmul(1) = sum(sx(idel:idel+1,  jdel-1,n,k,nn)*dmul(2:3))
          rmul(2) = sum(sx(idel-1:idel+2,jdel,  n,k,nn)*cmul(1:4))
          rmul(3) = sum(sx(idel-1:idel+2,jdel+1,n,k,nn)*cmul(1:4))
          rmul(4) = sum(sx(idel:idel+1,  jdel+2,n,k,nn)*dmul(2:3))
          sextra(ii)%a(nn+(iq-1)*ntr) = min( max( cmin, sum(rmul(1:4)*emul(1:4)) ), cmax ) ! Bermejo & Staniforth
        end do
       
      end do        ! iq loop
    end do          ! ii loop
  end if            ! (nfield<mh_bs)  .. else ..

  ! Send messages to other processors.  We then start the calculation for this processor while waiting for
  ! the messages to return, thereby overlapping computation with communication.
  call intssync_send(ntr)

  if ( nfield < mh_bs ) then
    do k = 1,kl
      do iq = 1,ifull    ! non Berm-Stan option
        idel = int(xg(iq,k))
        xxg = xg(iq,k) - idel
        jdel = int(yg(iq,k))
        yyg = yg(iq,k) - jdel
        idel = idel - ioff
        jdel = jdel - joff
        n = nface(iq,k) + noff ! Make this a local index

        if ( idel>=0 .and. idel<=ipan .and. jdel>=0 .and. jdel<=jpan .and. n>=1 .and. n<=npan ) then
          ! bi-cubic
          cmul(1) = (1.-xxg)*(2.-xxg)*(-xxg)/6.
          cmul(2) = (1.-xxg)*(2.-xxg)*(1.+xxg)/2.
          cmul(3) = xxg*(1.+xxg)*(2.-xxg)/2.
          cmul(4) = (1.-xxg)*(-xxg)*(1.+xxg)/6.
          dmul(2) = (1.-xxg)
          dmul(3) = xxg
          emul(1) = (1.-yyg)*(2.-yyg)*(-yyg)/6.
          emul(2) = (1.-yyg)*(2.-yyg)*(1.+yyg)/2.
          emul(3) = yyg*(1.+yyg)*(2.-yyg)/2.
          emul(4) = (1.-yyg)*(-yyg)*(1.+yyg)/6.
          do nn = 1,ntr
            rmul(1) = sum(sx(idel:idel+1,  jdel-1,n,k,nn)*dmul(2:3))
            rmul(2) = sum(sx(idel-1:idel+2,jdel,  n,k,nn)*cmul(1:4))
            rmul(3) = sum(sx(idel-1:idel+2,jdel+1,n,k,nn)*cmul(1:4))
            rmul(4) = sum(sx(idel:idel+1,  jdel+2,n,k,nn)*dmul(2:3))
            s(iq,k,nn) = sum(rmul(1:4)*emul(1:4))
          end do
        end if
        
      end do         ! iq loop
    end do           ! k loop
  else               ! (nfield<mh_bs)
    do k = 1,kl
      do iq = 1,ifull    ! Berm-Stan option here e.g. qg & gases
        idel = int(xg(iq,k))
        xxg = xg(iq,k) - idel
        jdel = int(yg(iq,k))
        yyg = yg(iq,k) - jdel
        ! Now make them proper indices in this processor's region
        idel = idel - ioff
        jdel = jdel - joff
        n = nface(iq,k) + noff ! Make this a local index

        if ( idel>=0 .and. idel<=ipan .and. jdel>=0 .and. jdel<=jpan .and. n>=1 .and. n<=npan ) then
          ! bi-cubic
          cmul(1) = (1.-xxg)*(2.-xxg)*(-xxg)/6.
          cmul(2) = (1.-xxg)*(2.-xxg)*(1.+xxg)/2.
          cmul(3) = xxg*(1.+xxg)*(2.-xxg)/2.
          cmul(4) = (1.-xxg)*(-xxg)*(1.+xxg)/6.
          dmul(2) = (1.-xxg)
          dmul(3) = xxg
          emul(1) = (1.-yyg)*(2.-yyg)*(-yyg)/6.
          emul(2) = (1.-yyg)*(2.-yyg)*(1.+yyg)/2.
          emul(3) = yyg*(1.+yyg)*(2.-yyg)/2.
          emul(4) = (1.-yyg)*(-yyg)*(1.+yyg)/6.
          do nn = 1,ntr
            cmin = minval(sx(idel:idel+1,jdel:jdel+1,n,k,nn))
            cmax = maxval(sx(idel:idel+1,jdel:jdel+1,n,k,nn))
            rmul(1) = sum(sx(idel:idel+1,  jdel-1,n,k,nn)*dmul(2:3))
            rmul(2) = sum(sx(idel-1:idel+2,jdel,  n,k,nn)*cmul(1:4))
            rmul(3) = sum(sx(idel-1:idel+2,jdel+1,n,k,nn)*cmul(1:4))
            rmul(4) = sum(sx(idel:idel+1,  jdel+2,n,k,nn)*dmul(2:3))
            s(iq,k,nn) = min( max( cmin, sum(rmul(1:4)*emul(1:4)) ), cmax ) ! Bermejo & Staniforth
          end do
        end if
      
      end do        ! iq loop
    end do          ! k loop
  end if            ! (nfield<mh_bs)  .. else ..
            
!========================   end of intsch=1 section ====================
else     ! if(intsch==1)then
!======================== start of intsch=2 section ====================
!       this is intsc           NS interps done first
!       first extend s arrays into sx - this one -1:il+2 & -1:il+2
  sx(1:ipan,1:jpan,1:npan,1:kl,1:ntr) = reshape(s(1:ipan*jpan*npan,1:kl,1:ntr), (/ipan,jpan,npan,kl,ntr/))
  do n = 1,npan
    do j = 1,jpan
      sx(0,j,n,:,:)      = s(iw(1+(j-1)*ipan+(n-1)*ipan*jpan),:,:)
      sx(-1,j,n,:,:)     = s(iww(1+(j-1)*ipan+(n-1)*ipan*jpan),:,:)
      sx(ipan+1,j,n,:,:) = s(ie(j*ipan+(n-1)*ipan*jpan),:,:)
      sx(ipan+2,j,n,:,:) = s(iee(j*ipan+(n-1)*ipan*jpan),:,:)
    end do            ! j loop
    do i = 1,ipan
      sx(i,0,n,:,:)      = s(is(i+(n-1)*ipan*jpan),:,:)
      sx(i,-1,n,:,:)     = s(iss(i+(n-1)*ipan*jpan),:,:)
      sx(i,jpan+1,n,:,:) = s(in(i-ipan+n*ipan*jpan),:,:)
      sx(i,jpan+2,n,:,:) = s(inn(i-ipan+n*ipan*jpan),:,:)
    end do            ! i loop
!   for ns interpolation, sometimes need (different from ew):
!       (-1,0),   (0,0),   (0,-1)   (-1,il+1),   (0,il+1),   (0,il+2)
!     (il+1,0),(il+2,0),(il+1,-1) (il+1,il+1),(il+2,il+1),(il+1,il+2)

    sx(-1,0,n,:,:)          = s(lsww(n),:,:)
    sx(0,0,n,:,:)           = s(isw(1+(n-1)*ipan*jpan),:,:)
    sx(0,-1,n,:,:)          = s(lssw(n),:,:)
    sx(ipan+2,0,n,:,:)      = s(lsee(n),:,:)
    sx(ipan+1,-1,n,:,:)     = s(lsse(n),:,:)
    sx(-1,jpan+1,n,:,:)     = s(lnww(n),:,:)
    sx(0,jpan+1,n,:,:)      = s(inw(1-ipan+n*ipan*jpan),:,:)
    sx(0,jpan+2,n,:,:)      = s(lnnw(n),:,:)
    sx(ipan+2,jpan+1,n,:,:) = s(lnee(n),:,:)
    sx(ipan+1,jpan+2,n,:,:) = s(lnne(n),:,:)
    sx(ipan+1,0,n,:,:)      = s(ise(ipan+(n-1)*ipan*jpan),:,:)
    sx(ipan+1,jpan+1,n,:,:) = s(ine(n*ipan*jpan),:,:)
  end do              ! n loop

  ! For other processes
  if ( nfield < mh_bs ) then
    do ii = neighnum,1,-1
      do iq = 1,drlen(ii)
        n = nint(dpoints(ii)%a(1,iq)) + noff ! Local index
        !  Need global face index in fproc call
        idel = int(dpoints(ii)%a(2,iq))
        xxg = dpoints(ii)%a(2,iq) - idel
        jdel = int(dpoints(ii)%a(3,iq))
        yyg = dpoints(ii)%a(3,iq) - jdel
        k = nint(dpoints(ii)%a(4,iq))
        idel = idel - ioff
        jdel = jdel - joff

        ! bi-cubic
        cmul(1) = (1.-yyg)*(2.-yyg)*(-yyg)/6.
        cmul(2) = (1.-yyg)*(2.-yyg)*(1.+yyg)/2.
        cmul(3) = yyg*(1.+yyg)*(2.-yyg)/2.
        cmul(4) = (1.-yyg)*(-yyg)*(1.+yyg)/6.
        dmul(2) = (1.-yyg)
        dmul(3) = yyg
        emul(1) = (1.-xxg)*(2.-xxg)*(-xxg)/6.
        emul(2) = (1.-xxg)*(2.-xxg)*(1.+xxg)/2.
        emul(3) = xxg*(1.+xxg)*(2.-xxg)/2.
        emul(4) = (1.-xxg)*(-xxg)*(1.+xxg)/6.
        do nn = 1,ntr
          rmul(1) = sum(sx(idel-1,jdel:jdel+1,  n,k,nn)*dmul(2:3))
          rmul(2) = sum(sx(idel,  jdel-1:jdel+2,n,k,nn)*cmul(1:4))
          rmul(3) = sum(sx(idel+1,jdel-1:jdel+2,n,k,nn)*cmul(1:4))
          rmul(4) = sum(sx(idel+2,jdel:jdel+1,  n,k,nn)*dmul(2:3))
          sextra(ii)%a(nn+(iq-1)*ntr) = sum(rmul(1:4)*emul(1:4))
        end do
        
      end do           ! iq loop
    end do             ! ii
  else                 ! (nfield<mh_bs)
    do ii = neighnum,1,-1
      do iq = 1,drlen(ii)
        n = nint(dpoints(ii)%a(1,iq)) + noff ! Local index
        !  Need global face index in fproc call
        idel = int(dpoints(ii)%a(2,iq))
        xxg = dpoints(ii)%a(2,iq) - idel
        jdel = int(dpoints(ii)%a(3,iq))
        yyg = dpoints(ii)%a(3,iq) - jdel
        k = nint(dpoints(ii)%a(4,iq))
        idel = idel - ioff
        jdel = jdel - joff

        ! bi-cubic
        cmul(1) = (1.-yyg)*(2.-yyg)*(-yyg)/6.
        cmul(2) = (1.-yyg)*(2.-yyg)*(1.+yyg)/2.
        cmul(3) = yyg*(1.+yyg)*(2.-yyg)/2.
        cmul(4) = (1.-yyg)*(-yyg)*(1.+yyg)/6.
        dmul(2) = (1.-yyg)
        dmul(3) = yyg
        emul(1) = (1.-xxg)*(2.-xxg)*(-xxg)/6.
        emul(2) = (1.-xxg)*(2.-xxg)*(1.+xxg)/2.
        emul(3) = xxg*(1.+xxg)*(2.-xxg)/2.
        emul(4) = (1.-xxg)*(-xxg)*(1.+xxg)/6.
        do nn = 1,ntr
          cmin = minval(sx(idel:idel+1,jdel:jdel+1,n,k,nn))
          cmax = maxval(sx(idel:idel+1,jdel:jdel+1,n,k,nn))
          rmul(1) = sum(sx(idel-1,jdel:jdel+1,  n,k,nn)*dmul(2:3))
          rmul(2) = sum(sx(idel,  jdel-1:jdel+2,n,k,nn)*cmul(1:4))
          rmul(3) = sum(sx(idel+1,jdel-1:jdel+2,n,k,nn)*cmul(1:4))
          rmul(4) = sum(sx(idel+2,jdel:jdel+1,  n,k,nn)*dmul(2:3))
          sextra(ii)%a(nn+(iq-1)*ntr) = min( max( cmin, sum(rmul(1:4)*emul(1:4)) ), cmax ) ! Bermejo & Staniforth
        end do

      end do        ! iq loop
    end do          ! ii loop
  end if            ! (nfield<mh_bs)  .. else ..

  call intssync_send(ntr)

  if ( nfield < mh_bs ) then
    do k = 1,kl
      do iq = 1,ifull    ! non Berm-Stan option
        ! Convert face index from 0:npanels to array indices
        idel = int(xg(iq,k))
        xxg = xg(iq,k) - idel
        jdel = int(yg(iq,k))
        yyg = yg(iq,k) - jdel
        ! Now make them proper indices in this processor's region
        idel = idel - ioff
        jdel = jdel - joff
        n = nface(iq,k) + noff ! Make this a local index

        if ( idel>=0 .and. idel<=ipan .and. jdel>=0 .and. jdel<=jpan .and. n>=1 .and. n<=npan ) then
          ! bi-cubic
          cmul(1) = (1.-yyg)*(2.-yyg)*(-yyg)/6.
          cmul(2) = (1.-yyg)*(2.-yyg)*(1.+yyg)/2.
          cmul(3) = yyg*(1.+yyg)*(2.-yyg)/2.
          cmul(4) = (1.-yyg)*(-yyg)*(1.+yyg)/6.
          dmul(2) = (1.-yyg)
          dmul(3) = yyg
          emul(1) = (1.-xxg)*(2.-xxg)*(-xxg)/6.
          emul(2) = (1.-xxg)*(2.-xxg)*(1.+xxg)/2.
          emul(3) = xxg*(1.+xxg)*(2.-xxg)/2.
          emul(4) = (1.-xxg)*(-xxg)*(1.+xxg)/6.
          do nn = 1,ntr
            rmul(1) = sum(sx(idel-1,jdel:jdel+1,  n,k,nn)*dmul(2:3))
            rmul(2) = sum(sx(idel,  jdel-1:jdel+2,n,k,nn)*cmul(1:4))
            rmul(3) = sum(sx(idel+1,jdel-1:jdel+2,n,k,nn)*cmul(1:4))
            rmul(4) = sum(sx(idel+2,jdel:jdel+1,  n,k,nn)*dmul(2:3))
            s(iq,k,nn) = sum(rmul(1:4)*emul(1:4))
          end do
        end if
        
      end do         ! iq loop
    end do           ! k loop
  else               ! (nfield<mh_bs)
    do k = 1,kl
      do iq = 1,ifull    ! Berm-Stan option here e.g. qg & gases
        idel = int(xg(iq,k))
        xxg = xg(iq,k) - idel
        jdel = int(yg(iq,k))
        yyg = yg(iq,k) - jdel
        ! Now make them proper indices in this processor's region
        idel = idel - ioff
        jdel = jdel - joff
        n = nface(iq,k) + noff ! Make this a local index

        if ( idel>=0 .and. idel<=ipan .and. jdel>=0 .and. jdel<=jpan .and. n>=1 .and. n<=npan ) then
          ! bi-cubic
          cmul(1) = (1.-yyg)*(2.-yyg)*(-yyg)/6.
          cmul(2) = (1.-yyg)*(2.-yyg)*(1.+yyg)/2.
          cmul(3) = yyg*(1.+yyg)*(2.-yyg)/2.
          cmul(4) = (1.-yyg)*(-yyg)*(1.+yyg)/6.
          dmul(2) = (1.-yyg)
          dmul(3) = yyg
          emul(1) = (1.-xxg)*(2.-xxg)*(-xxg)/6.
          emul(2) = (1.-xxg)*(2.-xxg)*(1.+xxg)/2.
          emul(3) = xxg*(1.+xxg)*(2.-xxg)/2.
          emul(4) = (1.-xxg)*(-xxg)*(1.+xxg)/6.
          do nn = 1,ntr
            cmin = minval(sx(idel:idel+1,jdel:jdel+1,n,k,nn))
            cmax = maxval(sx(idel:idel+1,jdel:jdel+1,n,k,nn))
            rmul(1) = sum(sx(idel-1,jdel:jdel+1,  n,k,nn)*dmul(2:3))
            rmul(2) = sum(sx(idel,  jdel-1:jdel+2,n,k,nn)*cmul(1:4))
            rmul(3) = sum(sx(idel+1,jdel-1:jdel+2,n,k,nn)*cmul(1:4))
            rmul(4) = sum(sx(idel+2,jdel:jdel+1,  n,k,nn)*dmul(2:3))
            s(iq,k,nn) = min( max( cmin, sum(rmul(1:4)*emul(1:4)) ), cmax ) ! Bermejo & Staniforth
          end do
        end if
        
      end do         ! iq loop
    end do           ! k loop
  end if             ! (nfield<mh_bs)  .. else ..

end if               ! (intsch==1) .. else ..
!========================   end of intsch=1 section ====================

call intssync_recv(s)
      
call END_LOG(ints_end)
return
end subroutine ints

subroutine ints_bl(s,intsch,nface,xg,yg)  ! not usually called

use cc_mpi             ! CC MPI routines
use indices_m          ! Grid index arrays

implicit none

include 'newmpar.h'    ! Grid parameters
include 'parm.h'       ! Model configuration
include 'parmhor.h'    ! Horizontal advection parameters

integer idel, iq, jdel
integer i, j, k, n
integer ii
integer, intent(in) :: intsch
integer, dimension(ifull,kl), intent(in) :: nface
real xxg, yyg
real, dimension(ifull,kl), intent(in) :: xg,yg      ! now passed through call
real, dimension(0:ipan+1,0:jpan+1,1:npan,kl) :: sx
real, dimension(ifull+iextra,kl,1), intent(inout) :: s

!     this one does bi-linear interpolation only
!     first extend s arrays into sx - this one -1:il+2 & -1:il+2
!                    but for bi-linear only need 0:il+1 &  0:il+1
call START_LOG(ints_begin)
call bounds(s,corner=.true.)
sx(1:ipan,1:jpan,1:npan,1:kl) = reshape(s(1:ipan*jpan*npan,1:kl,1), (/ipan,jpan,npan,kl/))
do n = 1,npan
  do j = 1,jpan
    sx(0,j,n,:)      = s(iw(1+(j-1)*ipan+(n-1)*ipan*jpan),:,1)
    sx(ipan+1,j,n,:) = s(ie(j*ipan+(n-1)*ipan*jpan),:,1)
  end do               ! j loop
  do i = 1,ipan
    sx(i,0,n,:)      = s(is(i+(n-1)*ipan*jpan),:,1)
    sx(i,jpan+1,n,:) = s(in(i-ipan+n*ipan*jpan),:,1)
  end do               ! i loop
  sx(0,0,n,:)           = s(iws(1+(n-1)*ipan*jpan),:,1)
  sx(ipan+1,0,n,:)      = s(ies(ipan+(n-1)*ipan*jpan),:,1)
  sx(0,jpan+1,n,:)      = s(iwn(1-ipan+n*ipan*jpan),:,1)
  sx(ipan+1,jpan+1,n,:) = s(ien(n*ipan*jpan),:,1)
end do                 ! n loop

! Loop over points that need to be calculated for other processes
do ii = neighnum,1,-1
  do iq = 1,drlen(ii)
    n = nint(dpoints(ii)%a(1,iq)) + noff ! Local index
    !  Need global face index in fproc call
    idel = int(dpoints(ii)%a(2,iq))
    xxg = dpoints(ii)%a(2,iq) - idel
    jdel = int(dpoints(ii)%a(3,iq))
    yyg = dpoints(ii)%a(3,iq) - jdel
    k = nint(dpoints(ii)%a(4,iq))
    idel = idel - ioff
    jdel = jdel - joff
    sextra(ii)%a(iq) =      yyg*(xxg*sx(idel+1,jdel+1,n,k)+(1.-xxg)*sx(idel,jdel+1,n,k)) &
                     + (1.-yyg)*(xxg*sx(idel+1,  jdel,n,k)+(1.-xxg)*sx(idel,  jdel,n,k))
  end do
end do

call intssync_send(1)

do k = 1,kl
  do iq = 1,ifull
    ! Convert face index from 0:npanels to array indices
    idel=int(xg(iq,k))
    xxg=xg(iq,k)-idel
    jdel=int(yg(iq,k))
    yyg=yg(iq,k)-jdel
    ! Now make them proper indices in this processor's region
    idel = idel - ioff
    jdel = jdel - joff
    n = nface(iq,k) + noff ! Make this a local index

    if ( idel>=0 .and. idel<=ipan .and. jdel>=0 .and. jdel<=jpan .and. n>=1 .and. n<=npan ) then
      s(iq,k,1) =      yyg*(xxg*sx(idel+1,jdel+1,n,k)+(1.-xxg)*sx(idel,jdel+1,n,k)) &
                + (1.-yyg)*(xxg*sx(idel+1,  jdel,n,k)+(1.-xxg)*sx(idel,  jdel,n,k))
    end if
    
  end do                  ! iq loop
end do                    ! k

call intssync_recv(s)

call END_LOG(ints_end)
return
end subroutine ints_bl
