! Conformal Cubic Atmospheric Model
    
! Copyright 2015-2020 Commonwealth Scientific Industrial Research Organisation (CSIRO)
    
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
    
subroutine ints(s,intsch,nface,xg,yg,nfield)

use cc_mpi             ! CC MPI routines
use indices_m          ! Grid index arrays
use newmpar_m          ! Grid parameters
use parm_m             ! Model configuration
use parmhor_m          ! Horizontal advection parameters

implicit none

integer, intent(in) :: intsch  ! method to interpolate panel corners
integer, intent(in) :: nfield  ! use B&S if nfield>=mh_bs
integer idel, iq, jdel
integer i, j, k, n, ii
integer, dimension(ifull,kl), intent(in) :: nface        ! interpolation coordinates
real, dimension(ifull,kl), intent(in) :: xg, yg          ! interpolation coordinates
real, dimension(ifull+iextra,kl), intent(inout) :: s ! array of tracers
real, dimension(-1:ipan+2,-1:jpan+2,1:npan,kl) :: sx ! unpacked tracer array
real xxg, yyg, cmin, cmax
real dmul_2, dmul_3, cmul_1, cmul_2, cmul_3, cmul_4
real emul_1, emul_2, emul_3, emul_4, rmul_1, rmul_2, rmul_3, rmul_4
!$acc declare present(xg,yg,nface)

call START_LOG(ints_begin)

call bounds(s,nrows=2) ! also includes corners

!$acc enter data create(sx,s)
!$acc update device(s)

!======================== start of intsch=1 section ====================
if ( intsch==1 ) then
  ! this is intsb           EW interps done first
  ! first extend s arrays into sx - this one -1:il+2 & -1:il+2

  !$acc parallel loop collapse(4)
  do concurrent (k = 1:kl)
    do concurrent (n = 1:npan)
      do concurrent (j = 1:jpan)
        do concurrent (i = 1:ipan)
          iq = i + (j-1)*ipan + (n-1)*ipan*jpan
          sx(i,j,n,k) = s(iq,k)
        end do
      end do
    end do
  end do
  !$acc end parallel loop
  !$acc parallel loop collapse(2)
  do concurrent (k = 1:kl)
    do concurrent (n = 1:npan)
      do concurrent (j = 1:jpan)
        iq = 1 + (j-1)*ipan + (n-1)*ipan*jpan
        sx(0,j,n,k)      = s(iw(iq),k)
        sx(-1,j,n,k)     = s(iww(iq),k)
        iq = j*ipan + (n-1)*ipan*jpan
        sx(ipan+1,j,n,k) = s(ie(iq),k)
        sx(ipan+2,j,n,k) = s(iee(iq),k)
      end do            ! j loop
      do concurrent (i = 1:ipan)
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
  !$acc end parallel loop

  !$acc update self(sx)

! Loop over points that need to be calculated for other processes
  if ( nfield<mh_bs ) then

    do concurrent (ii = 1:neighnum)
      do concurrent (iq = 1:drlen(ii))
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
        rmul_1 = sx(idel,  jdel-1,n,k)*dmul_2 + sx(idel+1,jdel-1,n,k)*dmul_3
        rmul_2 = sx(idel-1,jdel,  n,k)*cmul_1 + sx(idel,  jdel,  n,k)*cmul_2 + &
                 sx(idel+1,jdel,  n,k)*cmul_3 + sx(idel+2,jdel,  n,k)*cmul_4
        rmul_3 = sx(idel-1,jdel+1,n,k)*cmul_1 + sx(idel,  jdel+1,n,k)*cmul_2 + &
                 sx(idel+1,jdel+1,n,k)*cmul_3 + sx(idel+2,jdel+1,n,k)*cmul_4
        rmul_4 = sx(idel,  jdel+2,n,k)*dmul_2 + sx(idel+1,jdel+2,n,k)*dmul_3
        sextra(ii)%a(iq) = rmul_1*emul_1 + rmul_2*emul_2 + rmul_3*emul_3 + rmul_4*emul_4
      end do      ! iq loop
    end do          ! ii loop

    ! Send messages to other processors.  We then start the calculation for this processor while waiting for
    ! the messages to return, thereby overlapping computation with communication.
    call intssync_send(1)

    !$acc parallel loop collapse(2)
    do concurrent (k = 1:kl)
      do concurrent (iq = 1:ifull)    ! non Berm-Stan option
        idel = int(xg(iq,k))
        xxg = xg(iq,k) - idel
        jdel = int(yg(iq,k))
        yyg = yg(iq,k) - jdel
        idel = min( max( idel - ioff, 0), ipan )
        jdel = min( max( jdel - joff, 0), jpan )
        n = min( max( nface(iq,k) + noff, 1), npan )
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
        rmul_1 = sx(idel,  jdel-1,n,k)*dmul_2 + sx(idel+1,jdel-1,n,k)*dmul_3
        rmul_2 = sx(idel-1,jdel,  n,k)*cmul_1 + sx(idel,  jdel,  n,k)*cmul_2 + &
                 sx(idel+1,jdel,  n,k)*cmul_3 + sx(idel+2,jdel,  n,k)*cmul_4
        rmul_3 = sx(idel-1,jdel+1,n,k)*cmul_1 + sx(idel,  jdel+1,n,k)*cmul_2 + &
                 sx(idel+1,jdel+1,n,k)*cmul_3 + sx(idel+2,jdel+1,n,k)*cmul_4
        rmul_4 = sx(idel,  jdel+2,n,k)*dmul_2 + sx(idel+1,jdel+2,n,k)*dmul_3
        s(iq,k) = rmul_1*emul_1 + rmul_2*emul_2 + rmul_3*emul_3 + rmul_4*emul_4
      end do       ! iq loop
    end do         ! k loop
    !$acc end parallel loop
    
  else              ! (nfield<mh_bs)

    do concurrent (ii = 1:neighnum)
      do concurrent (iq = 1:drlen(ii))
        n = nint(dpoints(ii)%a(1,iq)) + noff ! Local index
        idel = int(dpoints(ii)%a(2,iq))
        xxg = dpoints(ii)%a(2,iq) - idel
        jdel = int(dpoints(ii)%a(3,iq))
        yyg = dpoints(ii)%a(3,iq) - jdel
        k = nint(dpoints(ii)%a(4,iq))
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
        cmin = min(sx(idel,  jdel,n,k),sx(idel+1,jdel,  n,k), &
                   sx(idel,jdel+1,n,k),sx(idel+1,jdel+1,n,k))
        cmax = max(sx(idel,  jdel,n,k),sx(idel+1,jdel,  n,k), &
                   sx(idel,jdel+1,n,k),sx(idel+1,jdel+1,n,k))
        rmul_1 = sx(idel,  jdel-1,n,k)*dmul_2 + sx(idel+1,jdel-1,n,k)*dmul_3
        rmul_2 = sx(idel-1,jdel,  n,k)*cmul_1 + sx(idel,  jdel,  n,k)*cmul_2 + &
                 sx(idel+1,jdel,  n,k)*cmul_3 + sx(idel+2,jdel,  n,k)*cmul_4
        rmul_3 = sx(idel-1,jdel+1,n,k)*cmul_1 + sx(idel,  jdel+1,n,k)*cmul_2 + &
                 sx(idel+1,jdel+1,n,k)*cmul_3 + sx(idel+2,jdel+1,n,k)*cmul_4
        rmul_4 = sx(idel,  jdel+2,n,k)*dmul_2 + sx(idel+1,jdel+2,n,k)*dmul_3
        sextra(ii)%a(iq) = min( max( cmin, &
            rmul_1*emul_1 + rmul_2*emul_2 + rmul_3*emul_3 + rmul_4*emul_4 ), cmax ) ! Bermejo & Staniforth
      end do      ! iq loop
    end do          ! ii loop

    ! Send messages to other processors.  We then start the calculation for this processor while waiting for
    ! the messages to return, thereby overlapping computation with communication.
    call intssync_send(1)

    !$acc parallel loop collapse(2) copyin(sx) copyout(s(1:ifull,1:kl))
    do concurrent (k = 1:kl)
      do concurrent (iq = 1:ifull)    ! Berm-Stan option here e.g. qg & gases
        idel = int(xg(iq,k))
        xxg = xg(iq,k) - idel
        jdel = int(yg(iq,k))
        yyg = yg(iq,k) - jdel
        idel = min( max( idel - ioff, 0), ipan )
        jdel = min( max( jdel - joff, 0), jpan )
        n = min( max( nface(iq,k) + noff, 1), npan )
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
        cmin = min(sx(idel,  jdel,n,k),sx(idel+1,jdel,  n,k), &
                   sx(idel,jdel+1,n,k),sx(idel+1,jdel+1,n,k))
        cmax = max(sx(idel,  jdel,n,k),sx(idel+1,jdel,  n,k), &
                   sx(idel,jdel+1,n,k),sx(idel+1,jdel+1,n,k))
        rmul_1 = sx(idel,  jdel-1,n,k)*dmul_2 + sx(idel+1,jdel-1,n,k)*dmul_3
        rmul_2 = sx(idel-1,jdel,  n,k)*cmul_1 + sx(idel,  jdel,  n,k)*cmul_2 + &
                 sx(idel+1,jdel,  n,k)*cmul_3 + sx(idel+2,jdel,  n,k)*cmul_4
        rmul_3 = sx(idel-1,jdel+1,n,k)*cmul_1 + sx(idel,  jdel+1,n,k)*cmul_2 + &
                 sx(idel+1,jdel+1,n,k)*cmul_3 + sx(idel+2,jdel+1,n,k)*cmul_4
        rmul_4 = sx(idel,  jdel+2,n,k)*dmul_2 + sx(idel+1,jdel+2,n,k)*dmul_3
        s(iq,k) = min( max( cmin, &
            rmul_1*emul_1 + rmul_2*emul_2 + rmul_3*emul_3 + rmul_4*emul_4 ), cmax ) ! Bermejo & Staniforth
      end do      ! iq loop
    end do        ! k loop
    !$acc end parallel loop
    
  end if            ! (nfield<mh_bs)  .. else ..
            
!========================   end of intsch=1 section ====================
else     ! if(intsch==1)then
!======================== start of intsch=2 section ====================
!       this is intsc           NS interps done first
!       first extend s arrays into sx - this one -1:il+2 & -1:il+2

  !$acc parallel loop collapse(4)
  do concurrent (k = 1:kl)
    do concurrent (n = 1:npan)
      do concurrent (j = 1:jpan)
        do concurrent (i = 1:ipan)
          iq = i + (j-1)*ipan + (n-1)*ipan*jpan
          sx(i,j,n,k) = s(iq,k)
        end do
      end do
    end do
  end do
  !$acc end parallel loop
  !$acc parallel loop collapse(2)
  do concurrent (k = 1:kl)
    do concurrent (n = 1:npan)
      do concurrent (j = 1:jpan)
        iq = 1+(j-1)*ipan+(n-1)*ipan*jpan
        sx(0,j,n,k)      = s(iw(iq),k)
        sx(-1,j,n,k)     = s(iww(iq),k)
        iq = j*ipan+(n-1)*ipan*jpan
        sx(ipan+1,j,n,k) = s(ie(iq),k)
        sx(ipan+2,j,n,k) = s(iee(iq),k)
      end do            ! j loop
      do concurrent (i = 1:ipan)
        iq = i+(n-1)*ipan*jpan
        sx(i,0,n,k)      = s(is(iq),k)
        sx(i,-1,n,k)     = s(iss(iq),k)
        iq = i-ipan+n*ipan*jpan
        sx(i,jpan+1,n,k) = s(in(iq),k)
        sx(i,jpan+2,n,k) = s(inn(iq),k)
      end do            ! i loop
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
  !$acc end parallel

  !$acc update self(sx)

  ! For other processes
  if ( nfield < mh_bs ) then

    do concurrent (ii = 1:neighnum)
      do concurrent (iq = 1:drlen(ii))
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
        rmul_1 = sx(idel-1,jdel,  n,k)*dmul_2 + sx(idel-1,jdel+1,n,k)*dmul_3
        rmul_2 = sx(idel,  jdel-1,n,k)*cmul_1 + sx(idel,  jdel,  n,k)*cmul_2 + &
                 sx(idel,  jdel+1,n,k)*cmul_3 + sx(idel,  jdel+2,n,k)*cmul_4
        rmul_3 = sx(idel+1,jdel-1,n,k)*cmul_1 + sx(idel+1,jdel,  n,k)*cmul_2 + &
                 sx(idel+1,jdel+1,n,k)*cmul_3 + sx(idel+1,jdel+2,n,k)*cmul_4
        rmul_4 = sx(idel+2,jdel,  n,k)*dmul_2 + sx(idel+2,jdel+1,n,k)*dmul_3
        sextra(ii)%a(iq) = rmul_1*emul_1 + rmul_2*emul_2 + rmul_3*emul_3 + rmul_4*emul_4
      end do         ! iq loop
    end do             ! ii
    
    call intssync_send(1)

    !$acc parallel loop collapse(2)
    do concurrent (k = 1:kl)
      do concurrent (iq = 1:ifull)    ! non Berm-Stan option
        ! Convert face index from 0:npanels to array indices
        idel = int(xg(iq,k))
        xxg = xg(iq,k) - idel
        jdel = int(yg(iq,k))
        yyg = yg(iq,k) - jdel
        idel = min( max( idel - ioff, 0), ipan )
        jdel = min( max( jdel - joff, 0), jpan )
        n = min( max( nface(iq,k) + noff, 1), npan )
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
        rmul_1 = sx(idel-1,jdel,  n,k)*dmul_2 + sx(idel-1,jdel+1,n,k)*dmul_3
        rmul_2 = sx(idel,  jdel-1,n,k)*cmul_1 + sx(idel,  jdel,  n,k)*cmul_2 + &
                 sx(idel,  jdel+1,n,k)*cmul_3 + sx(idel,  jdel+2,n,k)*cmul_4
        rmul_3 = sx(idel+1,jdel-1,n,k)*cmul_1 + sx(idel+1,jdel,  n,k)*cmul_2 + &
                 sx(idel+1,jdel+1,n,k)*cmul_3 + sx(idel+1,jdel+2,n,k)*cmul_4
        rmul_4 = sx(idel+2,jdel,  n,k)*dmul_2 + sx(idel+2,jdel+1,n,k)*dmul_3
        s(iq,k) = rmul_1*emul_1 + rmul_2*emul_2 + rmul_3*emul_3 + rmul_4*emul_4
      end do       ! iq loop
    end do         ! k loop
    !$acc end parallel loop
    
  else                 ! (nfield<mh_bs)

    do concurrent (ii = 1:neighnum)
      do concurrent (iq = 1:drlen(ii))
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
        cmin = min(sx(idel,  jdel,n,k),sx(idel+1,jdel,  n,k), &
                   sx(idel,jdel+1,n,k),sx(idel+1,jdel+1,n,k))
        cmax = max(sx(idel,  jdel,n,k),sx(idel+1,jdel,  n,k), &
                   sx(idel,jdel+1,n,k),sx(idel+1,jdel+1,n,k))
        rmul_1 = sx(idel-1,jdel,  n,k)*dmul_2 + sx(idel-1,jdel+1,n,k)*dmul_3
        rmul_2 = sx(idel,  jdel-1,n,k)*cmul_1 + sx(idel,  jdel,  n,k)*cmul_2 + &
                 sx(idel,  jdel+1,n,k)*cmul_3 + sx(idel,  jdel+2,n,k)*cmul_4
        rmul_3 = sx(idel+1,jdel-1,n,k)*cmul_1 + sx(idel+1,jdel,  n,k)*cmul_2 + &
                 sx(idel+1,jdel+1,n,k)*cmul_3 + sx(idel+1,jdel+2,n,k)*cmul_4
        rmul_4 = sx(idel+2,jdel,  n,k)*dmul_2 + sx(idel+2,jdel+1,n,k)*dmul_3
        sextra(ii)%a(iq) = min( max( cmin, &
            rmul_1*emul_1 + rmul_2*emul_2 + rmul_3*emul_3 + rmul_4*emul_4 ), cmax ) ! Bermejo & Staniforth
      end do      ! iq loop
    end do          ! ii loop
  
    call intssync_send(1)

    !$acc parallel loop collapse(2)
    do concurrent (k = 1:kl)
      do concurrent (iq = 1:ifull)    ! Berm-Stan option here e.g. qg & gases
        idel = int(xg(iq,k))
        xxg = xg(iq,k) - idel
        jdel = int(yg(iq,k))
        yyg = yg(iq,k) - jdel
        idel = min( max( idel - ioff, 0), ipan )
        jdel = min( max( jdel - joff, 0), jpan )
        n = min( max( nface(iq,k) + noff, 1), npan )
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
        cmin = min(sx(idel,  jdel,n,k),sx(idel+1,jdel,  n,k), &
                   sx(idel,jdel+1,n,k),sx(idel+1,jdel+1,n,k))
        cmax = max(sx(idel,  jdel,n,k),sx(idel+1,jdel,  n,k), &
                   sx(idel,jdel+1,n,k),sx(idel+1,jdel+1,n,k))
        rmul_1 = sx(idel-1,jdel,  n,k)*dmul_2 + sx(idel-1,jdel+1,n,k)*dmul_3
        rmul_2 = sx(idel,  jdel-1,n,k)*cmul_1 + sx(idel,  jdel,  n,k)*cmul_2 + &
                 sx(idel,  jdel+1,n,k)*cmul_3 + sx(idel,  jdel+2,n,k)*cmul_4
        rmul_3 = sx(idel+1,jdel-1,n,k)*cmul_1 + sx(idel+1,jdel,  n,k)*cmul_2 + &
                 sx(idel+1,jdel+1,n,k)*cmul_3 + sx(idel+1,jdel+2,n,k)*cmul_4
        rmul_4 = sx(idel+2,jdel,  n,k)*dmul_2 + sx(idel+2,jdel+1,n,k)*dmul_3
        s(iq,k) = min( max( cmin, &
            rmul_1*emul_1 + rmul_2*emul_2 + rmul_3*emul_3 + rmul_4*emul_4 ), cmax ) ! Bermejo & Staniforth
      end do       ! iq loop
    end do         ! k loop
    !$acc end parallel loop
    
  end if            ! (nfield<mh_bs)  .. else ..

end if               ! (intsch==1) .. else ..
!========================   end of intsch=1 section ====================

!$acc update self(s(1:ifull,1:kl))
!$acc exit data delete(sx,s)

call intssync_recv(s)
      
call END_LOG(ints_end)

return
end subroutine ints

subroutine ints_bl(s,intsch,nface,xg,yg)  ! not usually called

use cc_mpi             ! CC MPI routines
use indices_m          ! Grid index arrays
use newmpar_m          ! Grid parameters
use parm_m             ! Model configuration
use parmhor_m          ! Horizontal advection parameters

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
!$acc declare present(xg,yg,nface)

!     this one does bi-linear interpolation only
!     first extend s arrays into sx - this one -1:il+2 & -1:il+2
!                    but for bi-linear only need 0:il+1 &  0:il+1

call START_LOG(ints_begin)

call bounds(s,corner=.true.)

!$acc enter data create(sx,s)
!$acc update device(s)

!$acc parallel loop collapse(4)
do concurrent (k = 1:kl)
  do concurrent (n = 1:npan)
    do concurrent (j = 1:jpan)
      do concurrent (i = 1:ipan)
        iq = i + (j-1)*ipan + (n-1)*ipan*jpan
        sx(i,j,n,k) = s(iq,k)
      end do
    end do
  end do
end do
!$acc end parallel loop
!$acc parallel loop collapse(2)
do concurrent (k = 1:kl)
  do concurrent (n = 1:npan)
    do concurrent (j = 1:jpan)
      sx(0,j,n,k)      = s(iw(1+(j-1)*ipan+(n-1)*ipan*jpan),k)
      sx(ipan+1,j,n,k) = s(ie(j*ipan+(n-1)*ipan*jpan),k)
    end do               ! j loop
    do concurrent (i = 1:ipan)
      sx(i,0,n,k)      = s(is(i+(n-1)*ipan*jpan),k)
      sx(i,jpan+1,n,k) = s(in(i-ipan+n*ipan*jpan),k)
    end do               ! i loop
    sx(0,0,n,k)           = s(iws(1+(n-1)*ipan*jpan),k)
    sx(ipan+1,0,n,k)      = s(ies(ipan+(n-1)*ipan*jpan),k)
    sx(0,jpan+1,n,k)      = s(iwn(1-ipan+n*ipan*jpan),k)
    sx(ipan+1,jpan+1,n,k) = s(ien(n*ipan*jpan),k)
  end do                 ! n loop
end do                   ! k loop
!$acc end parallel loop

!$acc update self(sx)

! Loop over points that need to be calculated for other processes
do concurrent (ii = 1:neighnum)
  do concurrent (iq = 1:drlen(ii))
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

!$acc parallel loop collapse(2)
do concurrent (k = 1:kl)
  do concurrent (iq = 1:ifull)
    ! Convert face index from 0:npanels to array indices
    idel=int(xg(iq,k))
    xxg=xg(iq,k)-idel
    jdel=int(yg(iq,k))
    yyg=yg(iq,k)-jdel
    idel = min( max( idel - ioff, 0), ipan )
    jdel = min( max( jdel - joff, 0), jpan )
    n = min( max( nface(iq,k) + noff, 1), npan )
    s(iq,k) =      yyg*(xxg*sx(idel+1,jdel+1,n,k)+(1.-xxg)*sx(idel,jdel+1,n,k)) &
              + (1.-yyg)*(xxg*sx(idel+1,  jdel,n,k)+(1.-xxg)*sx(idel,  jdel,n,k))
  end do                  ! iq loop
end do                    ! k
!$acc end parallel loop

!$acc update self(s(1:ifull,1:kl))
!$acc exit data delete(s,sx)

call intssync_recv(s)

call END_LOG(ints_end)

return
end subroutine ints_bl
