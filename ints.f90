! Conformal Cubic Atmospheric Model
    
! Copyright 2015-2022 Commonwealth Scientific Industrial Research Organisation (CSIRO)
    
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
    
subroutine ints(s,ntr,intsch,nface,xg,yg,nfield)

use cc_acc             ! CC OpenACC routines
use cc_mpi             ! CC MPI routines
use indices_m          ! Grid index arrays
use newmpar_m          ! Grid parameters
use parm_m             ! Model configuration
use parmhor_m          ! Horizontal advection parameters

implicit none

integer, intent(in) :: intsch  ! method to interpolate panel corners
integer, intent(in) :: nfield  ! use B&S if nfield>=mh_bs
integer, intent(in) :: ntr     ! number of tracers to process
integer idel, iq, jdel
integer i, j, k, n, ii, nn, np
integer nstart, nend, nlen, async_counter
integer, dimension(ifull,kl), intent(in) :: nface        ! interpolation coordinates
real, dimension(ifull,kl), intent(in) :: xg, yg          ! interpolation coordinates
real, dimension(ifull+iextra,kl,ntr), intent(inout) :: s ! array of tracers
real, dimension(-1:ipan+2,-1:jpan+2,1:npan,kl,nagg) :: sx ! unpacked tracer array
real xxg, yyg, cmin, cmax
real dmul_2, dmul_3, cmul_1, cmul_2, cmul_3, cmul_4
real emul_1, emul_2, emul_3, emul_4, rmul_1, rmul_2, rmul_3, rmul_4

call START_LOG(ints_begin)

! now call bounds before calling ints

!======================== start of intsch=1 section ====================
if ( intsch==1 ) then
    
  do nstart = 1,ntr,nagg
    nend = min(nstart + nagg - 1, ntr )
    nlen = nend - nstart + 1
    
    ! this is intsb           EW interps done first
    ! first extend s arrays into sx - this one -1:il+2 & -1:il+2

    sx(1:ipan,1:jpan,1:npan,1:kl,1:nlen) = reshape(s(1:ipan*jpan*npan,1:kl,nstart:nend), (/ipan,jpan,npan,kl,nlen/))
    ! this is intsb           EW interps done first
    ! first extend s arrays into sx - this one -1:il+2 & -1:il+2
    do nn = 1,nlen
      np = nn - 1 + nstart  
      do k = 1,kl
        do n = 1,npan
          do j = 1,jpan
            iq = 1 + (j-1)*ipan + (n-1)*ipan*jpan
            sx(0,j,n,k,nn)      = s(iw(iq),k,np)
            sx(-1,j,n,k,nn)     = s(iww(iq),k,np)
            iq = j*ipan + (n-1)*ipan*jpan
            sx(ipan+1,j,n,k,nn) = s(ie(iq),k,np)
            sx(ipan+2,j,n,k,nn) = s(iee(iq),k,np)
          end do            ! j loop
          do i = 1,ipan
            iq = i + (n-1)*ipan*jpan
            sx(i,0,n,k,nn)      = s(is(iq),k,np)
            sx(i,-1,n,k,nn)     = s(iss(iq),k,np)
            iq = i - ipan + n*ipan*jpan
            sx(i,jpan+1,n,k,nn) = s(in(iq),k,np)
            sx(i,jpan+2,n,k,nn) = s(inn(iq),k,np)
          end do            ! i loop
        end do
        do n = 1,npan
          sx(-1,0,n,k,nn)          = s(lwws(n),k,np)
          sx(0,0,n,k,nn)           = s(iws(1+(n-1)*ipan*jpan),k,np)
          sx(0,-1,n,k,nn)          = s(lwss(n),k,np)
          sx(ipan+1,0,n,k,nn)      = s(ies(ipan+(n-1)*ipan*jpan),k,np)
          sx(ipan+2,0,n,k,nn)      = s(lees(n),k,np)
          sx(ipan+1,-1,n,k,nn)     = s(less(n),k,np)
          sx(-1,jpan+1,n,k,nn)     = s(lwwn(n),k,np)
          sx(0,jpan+2,n,k,nn)      = s(lwnn(n),k,np)
          sx(ipan+2,jpan+1,n,k,nn) = s(leen(n),k,np)
          sx(ipan+1,jpan+2,n,k,nn) = s(lenn(n),k,np)
          sx(0,jpan+1,n,k,nn)      = s(iwn(1-ipan+n*ipan*jpan),k,np)
          sx(ipan+1,jpan+1,n,k,nn) = s(ien(n*ipan*jpan),k,np)
        end do          ! n loop
      end do            ! k loop
    end do              ! nn loop  

    ! Loop over points that need to be calculated for other processes
    if ( nfield<mh_bs ) then

      do nn = 1,nlen 
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
            rmul_1 = sx(idel,  jdel-1,n,k,nn)*dmul_2 + sx(idel+1,jdel-1,n,k,nn)*dmul_3
            rmul_2 = sx(idel-1,jdel,  n,k,nn)*cmul_1 + sx(idel,  jdel,  n,k,nn)*cmul_2 + &
                     sx(idel+1,jdel,  n,k,nn)*cmul_3 + sx(idel+2,jdel,  n,k,nn)*cmul_4
            rmul_3 = sx(idel-1,jdel+1,n,k,nn)*cmul_1 + sx(idel,  jdel+1,n,k,nn)*cmul_2 + &
                     sx(idel+1,jdel+1,n,k,nn)*cmul_3 + sx(idel+2,jdel+1,n,k,nn)*cmul_4
            rmul_4 = sx(idel,  jdel+2,n,k,nn)*dmul_2 + sx(idel+1,jdel+2,n,k,nn)*dmul_3
            sextra(ii)%a(iq+(nn-1)*drlen(ii)) = rmul_1*emul_1 + rmul_2*emul_2 + rmul_3*emul_3 + rmul_4*emul_4
          end do      ! iq loop
        end do        ! ii loop
      end do          ! nn loop  

      ! Send messages to other processors.  We then start the calculation for this processor while waiting for
      ! the messages to return, thereby overlapping computation with communication.
      call intssync_send(nlen)

#ifndef GPU
      !$omp parallel do collapse(2) schedule(static) private(k,iq,idel,xxg,jdel,yyg),  &
      !$omp private(n,cmul_1,cmul_2,cmul_3,cmul_4,dmul_2,dmul_3,emul_1,emul_2,emul_3), &
      !$omp private(emul_4,rmul_1,rmul_2,rmul_3,rmul_4,nn,np)
#endif
      do nn = 1,nlen
#ifdef _OPENMP
#ifdef GPU
        !$omp target teams distribute parallel do collapse(2) schedule(static),                   &
        !$omp map(to:sx(:,:,:,:,nn)) map(from:s(:,:,nn-1+nstart)) private(k,iq,idel,xxg,jdel,yyg) &
        !$omp private(n,cmul_1,cmul_2,cmul_3,cmul_4,dmul_2,dmul_3,emul_1,emul_2,emul_3)           &
        !$omp private(emul_4,rmul_1,rmul_2,rmul_3,rmul_4)
#endif
#else
        async_counter = mod(nn-1, async_length)+1
        !$acc parallel loop collapse(2) copyin(sx(:,:,:,:,nn)) copyout(s(:,:,nn-1+nstart))        &
        !$acc   present(xg,yg,nface) async(async_counter)
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
            rmul_1 = sx(idel,  jdel-1,n,k,nn)*dmul_2 + sx(idel+1,jdel-1,n,k,nn)*dmul_3
            rmul_2 = sx(idel-1,jdel,  n,k,nn)*cmul_1 + sx(idel,  jdel,  n,k,nn)*cmul_2 + &
                     sx(idel+1,jdel,  n,k,nn)*cmul_3 + sx(idel+2,jdel,  n,k,nn)*cmul_4
            rmul_3 = sx(idel-1,jdel+1,n,k,nn)*cmul_1 + sx(idel,  jdel+1,n,k,nn)*cmul_2 + &
                     sx(idel+1,jdel+1,n,k,nn)*cmul_3 + sx(idel+2,jdel+1,n,k,nn)*cmul_4
            rmul_4 = sx(idel,  jdel+2,n,k,nn)*dmul_2 + sx(idel+1,jdel+2,n,k,nn)*dmul_3
            s(iq,k,nn-1+nstart) = rmul_1*emul_1 + rmul_2*emul_2 + rmul_3*emul_3 + rmul_4*emul_4
          end do       ! iq loop
        end do         ! k loop
#ifdef _OPENMP
#ifdef GPU
        !$omp end target teams distribute parallel do
#endif
#else
        !$acc end parallel loop
#endif
      end do           ! nn loop
#ifndef GPU
      !$omp end parallel do
#endif  
      !$acc wait
    
    else              ! (nfield<mh_bs)

      do nn = 1,nlen  
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
            cmin = min(sx(idel,  jdel,n,k,nn),sx(idel+1,jdel,  n,k,nn), &
                       sx(idel,jdel+1,n,k,nn),sx(idel+1,jdel+1,n,k,nn))
            cmax = max(sx(idel,  jdel,n,k,nn),sx(idel+1,jdel,  n,k,nn), &
                       sx(idel,jdel+1,n,k,nn),sx(idel+1,jdel+1,n,k,nn))
            rmul_1 = sx(idel,  jdel-1,n,k,nn)*dmul_2 + sx(idel+1,jdel-1,n,k,nn)*dmul_3
            rmul_2 = sx(idel-1,jdel,  n,k,nn)*cmul_1 + sx(idel,  jdel,  n,k,nn)*cmul_2 + &
                     sx(idel+1,jdel,  n,k,nn)*cmul_3 + sx(idel+2,jdel,  n,k,nn)*cmul_4
            rmul_3 = sx(idel-1,jdel+1,n,k,nn)*cmul_1 + sx(idel,  jdel+1,n,k,nn)*cmul_2 + &
                     sx(idel+1,jdel+1,n,k,nn)*cmul_3 + sx(idel+2,jdel+1,n,k,nn)*cmul_4
            rmul_4 = sx(idel,  jdel+2,n,k,nn)*dmul_2 + sx(idel+1,jdel+2,n,k,nn)*dmul_3
            sextra(ii)%a(iq+(nn-1)*drlen(ii)) = min( max( cmin, &
                rmul_1*emul_1 + rmul_2*emul_2 + rmul_3*emul_3 + rmul_4*emul_4 ), cmax ) ! Bermejo & Staniforth
          end do      ! iq loop
        end do        ! ii loop
      end do          ! nn loop
      ! Send messages to other processors.  We then start the calculation for this processor while waiting for
      ! the messages to return, thereby overlapping computation with communication.
      call intssync_send(nlen)

#ifndef GPU
      !$omp parallel do collapse(2) schedule(static) private(k,iq,idel,xxg,jdel,yyg),  &
      !$omp private(n,cmul_1,cmul_2,cmul_3,cmul_4,dmul_2,dmul_3,emul_1,emul_2,emul_3), &
      !$omp private(emul_4,rmul_1,rmul_2,rmul_3,rmul_4,nn,np)
#endif
      do nn = 1,nlen
#ifdef _OPENMP
#ifdef GPU
        !$omp target teams distribute parallel do collapse(2) schedule(static)                    &
        !$omp map(to:sx(:,:,:,:,nn)) map(from:s(:,:,nn-1+nstart)) private(k,iq,idel,xxg,jdel,yyg) &
        !$omp private(n,cmul_1,cmul_2,cmul_3,cmul_4,dmul_2,dmul_3,emul_1,emul_2,emul_3)           &
        !$omp private(emul_4,rmul_1,rmul_2,rmul_3,rmul_4)
#endif
#else
        async_counter = mod(nn-1, async_length)+1
        !$acc parallel loop collapse(2) copyin(sx(:,:,:,:,nn)) copyout(s(:,:,nn-1+nstart))        &
        !$acc   present(xg,yg,nface) async(async_counter)
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
            cmin = min(sx(idel,  jdel,n,k,nn),sx(idel+1,jdel,  n,k,nn), &
                       sx(idel,jdel+1,n,k,nn),sx(idel+1,jdel+1,n,k,nn))
            cmax = max(sx(idel,  jdel,n,k,nn),sx(idel+1,jdel,  n,k,nn), &
                       sx(idel,jdel+1,n,k,nn),sx(idel+1,jdel+1,n,k,nn))
            rmul_1 = sx(idel,  jdel-1,n,k,nn)*dmul_2 + sx(idel+1,jdel-1,n,k,nn)*dmul_3
            rmul_2 = sx(idel-1,jdel,  n,k,nn)*cmul_1 + sx(idel,  jdel,  n,k,nn)*cmul_2 + &
                     sx(idel+1,jdel,  n,k,nn)*cmul_3 + sx(idel+2,jdel,  n,k,nn)*cmul_4
            rmul_3 = sx(idel-1,jdel+1,n,k,nn)*cmul_1 + sx(idel,  jdel+1,n,k,nn)*cmul_2 + &
                     sx(idel+1,jdel+1,n,k,nn)*cmul_3 + sx(idel+2,jdel+1,n,k,nn)*cmul_4
            rmul_4 = sx(idel,  jdel+2,n,k,nn)*dmul_2 + sx(idel+1,jdel+2,n,k,nn)*dmul_3
            s(iq,k,nn-1+nstart) = min( max( cmin, &
                rmul_1*emul_1 + rmul_2*emul_2 + rmul_3*emul_3 + rmul_4*emul_4 ), cmax ) ! Bermejo & Staniforth
          end do      ! iq loop
        end do        ! k loop
#ifdef _OPENMP
#ifdef GPU
        !$omp end target teams distribute parallel do
#endif
#else
        !$acc end parallel loop
#endif
      end do          ! nn loop  
#ifndef GPU
      !$omp end parallel do
#endif
      !$acc wait

    end if            ! (nfield<mh_bs)  .. else ..
  
    call intssync_recv(s(:,:,nstart:nend))  
    
  end do ! ntr
  
!========================   end of intsch=1 section ====================
else     ! if(intsch==1)then
!======================== start of intsch=2 section ====================
!       this is intsc           NS interps done first
!       first extend s arrays into sx - this one -1:il+2 & -1:il+2

  do nstart = 1,ntr,nagg
    nend = min(nstart + nagg - 1, ntr )
    nlen = nend - nstart + 1    
    
    sx(1:ipan,1:jpan,1:npan,1:kl,1:nlen) = reshape(s(1:ipan*jpan*npan,1:kl,nstart:nend), (/ipan,jpan,npan,kl,nlen/))
    do nn = 1,nlen
      np = nn - 1 + nstart  
      do k = 1,kl
        do n = 1,npan
          do j = 1,jpan
            iq = 1+(j-1)*ipan+(n-1)*ipan*jpan
            sx(0,j,n,k,nn)      = s(iw(iq),k,np)
            sx(-1,j,n,k,nn)     = s(iww(iq),k,np)
            iq = j*ipan+(n-1)*ipan*jpan
            sx(ipan+1,j,n,k,nn) = s(ie(iq),k,np)
            sx(ipan+2,j,n,k,nn) = s(iee(iq),k,np)
          end do            ! j loop
          do i = 1,ipan
            iq = i+(n-1)*ipan*jpan
            sx(i,0,n,k,nn)      = s(is(iq),k,np)
            sx(i,-1,n,k,nn)     = s(iss(iq),k,np)
            iq = i-ipan+n*ipan*jpan
            sx(i,jpan+1,n,k,nn) = s(in(iq),k,np)
            sx(i,jpan+2,n,k,nn) = s(inn(iq),k,np)
          end do            ! i loop
        end do
        do n = 1,npan
          sx(-1,0,n,k,nn)          = s(lsww(n),k,np)
          sx(0,0,n,k,nn)           = s(isw(1+(n-1)*ipan*jpan),k,np)
          sx(0,-1,n,k,nn)          = s(lssw(n),k,np)
          sx(ipan+2,0,n,k,nn)      = s(lsee(n),k,np)
          sx(ipan+1,-1,n,k,nn)     = s(lsse(n),k,np)
          sx(-1,jpan+1,n,k,nn)     = s(lnww(n),k,np)
          sx(0,jpan+1,n,k,nn)      = s(inw(1-ipan+n*ipan*jpan),k,np)
          sx(0,jpan+2,n,k,nn)      = s(lnnw(n),k,np)
          sx(ipan+2,jpan+1,n,k,nn) = s(lnee(n),k,np)
          sx(ipan+1,jpan+2,n,k,nn) = s(lnne(n),k,np)
          sx(ipan+1,0,n,k,nn)      = s(ise(ipan+(n-1)*ipan*jpan),k,np)
          sx(ipan+1,jpan+1,n,k,nn) = s(ine(n*ipan*jpan),k,np)
        end do              ! n loop
      end do                ! k loop
    end do                  ! nn loop

    ! For other processes
    if ( nfield < mh_bs ) then

      do nn = 1,nlen  
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
            rmul_1 = sx(idel-1,jdel,  n,k,nn)*dmul_2 + sx(idel-1,jdel+1,n,k,nn)*dmul_3
            rmul_2 = sx(idel,  jdel-1,n,k,nn)*cmul_1 + sx(idel,  jdel,  n,k,nn)*cmul_2 + &
                     sx(idel,  jdel+1,n,k,nn)*cmul_3 + sx(idel,  jdel+2,n,k,nn)*cmul_4
            rmul_3 = sx(idel+1,jdel-1,n,k,nn)*cmul_1 + sx(idel+1,jdel,  n,k,nn)*cmul_2 + &
                     sx(idel+1,jdel+1,n,k,nn)*cmul_3 + sx(idel+1,jdel+2,n,k,nn)*cmul_4
            rmul_4 = sx(idel+2,jdel,  n,k,nn)*dmul_2 + sx(idel+2,jdel+1,n,k,nn)*dmul_3
            sextra(ii)%a(iq+(nn-1)*drlen(ii)) = rmul_1*emul_1 + rmul_2*emul_2 + rmul_3*emul_3 + rmul_4*emul_4
          end do         ! iq loop
        end do           ! ii loop
      end do             ! nn loop  
    
      call intssync_send(nlen)

#ifndef GPU
      !$omp parallel do collapse(2) schedule(static) private(k,iq,idel,xxg,jdel,yyg),  &
      !$omp private(n,cmul_1,cmul_2,cmul_3,cmul_4,dmul_2,dmul_3,emul_1,emul_2,emul_3), &
      !$omp private(emul_4,rmul_1,rmul_2,rmul_3,rmul_4,nn,np)
#endif
      do nn = 1,nlen
#ifdef _OPENMP
#ifdef GPU
        !$omp target teams distribute parallel do collapse(2) schedule(static)                    &
        !$omp map(to:sx(:,:,:,:,nn)) map(from:s(:,:,nn-1+nstart)) private(k,iq,idel,xxg,jdel,yyg) &
        !$omp private(n,cmul_1,cmul_2,cmul_3,cmul_4,dmul_2,dmul_3,emul_1,emul_2,emul_3)           &
        !$omp private(emul_4,rmul_1,rmul_2,rmul_3,rmul_4)
#endif
#else
        async_counter = mod(nn-1, async_length)+1
        !$acc parallel loop collapse(2) copyin(sx(:,:,:,:,nn)) copyout(s(:,:,nn-1+nstart))        &
        !$acc   present(xg,yg,nface) async(async_counter)
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
            rmul_1 = sx(idel-1,jdel,  n,k,nn)*dmul_2 + sx(idel-1,jdel+1,n,k,nn)*dmul_3
            rmul_2 = sx(idel,  jdel-1,n,k,nn)*cmul_1 + sx(idel,  jdel,  n,k,nn)*cmul_2 + &
                     sx(idel,  jdel+1,n,k,nn)*cmul_3 + sx(idel,  jdel+2,n,k,nn)*cmul_4
            rmul_3 = sx(idel+1,jdel-1,n,k,nn)*cmul_1 + sx(idel+1,jdel,  n,k,nn)*cmul_2 + &
                     sx(idel+1,jdel+1,n,k,nn)*cmul_3 + sx(idel+1,jdel+2,n,k,nn)*cmul_4
            rmul_4 = sx(idel+2,jdel,  n,k,nn)*dmul_2 + sx(idel+2,jdel+1,n,k,nn)*dmul_3
            s(iq,k,nn-1+nstart) = rmul_1*emul_1 + rmul_2*emul_2 + rmul_3*emul_3 + rmul_4*emul_4
          end do       ! iq loop
        end do         ! k loop
#ifdef _OPENMP
#ifdef GPU
        !$omp end target teams distribute parallel do
#endif
#else
        !$acc end parallel loop
#endif
      end do           ! nn loop  
#ifndef GPU
      !$omp end parallel do
#endif
      !$acc wait
    
    else                 ! (nfield<mh_bs)

      do nn = 1,nlen  
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
            cmin = min(sx(idel,  jdel,n,k,nn),sx(idel+1,jdel,  n,k,nn), &
                       sx(idel,jdel+1,n,k,nn),sx(idel+1,jdel+1,n,k,nn))
            cmax = max(sx(idel,  jdel,n,k,nn),sx(idel+1,jdel,  n,k,nn), &
                       sx(idel,jdel+1,n,k,nn),sx(idel+1,jdel+1,n,k,nn))
            rmul_1 = sx(idel-1,jdel,  n,k,nn)*dmul_2 + sx(idel-1,jdel+1,n,k,nn)*dmul_3
            rmul_2 = sx(idel,  jdel-1,n,k,nn)*cmul_1 + sx(idel,  jdel,  n,k,nn)*cmul_2 + &
                     sx(idel,  jdel+1,n,k,nn)*cmul_3 + sx(idel,  jdel+2,n,k,nn)*cmul_4
            rmul_3 = sx(idel+1,jdel-1,n,k,nn)*cmul_1 + sx(idel+1,jdel,  n,k,nn)*cmul_2 + &
                     sx(idel+1,jdel+1,n,k,nn)*cmul_3 + sx(idel+1,jdel+2,n,k,nn)*cmul_4
            rmul_4 = sx(idel+2,jdel,  n,k,nn)*dmul_2 + sx(idel+2,jdel+1,n,k,nn)*dmul_3
            sextra(ii)%a(iq+(nn-1)*drlen(ii)) = min( max( cmin, &
                rmul_1*emul_1 + rmul_2*emul_2 + rmul_3*emul_3 + rmul_4*emul_4 ), cmax ) ! Bermejo & Staniforth
          end do      ! iq loop
        end do        ! ii loop
      end do          ! nn loop  
  
      call intssync_send(nlen)

#ifndef GPU
      !$omp parallel do collapse(2) schedule(static) private(k,iq,idel,xxg,jdel,yyg),  &
      !$omp private(n,cmul_1,cmul_2,cmul_3,cmul_4,dmul_2,dmul_3,emul_1,emul_2,emul_3), &
      !$omp private(emul_4,rmul_1,rmul_2,rmul_3,rmul_4,nn,np)
#endif
      do nn = 1,nlen
#ifdef _OPENMP
#ifdef GPU
        !$omp target teams distribute parallel do collapse(2) schedule(static)                    &
        !$omp map(to:sx(:,:,:,:,nn)) map(from:s(:,:,nn-1+nstart)) private(k,iq,idel,xxg,jdel,yyg) &
        !$omp private(n,cmul_1,cmul_2,cmul_3,cmul_4,dmul_2,dmul_3,emul_1,emul_2,emul_3)           &
        !$omp private(emul_4,rmul_1,rmul_2,rmul_3,rmul_4)
#endif
#else
        async_counter = mod(nn-1, async_length)+1
        !$acc parallel loop collapse(2) copyin(sx(:,:,:,:,nn)) copyout(s(:,:,nn-1+nstart))        &
        !$acc   present(xg,yg,nface) async(async_counter)
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
            cmin = min(sx(idel,  jdel,n,k,nn),sx(idel+1,jdel,  n,k,nn), &
                       sx(idel,jdel+1,n,k,nn),sx(idel+1,jdel+1,n,k,nn))
            cmax = max(sx(idel,  jdel,n,k,nn),sx(idel+1,jdel,  n,k,nn), &
                       sx(idel,jdel+1,n,k,nn),sx(idel+1,jdel+1,n,k,nn))
            rmul_1 = sx(idel-1,jdel,  n,k,nn)*dmul_2 + sx(idel-1,jdel+1,n,k,nn)*dmul_3
            rmul_2 = sx(idel,  jdel-1,n,k,nn)*cmul_1 + sx(idel,  jdel,  n,k,nn)*cmul_2 + &
                     sx(idel,  jdel+1,n,k,nn)*cmul_3 + sx(idel,  jdel+2,n,k,nn)*cmul_4
            rmul_3 = sx(idel+1,jdel-1,n,k,nn)*cmul_1 + sx(idel+1,jdel,  n,k,nn)*cmul_2 + &
                     sx(idel+1,jdel+1,n,k,nn)*cmul_3 + sx(idel+1,jdel+2,n,k,nn)*cmul_4
            rmul_4 = sx(idel+2,jdel,  n,k,nn)*dmul_2 + sx(idel+2,jdel+1,n,k,nn)*dmul_3
            s(iq,k,nn-1+nstart) = min( max( cmin, &
                rmul_1*emul_1 + rmul_2*emul_2 + rmul_3*emul_3 + rmul_4*emul_4 ), cmax ) ! Bermejo & Staniforth
          end do       ! iq loop
        end do         ! k loop
#ifdef _OPENMP
#ifdef GPU
        !$omp end target teams distribute parallel do
#endif
#else
        !$acc end parallel loop
#endif
      end do           ! nn loop  
#ifndef GPU
      !$omp end parallel do
#endif
      !$acc wait
    
    end if            ! (nfield<mh_bs)  .. else ..

    call intssync_recv(s(:,:,nstart:nend))  
    
  end do ! ntr
  
end if               ! (intsch==1) .. else ..
!========================   end of intsch=1 section ====================

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

!     this one does bi-linear interpolation only
!     first extend s arrays into sx - this one -1:il+2 & -1:il+2
!                    but for bi-linear only need 0:il+1 &  0:il+1

call START_LOG(ints_begin)

! now call bounds before calling ints_bl
!call bounds(s,corner=.true.)

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

!$omp parallel do schedule(static) private(k,iq,idel,xxg,jdel,yyg,n)
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
!$omp end parallel do

call intssync_recv(s)

call END_LOG(ints_end)

return
end subroutine ints_bl
