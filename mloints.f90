! Conformal Cubic Atmospheric Model
    
! Copyright 2015-2023 Commonwealth Scientific Industrial Research Organisation (CSIRO)
    
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
public mlob2ints_bs

interface mlob2ints_bs
  module procedure mlob2ints_bs_2, mlob2ints_bs_3
end interface

contains


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Calculate depature points for MLO semi-Lagrangian advection

subroutine mlob2ints_bs_3(s,nface,xg,yg,wtr,sal_test,bs_test)

use cc_acc
use cc_mpi
use indices_m
use mlo
use newmpar_m
use parm_m
use parmhor_m

implicit none

integer idel, iq, jdel
integer i, j, k, n, intsch, nn, np
integer ii, ntr, nstart, nend, nlen, async_counter
integer, dimension(ifull,wlev), intent(in) :: nface
real, dimension(ifull,wlev), intent(in) :: xg, yg
real, dimension(:,:,:), intent(inout) :: s
real, dimension(-1:ipan+2,-1:jpan+2,1:npan,wlev,size(s,3)) :: sx
real, dimension(ifull+iextra,wlev,size(s,3)) :: s_old
real, dimension(ifull,wlev,size(s,3)) :: s_store
real s_tot, s_count
real cmax, cmin, xxg, yyg
real dmul_2, dmul_3, cmul_1, cmul_2, cmul_3, cmul_4
real emul_1, emul_2, emul_3, emul_4, rmul_1, rmul_2, rmul_3, rmul_4
real, parameter :: cxx = -9999. ! missing value flag
logical, dimension(ifull+iextra,wlev), intent(in) :: wtr
logical, dimension(size(s,3)), intent(in) :: sal_test, bs_test

call START_LOG(waterints_begin)

ntr = size(s,3)
intsch = mod(ktau, 2)

do nn = 1,ntr
  do k = 1,wlev
    s_store(1:ifull,k,nn) = s(1:ifull,k,nn)
    where (.not.wtr(1:ifull,k))
      s(1:ifull,k,nn) = cxx - 1. ! missing value flag
    end where
  end do
  if ( sal_test(nn) ) then
    do k = 1,wlev  
      where ( s(1:ifull,k,nn)<2. )
        s(1:ifull,k,nn) = cxx - 1. ! missing value flag
      end where
    end do
  end if  
end do

! fill
do ii = 1,6 ! 6 iterations of fill should be enough
  s_old(1:ifull,:,:) = s(1:ifull,:,:)
  call bounds(s_old)
  !$omp parallel do collapse(2) private(nn,k,iq,s_tot,s_count)
  do nn = 1,ntr
    do k = 1,wlev
      do iq = 1,ifull
        if ( s(iq,k,nn)<cxx ) then
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
            s(iq,k,nn) = s_tot/s_count
          end if  
        end if
      end do  
    end do
  end do
  !$omp end parallel do
end do

call bounds(s,nrows=2)

!======================== start of intsch=1 section ====================
if ( intsch==1 ) then
    
  do nstart = 1,ntr,nagg
    nend = min(nstart + nagg - 1, ntr)
    nlen = nend - nstart + 1

    do nn = 1,nlen
      np = nn - 1 + nstart  
      do k = 1,wlev
        sx(1:ipan,1:jpan,1:npan,k,nn) = &
          reshape( s(1:ipan*jpan*npan,k,np), (/ ipan, jpan, npan /) )
        do n = 1,npan
          do j = 1,jpan
            iq = 1+(j-1)*ipan+(n-1)*ipan*jpan
            sx(0,j,n,k,nn)      = s( iw(iq),k,np)
            sx(-1,j,n,k,nn)     = s(iww(iq),k,np)
            iq = j*ipan+(n-1)*ipan*jpan
            sx(ipan+1,j,n,k,nn) = s( ie(iq),k,np)
            sx(ipan+2,j,n,k,nn) = s(iee(iq),k,np)
          end do            ! j loop
          do i = 1,ipan
            iq = i+(n-1)*ipan*jpan
            sx(i,0,n,k,nn)      = s( is(iq),k,np)
            sx(i,-1,n,k,nn)     = s(iss(iq),k,np)
            iq = i-ipan+n*ipan*jpan
            sx(i,jpan+1,n,k,nn) = s( in(iq),k,np)
            sx(i,jpan+2,n,k,nn) = s(inn(iq),k,np)
          end do            ! i loop
        end do
!   for ew interpolation, sometimes need (different from ns):
!       (-1,0),   (0,0),   (0,-1)   (-1,il+1),   (0,il+1),   (0,il+2)
!     (il+1,0),(il+2,0),(il+1,-1) (il+1,il+1),(il+2,il+1),(il+1,il+2)
        do n = 1,npan
          sx(-1,0,n,k,nn)          = s(lwws(n),                  k,np)
          sx(0,0,n,k,nn)           = s(iws(1+(n-1)*ipan*jpan),   k,np)
          sx(0,-1,n,k,nn)          = s(lwss(n),                  k,np)
          sx(ipan+1,0,n,k,nn)      = s(ies(ipan+(n-1)*ipan*jpan),k,np)
          sx(ipan+2,0,n,k,nn)      = s(lees(n),                  k,np)
          sx(ipan+1,-1,n,k,nn)     = s(less(n),                  k,np)
          sx(-1,jpan+1,n,k,nn)     = s(lwwn(n),                  k,np)
          sx(0,jpan+2,n,k,nn)      = s(lwnn(n),                  k,np)
          sx(ipan+2,jpan+1,n,k,nn) = s(leen(n),                  k,np)
          sx(ipan+1,jpan+2,n,k,nn) = s(lenn(n),                  k,np)
          sx(0,jpan+1,n,k,nn)      = s(iwn(1-ipan+n*ipan*jpan),  k,np)
          sx(ipan+1,jpan+1,n,k,nn) = s(ien(n*ipan*jpan),         k,np)
        end do           ! n loop
      end do             ! k loop
    end do               ! nn loop  

    ! Loop over points that need to be calculated for other processes
    do nn = 1,nlen
      if ( bs_test(nn-1+nstart) ) then
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
            ! bi-cubic interpolation
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
            ! bi-cubic interpolation
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
            sextra(ii)%a(iq+(nn-1)*drlen(ii)) = &
                rmul_1*emul_1 + rmul_2*emul_2 + rmul_3*emul_3 + rmul_4*emul_4
          end do          ! iq loop
        end do            ! ii loop
      end if              ! bs_test ..else.. 
    end do                ! nn loop  

    call intssync_send(nlen)

#ifndef GPU
    !$omp parallel do schedule(static) private(nn,async_counter,k,iq,idel,jdel,n,xxg,yyg)  &
    !$omp   private(cmul_1,cmul_2,cmul_3,cmul_4,dmul_2,dmul_3,emul_1,emul_2,emul_3,emul_4) &
    !$omp   private(rmul_1,rmul_2,rmul_3,rmul_4,cmin,cmax)
#endif    
    do nn = 1,nlen
      async_counter = mod(nn-1, async_length)
      if ( bs_test(nn-1+nstart) ) then
        !$acc parallel loop collapse(2) copyin(sx(:,:,:,:,nn)) copyout(s(:,:,nn-1+nstart))        &
        !$acc   present(xg,yg,nface) async(async_counter)
        do k = 1,wlev      
          do iq = 1,ifull
            idel = int(xg(iq,k))
            xxg = xg(iq,k) - real(idel)
            jdel = int(yg(iq,k))
            yyg = yg(iq,k) - real(jdel)
            idel = min( max(idel - ioff, 0), ipan)
            jdel = min( max(jdel - joff, 0), jpan)
            n = min( max(nface(iq,k) + noff, 1), npan)
            ! bi-cubic interpolation
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
          end do       ! iq loop
        end do         ! k loop
        !$acc end parallel loop
      else
        !$acc parallel loop collapse(2) copyin(sx(:,:,:,:,nn)) copyout(s(:,:,nn-1+nstart))        &
        !$acc   present(xg,yg,nface) async(async_counter)
        do k = 1,wlev      
          do iq = 1,ifull
            idel = int(xg(iq,k))
            xxg = xg(iq,k) - real(idel)
            jdel = int(yg(iq,k))
            yyg = yg(iq,k) - real(jdel)
            idel = min( max(idel - ioff, 0), ipan)
            jdel = min( max(jdel - joff, 0), jpan)
            n = min( max(nface(iq,k) + noff, 1), npan)
            ! bi-cubic interpolation
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
            s(iq,k,nn-1+nstart) = &
                rmul_1*emul_1 + rmul_2*emul_2 + rmul_3*emul_3 + rmul_4*emul_4
          end do       ! iq loop
        end do         ! k loop
        !$acc end parallel loop
      end if           ! bs_test ..else..
    end do             ! nn loop
#ifndef GPU
    !$omp end parallel do
#endif
    !$acc wait

    call intssync_recv(s(:,:,nstart:nend))  
    
  end do ! ntr

!========================   end of intsch=1 section ====================
else     ! if(intsch==1)then
!======================== start of intsch=2 section ====================
!       this is intsc           NS interps done first
!       first extend s arrays into sx - this one -1:il+2 & -1:il+2

  do nstart = 1,ntr,nagg
    nend = min(nstart + nagg - 1, ntr)
    nlen = nend - nstart + 1
    
    do nn = 1,nlen
      np = nn - 1 + nstart      
      do k = 1,wlev
        sx(1:ipan,1:jpan,1:npan,k,nn) = &
          reshape( s(1:ipan*jpan*npan,k,np), (/ ipan, jpan, npan /) )
        do n = 1,npan
          do j = 1,jpan
            iq = 1+(j-1)*ipan+(n-1)*ipan*jpan
            sx(0,j,n,k,nn)      = s( iw(iq),k,np)
            sx(-1,j,n,k,nn)     = s(iww(iq),k,np)
            iq = j*ipan+(n-1)*ipan*jpan
            sx(ipan+1,j,n,k,nn) = s( ie(iq),k,np)
            sx(ipan+2,j,n,k,nn) = s(iee(iq),k,np)
          end do            ! j loop
          do i = 1,ipan
            iq = i+(n-1)*ipan*jpan
            sx(i,0,n,k,nn)      = s( is(iq),k,np)
            sx(i,-1,n,k,nn)     = s(iss(iq),k,np)
            iq = i-ipan+n*ipan*jpan
            sx(i,jpan+1,n,k,nn) = s( in(iq),k,np)
            sx(i,jpan+2,n,k,nn) = s(inn(iq),k,np)
          end do            ! i loop
        end do
!   for ns interpolation, sometimes need (different from ew):
!        (-1,0),   (0,0),   (0,-1)   (-1,il+1),   (0,il+1),   (0,il+2)
!      (il+1,0),(il+2,0),(il+1,-1) (il+1,il+1),(il+2,il+1),(il+1,il+2)
        do n = 1,npan
          sx(-1,0,n,k,nn)          = s(lsww(n),k,np)
          sx(0,0,n,k,nn)           = s(isw(1+(n-1)*ipan*jpan),   k,np)
          sx(0,-1,n,k,nn)          = s(lssw(n),k,np)
          sx(ipan+2,0,n,k,nn)      = s(lsee(n),k,np)
          sx(ipan+1,-1,n,k,nn)     = s(lsse(n),k,np)
          sx(-1,jpan+1,n,k,nn)     = s(lnww(n),k,np)
          sx(0,jpan+1,n,k,nn)      = s(inw(1-ipan+n*ipan*jpan),  k,np)
          sx(0,jpan+2,n,k,nn)      = s(lnnw(n),k,np)
          sx(ipan+2,jpan+1,n,k,nn) = s(lnee(n),k,np)
          sx(ipan+1,jpan+2,n,k,nn) = s(lnne(n),k,np)
          sx(ipan+1,0,n,k,nn)      = s(ise(ipan+(n-1)*ipan*jpan),k,np)
          sx(ipan+1,jpan+1,n,k,nn) = s(ine(n*ipan*jpan),         k,np)
        end do           ! n loop
      end do             ! k loop
    end do               ! nn loop  

    ! For other processes
    do nn = 1,nlen
      if ( bs_test(nn-1+nstart) ) then
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
            ! bi-cubic interpolation
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
            ! bi-cubic interpolation
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
            sextra(ii)%a(iq+(nn-1)*drlen(ii)) = &
                rmul_1*emul_1 + rmul_2*emul_2 + rmul_3*emul_3 + rmul_4*emul_4
          end do          ! iq loop
        end do            ! ii loop
      end if              ! bs_test ..else.. 
    end do                ! nn loop  

    call intssync_send(nlen)

#ifndef GPU
    !$omp parallel do schedule(static) private(nn,async_counter,k,iq,idel,jdel,n,xxg,yyg)  &
    !$omp   private(cmul_1,cmul_2,cmul_3,cmul_4,dmul_2,dmul_3,emul_1,emul_2,emul_3,emul_4) &
    !$omp   private(rmul_1,rmul_2,rmul_3,rmul_4,cmin,cmax)
#endif  
    do nn = 1,nlen
      async_counter = mod(nn-1, async_length)
      if ( bs_test(nn-1+nstart) ) then
        !$acc parallel loop collapse(2) copyin(sx(:,:,:,:,nn)) copyout(s(:,:,nn-1+nstart))        &
        !$acc   present(xg,yg,nface) async(async_counter)
        do k = 1,wlev
          do iq = 1,ifull
            idel = int(xg(iq,k))
            xxg = xg(iq,k) - real(idel)
            jdel = int(yg(iq,k))
            yyg = yg(iq,k) - real(jdel)
            idel = min( max(idel - ioff, 0), ipan)
            jdel = min( max(jdel - joff, 0), jpan)
            n = min( max(nface(iq,k) + noff, 1), npan)
            ! bi-cubic interpolation
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
          end do
        end do
        !$acc end parallel loop
      else
        !$acc parallel loop collapse(2) copyin(sx(:,:,:,:,nn)) copyout(s(:,:,nn-1+nstart))        &
        !$acc   present(xg,yg,nface) async(async_counter)
        do k = 1,wlev
          do iq = 1,ifull
            idel = int(xg(iq,k))
            xxg = xg(iq,k) - real(idel)
            jdel = int(yg(iq,k))
            yyg = yg(iq,k) - real(jdel)
            idel = min( max(idel - ioff, 0), ipan)
            jdel = min( max(jdel - joff, 0), jpan)
            n = min( max(nface(iq,k) + noff, 1), npan)
            ! bi-cubic interpolation
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
            s(iq,k,nn-1+nstart) = &
                rmul_1*emul_1 + rmul_2*emul_2 + rmul_3*emul_3 + rmul_4*emul_4
          end do
        end do
        !$acc end parallel loop
      end if
    end do           ! nn loop
#ifndef GPU
    !$omp end parallel do
#endif
    !$acc wait

    call intssync_recv(s(:,:,nstart:nend))  

  end do ! ntr    
    
end if                     ! (intsch==1) .. else ..
!========================   end of intsch=1 section ====================

do nn = 1,ntr
  do k = 1,wlev
    where ( .not.wtr(1:ifull,k) .or. s(1:ifull,k,nn)<cxx+10. )
      s(1:ifull,k,nn) = s_store(1:ifull,k,nn)
    end where
  end do
  if ( sal_test(nn) ) then
    do k = 1,wlev
      where ( s_store(1:ifull,k,nn)<2. )
        s(1:ifull,k,nn) = s_store(1:ifull,k,nn)  
      end where
    end do
  end if
end do

call END_LOG(waterints_end)

return
end subroutine mlob2ints_bs_3

subroutine mlob2ints_bs_2(s,nface,xg,yg,wtr,sal_test,bs_test)

use cc_acc
use cc_mpi
use indices_m
use mlo
use newmpar_m
use parm_m
use parmhor_m

implicit none

integer idel, iq, jdel
integer i, j, k, n, intsch
integer ii
integer, dimension(ifull,wlev), intent(in) :: nface
real, dimension(ifull,wlev), intent(in) :: xg, yg
real, dimension(:,:), intent(inout) :: s
real, dimension(-1:ipan+2,-1:jpan+2,1:npan,wlev) :: sx
real, dimension(ifull+iextra,wlev) :: s_old
real, dimension(ifull,wlev) :: s_store
real s_tot, s_count
real cmax, cmin, xxg, yyg
real dmul_2, dmul_3, cmul_1, cmul_2, cmul_3, cmul_4
real emul_1, emul_2, emul_3, emul_4, rmul_1, rmul_2, rmul_3, rmul_4
real, parameter :: cxx = -9999. ! missing value flag
logical, dimension(ifull+iextra,wlev), intent(in) :: wtr
logical, intent(in) :: sal_test, bs_test

call START_LOG(waterints_begin)

intsch = mod(ktau,2)

do k = 1,wlev
  s_store(1:ifull,k) = s(1:ifull,k)
  where (.not.wtr(1:ifull,k))
    s(1:ifull,k) = cxx - 1. ! missing value flag
  end where
end do
if ( sal_test ) then
  do k = 1,wlev  
    where ( s(1:ifull,k)<2. )
      s(1:ifull,k) = cxx - 1. ! missing value flag
    end where
  end do
end if

! fill
do ii = 1,6 ! 6 iterations of fill should be enough
  s_old(1:ifull,:) = s(1:ifull,:)
  call bounds(s_old)
  !$omp parallel do schedule(static) private(k,iq,s_tot,s_count)
  do k = 1,wlev
    do iq = 1,ifull
      if ( s(iq,k)<cxx ) then
        s_tot = 0.
        s_count = 0.
        if ( s_old(is(iq),k)>=cxx ) then
          s_tot = s_tot + s_old(is(iq),k)
          s_count = s_count + 1.
        end if
        if ( s_old(in(iq),k)>=cxx ) then
          s_tot = s_tot + s_old(in(iq),k)
          s_count = s_count + 1.
        end if
        if ( s_old(iw(iq),k)>=cxx ) then
          s_tot = s_tot + s_old(iw(iq),k)
          s_count = s_count + 1.
        end if
        if ( s_old(ie(iq),k)>=cxx ) then
          s_tot = s_tot + s_old(ie(iq),k)
          s_count = s_count + 1.
        end if
        if ( s_count>0. ) then
          s(iq,k) = s_tot/s_count
        end if  
      end if
    end do  
  end do
  !$omp end parallel do
end do

call bounds(s,nrows=2)

!======================== start of intsch=1 section ====================
if ( intsch==1 ) then
    
  do k = 1,wlev
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
        ! bi-cubic interpolation
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
        ! bi-cubic interpolation
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
        sextra(ii)%a(iq) = &
            rmul_1*emul_1 + rmul_2*emul_2 + rmul_3*emul_3 + rmul_4*emul_4
      end do          ! iq loop
    end do            ! ii loop
  end if

  call intssync_send(1)

  if ( bs_test ) then
    !$acc parallel loop collapse(2) copyin(sx) copyout(s)        &
    !$acc   present(xg,yg,nface)
    do k = 1,wlev      
      do iq = 1,ifull
        idel = int(xg(iq,k))
        xxg = xg(iq,k) - real(idel)
        jdel = int(yg(iq,k))
        yyg = yg(iq,k) - real(jdel)
        idel = min( max(idel - ioff, 0), ipan)
        jdel = min( max(jdel - joff, 0), jpan)
        n = min( max(nface(iq,k) + noff, 1), npan)
        ! bi-cubic interpolation
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
      end do       ! iq loop
    end do         ! k loop
    !$acc end parallel loop
  else
    !$acc parallel loop collapse(2) copyin(sx) copyout(s)        &
    !$acc   present(xg,yg,nface)
    do k = 1,wlev      
      do iq = 1,ifull
        idel = int(xg(iq,k))
        xxg = xg(iq,k) - real(idel)
        jdel = int(yg(iq,k))
        yyg = yg(iq,k) - real(jdel)
        idel = min( max(idel - ioff, 0), ipan)
        jdel = min( max(jdel - joff, 0), jpan)
        n = min( max(nface(iq,k) + noff, 1), npan)
        ! bi-cubic interpolation
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
        s(iq,k) = &
            rmul_1*emul_1 + rmul_2*emul_2 + rmul_3*emul_3 + rmul_4*emul_4
      end do       ! iq loop
    end do         ! k loop
    !$acc end parallel loop
  end if

  call intssync_recv(s)  

!========================   end of intsch=1 section ====================
else     ! if(intsch==1)then
!======================== start of intsch=2 section ====================
!       this is intsc           NS interps done first
!       first extend s arrays into sx - this one -1:il+2 & -1:il+2
    
  do k = 1,wlev
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

  ! For other processes
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
        ! bi-cubic interpolation
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
        ! bi-cubic interpolation
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
        sextra(ii)%a(iq) = &
            rmul_1*emul_1 + rmul_2*emul_2 + rmul_3*emul_3 + rmul_4*emul_4
      end do          ! iq loop
    end do            ! ii loop
  end if

  call intssync_send(1)

  if ( bs_test ) then
    !$acc parallel loop collapse(2) copyin(sx) copyout(s)        &
    !$acc   present(xg,yg,nface)
    do k = 1,wlev
      do iq = 1,ifull
        idel = int(xg(iq,k))
        xxg = xg(iq,k) - real(idel)
        jdel = int(yg(iq,k))
        yyg = yg(iq,k) - real(jdel)
        idel = min( max(idel - ioff, 0), ipan)
        jdel = min( max(jdel - joff, 0), jpan)
        n = min( max(nface(iq,k) + noff, 1), npan)
        ! bi-cubic interpolation
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
      end do
    end do
    !$acc end parallel loop
  else
    !$acc parallel loop collapse(2) copyin(sx) copyout(s)        &
    !$acc   present(xg,yg,nface)
    do k = 1,wlev
      do iq = 1,ifull
        idel = int(xg(iq,k))
        xxg = xg(iq,k) - real(idel)
        jdel = int(yg(iq,k))
        yyg = yg(iq,k) - real(jdel)
        idel = min( max(idel - ioff, 0), ipan)
        jdel = min( max(jdel - joff, 0), jpan)
        n = min( max(nface(iq,k) + noff, 1), npan)
        ! bi-cubic interpolation
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
        s(iq,k) = &
            rmul_1*emul_1 + rmul_2*emul_2 + rmul_3*emul_3 + rmul_4*emul_4
      end do
    end do
    !$acc end parallel loop
  end if

  call intssync_recv(s)  
    
end if                     ! (intsch==1) .. else ..
!========================   end of intsch=1 section ====================

do k = 1,wlev
  where ( .not.wtr(1:ifull,k) .or. s(1:ifull,k)<cxx+10. )
    s(1:ifull,k) = s_store(1:ifull,k)
  end where
end do
if ( sal_test ) then
  do k = 1,wlev
    where ( s_store(1:ifull,k)<2. )
      s(1:ifull,k) = s_store(1:ifull,k)  
    end where
  end do
end if

call END_LOG(waterints_end)

end subroutine mlob2ints_bs_2

end module mloints
