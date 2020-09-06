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
    
module mlodepts

implicit none

private
public mlodeps

contains


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Calculate depature points for MLO semi-Lagrangian advection
! (This subroutine is based on depts.f)

subroutine mlodeps(ubar,vbar,nface,xg,yg,x3d,y3d,z3d,wtr)

use cc_mpi
use const_phys
use indices_m
use mlo
use newmpar_m
use parm_m
use parmhor_m
use vecsuv_m
use xyzinfo_m

implicit none

integer iq,i,j,k,n,nn,idel,jdel,intsch,ii
integer, dimension(ifull,wlev), intent(out) :: nface
integer s_count
real, dimension(ifull,wlev), intent(in) :: ubar,vbar
real, dimension(ifull,wlev), intent(out) :: xg,yg
real(kind=8), dimension(ifull,wlev), intent(out) :: x3d,y3d,z3d
real, dimension(ifull,wlev) :: uc,vc,wc
real, dimension(ifull+iextra,wlev,3) :: s, s_old
real, dimension(-1:ipan+2,-1:jpan+2,1:npan,wlev,3) :: sx
real s_tot, s_test
real dmul_2, dmul_3, cmul_1, cmul_2, cmul_3, cmul_4
real emul_1, emul_2, emul_3, emul_4, rmul_1, rmul_2, rmul_3, rmul_4
real xxg,yyg
real, parameter :: cxx = -9999. ! missing value flag
logical, dimension(ifull+iextra,wlev), intent(in) :: wtr

call START_LOG(waterdeps_begin)

! departure point x, y, z is called x3d, y3d, z3d
! first find corresponding cartesian vels
do k = 1,wlev
  uc(:,k) = (ax(1:ifull)*ubar(:,k)+bx(1:ifull)*vbar(:,k))*dt/rearth ! unit sphere 
  vc(:,k) = (ay(1:ifull)*ubar(:,k)+by(1:ifull)*vbar(:,k))*dt/rearth ! unit sphere 
  wc(:,k) = (az(1:ifull)*ubar(:,k)+bz(1:ifull)*vbar(:,k))*dt/rearth ! unit sphere 
  x3d(:,k) = x - uc(:,k) ! 1st guess
  y3d(:,k) = y - vc(:,k)
  z3d(:,k) = z - wc(:,k)
end do

!$acc data create(xg,yg,nface,sx,s)

! convert to grid point numbering
call mlotoij5(x3d,y3d,z3d,nface,xg,yg)

! Share off processor departure points.
call deptsync(nface,xg,yg)

intsch = mod(ktau,2)

do k = 1,wlev
  where (wtr(1:ifull,k))
    s(1:ifull,k,1) = uc(1:ifull,k) 
    s(1:ifull,k,2) = vc(1:ifull,k)
    s(1:ifull,k,3) = wc(1:ifull,k)
  else where    
    s(1:ifull,k,1) = cxx - 1. ! missing value flag
    s(1:ifull,k,2) = cxx - 1.
    s(1:ifull,k,3) = cxx - 1.
  end where
end do

! fill
do ii = 1,3 ! 3 iterations of fill should be enough
  s_old(1:ifull,:,:) = s(1:ifull,:,:)
  call bounds(s_old)
  do concurrent (nn = 1:3)
    do concurrent (k = 1:wlev)
      do concurrent (iq = 1:ifull)
        if ( s(iq,k,nn)<cxx ) then
          s_tot = 0.
          s_count = 0
          s_test = s_old(is(iq),k,nn)
          if ( s_test>cxx ) then
            s_tot = s_tot + s_test
            s_count = s_count +1
          end if
          s_test = s_old(in(iq),k,nn)
          if ( s_test>cxx ) then
            s_tot = s_tot + s_test
            s_count = s_count +1
          end if 
          s_test = s_old(ie(iq),k,nn)
          if ( s_test>cxx ) then
            s_tot = s_tot + s_test
            s_count = s_count +1
          end if
          s_test = s_old(iw(iq),k,nn)
          if ( s_test>cxx ) then
            s_tot = s_tot + s_test
            s_count = s_count +1
          end if
          if ( s_count>0 ) then
            s(iq,k,nn) = s_tot/real(s_count)
          end if
        end if
      end do
    end do
  end do
end do

call bounds(s,nrows=2)

!$acc update device(s)

!======================== start of intsch=1 section ====================
if ( intsch==1 ) then

  !$acc parallel loop collapse(5) present(sx,s)
  do concurrent (nn = 1:3)
    do concurrent (k = 1:wlev)
      do concurrent (n = 1:npan)
        do concurrent (j = 1:jpan)
          do concurrent (i = 1:ipan)
            iq = i + (j-1)*ipan + (n-1)*ipan*jpan
            sx(i,j,n,k,nn) = s(iq,k,nn)
          end do
        end do
      end do
    end do
  end do
  !$acc end parallel loop
  !$acc parallel loop collapse(3) present(s,sx)
  do concurrent (nn = 1:3)
    do concurrent (k = 1:wlev)
      do concurrent (n = 1:npan)
        do concurrent (j = 1:jpan)
          iq = 1+(j-1)*ipan+(n-1)*ipan*jpan
          sx(0,j,n,k,nn)      = s( iw(iq),k,nn)
          sx(-1,j,n,k,nn)     = s(iww(iq),k,nn)
          iq = j*ipan+(n-1)*ipan*jpan
          sx(ipan+1,j,n,k,nn) = s( ie(iq),k,nn)
          sx(ipan+2,j,n,k,nn) = s(iee(iq),k,nn)
        end do
        do concurrent (i = 1:ipan)
          iq = i+(n-1)*ipan*jpan
          sx(i,0,n,k,nn)      = s( is(iq),k,nn)
          sx(i,-1,n,k,nn)     = s(iss(iq),k,nn)
          iq = i-ipan+n*ipan*jpan
          sx(i,jpan+1,n,k,nn) = s( in(iq),k,nn)
          sx(i,jpan+2,n,k,nn) = s(inn(iq),k,nn)
        end do
        sx(-1,0,n,k,nn)          = s(lwws(n),k,nn)
        sx(0,0,n,k,nn)           = s(iws(1+(n-1)*ipan*jpan),   k,nn)
        sx(0,-1,n,k,nn)          = s(lwss(n),k,nn)
        sx(ipan+1,0,n,k,nn)      = s(ies(ipan+(n-1)*ipan*jpan),k,nn)
        sx(ipan+2,0,n,k,nn)      = s(lees(n),k,nn)
        sx(ipan+1,-1,n,k,nn)     = s(less(n),k,nn)
        sx(-1,jpan+1,n,k,nn)     = s(lwwn(n),k,nn)
        sx(0,jpan+2,n,k,nn)      = s(lwnn(n),k,nn)
        sx(ipan+2,jpan+1,n,k,nn) = s(leen(n),k,nn)
        sx(ipan+1,jpan+2,n,k,nn) = s(lenn(n),k,nn)
        sx(0,jpan+1,n,k,nn)      = s(iwn(1-ipan+n*ipan*jpan),k,nn)
        sx(ipan+1,jpan+1,n,k,nn) = s(ien(n*ipan*jpan),       k,nn)
      end do               ! n loop
    end do                 ! k loop
  end do                   ! nn loop
  !$acc end parallel loop

  !$acc update self(sx)
 
! Loop over points that need to be calculated for other processes
  do ii = 1,neighnum
    do nn = 1,3
      do iq = 1,drlen(ii)
        n = nint(dpoints(ii)%a(1,iq)) + noff ! Local index
        !  Need global face index in fproc call
        idel = int(dpoints(ii)%a(2,iq))
        xxg = dpoints(ii)%a(2,iq) - real(idel)
        jdel = int(dpoints(ii)%a(3,iq))
        yyg = dpoints(ii)%a(3,iq) - real(jdel)
        k = nint(dpoints(ii)%a(4,iq))
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
        sextra(ii)%a(iq+(nn-1)*drlen(ii)) = rmul_1*emul_1 + rmul_2*emul_2 + rmul_3*emul_3 + rmul_4*emul_4
      end do          ! iq loop
    end do            ! nn loop
  end do              ! ii loop
  
  call intssync_send(3)

  !$acc parallel loop collapse(3) present(xg,yg,nface,sx,s)
  do concurrent (nn = 1:3)
    do concurrent (k = 1:wlev)
      do concurrent (iq = 1:ifull)
        idel = int(xg(iq,k))
        xxg  = xg(iq,k) - real(idel)
        jdel = int(yg(iq,k))
        yyg  = yg(iq,k) - real(jdel)
        idel = min( max( idel - ioff, 0), ipan )
        jdel = min( max( jdel - joff, 0), jpan )
        n = min( max( nface(iq,k) + noff, 1), npan )
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
        s(iq,k,nn) = rmul_1*emul_1 + rmul_2*emul_2 + rmul_3*emul_3 + rmul_4*emul_4
      end do     ! iq loop
    end do       ! k loop
  end do         ! nn loop
  !$acc end parallel loop
  !$acc update self(s)
       
!========================   end of intsch=1 section ====================
else     ! if(intsch==1)then
!======================== start of intsch=2 section ====================

  !$acc parallel loop collapse(5) present(sx,s)
  do concurrent (nn = 1:3)
    do concurrent (k = 1:wlev)
      do concurrent (n = 1:npan)
        do concurrent (j = 1:jpan)
          do concurrent (i = 1:ipan)
            iq = i + (j-1)*ipan + (n-1)*ipan*jpan
            sx(i,j,n,k,nn) = s(iq,k,nn)
          end do
        end do
      end do
    end do
  end do
  !$acc end parallel loop
  !$acc parallel loop collapse(3) present(s,sx)
  do concurrent (nn = 1:3)
    do concurrent (k = 1:wlev)
      do concurrent (n = 1:npan)
        do concurrent (j = 1:jpan)
          iq = 1+(j-1)*ipan+(n-1)*ipan*jpan
          sx(0,j,n,k,nn)      = s( iw(iq),k,nn)
          sx(-1,j,n,k,nn)     = s(iww(iq),k,nn)
          iq = j*ipan+(n-1)*ipan*jpan
          sx(ipan+1,j,n,k,nn) = s( ie(iq),k,nn)
          sx(ipan+2,j,n,k,nn) = s(iee(iq),k,nn)
        end do            ! j loop
        do concurrent (i = 1:ipan)
          iq = i+(n-1)*ipan*jpan
          sx(i,0,n,k,nn)      = s( is(iq),k,nn)
          sx(i,-1,n,k,nn)     = s(iss(iq),k,nn)
          iq = i-ipan+n*ipan*jpan
          sx(i,jpan+1,n,k,nn) = s( in(iq),k,nn)
          sx(i,jpan+2,n,k,nn) = s(inn(iq),k,nn)
        end do            ! i loop
        sx(-1,0,n,k,nn)          = s(lsww(n),k,nn)
        sx(0,0,n,k,nn)           = s(isw(1+(n-1)*ipan*jpan),k,nn)
        sx(0,-1,n,k,nn)          = s(lssw(n),k,nn)
        sx(ipan+2,0,n,k,nn)      = s(lsee(n),k,nn)
        sx(ipan+1,-1,n,k,nn)     = s(lsse(n),k,nn)
        sx(-1,jpan+1,n,k,nn)     = s(lnww(n),k,nn)
        sx(0,jpan+1,n,k,nn)      = s(inw(1-ipan+n*ipan*jpan),k,nn)
        sx(0,jpan+2,n,k,nn)      = s(lnnw(n),k,nn)
        sx(ipan+2,jpan+1,n,k,nn) = s(lnee(n),k,nn)
        sx(ipan+1,jpan+2,n,k,nn) = s(lnne(n),k,nn)
        sx(ipan+1,0,n,k,nn)      = s(ise(ipan+(n-1)*ipan*jpan),k,nn)
        sx(ipan+1,jpan+1,n,k,nn) = s(ine(n*ipan*jpan),         k,nn)
      end do           ! n loop
    end do             ! k loop
  end do               ! nn loop
  !$acc end parallel loop

  !$acc update self(sx)

! For other processes
  do ii = 1,neighnum
    do nn = 1,3
      do iq = 1,drlen(ii)
        n = nint(dpoints(ii)%a(1,iq)) + noff ! Local index
        idel = int(dpoints(ii)%a(2,iq))
        xxg = dpoints(ii)%a(2,iq) - real(idel)
        jdel = int(dpoints(ii)%a(3,iq))
        yyg = dpoints(ii)%a(3,iq) - real(jdel)
        k = nint(dpoints(ii)%a(4,iq))
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
        sextra(ii)%a(iq+(nn-1)*drlen(ii)) = rmul_1*emul_1 + rmul_2*emul_2 + rmul_3*emul_3 + rmul_4*emul_4
      end do          ! iq loop
    end do            ! nn loop
  end do              ! ii loop

  call intssync_send(3)

  !$acc parallel loop collapse(3) present(xg,yg,nface,sx,s)
  do concurrent (nn = 1:3)
    do concurrent (k = 1:wlev)
      do concurrent (iq = 1:ifull)
        idel = int(xg(iq,k))
        xxg = xg(iq,k) - real(idel)
        jdel = int(yg(iq,k))
        yyg = yg(iq,k) - real(jdel)
        idel = min( max( idel - ioff, 0), ipan )
        jdel = min( max( jdel - joff, 0), jpan )
        n = min( max( nface(iq,k) + noff, 1), npan )
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
        s(iq,k,nn) = rmul_1*emul_1 + rmul_2*emul_2 + rmul_3*emul_3 + rmul_4*emul_4
      end do
    end do
  end do
  !$acc end parallel loop
  !$acc update self(s)

end if                     ! (intsch==1) .. else ..
!========================   end of intsch=1 section ====================

call intssync_recv(s)

do k = 1,wlev
  where ( wtr(1:ifull,k) )
    x3d(:,k) = x - 0.5*(uc(:,k)+s(1:ifull,k,1)) ! n+1 guess
    y3d(:,k) = y - 0.5*(vc(:,k)+s(1:ifull,k,2)) ! n+1 guess
    z3d(:,k) = z - 0.5*(wc(:,k)+s(1:ifull,k,3)) ! n+1 guess
  elsewhere
    x3d(:,k) = x
    y3d(:,k) = y
    z3d(:,k) = z
  end where
end do

call mlotoij5(x3d,y3d,z3d,nface,xg,yg)
!     Share off processor departure points.
call deptsync(nface,xg,yg)

!======================== start of intsch=1 section ====================
if ( intsch==1 ) then
 
  ! Loop over points that need to be calculated for other processes
  do ii = 1,neighnum
    do nn = 1,3
      do iq = 1,drlen(ii)
        n = nint(dpoints(ii)%a(1,iq)) + noff ! Local index
        idel = int(dpoints(ii)%a(2,iq))
        xxg = dpoints(ii)%a(2,iq) - real(idel)
        jdel = int(dpoints(ii)%a(3,iq))
        yyg = dpoints(ii)%a(3,iq) - real(jdel)
        k = nint(dpoints(ii)%a(4,iq))
        idel = idel - ioff
        jdel = jdel - joff
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
      end do          ! iq loop
    end do            ! nn loop
  end do              ! ii loop
  
  call intssync_send(3)

  !$acc parallel loop collapse(3) present(xg,yg,nface,sx,s)
  do concurrent (nn = 1:3)
    do concurrent (k = 1:wlev)
      do concurrent (iq = 1:ifull)
        idel = int(xg(iq,k))
        xxg  = xg(iq,k) - real(idel)
        jdel = int(yg(iq,k))
        yyg  = yg(iq,k) - real(jdel)
        idel = min( max( idel - ioff, 0), ipan )
        jdel = min( max( jdel - joff, 0), jpan )
        n = min( max( nface(iq,k) + noff, 1), npan )
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
        s(iq,k,nn) = rmul_1*emul_1 + rmul_2*emul_2 + rmul_3*emul_3 + rmul_4*emul_4
      end do     ! iq loop
    end do       ! k loop
  end do         ! nn loop
  !$acc end parallel loop
  !$acc update self(s)
       
!========================   end of intsch=1 section ====================
else     ! if(intsch==1)then
!======================== start of intsch=2 section ====================

  ! For other processes
  do ii = 1,neighnum
    do nn = 1,3
      do iq = 1,drlen(ii)
        n = nint(dpoints(ii)%a(1,iq)) + noff ! Local index
        idel = int(dpoints(ii)%a(2,iq))
        xxg = dpoints(ii)%a(2,iq) - real(idel)
        jdel = int(dpoints(ii)%a(3,iq))
        yyg = dpoints(ii)%a(3,iq) - real(jdel)
        k = nint(dpoints(ii)%a(4,iq))
        idel = idel - ioff
        jdel = jdel - joff
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
      end do          ! iq loop
    end do            ! nn loop
  end do              ! ii loop

  call intssync_send(3)

  !$acc parallel loop collapse(3) present(xg,yg,nface,sx,s)
  do concurrent (nn = 1:3)
    do concurrent (k = 1:wlev)
      do concurrent (iq = 1:ifull)
        idel = int(xg(iq,k))
        xxg = xg(iq,k) - real(idel)
        jdel = int(yg(iq,k))
        yyg = yg(iq,k) - real(jdel)
        idel = min( max( idel - ioff, 0), ipan )
        jdel = min( max( jdel - joff, 0), jpan )
        n = min( max( nface(iq,k) + noff, 1), npan )
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
        s(iq,k,nn) = rmul_1*emul_1 + rmul_2*emul_2 + rmul_3*emul_3 + rmul_4*emul_4
      end do
    end do
  end do
  !$acc end parallel loop
  !$acc update self(s)

end if                     ! (intsch==1) .. else ..
!========================   end of intsch=1 section ====================

call intssync_recv(s)

do k = 1,wlev
  where (wtr(1:ifull,k))
    x3d(:,k) = x - 0.5*(uc(:,k)+s(1:ifull,k,1)) ! n+1 guess
    y3d(:,k) = y - 0.5*(vc(:,k)+s(1:ifull,k,2)) ! n+1 guess
    z3d(:,k) = z - 0.5*(wc(:,k)+s(1:ifull,k,3)) ! n+1 guess
  elsewhere
    x3d(:,k) = x
    y3d(:,k) = y
    z3d(:,k) = z
  end where
end do

call mlotoij5(x3d,y3d,z3d,nface,xg,yg)
!     Share off processor departure points.
call deptsync(nface,xg,yg)

!$acc end data

call END_LOG(waterdeps_end)

return
end subroutine mlodeps

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Calculate indices
! This code is from depts.f90

pure subroutine mlotoij5(x3d,y3d,z3d,nface,xg,yg)

use bigxy4_m
use cc_mpi
use mlo
use newmpar_m
use parm_m
use parmgeom_m
use xyzinfo_m

implicit none

integer loop,iq,i,j,is,js
integer ii
integer, dimension(ifull,wlev), intent(out) :: nface
real, dimension(ifull,wlev), intent(out) :: xg,yg
real xstr,ystr,zstr
real denxyz,xd,yd,zd
real ri,rj
real(kind=8), dimension(ifull,wlev), intent(inout) :: x3d,y3d,z3d
real(kind=8) den
real(kind=8) alf,alfonsch
real(kind=8) dxx,dxy,dyx,dyy
integer, parameter :: nmaploop = 3

alf = (1._8-schmidt**2)/(1._8+schmidt**2)
alfonsch = 2._8*schmidt/(1._8+schmidt**2)

!$acc parallel loop collapse(2) copyin(x3d,y3d,z3d,xx4,yy4) present(xg,yg,nface)
do concurrent (ii = 1:wlev)
  do concurrent (iq = 1:ifull)

    !     if necessary, transform (x3d, y3d, z3d) to equivalent
    !     coordinates (xstr, ystr, zstr) on regular gnomonic panels
    den=1._8-alf*z3d(iq,ii) ! to force real*8
    xstr=real(x3d(iq,ii)*(alfonsch/den))
    ystr=real(y3d(iq,ii)*(alfonsch/den))
    zstr=real((z3d(iq,ii)-alf)/den)

    !      first deduce departure faces
    !      instead calculate cubic coordinates
    !      The faces are:
    !      0: X=1   1: Z=1   2: Y=1   3: X=-1   4: Z=-1   5: Y=-1
    denxyz=max( abs(xstr),abs(ystr),abs(zstr) )
    xd=xstr/denxyz
    yd=ystr/denxyz
    zd=zstr/denxyz

    if (abs(xstr-denxyz)<1.E-8) then
      nface(iq,ii)    =0
      xg(iq,ii) =      yd
      yg(iq,ii) =      zd
    else if (abs(xstr+denxyz)<1.E-8) then
      nface(iq,ii)    =3
      xg(iq,ii) =     -zd
      yg(iq,ii) =     -yd
    else if (abs(zstr-denxyz)<1.E-8) then
      nface(iq,ii)    =1
      xg(iq,ii) =      yd
      yg(iq,ii) =     -xd
    else if (abs(zstr+denxyz)<1.E-8) then
      nface(iq,ii)    =4
      xg(iq,ii) =      xd
      yg(iq,ii) =     -yd
    else if (abs(ystr-denxyz)<1.E-8) then
      nface(iq,ii)    =2
      xg(iq,ii) =     -zd
      yg(iq,ii) =     -xd
    else
      nface(iq,ii)    =5
      xg(iq,ii) =      xd
      yg(iq,ii) =      zd
    end if

    !     use 4* resolution grid il --> 4*il
    xg(iq,ii)=min(max(-.999999,xg(iq,ii)),.999999)
    yg(iq,ii)=min(max(-.999999,yg(iq,ii)),.999999)
    !      first guess for ri, rj and nearest i,j
    ri=1.+(1.+xg(iq,ii))*real(2*il_g)
    rj=1.+(1.+yg(iq,ii))*real(2*il_g)
    do loop=1,nmaploop
      i=nint(ri)
      j=nint(rj)
      is=nint(sign(1.,ri-real(i)))
      js=nint(sign(1.,rj-real(j)))
      ! predict new value for ri, rj
      dxx=xx4(i+is,j)-xx4(i,j)
      dyx=xx4(i,j+js)-xx4(i,j)
      dxy=yy4(i+is,j)-yy4(i,j)
      dyy=yy4(i,j+js)-yy4(i,j)       
      den=dxx*dyy-dyx*dxy
      ri = real(i) + real(is)*real(((xg(iq,ii)-xx4(i,j))*dyy-(yg(iq,ii)-yy4(i,j))*dyx)/den)
      rj = real(j) + real(js)*real(((yg(iq,ii)-yy4(i,j))*dxx-(xg(iq,ii)-xx4(i,j))*dxy)/den)
      ri = min(ri,1.0+1.999999*real(2*il_g))
      ri = max(ri,1.0+0.000001*real(2*il_g))
      rj = min(rj,1.0+1.999999*real(2*il_g))
      rj = max(rj,1.0+0.000001*real(2*il_g))
    end do  ! loop loop
    !      expect xg, yg to range between .5 and il+.5
    xg(iq,ii)=0.25*(ri+3.)-0.5  ! -.5 for stag; back to normal ri, rj defn
    yg(iq,ii)=0.25*(rj+3.)-0.5  ! -.5 for stag

  end do
end do
!$acc end parallel loop

!$acc update self(xg,yg,nface)

return
end subroutine mlotoij5

end module mlodepts
