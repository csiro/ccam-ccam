! Conformal Cubic Atmospheric Model
    
! Copyright 2015-2016 Commonwealth Scientific Industrial Research Organisation (CSIRO)
    
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
    
subroutine depts1(x3d,y3d,z3d,intsch)  ! input ubar,vbar are unstaggered vels for level k

use cc_mpi
use const_phys
use indices_m
use map_m
use newmpar_m
use parm_m
use parmhor_m
use uvbar_m
use vecsuv_m
use work3f_m
use xyzinfo_m

implicit none

integer, intent(in) :: intsch
integer iq, k, idel, jdel, nn
integer i, j, n, ii
real xxg, yyg
real, dimension(ifull,kl) :: uc, vc, wc
real, dimension(ifull+iextra,kl,3) :: s
real, dimension(-1:ipan+2,-1:jpan+2,1:npan,kl,3) :: sx
real(kind=8), dimension(ifull,kl) :: x3d, y3d, z3d   ! upglobal depts 
real dmul_2, dmul_3, cmul_1, cmul_2, cmul_3, cmul_4
real emul_1, emul_2, emul_3, emul_4, rmul_1, rmul_2, rmul_3, rmul_4
      
call START_LOG(depts_begin)

do k = 1,kl
  ! departure point x, y, z is called x3d, y3d, z3d
  ! first find corresponding cartesian vels
  uc(1:ifull,k) = (ax(1:ifull)*ubar(1:ifull,k) + bx(1:ifull)*vbar(1:ifull,k))*dt/rearth ! unit sphere 
  vc(1:ifull,k) = (ay(1:ifull)*ubar(1:ifull,k) + by(1:ifull)*vbar(1:ifull,k))*dt/rearth ! unit sphere 
  wc(1:ifull,k) = (az(1:ifull)*ubar(1:ifull,k) + bz(1:ifull)*vbar(1:ifull,k))*dt/rearth ! unit sphere 
  x3d(1:ifull,k) = x(1:ifull) - real(uc(1:ifull,k),8) ! 1st guess
  y3d(1:ifull,k) = y(1:ifull) - real(vc(1:ifull,k),8)
  z3d(1:ifull,k) = z(1:ifull) - real(wc(1:ifull,k),8)
end do

! convert to grid point numbering
call toij5(x3d,y3d,z3d)

! Share off processor departure points.
call deptsync(nface,xg,yg)

if ( diag .and. mydiag ) then
  write(6,*) 'ubar,vbar ',ubar(idjd,nlv),vbar(idjd,nlv)
  write(6,*) 'uc,vc,wc ',uc(idjd,nlv),vc(idjd,nlv),wc(idjd,nlv)
  write(6,*) 'x,y,z ',x(idjd),y(idjd),z(idjd)
  write(6,*) '1st guess for k = ',nlv
  write(6,*) 'x3d,y3d,z3d ',x3d(idjd,nlv),y3d(idjd,nlv),z3d(idjd,nlv)
  write(6,*) 'xg,yg,nface ',xg(idjd,nlv),yg(idjd,nlv),nface(idjd,nlv)
endif

s(1:ifull,:,1) = uc(1:ifull,:)
s(1:ifull,:,2) = vc(1:ifull,:)
s(1:ifull,:,3) = wc(1:ifull,:)

call bounds(s,nrows=2)

!$acc data create(xg,yg,nface,sx,s)
!$acc update device(xg,yg,nface)

!======================== start of intsch=1 section ====================
if ( intsch==1 ) then
  sx(1:ipan,1:jpan,1:npan,1:kl,1:3) = reshape(s(1:ipan*jpan*npan,1:kl,1:3), (/ipan,jpan,npan,kl,3/))
  do nn = 1,3
    do k = 1,kl
      do n = 1,npan
        do j = 1,jpan
          iq = 1+(j-1)*ipan+(n-1)*ipan*jpan
          sx(0,j,n,k,nn)      = s(iw(iq),k,nn)
          sx(-1,j,n,k,nn)     = s(iww(iq),k,nn)
          iq = j*ipan+(n-1)*ipan*jpan
          sx(ipan+1,j,n,k,nn) = s(ie(iq),k,nn)
          sx(ipan+2,j,n,k,nn) = s(iee(iq),k,nn)
        end do            ! j loop
        do i = 1,ipan
          iq = i+(n-1)*ipan*jpan
          sx(i,0,n,k,nn)      = s(is(iq),k,nn)
          sx(i,-1,n,k,nn)     = s(iss(iq),k,nn)
          iq = i-ipan+n*ipan*jpan
          sx(i,jpan+1,n,k,nn) = s(in(iq),k,nn)
          sx(i,jpan+2,n,k,nn) = s(inn(iq),k,nn)
        end do            ! i loop
        sx(-1,0,n,k,nn)          = s(lwws(n),k,nn)
        sx(0,0,n,k,nn)           = s(iws(1+(n-1)*ipan*jpan),k,nn)
        sx(0,-1,n,k,nn)          = s(lwss(n),k,nn)
        sx(ipan+1,0,n,k,nn)      = s(ies(ipan+(n-1)*ipan*jpan),k,nn)
        sx(ipan+2,0,n,k,nn)      = s(lees(n),k,nn)
        sx(ipan+1,-1,n,k,nn)     = s(less(n),k,nn)
        sx(-1,jpan+1,n,k,nn)     = s(lwwn(n),k,nn)
        sx(0,jpan+2,n,k,nn)      = s(lwnn(n),k,nn)
        sx(ipan+2,jpan+1,n,k,nn) = s(leen(n),k,nn)
        sx(ipan+1,jpan+2,n,k,nn) = s(lenn(n),k,nn)
        sx(0,jpan+1,n,k,nn)      = s(iwn(1-ipan+n*ipan*jpan),k,nn)
        sx(ipan+1,jpan+1,n,k,nn) = s(ien(n*ipan*jpan),k,nn)
      end do          ! n loop
    end do            ! k loop
  end do              ! nn loop

  !$acc update device(sx)

  ! Loop over points that need to be calculated for other processes
  do ii = neighnum,1,-1
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
      end do        ! iq loop
    end do          ! nn loop
  end do            ! ii loop

  call intssync_send(3)
  
  !$acc parallel loop collapse(3) present(xg,yg,nface,sx,s) 
  do concurrent (nn = 1:3)
    do concurrent (k = 1:kl)
      do concurrent (iq = 1:ifull)    ! non Berm-Stan option
        idel = int(xg(iq,k))
        xxg = xg(iq,k) - real(idel)
        jdel = int(yg(iq,k))
        yyg = yg(iq,k) - real(jdel)
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
        rmul_1 = sx(idel,  jdel-1,n,k,nn)*dmul_2 + sx(idel+1,jdel-1,n,k,nn)*dmul_3
        rmul_2 = sx(idel-1,jdel,  n,k,nn)*cmul_1 + sx(idel,  jdel,  n,k,nn)*cmul_2 + &
                 sx(idel+1,jdel,  n,k,nn)*cmul_3 + sx(idel+2,jdel,  n,k,nn)*cmul_4
        rmul_3 = sx(idel-1,jdel+1,n,k,nn)*cmul_1 + sx(idel,  jdel+1,n,k,nn)*cmul_2 + &
                 sx(idel+1,jdel+1,n,k,nn)*cmul_3 + sx(idel+2,jdel+1,n,k,nn)*cmul_4
        rmul_4 = sx(idel,  jdel+2,n,k,nn)*dmul_2 + sx(idel+1,jdel+2,n,k,nn)*dmul_3
        s(iq,k,nn) = rmul_1*emul_1 + rmul_2*emul_2 + rmul_3*emul_3 + rmul_4*emul_4
      end do   ! iq loop
    end do     ! k loop
  end do       ! nn loop
  !$acc end parallel loop
  !$acc update self(s)
            
!========================   end of intsch=1 section ====================
else     ! if(intsch==1)then
!======================== start of intsch=2 section ====================

  sx(1:ipan,1:jpan,1:npan,1:kl,1:3) = reshape(s(1:ipan*jpan*npan,1:kl,1:3), (/ipan,jpan,npan,kl,3/))
  do nn = 1,3
    do k = 1,kl
      do n = 1,npan
        do j = 1,jpan
          iq = 1+(j-1)*ipan+(n-1)*ipan*jpan
          sx(0,j,n,k,nn)      = s(iw(iq),k,nn)
          sx(-1,j,n,k,nn)     = s(iww(iq),k,nn)
          iq = j*ipan+(n-1)*ipan*jpan
          sx(ipan+1,j,n,k,nn) = s(ie(iq),k,nn)
          sx(ipan+2,j,n,k,nn) = s(iee(iq),k,nn)
        end do            ! j loop
        do i = 1,ipan
          iq = i+(n-1)*ipan*jpan
          sx(i,0,n,k,nn)      = s(is(iq),k,nn)
          sx(i,-1,n,k,nn)     = s(iss(iq),k,nn)
          iq = i-ipan+n*ipan*jpan
          sx(i,jpan+1,n,k,nn) = s(in(iq),k,nn)
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
        sx(ipan+1,jpan+1,n,k,nn) = s(ine(n*ipan*jpan),k,nn)
      end do              ! n loop
    end do                ! k loop
  end do                  ! nn loop

  !$acc update device(sx)

  ! For other processes
  do ii = neighnum,1,-1
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
      end do          ! iq loop
    end do            ! nn loop
  end do              ! ii

  call intssync_send(3)

  !$acc parallel loop collapse(3) present(xg,yg,nface,sx,s) 
  do concurrent (nn = 1:3)
    do concurrent (k = 1:kl)
      do concurrent (iq = 1:ifull)    ! non Berm-Stan option
        idel = int(xg(iq,k))
        xxg = xg(iq,k) - real(idel)
        jdel = int(yg(iq,k))
        yyg = yg(iq,k) - real(jdel)
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
        rmul_1 = sx(idel-1,jdel,  n,k,nn)*dmul_2 + sx(idel-1,jdel+1,n,k,nn)*dmul_3
        rmul_2 = sx(idel,  jdel-1,n,k,nn)*cmul_1 + sx(idel,  jdel,  n,k,nn)*cmul_2 + &
                 sx(idel,  jdel+1,n,k,nn)*cmul_3 + sx(idel,  jdel+2,n,k,nn)*cmul_4
        rmul_3 = sx(idel+1,jdel-1,n,k,nn)*cmul_1 + sx(idel+1,jdel,  n,k,nn)*cmul_2 + &
                 sx(idel+1,jdel+1,n,k,nn)*cmul_3 + sx(idel+1,jdel+2,n,k,nn)*cmul_4
        rmul_4 = sx(idel+2,jdel,  n,k,nn)*dmul_2 + sx(idel+2,jdel+1,n,k,nn)*dmul_3
        s(iq,k,nn) = rmul_1*emul_1 + rmul_2*emul_2 + rmul_3*emul_3 + rmul_4*emul_4
      end do          ! iq loop
    end do            ! k loop
  end do              ! nn loop
  !$acc end parallel loop
  !$acc update self(s)

endif                     ! (intsch==1) .. else ..
!========================   end of intsch=1 section ====================

call intssync_recv(s)

do k = 1,kl
  x3d(1:ifull,k) = x(1:ifull) - 0.5_8*(real(uc(1:ifull,k),8)+real(s(1:ifull,k,1),8)) ! 2nd guess
  y3d(1:ifull,k) = y(1:ifull) - 0.5_8*(real(vc(1:ifull,k),8)+real(s(1:ifull,k,2),8)) ! 2nd guess
  z3d(1:ifull,k) = z(1:ifull) - 0.5_8*(real(wc(1:ifull,k),8)+real(s(1:ifull,k,3),8)) ! 2nd guess
end do

call toij5 (x3d,y3d,z3d)
!     Share off processor departure points.
call deptsync(nface,xg,yg)
!$acc update device(nface,xg,yg)

if ( diag .and. mydiag ) then
  write(6,*) '2nd guess for k = ',nlv
  write(6,*) 'x3d,y3d,z3d ',x3d(idjd,nlv),y3d(idjd,nlv),z3d(idjd,nlv)
  write(6,*) 'xg,yg,nface ',xg(idjd,nlv),yg(idjd,nlv),nface(idjd,nlv)
end if

!======================== start of intsch=1 section ====================
if ( intsch==1 ) then

  ! Loop over points that need to be calculated for other processes
  do ii = neighnum,1,-1
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
      end do        ! iq loop
    end do          ! nn loop
  end do            ! ii loop

  call intssync_send(3)

  !$acc parallel loop collapse(3) present(xg,yg,nface,sx,s) 
  do concurrent (nn = 1:3)
    do concurrent (k = 1:kl)
      do concurrent (iq = 1:ifull)    ! non Berm-Stan option
        idel = int(xg(iq,k))
        xxg = xg(iq,k) - real(idel)
        jdel = int(yg(iq,k))
        yyg = yg(iq,k) - real(jdel)
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
        rmul_1 = sx(idel,  jdel-1,n,k,nn)*dmul_2 + sx(idel+1,jdel-1,n,k,nn)*dmul_3
        rmul_2 = sx(idel-1,jdel,  n,k,nn)*cmul_1 + sx(idel,  jdel,  n,k,nn)*cmul_2 + &
                 sx(idel+1,jdel,  n,k,nn)*cmul_3 + sx(idel+2,jdel,  n,k,nn)*cmul_4
        rmul_3 = sx(idel-1,jdel+1,n,k,nn)*cmul_1 + sx(idel,  jdel+1,n,k,nn)*cmul_2 + &
                 sx(idel+1,jdel+1,n,k,nn)*cmul_3 + sx(idel+2,jdel+1,n,k,nn)*cmul_4
        rmul_4 = sx(idel,  jdel+2,n,k,nn)*dmul_2 + sx(idel+1,jdel+2,n,k,nn)*dmul_3
        s(iq,k,nn) = rmul_1*emul_1 + rmul_2*emul_2 + rmul_3*emul_3 + rmul_4*emul_4
      end do   ! iq loop
    end do     ! k loop
  end do       ! nn loop
  !$acc end parallel loop
  !$acc update self(s)
            
!========================   end of intsch=1 section ====================
else     ! if(intsch==1)then
!======================== start of intsch=2 section ====================

  ! For other processes
  do ii = neighnum,1,-1
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
      end do          ! iq loop
    end do            ! nn loop
  end do              ! ii

  call intssync_send(3)

  !$acc parallel loop collapse(3) present(xg,yg,nface,sx,s) 
  do concurrent (nn = 1:3)
    do concurrent (k = 1:kl)
      do concurrent (iq = 1:ifull)    ! non Berm-Stan option
        ! Convert face index from 0:npanels to array indices
        idel = int(xg(iq,k))
        xxg = xg(iq,k) - real(idel)
        jdel = int(yg(iq,k))
        yyg = yg(iq,k) - real(jdel)
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
        rmul_1 = sx(idel-1,jdel,  n,k,nn)*dmul_2 + sx(idel-1,jdel+1,n,k,nn)*dmul_3
        rmul_2 = sx(idel,  jdel-1,n,k,nn)*cmul_1 + sx(idel,  jdel,  n,k,nn)*cmul_2 + &
                 sx(idel,  jdel+1,n,k,nn)*cmul_3 + sx(idel,  jdel+2,n,k,nn)*cmul_4
        rmul_3 = sx(idel+1,jdel-1,n,k,nn)*cmul_1 + sx(idel+1,jdel,  n,k,nn)*cmul_2 + &
                 sx(idel+1,jdel+1,n,k,nn)*cmul_3 + sx(idel+1,jdel+2,n,k,nn)*cmul_4
        rmul_4 = sx(idel+2,jdel,  n,k,nn)*dmul_2 + sx(idel+2,jdel+1,n,k,nn)*dmul_3
        s(iq,k,nn) = rmul_1*emul_1 + rmul_2*emul_2 + rmul_3*emul_3 + rmul_4*emul_4
      end do          ! iq loop
    end do            ! nn loop
  end do              ! k loop
  !$acc end parallel loop
  !$acc update self(s)

endif                     ! (intsch==1) .. else ..
!========================   end of intsch=1 section ====================

!$acc end data

call intssync_recv(s)

do k = 1,kl
  x3d(1:ifull,k) = x(1:ifull) - 0.5_8*(real(uc(1:ifull,k),8)+real(s(1:ifull,k,1),8)) ! 3rd guess
  y3d(1:ifull,k) = y(1:ifull) - 0.5_8*(real(vc(1:ifull,k),8)+real(s(1:ifull,k,2),8)) ! 3rd guess
  z3d(1:ifull,k) = z(1:ifull) - 0.5_8*(real(wc(1:ifull,k),8)+real(s(1:ifull,k,3),8)) ! 3rd guess
end do

call toij5(x3d,y3d,z3d)

if ( diag .and. mydiag ) then
  write(6,*) '3rd guess for k = ',nlv
  write(6,*) 'x3d,y3d,z3d ',x3d(idjd,nlv),y3d(idjd,nlv),z3d(idjd,nlv)
  write(6,*) 'xg,yg,nface ',xg(idjd,nlv),yg(idjd,nlv),nface(idjd,nlv)
end if

! Share off processor departure points.
call deptsync(nface,xg,yg)

call END_LOG(depts_end)
      
return
end subroutine depts1


subroutine toij5(x3d,y3d,z3d)

use bigxy4_m
use cc_mpi
use newmpar_m
use parm_m
use parmgeom_m
use work3f_m
use xyzinfo_m

implicit none

#ifdef debug
integer, parameter :: ntest = 0
integer, parameter :: ndiag = 0
#endif

#ifdef cray
integer, save :: num = 0
integer, dimension(ifull) :: nf
real, dimension(0:5), parameter :: xgx = (/ 0., 0., 0., 0., 1., 1. /)
real, dimension(0:5), parameter :: xgy = (/ 1., 1., 0., 0., 0., 0. /)
real, dimension(0:5), parameter :: xgz = (/ 0., 0.,-1.,-1., 0., 0. /)
real, dimension(0:5), parameter :: ygx = (/ 0.,-1.,-1., 0., 0., 0. /)
real, dimension(0:5), parameter :: ygy = (/ 0., 0., 0.,-1.,-1., 0. /)
real, dimension(0:5), parameter :: ygz = (/ 1., 0., 0., 0., 0., 1. /)
#endif

integer, parameter :: nmaploop = 3
integer k, iq, loop, i, j, is, js
real, dimension(ifull) :: ri,rj
real, dimension(ifull) :: xstr,ystr,zstr
real, dimension(ifull) :: denxyz,xd,yd,zd
real(kind=8) alf,alfonsch  ! 6/11/07 esp for 200m
real(kind=8) dxx,dxy,dyx,dyy
real(kind=8), dimension(ifull,kl), intent(in) :: x3d,y3d,z3d
real(kind=8), dimension(ifull) :: den
logical, dimension(ifull) :: xytest, xztest, yztest

#ifdef cray
! check if divide by itself is working
if ( num==0 ) then
  if ( myid==0 ) write(6,*)'checking for ncray'
  call checkdiv(xstr,ystr,zstr)
end if
num = 1
#endif

! if necessary, transform (x3d, y3d, z3d) to equivalent
! coordinates (xstr, ystr, zstr) on regular gnomonic panels
alf           = (1._8-real(schmidt,8)*real(schmidt,8))/(1._8+real(schmidt,8)*real(schmidt,8))
alfonsch      = 2._8*real(schmidt,8)/(1._8+real(schmidt,8)*real(schmidt,8))  ! same but bit more accurate

do k = 1,kl
    
  den(1:ifull)  = 1._8 - alf*z3d(1:ifull,k)
  xstr(1:ifull) = real(x3d(1:ifull,k)*(alfonsch/den(1:ifull)))
  ystr(1:ifull) = real(y3d(1:ifull,k)*(alfonsch/den(1:ifull)))
  zstr(1:ifull) = real(   (z3d(1:ifull,k)-alf)/den(1:ifull))

!      first deduce departure faces
!      instead calculate cubic coordinates
!      The faces are:
!      0: X=1   1: Z=1   2: Y=1   3: X=-1   4: Z=-1   5: Y=-1

  denxyz(1:ifull) = max(abs(xstr(1:ifull)),abs(ystr(1:ifull)),abs(zstr(1:ifull)) )
  xd(1:ifull) = xstr(1:ifull)/denxyz(1:ifull)
  yd(1:ifull) = ystr(1:ifull)/denxyz(1:ifull)
  zd(1:ifull) = zstr(1:ifull)/denxyz(1:ifull)

#ifndef cray
  xytest = abs(xd(1:ifull))>abs(yd(1:ifull))
  xztest = abs(xd(1:ifull))>abs(zd(1:ifull))
  yztest = abs(yd(1:ifull))>abs(zd(1:ifull))
  ! all these if statements are replaced by the subsequent cunning code
  where ( xytest .and. xztest .and. xd(1:ifull)>0. )    ! Cray
    nface(1:ifull,k)    =0                              ! Cray
    xg(1:ifull,k) =       yd(1:ifull)                   ! Cray
    yg(1:ifull,k) =       zd(1:ifull)                   ! Cray
  elsewhere ( xytest .and. xztest )                     ! Cray
    nface(1:ifull,k)    =3                              ! Cray
    xg(1:ifull,k) =     -zd(1:ifull)                    ! Cray
    yg(1:ifull,k) =     -yd(1:ifull)                    ! Cray
  elsewhere ( yztest .and. yd(1:ifull)>0. )             ! Cray
    nface(1:ifull,k)    =2                              ! Cray
    xg(1:ifull,k) =     -zd(1:ifull)                    ! Cray
    yg(1:ifull,k) =     -xd(1:ifull)                    ! Cray
  elsewhere ( yztest )                                  ! Cray
    nface(1:ifull,k)    =5                              ! Cray
    xg(1:ifull,k) =      xd(1:ifull)                    ! Cray
    yg(1:ifull,k) =      zd(1:ifull)                    ! Cray
  elsewhere ( zd(1:ifull)>0. )                          ! Cray
    nface(1:ifull,k)    =1                              ! Cray
    xg(1:ifull,k) =      yd(1:ifull)                    ! Cray
    yg(1:ifull,k) =     -xd(1:ifull)                    ! Cray
  elsewhere                                             ! Cray
    nface(1:ifull,k)    =4                              ! Cray
    xg(1:ifull,k) =      xd(1:ifull)                    ! Cray
    yg(1:ifull,k) =     -yd(1:ifull)                    ! Cray
  end where                                             ! Cray
#else
  ! N.B. the Cray copes poorly with the following (sometimes .ne.1),
  ! with e.g. division of  .978 by itself giving  .99999.....53453
  ! max() allows for 2 of x,y,z being 1.  This is the cunning code:
  nf(1:ifull)=max( int(xd(1:ifull))*(3*int(xd(1:ifull))-3) , int(zd(1:ifull))*(5*int(zd(1:ifull))-3) , &
                   int(yd(1:ifull))*(7*int(yd(1:ifull))-3) )/2.
  nface(1:ifull,k)=nf(1:ifull)
  xg(1:ifull,k)=xgx(nf(1:ifull))*xd(1:ifull)+xgy(nf(1:ifull))*yd(1:ifull)+xgz(nf(1:ifull))*zd(1:ifull)  ! -1 to 1
  yg(1:ifull,k)=ygx(nf(1:ifull))*xd(1:ifull)+ygy(nf(1:ifull))*yd(1:ifull)+ygz(nf(1:ifull))*zd(1:ifull)
#endif

#ifdef debug
  if(ntest==1.and.k==nlv)then
    iq=idjd
    write(6,*) 'x3d,y3d,z3d ',x3d(iq),y3d(iq),z3d(iq)
    den(iq)=1._8-alf*z3d(iq) ! to force real*8
    write(6,*) 'den ',den(iq)
    denxyz(iq)=max( abs(xstr(iq)),abs(ystr(iq)),abs(zstr(iq)) )
    xd(iq)=xstr(iq)/denxyz(iq)
    yd(iq)=ystr(iq)/denxyz(iq)
    zd(iq)=zstr(iq)/denxyz(iq)
    write(6,*) 'k,xstr,ystr,zstr,denxyz ',k,xstr(iq),ystr(iq),zstr(iq),denxyz
    write(6,*) 'abs(xstr,ystr,zstr) ',abs(xstr(iq)),abs(ystr(iq)),abs(zstr(iq))
    write(6,*) 'xd,yd,zd,nface ',xd(iq),yd(iq),zd(iq),nface(iq,k)
    write(6,*) 'alf,alfonsch ',alf,alfonsch
  endif
  if(ndiag==2)then
    call printp('xcub',xd)  ! need to reinstate as arrays for this diag
    call printp('ycub',yd)
    call printp('zcub',zd)
    write(6,*) 'before xytoiq'
    call printp('xg  ',xg)
    call printp('yg  ',yg)
  endif
#endif

  ! use 4* resolution grid il --> 4*il
  xg(1:ifull,k) = min(max(-.99999,xg(1:ifull,k)),.99999)
  yg(1:ifull,k) = min(max(-.99999,yg(1:ifull,k)),.99999)
  ! first guess for ri, rj and nearest i,j
  ri(1:ifull) = 1. + (1.+xg(1:ifull,k))*real(2*il_g)
  rj(1:ifull) = 1. + (1.+yg(1:ifull,k))*real(2*il_g)
  do loop = 1,nmaploop
    do iq = 1,ifull
      i = nint(ri(iq))
      j = nint(rj(iq))
      is = nint(sign(1.,ri(iq)-real(i)))
      js = nint(sign(1.,rj(iq)-real(j)))
      ! predict new value for ri, rj
      dxx = xx4(i+is,j) - xx4(i,j)
      dyx = xx4(i,j+js) - xx4(i,j)
      dxy = yy4(i+is,j) - yy4(i,j)
      dyy = yy4(i,j+js) - yy4(i,j)       
      den(iq) = dxx*dyy - dyx*dxy
      ri(iq) = real(i) + real(is)*real(((xg(iq,k)-xx4(i,j))*dyy-(yg(iq,k)-yy4(i,j))*dyx)/real(den(iq),8))
      rj(iq) = real(j) + real(js)*real(((yg(iq,k)-yy4(i,j))*dxx-(xg(iq,k)-xx4(i,j))*dxy)/real(den(iq),8))
    end do        
    ri(1:ifull) = min( ri(1:ifull), 1.+1.99999*real(2*il_g) )
    ri(1:ifull) = max( ri(1:ifull), 1.+0.00001*real(2*il_g) )
    rj(1:ifull) = min( rj(1:ifull), 1.+1.99999*real(2*il_g) )
    rj(1:ifull) = max( rj(1:ifull), 1.+0.00001*real(2*il_g) )
  end do  ! loop loop
  ! expect xg, yg to range between .5 and il+.5
  xg(1:ifull,k) = 0.25*(ri(1:ifull)+3.) - 0.5  ! -.5 for stag; back to normal ri, rj defn
  yg(1:ifull,k) = 0.25*(rj(1:ifull)+3.) - 0.5  ! -.5 for stag
  
end do

return
end subroutine toij5

#ifdef cray    
subroutine checkdiv(xstr,ystr,zstr)
! Check whether optimisation uses multiplication by reciprocal so
! that x/x /= 1.

use newmpar_m

implicit none

real, dimension(1:ifull) :: xstr,ystr,zstr
real, dimension(1:ifull) :: denxyz
integer, parameter :: n=100

call random_number(xstr(1:n))
ystr(1:n) = 0.9*xstr(1:n)
zstr(1:n) = 0.8*xstr(1:n)
! By construction here, xstr is largest, so xstr(iq)/denxyz should be
! 1.
denxyz(1:ifull)=max( abs(xstr(1:ifull)),abs(ystr(1:ifull)),abs(zstr(1:ifull)) )
xstr(1:ifull) = xstr(1:ifull)/denxyz(1:ifull)
ystr(1:ifull) = ystr(1:ifull)/denxyz(1:ifull)
zstr(1:ifull) = zstr(1:ifull)/denxyz(1:ifull)
if ( any(xstr(1:n)/=1.0) ) then
  write(6,*) "Error, must not use -Dcray on this machine"
  stop
end if
return
end subroutine checkdiv
#endif
