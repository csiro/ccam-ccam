! Conformal Cubic Atmospheric Model
    
! Copyright 2015-2021 Commonwealth Scientific Industrial Research Organisation (CSIRO)
    
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

use bigxy4_m
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
real(kind=8), dimension(ifull,kl), intent(out) :: x3d, y3d, z3d   ! upglobal depts 
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
  s(1:ifull,k,1) = uc(1:ifull,k)
  s(1:ifull,k,2) = vc(1:ifull,k)
  s(1:ifull,k,3) = wc(1:ifull,k)
end do  

call bounds(s,nrows=2)

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
else
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

end if

#ifdef GPU
!$omp target data map(to:sx,xx4,yy4)
#endif
!$acc data create(sx,xx4,yy4)
!$acc update device(sx,xx4,yy4)

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

!======================== start of intsch=1 section ====================
if ( intsch==1 ) then  
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
  
#ifdef _OPENMP
#ifdef GPU
  !$omp target teams distribute parallel do collapse(3) schedule(static),               &
  !$omp map(to:xg,yg,nface) map(from:s) private(nn,k,iq,idel,xxg,jdel,yyg,n),           &
  !$omp private(cmul_1,cmul_2,cmul_3,cmul_4,dmul_2,dmul_3,emul_1,emul_2,emul_3,emul_4), &
  !$omp private(rmul_1,rmul_2,rmul_3,rmul_4)
#else
  !$omp parallel do collapse(2) schedule(static) private(nn,k,iq,idel,xxg,jdel,yyg,n),  &
  !$omp private(cmul_1,cmul_2,cmul_3,cmul_4,dmul_2,dmul_3,emul_1,emul_2,emul_3,emul_4), &
  !$omp private(rmul_1,rmul_2,rmul_3,rmul_4)
#endif
#else
  !$acc parallel loop collapse(3) copyin(xg,yg,nface) copyout(s) present(sx)
#endif
  do nn = 1,3
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
        s(iq,k,nn) = rmul_1*emul_1 + rmul_2*emul_2 + rmul_3*emul_3 + rmul_4*emul_4
      end do   ! iq loop
    end do     ! k loop
  end do       ! nn loop
#ifdef _OPENMP
#ifdef GPU
  !$omp end target teams distribute parallel do
#else
  !$omp end parallel do
#endif
#else
  !$acc end parallel loop
#endif
            
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

#ifdef _OPENMP
#ifdef GPU
  !$omp target teams distribute parallel do collapse(3) schedule(static),               &
  !$omp map(to:xg,yg,nface) map(from:s) private(nn,k,iq,idel,xxg,jdel,yyg,n),           &
  !$omp private(cmul_1,cmul_2,cmul_3,cmul_4,dmul_2,dmul_3,emul_1,emul_2,emul_3,emul_4), &
  !$omp private(rmul_1,rmul_2,rmul_3,rmul_4)
#else
  !$omp parallel do collapse(2) schedule(static) private(nn,k,iq,idel,xxg,jdel,yyg,n),  &
  !$omp private(cmul_1,cmul_2,cmul_3,cmul_4,dmul_2,dmul_3,emul_1,emul_2,emul_3,emul_4), &
  !$omp private(rmul_1,rmul_2,rmul_3,rmul_4)
#endif
#else
  !$acc parallel loop collapse(3) copyin(xg,yg,nface) copyout(s) present(sx)
#endif
  do nn = 1,3
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
#ifdef _OPENMP
#ifdef GPU
  !$omp end target teams distribute parallel do
#else
  !$omp end parallel do
#endif
#else
  !$acc end parallel loop
#endif
  
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

if ( diag .and. mydiag ) then
  write(6,*) '2nd guess for k = ',nlv
  write(6,*) 'x3d,y3d,z3d ',x3d(idjd,nlv),y3d(idjd,nlv),z3d(idjd,nlv)
  write(6,*) 'xg,yg,nface ',xg(idjd,nlv),yg(idjd,nlv),nface(idjd,nlv)
end if

!======================== start of intsch=1 section ====================
if ( intsch==1 ) then

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

#ifdef _OPENMP
#ifdef GPU
  !$omp target teams distribute parallel do collapse(3) schedule(static),               &
  !$omp map(to:xg,yg,nface) map(from:s) private(nn,k,iq,idel,xxg,jdel,yyg,n),           &
  !$omp private(cmul_1,cmul_2,cmul_3,cmul_4,dmul_2,dmul_3,emul_1,emul_2,emul_3,emul_4), &
  !$omp private(rmul_1,rmul_2,rmul_3,rmul_4)
#else
  !$omp parallel do collapse(2) schedule(static) private(nn,k,iq,idel,xxg,jdel,yyg,n),  &
  !$omp private(cmul_1,cmul_2,cmul_3,cmul_4,dmul_2,dmul_3,emul_1,emul_2,emul_3,emul_4), &
  !$omp private(rmul_1,rmul_2,rmul_3,rmul_4)
#endif
#else
  !$acc parallel loop collapse(3) copyin(xg,yg,nface) copyout(s) present(sx)
#endif
  do nn = 1,3
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
        s(iq,k,nn) = rmul_1*emul_1 + rmul_2*emul_2 + rmul_3*emul_3 + rmul_4*emul_4
      end do   ! iq loop
    end do     ! k loop
  end do       ! nn loop
#ifdef _OPENMP
#ifdef GPU
  !$omp end target teams distribute parallel do
#else
  !$omp end parallel do
#endif
#else
  !$acc end parallel loop
#endif
            
!========================   end of intsch=1 section ====================
else     ! if(intsch==1)then
!======================== start of intsch=2 section ====================

  ! For other processes
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

#ifdef _OPENMP
#ifdef GPU
  !$omp target teams distribute parallel do collapse(3) schedule(static),               &
  !$omp map(to:xg,yg,nface) map(from:s) private(nn,k,iq,idel,xxg,jdel,yyg,n),           &
  !$omp private(cmul_1,cmul_2,cmul_3,cmul_4,dmul_2,dmul_3,emul_1,emul_2,emul_3,emul_4), &
  !$omp private(rmul_1,rmul_2,rmul_3,rmul_4)
#else
  !$omp parallel do collapse(2) schedule(static) private(nn,k,iq,idel,xxg,jdel,yyg,n),  &
  !$omp private(cmul_1,cmul_2,cmul_3,cmul_4,dmul_2,dmul_3,emul_1,emul_2,emul_3,emul_4), &
  !$omp private(rmul_1,rmul_2,rmul_3,rmul_4)
#endif
#else
  !$acc parallel loop collapse(3) copyin(xg,yg,nface) copyout(s) present(sx)
#endif
  do nn = 1,3
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
        s(iq,k,nn) = rmul_1*emul_1 + rmul_2*emul_2 + rmul_3*emul_3 + rmul_4*emul_4
      end do          ! iq loop
    end do            ! nn loop
  end do              ! k loop
#ifdef _OPENMP
#ifdef GPU
  !$omp end target teams distribute parallel do
#else
  !$omp end parallel do
#endif
#else
  !$acc end parallel loop
#endif

endif                     ! (intsch==1) .. else ..
!========================   end of intsch=1 section ====================

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

#ifdef GPU
!$omp end target data
#endif
!$acc end data

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

integer, parameter :: nmaploop = 3
integer k, iq, loop, i, j, is, js
real ri,rj
real xstr,ystr,zstr
real denxyz,xd,yd,zd
real(kind=8) alf,alfonsch  ! 6/11/07 esp for 200m
real(kind=8) dxx,dxy,dyx,dyy
real(kind=8), dimension(ifull,kl), intent(in) :: x3d,y3d,z3d
real(kind=8) den
logical xytest, xztest, yztest

! if necessary, transform (x3d, y3d, z3d) to equivalent
! coordinates (xstr, ystr, zstr) on regular gnomonic panels
alf      = (1._8-real(schmidt,8)**2)/(1._8+real(schmidt,8)**2)
alfonsch = 2._8*real(schmidt,8)/(1._8+real(schmidt,8)**2)  ! same but bit more accurate

#ifdef _OPENMP
#ifdef GPU
!$omp target teams distribute parallel do collapse(2) schedule(static),             &
!$omp map(to:x3d,y3d,z3d) map(from:xg,yg,nface) private(k,iq,den),                  &
!$omp private(xstr,ystr,zstr,denxyz,xd,yd,zd,xytest,xztest,yztest,ri,rj,loop,i,j),  &
!$omp private(is,js,dxx,dyx,dxy,dyy)
#else
!$omp parallel do collapse(2) schedule(static) private(k,iq,den),                   &
!$omp private(xstr,ystr,zstr,denxyz,xd,yd,zd,xytest,xztest,yztest,ri,rj,loop,i,j),  &
!$omp private(is,js,dxx,dyx,dxy,dyy)
#endif
#else
!$acc parallel loop collapse(2) copyin(x3d,y3d,z3d) copyout(xg,yg,nface) present(xx4,yy4)
#endif
do k = 1,kl
  do iq = 1,ifull    
    den  = 1._8 - alf*z3d(iq,k)
    xstr = real(x3d(iq,k)*(alfonsch/den))
    ystr = real(y3d(iq,k)*(alfonsch/den))
    zstr = real(   (z3d(iq,k)-alf)/den)

!      first deduce departure faces
!      instead calculate cubic coordinates
!      The faces are:
!      0: X=1   1: Z=1   2: Y=1   3: X=-1   4: Z=-1   5: Y=-1

    denxyz = max(abs(xstr),abs(ystr),abs(zstr) )
    xd = xstr/denxyz
    yd = ystr/denxyz
    zd = zstr/denxyz

    xytest = abs(xd)>abs(yd)
    xztest = abs(xd)>abs(zd)
    yztest = abs(yd)>abs(zd)
    if ( xytest .and. xztest .and. xd>0. ) then
      nface(iq,k)    =0
      xg(iq,k) =       yd
      yg(iq,k) =       zd
    else if ( xytest .and. xztest ) then
      nface(iq,k)    =3
      xg(iq,k) =     -zd
      yg(iq,k) =     -yd
    else if ( yztest .and. yd>0. ) then
      nface(iq,k)    =2
      xg(iq,k) =     -zd
      yg(iq,k) =     -xd
    else if ( yztest ) then
      nface(iq,k)    =5
      xg(iq,k) =      xd
      yg(iq,k) =      zd
    else if ( zd>0. ) then
      nface(iq,k)    =1
      xg(iq,k) =      yd
      yg(iq,k) =     -xd
    else
      nface(iq,k)    =4
      xg(iq,k) =      xd
      yg(iq,k) =     -yd
    end if

#ifdef debug
    if(ntest==1.and.k==nlv.and.iq=idjd)then
      iq=idjd
      write(6,*) 'x3d,y3d,z3d ',x3d(iq),y3d(iq),z3d(iq)
      den(iq)=1._8-alf*z3d(iq) ! to force real*8
      write(6,*) 'den ',den
      denxyz=max( abs(xstr),abs(ystr),abs(zstr) )
      xd=xstr/denxyz
      yd=ystr/denxyz
      zd=zstr/denxyz
      write(6,*) 'k,xstr,ystr,zstr,denxyz ',k,xstr,ystr,zstr,denxyz
      write(6,*) 'abs(xstr,ystr,zstr) ',abs(xstr),abs(ystr),abs(zstr)
      write(6,*) 'xd,yd,zd,nface ',xd,yd,zd,nface(iq,k)
      write(6,*) 'alf,alfonsch ',alf,alfonsch
    endif
    if(ndiag==2)then
      write(6,*) 'before xytoiq'
      call printp('xg  ',xg)
      call printp('yg  ',yg)
    endif
#endif

    ! use 4* resolution grid il --> 4*il
    xg(iq,k) = min(max(-.99999,xg(iq,k)),.99999)
    yg(iq,k) = min(max(-.99999,yg(iq,k)),.99999)
    ! first guess for ri, rj and nearest i,j
    ri = 1. + (1.+xg(iq,k))*real(2*il_g)
    rj = 1. + (1.+yg(iq,k))*real(2*il_g)
    do loop = 1,nmaploop
      i = nint(ri)
      j = nint(rj)
      is = nint(sign(1.,ri-real(i)))
      js = nint(sign(1.,rj-real(j)))
      ! predict new value for ri, rj
      dxx = xx4(i+is,j) - xx4(i,j)
      dyx = xx4(i,j+js) - xx4(i,j)
      dxy = yy4(i+is,j) - yy4(i,j)
      dyy = yy4(i,j+js) - yy4(i,j)       
      den = dxx*dyy - dyx*dxy
      ri = real(i) + real(is)*real(((xg(iq,k)-xx4(i,j))*dyy-(yg(iq,k)-yy4(i,j))*dyx)/real(den,8))
      rj = real(j) + real(js)*real(((yg(iq,k)-yy4(i,j))*dxx-(xg(iq,k)-xx4(i,j))*dxy)/real(den,8))        
      ri = min( ri, 1.+1.99999*real(2*il_g) )
      ri = max( ri, 1.+0.00001*real(2*il_g) )
      rj = min( rj, 1.+1.99999*real(2*il_g) )
      rj = max( rj, 1.+0.00001*real(2*il_g) )
    end do  ! loop loop
    ! expect xg, yg to range between .5 and il+.5
    xg(iq,k) = 0.25*(ri+3.) - 0.5  ! -.5 for stag; back to normal ri, rj defn
    yg(iq,k) = 0.25*(rj+3.) - 0.5  ! -.5 for stag  
  end do
end do
#ifdef _OPENMP
#ifdef GPU
!$omp end target teams distribute parallel do
#else
!$omp end parallel do
#endif
#else
!$acc end parallel loop
#endif

return
end subroutine toij5
