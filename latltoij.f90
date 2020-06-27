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
    
module latltoij_m

private
public latltoij

interface latltoij
  module procedure latltoij_s, latltoij_v
end interface

contains
    
subroutine latltoij_s(rlongin,rlatin,rlong0,rlat0,schmidt,xout,yout,nf,xx4,yy4,ik)
!     given a pair of latitudes and longitudes (in degrees),
!     returns i and j values on the conformal-cubic grid as
!     xout ranging between .5 and   il +.5, and
!     yout ranging between .5 and   il +.5
!     Note that the parallel version still returns xout, yout on the 
!     global grid.

!     modify for Cray; used by plotg.f and topgencc.f

use const_phys
use newmpar_m
use parm_m
use parmdyn_m
use utilities

implicit none

integer, intent(in) :: ik
integer, intent(out) :: nf
integer ig,jg,is,js,loop
integer, parameter :: nmaploop=3
#ifdef debug
integer, parameter :: ntest=0, numtst=-1
integer, save :: num=0
#endif
real, intent(in) :: rlongin, rlatin, rlong0, rlat0, schmidt
real, intent(out) :: xout, yout
real, dimension(3,3) :: rotpolei
real ri,rj,xa,ya,za,xgrid,ygrid,xx,yy,zz
real(kind=8), dimension(:,:), pointer :: xx4, yy4 ! avoid intent for pointers
real(kind=8) dxx,dyy,dxy,dyx,denxyz,x,y,z 
real(kind=8) alf, den

alf = (1._8-schmidt**2)/(1._8+schmidt**2)
rotpolei = transpose(calc_rotpole(rlong0,rlat0))

#ifdef debug
num=num+1
if (num<numtst) write(6,*) 'a rlongin,rlatin ',rlongin,rlatin
#endif

xa=cos(rlongin*pi/180.)*cos(rlatin*pi/180.)
ya=sin(rlongin*pi/180.)*cos(rlatin*pi/180.)
za=sin(rlatin*pi/180.)

#ifdef debug
if (num<numtst) write(6,*) 'b xa,ya,za ',xa,ya,za
#endif

x=rotpolei(1,1)*xa+rotpolei(1,2)*ya+rotpolei(1,3)*za
y=rotpolei(2,1)*xa+rotpolei(2,2)*ya+rotpolei(2,3)*za
z=rotpolei(3,1)*xa+rotpolei(3,2)*ya+rotpolei(3,3)*za

#ifdef debug
if (num<numtst) write(6,*) 'c x,y,z ',x,y,z
#endif

!     if necessary, transform physical (x, y, z) to equivalent coordinates
!     on regular gnomonic panels
x=x*(1._8-alf)/(real(schmidt,8)*(1._8-alf*z))
y=y*(1._8-alf)/(real(schmidt,8)*(1._8-alf*z))
z=(z-alf)/(1._8-alf*z)

denxyz=max( abs(x),abs(y),abs(z) )
xx=real(x/denxyz)
yy=real(y/denxyz)
zz=real(z/denxyz)
!       deduce corresponding face
!        if(ncray.eq.1)then
!         all these if statements are replaced by the subsequent cunning code
if(abs(abs(x )-denxyz)<1.e-30_8)then     ! Cray
  if(abs(x-denxyz)<1.e-30_8)then         ! Cray
    nf    =0                             ! Cray
    xgrid =       yy                     ! Cray
    ygrid =       zz                     ! Cray
  else                                   ! Cray
    nf    =3                             ! Cray
    xgrid =     -zz                      ! Cray
    ygrid =     -yy                      ! Cray
  endif                                  ! Cray
elseif(abs(abs(z )-denxyz)<1.e-30_8)then ! Cray
  if(abs(z-denxyz)<1.e-30_8)then         ! Cray
    nf    =1                             ! Cray
    xgrid =      yy                      ! Cray
    ygrid =     -xx                      ! Cray
  else                                   ! Cray
    nf    =4                             ! Cray
    xgrid =      xx                      ! Cray
    ygrid =     -yy                      ! Cray
  endif                                  ! Cray
else                                     ! Cray
  if(abs(y-denxyz)<1.e-30_8)then         ! Cray
    nf    =2                             ! Cray
    xgrid =     -zz                      ! Cray
    ygrid =     -xx                      ! Cray
  else                                   ! Cray
    nf    =5                             ! Cray
    xgrid =      xx                      ! Cray
    ygrid =      zz                      ! Cray
  endif                                  ! Cray
endif                                    ! Cray
!        else  ! e.g. ncray=0
!          nf=max( int(xx)*(3*int(xx)-3) ,   ! ** n=0,5 version
!     .            int(zz)*(5*int(zz)-3) ,
!     .            int(yy)*(7*int(yy)-3) )/2
!          xgrid=xgx(nf)*xx+xgy(nf)*yy+xgz(nf)*zz  ! -1 to 1
!          ygrid=ygx(nf)*xx+ygy(nf)*yy+ygz(nf)*zz
!        endif    ! (ncray.eq.1)

!       convert to grid point numbering
!       the xytoij routine follows

!       use 4* resolution grid il --> 4*il
xgrid=min(max(-.99999,xgrid),.99999)
ygrid=min(max(-.99999,ygrid),.99999)
!       first guess for ri, rj and nearest ig,jg
ri=1.+(1.+xgrid)*real(2*ik)
rj=1.+(1.+ygrid)*real(2*ik)
do loop = 1,nmaploop
  ig=nint(ri)
  jg=nint(rj)
  is=nint(sign(1.,ri-real(ig)))
  js=nint(sign(1.,rj-real(jg)))
  ! predict new value for ri, rj
  dxx=xx4(ig+is,jg)-xx4(ig,jg)
  dyx=xx4(ig,jg+js)-xx4(ig,jg)
  dxy=yy4(ig+is,jg)-yy4(ig,jg)
  dyy=yy4(ig,jg+js)-yy4(ig,jg)
  den = dxx*dyy - dyx*dxy
  ri=real(ig)+real(is)*real(((xgrid-xx4(ig,jg))*dyy-(ygrid-yy4(ig,jg))*dyx)/den)
  rj=real(jg)+real(js)*real(((ygrid-yy4(ig,jg))*dxx-(xgrid-xx4(ig,jg))*dxy)/den)
  ri=min(ri,1.0+1.99999*real(2*ik))
  ri=max(ri,1.0+0.00001*real(2*ik))
  rj=min(rj,1.0+1.99999*real(2*ik))
  rj=max(rj,1.0+0.00001*real(2*ik))
end do  ! loop loop
xout=.25*(ri+3.) -.5  ! -.5 for stag; back to normal ri, rj defn
yout=.25*(rj+3.) -.5  ! -.5 for stag
!       expect xout, yout (at this point) to range between .5 and il+.5

#ifdef debug
if(ntest.eq.1.and.rlongin.gt.43.1.and.rlongin.lt.49.9)then
  if(rlatin.gt.-24.2.and.rlatin.lt.-23.8)then
    write(6,*) 'lat,long,x,y,z,den ',rlatin,rlongin,x,y,z,denxyz
    write(6,*) 'xx,yy,zz ',xx,yy,zz
    write(6,*) 'nf,xout,yout,ri,rj ',nf,xout,yout,ri,rj
    write(6,*) 'youtb ',yout+nf*ik
  endif
endif
#endif

return
end subroutine latltoij_s

subroutine latltoij_v(rlongin,rlatin,rlong0,rlat0,schmidt,xout,yout,nf,xx4,yy4,ik)
!     given a pair of latitudes and longitudes (in degrees),
!     returns i and j values on the conformal-cubic grid as
!     xout ranging between .5 and   il +.5, and
!     yout ranging between .5 and   il +.5
!     Note that the parallel version still returns xout, yout on the 
!     global grid.

!     modify for Cray; used by plotg.f and topgencc.f

use const_phys
use newmpar_m
use parm_m
use parmdyn_m
use utilities

implicit none

integer, intent(in) :: ik
integer, dimension(:), intent(out) :: nf
integer loop, ig, jg, is, js, iq
integer, parameter :: nmaploop=3
#ifdef debug
integer, parameter :: ntest=0, numtst=-1
integer, save :: num=0
#endif
real, dimension(:), intent(in) :: rlongin, rlatin
real, intent(in) :: rlong0, rlat0, schmidt
real, dimension(:), intent(out) :: xout, yout
real, dimension(3,3) :: rotpolei
real, dimension(size(nf)) :: xa, ya, za, xx, yy, zz, xgrid, ygrid
real, dimension(size(nf)) :: ri, rj
real(kind=8), dimension(:,:), pointer :: xx4, yy4 ! avoid intent for pointers
real(kind=8) dxx,dyy,dxy,dyx
real(kind=8), dimension(size(nf)) :: x, y, z, denxyz
real(kind=8) alf, den

alf = (1._8-schmidt**2)/(1._8+schmidt**2)
rotpolei = transpose(calc_rotpole(rlong0,rlat0))

#ifdef debug
num=num+1
if (num<numtst) write(6,*) 'a rlongin,rlatin ',rlongin(1),rlatin(1)
#endif

xa=cos(rlongin*pi/180.)*cos(rlatin*pi/180.)
ya=sin(rlongin*pi/180.)*cos(rlatin*pi/180.)
za=sin(rlatin*pi/180.)

#ifdef debug
if (num<numtst) write(6,*) 'b xa,ya,za ',xa(1),ya(1),za(1)
#endif

x=rotpolei(1,1)*xa+rotpolei(1,2)*ya+rotpolei(1,3)*za
y=rotpolei(2,1)*xa+rotpolei(2,2)*ya+rotpolei(2,3)*za
z=rotpolei(3,1)*xa+rotpolei(3,2)*ya+rotpolei(3,3)*za

#ifdef debug
if (num<numtst) write(6,*) 'c x,y,z ',x(1),y(1),z(1)
#endif

!     if necessary, transform physical (x, y, z) to equivalent coordinates
!     on regular gnomonic panels
x=x*(1._8-alf)/(real(schmidt,8)*(1._8-alf*z))
y=y*(1._8-alf)/(real(schmidt,8)*(1._8-alf*z))
z=(z-alf)/(1._8-alf*z)

denxyz=max( abs(x),abs(y),abs(z) )
xx=real(x/denxyz)
yy=real(y/denxyz)
zz=real(z/denxyz)
!       deduce corresponding face
!        if(ncray.eq.1)then
!         all these if statements are replaced by the subsequent cunning code
where ( abs(abs(x)-denxyz)<1.e-30_8 .and. abs(x-denxyz)<1.e-30_8 )     ! Cray
  nf    =0                                                             ! Cray
  xgrid =       yy                                                     ! Cray
  ygrid =       zz                                                     ! Cray
elsewhere ( abs(abs(x)-denxyz)<1.e-30_8 )                              ! Cray  
  nf    =3                                                             ! Cray
  xgrid =     -zz                                                      ! Cray
  ygrid =     -yy                                                      ! Cray
elsewhere ( abs(abs(z)-denxyz)<1.e-30_8 .and. abs(z-denxyz)<1.e-30_8 ) ! Cray
  nf    =1                                                             ! Cray
  xgrid =      yy                                                      ! Cray
  ygrid =     -xx                                                      ! Cray
elsewhere ( abs(abs(z)-denxyz)<1.e-30_8 )                              ! Cray
  nf    =4                                                             ! Cray
  xgrid =      xx                                                      ! Cray
  ygrid =     -yy                                                      ! Cray
elsewhere ( abs(y-denxyz)<1.e-30_8 )                                   ! Cray
  nf    =2                                                             ! Cray
  xgrid =     -zz                                                      ! Cray
  ygrid =     -xx                                                      ! Cray
elsewhere                                                              ! Cray
  nf    =5                                                             ! Cray
  xgrid =      xx                                                      ! Cray
  ygrid =      zz                                                      ! Cray
end where                                                              ! Cray

!       convert to grid point numbering
!       the xytoij routine follows

!       use 4* resolution grid il --> 4*il
xgrid = min(max(-.99999,xgrid),.99999)
ygrid = min(max(-.99999,ygrid),.99999)
!       first guess for ri, rj and nearest ig,jg
ri = 1. + (1.+xgrid)*real(2*ik)
rj = 1. + (1.+ygrid)*real(2*ik)
do loop = 1,nmaploop
  do iq = 1,size(nf)
    ig = nint(ri(iq))
    jg = nint(rj(iq))
    is = nint(sign(1.,ri(iq)-real(ig)))
    js = nint(sign(1.,rj(iq)-real(jg)))
    ! predict new value for ri, rj
    dxx = xx4(ig+is,jg) - xx4(ig,jg)
    dyx = xx4(ig,jg+js) - xx4(ig,jg)
    dxy = yy4(ig+is,jg) - yy4(ig,jg)
    dyy = yy4(ig,jg+js) - yy4(ig,jg)
    den = dxx*dyy - dyx*dxy
    ri(iq) = real(ig) + real(is)*real(((xgrid(iq)-xx4(ig,jg))*dyy-      &
                                       (ygrid(iq)-yy4(ig,jg))*dyx)/den)
    rj(iq) = real(jg) + real(js)*real(((ygrid(iq)-yy4(ig,jg))*dxx-      &
                                       (xgrid(iq)-xx4(ig,jg))*dxy)/den)
  end do  ! loop iq
  ri = min(ri, 1.0+1.99999*real(2*ik))
  ri = max(ri, 1.0+0.00001*real(2*ik))
  rj = min(rj, 1.0+1.99999*real(2*ik))
  rj = max(rj, 1.0+0.00001*real(2*ik))
end do   ! loop loop
xout = .25*(ri+3.) - .5  ! -.5 for stag; back to normal ri, rj defn
yout = .25*(rj+3.) - .5  ! -.5 for stag
!       expect xout, yout (at this point) to range between .5 and il+.5

#ifdef debug
if(ntest.eq.1.and.rlongin(1).gt.43.1.and.rlongin(1).lt.49.9)then
  if(rlatin(1).gt.-24.2.and.rlatin(1).lt.-23.8)then
    write(6,*) 'lat,long,x,y,z,den ',rlatin(1),rlongin(1),x(1),y(1),z(1),denxyz(1)
    write(6,*) 'xx,yy,zz ',xx(1),yy(1),zz(1)
    write(6,*) 'nf,xout,yout,ri,rj ',nf(1),xout(1),yout(1),ri(1),rj(1)
    write(6,*) 'youtb ',yout(1)+nf(1)*ik
  endif
endif
#endif

return
end subroutine latltoij_v

end module latltoij_m