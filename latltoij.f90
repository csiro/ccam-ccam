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

contains
    
subroutine latltoij(rlongin,rlatin,rlong0,rlat0,schmidt,xout,yout,nf,xx4,yy4,ik)
!     given a pair of latitudes and longitudes (in degrees),
!     returns i and j values on the conformal-cubic grid as
!     xout ranging between .5 and   il +.5, and
!     yout ranging between .5 and   il +.5
!     Note that the parallel version still returns xout, yout on the 
!     global grid.

!     modify for Cray; used by plotg.f and topgencc.f

use utilities

implicit none

include 'newmpar.h'
include 'const_phys.h'
include 'parm.h'
include 'parmdyn.h'

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
real ri,rj,xa,ya,za,xgrid,ygrid,xx,yy,zz,x1,z1
real(kind=8), dimension(:,:), pointer :: xx4, yy4 ! avoid intent for pointers
real(kind=8) dxx,dyy,dxy,dyx,denxyz,x,y,z 
real(kind=8) alf, den
real(kind=8), parameter :: one=1._8

alf=(one-schmidt**2)/(one+schmidt**2)
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
x1=real(x)
z1=real(z)
x=x*(1.-alf)/(schmidt*(1.-alf*z))
y=y*(1.-alf)/(schmidt*(1.-alf*z))
z=(z-alf)/(1.-alf*z)

#ifdef debug
if(ntest.eq.1.and.z1.gt..82.and.z1.lt..821)then
  write(6,*) 'latltoij: rlongin, rlatin ',rlongin, rlatin
  write(6,*) 'latltoij: xa,ya,za ',xa,ya,za
  write(6,*) 'latltoij: x1,z1 ',x1,z1
  write(6,*) 'latltoij: x,y,z ',x,y,z
endif
#endif

denxyz=max( abs(x),abs(y),abs(z) )
xx=real(x/denxyz)
yy=real(y/denxyz)
zz=real(z/denxyz)
!       deduce corresponding face
!        if(ncray.eq.1)then
!         all these if statements are replaced by the subsequent cunning code
if(abs(x ).eq.denxyz)then             ! Cray
  if(x .eq.denxyz)then                ! Cray
    nf    =0                          ! Cray
    xgrid =       yy                  ! Cray
    ygrid =       zz                  ! Cray
  else                                ! Cray
    nf    =3                          ! Cray
    xgrid =     -zz                   ! Cray
    ygrid =     -yy                   ! Cray
  endif                               ! Cray
elseif(abs(z ).eq.denxyz)then         ! Cray
  if(z .eq.denxyz)then                ! Cray
    nf    =1                          ! Cray
    xgrid =      yy                   ! Cray
    ygrid =     -xx                   ! Cray
  else                                ! Cray
    nf    =4                          ! Cray
    xgrid =      xx                   ! Cray
    ygrid =     -yy                   ! Cray
  endif                               ! Cray
else                                  ! Cray
  if(y .eq.denxyz)then                ! Cray
    nf    =2                          ! Cray
    xgrid =     -zz                   ! Cray
    ygrid =     -xx                   ! Cray
  else                                ! Cray
    nf    =5                          ! Cray
    xgrid =      xx                   ! Cray
    ygrid =      zz                   ! Cray
  endif                               ! Cray
endif                                 ! Cray
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
end subroutine latltoij

end module latltoij_m
